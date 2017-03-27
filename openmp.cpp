#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "omp.h"
#include "common.h"
#include "Quadtree.cpp"

using namespace std;

double sizeOfBin;
int binsPerRow;
typedef vector< particle_t > gridBin;
typedef vector< gridBin > particleGrid;
int navg;
double dmin;
double davg;



void buildBins( particleGrid& bins, particle_t* particles, int n )
{
    sizeOfBin = cutoff * 2;
    binsPerRow = floor( sqrt( n*density ) / sizeOfBin ) + 1; // Should be around sqrt(N/2)

    bins.resize( binsPerRow * binsPerRow );

    for ( int i = 0; i < n; i++ )
    {
        int x = int( particles[ i ].x / sizeOfBin );
        int y = int( particles[ i ].y / sizeOfBin );
        bins[ x * binsPerRow + y].push_back( particles[ i ] );
    }
}

void relocate( int numthreads, particleGrid& displacedParticles, particleGrid& binsOfParticles )
{
	for( int j = 0; j < numthreads; j++ )
	{
		gridBin& bin = displacedParticles[ j ];

		for ( int i = 0; i < bin.size( ); i++ )  // Put them into the new bin
		{
			int x = ( int ) floor( bin[ i ].x / sizeOfBin );
			int y = ( int ) floor( bin[ i ].y / sizeOfBin );

			//If using multiple threads, below is critical area
			binsOfParticles[ x * binsPerRow + y ].push_back( bin[ i ] );
		}
	}
}

void compute_forces( particleGrid& binsOfParticles )//, int navg, double davg, double dmin)
{
	#pragma omp for reduction (+:navg) reduction(+:davg)
	for ( int i = 0; i < binsPerRow; i++ )
	{
		for ( int j = 0; j < binsPerRow; j++ )
		{
			gridBin& bin1 = binsOfParticles[ i * binsPerRow + j ];

			for ( int k = 0; k < bin1.size( ); k++ )
			{
				bin1[ k ].ax = bin1[ k ].ay = 0;
			}
			
			for ( int ii = -1; ii <= 1; ii++ )  
			{
				for ( int jj = -1; jj <= 1; jj++ )
				{
					if( ( i + ii >= 0 ) && ( i + ii < binsPerRow ) && ( j + jj >= 0 ) && ( j + jj < binsPerRow ) )
					{
						gridBin& bin2 = binsOfParticles[ ( i + ii ) * binsPerRow + j + jj ];

						for ( int k = 0; k < bin1.size( ); k++ )
						{
							for ( int l = 0; l < bin2.size( ); l++ )
							{
								apply_force( bin1[ k ], bin2[ l ], &dmin, &davg, &navg );
							}
						}
					}

				}
			}
		}
	}
}

void move_particles( particleGrid& binsOfParticles, gridBin bin )
{
	#pragma omp for
	for ( int i = 0; i < binsPerRow; i++ )
	{
		for( int j = 0; j < binsPerRow; j++ )
		{
			gridBin& vec = binsOfParticles[ i * binsPerRow + j ];
			int tail = vec.size( ), k = 0;
			
			for(; k < tail; )
			{
				move( vec[ k ] );
				
				int x = ( int ) floor( vec[ k ].x / sizeOfBin );  //Check the position
				int y = ( int ) floor( vec[ k ].y / sizeOfBin );
				
				if ( x == i && y == j )  // Still inside original bin
				{
					k++;
				}
				else
				{
					bin.push_back( vec[ k ] );  // Store paricles that have changed bin.
					vec[ k ] = vec[ --tail ]; //Remove it from the current bin.
				}
			}
			
			vec.resize( k );
		}
	}
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int nabsavg = 0;
    double absmin = 1.0, absavg = 0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n" );
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles; 

    if ( ( particles = ( particle_t* ) malloc( n * sizeof( particle_t ) ) ) == NULL )
    {
	fprintf( stderr, "particles malloc NULL\n");
	return 0;	
    }

    particleGrid binsOfParticles;
    particleGrid displacedParticles;

    set_size( n );
    init_particles( n, particles );

    buildBins(binsOfParticles, particles, n);
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    
    int numthreads;
    
    #pragma omp parallel private(dmin)
    {
        #pragma omp master
        {
            numthreads = omp_get_num_threads( );
            displacedParticles.resize( numthreads );
        }

        for ( int step = 0; step < NSTEPS; step++ )
        {
		navg = 0;
		davg = 0.0;
		dmin = 1.0;

		compute_forces( binsOfParticles );//, navg, davg, dmin);

		int id = omp_get_thread_num( );  // Each thread has a seperate bin vector
		gridBin& bin = displacedParticles[ id ];
		bin.clear( );

		move_particles( binsOfParticles, bin );


		//Scan over all bin vectors using one threads...
		//Using multiple threads with critical area seems to decrease performance
		#pragma omp master
		{
			relocate( numthreads, displacedParticles, binsOfParticles );
		}


		if( find_option( argc, argv, "-no" ) == -1 )
		{
			#pragma omp master
			if ( navg )
			{
				absavg +=  davg / navg;
				nabsavg++;
			}

			#pragma omp critical
			if ( dmin < absmin )
			{
				absmin = dmin;
			}

			#pragma omp master
			if( fsave && ( step % SAVEFREQ ) == 0 )
			{
				save( fsave, n, particles );
			}
		}
	}

	#pragma omp barrier // Wait for all threads (mostly the master) to finish
    }
    
    simulation_time = read_timer( ) - simulation_time;

    printf( "NThreads: %d, n = %d, simulation time = %g seconds", numthreads,n, simulation_time );

    if( find_option( argc, argv, "-no" ) == -1 )
    {
        if ( nabsavg ) absavg /= nabsavg;
        //
        //  -The minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf( ", absmin = %lf, absavg = %lf", absmin, absavg );
        if ( absmin < 0.4 ) printf ( "\nThe minimum distance is below 0.4 meaning that some particle is not interacting" );
        if ( absavg < 0.8 ) printf ( "\nThe average distance is below 0.8 meaning that most particles are not interacting" );
    }
    printf( "\n" );

    //
    // Printing summary data
    //
    if( fsum )
    {
        fprintf( fsum, "%d %g\n", n, simulation_time);
    }

    //
    // Clearing space
    //
    if( fsum )
    {
        fclose( fsum );
    }

    free( particles );

    if( fsave )
    {
        fclose( fsave );
    }
 
    return 0;
}
