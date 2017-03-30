#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <iostream>
#include "Quadtree.cpp"

using namespace std;

typedef vector< particle_t* > gridBin;
typedef vector< gridBin > particleGrid;

int navg, nabsavg = 0;
double davg, dmin, absmin = 1.0, absavg = 0.0;


//
//create bins with length of cutoff
//
int create_bins( particleGrid &bins, particle_t* particles, int n ) 
{
	int binsPerRow = ( int ) ceil( sqrt( density * n ) / ( 2 * cutoff ));

	//resize the vector to the exact size of bins.
	bins.resize( binsPerRow * binsPerRow );

	//put particles in bins according to their locations
	for ( int j = 0; j < n; j++ ) 
	{
		int x = ( int ) floor( particles[ j ].x / binsize );
		int y = ( int ) floor( particles[ j ].y / binsize );

		bins[ x + y * binsPerRow ].push_back( &particles[ j ] );
	}

	return binsPerRow;
}

void compute_forces( particleGrid &bins, int numOfBins, int binsPerRow )
{
	for ( int i = 0; i < numOfBins; i++ )
	{
		gridBin bin1 = bins[ i ];
		int numOfParticles = bin1.size( );

		for ( int j = 0; j < numOfParticles; j++ ) 
		{
			bin1[ j ]->ax = bin1[ j ]->ay = 0;
		}


		//Search within current bin and a 'halo region'
		//Halo region may be defined as the region around the bin of size 'cutoff'
		//For ease in computation, the neighboring bins may be used as the halo.

		int location = i;

		vector< int > x;
		vector< int > y;

		x.push_back( 0 );
		y.push_back( 0 );

		if ( location >= binsPerRow ) 
		{
			y.push_back( -1 );
		}
		if ( location < binsPerRow*( binsPerRow - 1 ) ) 
		{
			y.push_back( 1 );
		}
		if ( location % binsPerRow != 0 ) 
		{
			x.push_back( -1 );
		}
		if ( location % binsPerRow != binsPerRow - 1 ) 
		{
			x.push_back( 1 );
		}

		for ( int a = 0; a < x.size( ); a++ ) 
		{
			for ( int b = 0; b < y.size( ); b++ ) 
			{
				int bin_num = i + x[ a ] + binsPerRow * y[ b ];


				for ( int c = 0; c < bins[ i ].size( ); c++ ) 
				{
					for ( int d = 0; d < bins[ bin_num ].size( ); d++ ) 
					{

						apply_force( *bins[ i ][ c ], *bins[ bin_num ][ d ], &dmin, &davg, &navg );

					}
				}
			}
		}

	}
}

void relocate( particleGrid &bins, int binsPerRow, particle_t* particles, int n )
{
	//put particles in bins according to their locations
	for( int j = 0; j < n; j++ )
	{
		int x = floor( particles[ j ].x/cutoff );
		int y = floor( particles[ j ].y/cutoff );
		bins[ x + binsPerRow * y ].push_back( particles + j );
	}
}

void clean_bins( particleGrid &bins, int numOfBins )
{
	//clean the bins
	for( int i = 0; i < numOfBins; i++ )
	{
		bins[ i ].clear( );
	}
}

void move_particles( particle_t* particles, particleGrid &bins, gridBin &displacedParticles, int n, int numOfBins, int binsPerRow )
{
	double binsizes = cutoff * 2;

	for ( int b = 0; b < numOfBins; b++ )
	{
		int size = bins[ b ].size( );
		for ( int p = 0; p < size; ) 
		{
			move( *bins[ b ][ p ] );

			int x = ( int ) floor( bins[ b ][ p ]->x / binsizes );
			int y = ( int ) floor( bins[ b ][ p ]->y / binsizes );
			
			if ( y * binsPerRow + x != b )
			{
				displacedParticles.push_back( bins[ b ][ p ] );
				bins[ b ].erase( bins[ b ].begin( ) + p );
				size--;
			}
			else
			{
				p++;
			}
		}
	}

	for( int i = 0; i < displacedParticles.size( ); i++ )
	{
		int x = ( int ) floor( displacedParticles[ i ]->x / binsize );
		int y = ( int ) floor(displacedParticles[ i ]->y / binsize );

		bins[ x + y * binsPerRow ].push_back( displacedParticles[ i ] );
	}
	displacedParticles.clear( );
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{
	if ( find_option( argc, argv, "-h" ) >= 0 )
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
	FILE *fsum = sumname ? fopen( sumname, "a" ) : NULL;

	particle_t *particles;

	if ( ( particles = ( particle_t* ) malloc( n * sizeof( particle_t ) ) ) == NULL ) 
	{
		fprintf( stderr,"particles malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
	
		return 0;
	}

	set_size( n );
	init_particles( n, particles );

	//create the bins to contain the particles
	particleGrid bins;
	gridBin displacedParticles;

	int binsPerRow = create_bins( bins, particles, n );
	int numOfBins = binsPerRow * binsPerRow;

	//
	//  simulate a number of time steps
	//
	double simulation_time = read_timer( );

	for ( int step = 0; step < NSTEPS; step++ )
	{
		navg = 0;
		davg = 0.0;
		dmin = 1.0;
		
		//
		//  compute forces
		//
		compute_forces( bins, numOfBins, binsPerRow );

		//
		//  move particles
		//  The particles must also be moved between bins as necessary.
		//
		move_particles( particles, bins, displacedParticles, n, numOfBins, binsPerRow );

		if ( find_option( argc, argv, "-no" ) == -1 )
		{
			//
			// Computing statistical data
			//
			if ( navg ) 
			{
				absavg += davg / navg;
				nabsavg++;
			}
			
			if ( dmin < absmin ) absmin = dmin;

			//
			//  save if necessary
			//
			if ( fsave && ( step%SAVEFREQ ) == 0 )
			{
				save( fsave, n, particles );
			}
		}

	}
	simulation_time = read_timer( ) - simulation_time;

	printf( "n = %d, simulation time = %g seconds", n, simulation_time );

	if ( find_option( argc, argv, "-no" ) == -1 )
	{
		if ( nabsavg ) absavg /= nabsavg;
		
		//
		//  -the minimum distance absmin between 2 particles during the run of the simulation
		//  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
		//  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
		//
		//  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
		//
		printf( ", absmin = %lf, absavg = %lf", absmin, absavg );
		if ( absmin < 0.4 ) printf( "\nThe minimum distance is below 0.4 meaning that some particle is not interacting" );
		if ( absavg < 0.8 ) printf( "\nThe average distance is below 0.8 meaning that most particles are not interacting" );
	}

	printf( "\n" );

	//
	// Printing summary data
	//
	if ( fsum )
	{
		fprintf( fsum, "%d %g\n", n, simulation_time );
	}
	//
	// Clearing space
	//
	if ( fsum )
	{
		fclose( fsum );
	}
	
	free( particles );
	
	if ( fsave )
	{
		fclose( fsave );
	}

	return 0;
}
