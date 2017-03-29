#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <iostream>
#include "Quadtree.cpp"

using namespace std;

typedef vector< vector< particle_t* > > particleGrid;
typedef vector< particle_t* > gridBin;

//
//create bins with length of cutoff
//
int create_bins( particleGrid &bins, particle_t* particles, int n ) 
{
	int binsPerRow = ceil( sqrt( density * n ) / cutoff );

	//resize the vector to the exact size of bins.
	bins.resize( binsPerRow * binsPerRow );

	return binsPerRow;
}

void move_particles( particle_t* particles, int n )
{
	//  move particles
	for( int p = 0; p < n; p++ )
	{
		move( particles[ p ] );
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

void compute_forces( int numOfBins, int binsPerRow, particleGrid &bins, particle_t* particles, double davg, double dmin, int navg )
{
	//
	//  compute forces
	//
	for ( int i = 0; i < numOfBins; i++ )
	{
		gridBin binQ = bins[ i ];

		int numParticlesPerBin = binQ.size( );

		for ( int j = 0; j < numParticlesPerBin; j++ ) 
		{
			particles[ j ].ax = particles[ j ].ay = 0;
		}

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


		//This should manage the ranges such that the halo region is searched.
		int xSize = x.size( );
		int ySize = y.size( );
		int binSize = bins[ i ].size( );
		for ( int a = 0; a < xSize; a++ ) 
		{
			for ( int b = 0; b < ySize; b++ ) 
			{
				int bin_num = i + x[ a ] + binsPerRow*y[ b ];
				int binSize2 = bins[ bin_num ].size( );
				for ( int c = 0; c < binSize; c++ ) 
				{
					for ( int d = 0; d < binSize2; d++ ) 
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

//
//  benchmarking program
//
int main( int argc, char **argv )
{
	int navg, nabsavg = 0;
	double davg, dmin, absmin = 1.0, absavg = 0.0;

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
		fprintf( stderr, "particles malloc NULL\n" );
		return 0;
	}

	set_size( n );
	init_particles( n, particles );

	//create the bins to contain the particles
	particleGrid bins;
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

		clean_bins( bins, numOfBins );

		relocate( bins, binsPerRow, particles, n );

		compute_forces( numOfBins, binsPerRow, bins, particles, davg, dmin, navg );

		move_particles( particles, n );

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
		fprintf( fsum, "%d %g\n", n, simulation_time );

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
