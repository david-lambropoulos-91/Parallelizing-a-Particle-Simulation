#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"
#include <list>
#include "Quadtree.cpp"

//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    double sizeOfBin;         // Size of the each bin
    double sizeOfGrid;        // Size of the grid of bins
    int binsPerRow;           // Number of bins in row of (binsPerRow)x(binsPerRow) grid
    int numOfBins;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    // Vector containing vectors (bins) of particles
    //particleBins = <b1,b2,...,bn>
    //b1 = <p1,...,pn>
    //b2 = <p1,...,pn>
    //...
    //bn = <p1,...,pn>
    std::vector< std::vector< particle_t > > particleBins;

    //
    //  Construct grid of bins
    //

    sizeOfGrid = sqrt( n * density );
    sizeOfBin = cutoff * 2;
    binsPerRow = ( int ) (sizeOfGrid / sizeOfBin) + 1;
    numOfBins = binsPerRow * binsPerRow;

    // Initialize the root node
    Node * rootNode = new Node;    //Initialised node
    
    setnode(rootNode, 0, 0, 500, 500);
    
    for(int i = 0; i < n; i++)
    {
    	
    }


    // Set size of the particle bins vector
    particleBins.resize( numOfBins );

    // Populate particle vector with particles
    for( int i = 0; i < n; i++ )
    {	
	// Obtain x and y coordinates
	int x = ( int ) particles[ i ].x / sizeOfBin;
	int y = ( int ) particles[ i ].y / sizeOfBin;

	// Place particles in respective bin based on location
	particleBins[ x * binsPerRow + y ].push_back( particles[ i ] );
    }

        // List of particles no longer in their original bin
        std::list< particle_t > displacedParticles;
     	std::list< particle_t >::iterator it;

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;
	dmin = 1.0;
        
	//
        //  compute forces
        //
        for( int i = 0; i < binsPerRow; i++ )
        {
		for( int j = 0; j < binsPerRow; j++ )
		{
			// Grab first bin
			std::vector< particle_t >& bin1 = particleBins[ i * binsPerRow + j ];
			
			// Obtain size of that bin
			int bin1_size = bin1.size( );

			// Set acceleration in both x and y direction to 0
			for( int k = 0; k < bin1_size; k++ )
			{
				bin1[ k ].ax = bin1[ k ].ay = 0; 
			}
			
			// Go through the Moore's neighborhood
			// https://en.wikipedia.org/wiki/Moore_neighborhood 
			//    +---+---+---+
			//    |   |   |   |
			//    +---+---+---+
			//    |   | x |   |
			//    +---+---+---+
			//    |   |   |   |
			//    +---+---+---+
			//
			//    *** Possible Improvement ***
			//    As suggested in http://www.cs.cornell.edu/~bindel/class/cs5220-f11/notes/spatial.pdf
			//    Given binsize >= 2*radius only need to check at most 3 neighbors and not the entire
			//    neighborhood.
			for( int ii = -1; ii <= 1; ii++ )
			{
				for( int jj = -1; jj <= 1; jj++ )
				{
					// Get (x,y) coordinate of neighbor
					int neighborX = i + ii;
					int neighborY = j + jj;
	
					// Check if neighbor coordinate isn't out of the grid
					if( neighborX >= 0 && neighborX < binsPerRow  &&  neighborY >= 0 && neighborY < binsPerRow  )
					{
						// Grab neighbor bin
						std::vector< particle_t >& bin2 = particleBins[ neighborX * binsPerRow + neighborY ];
						
						// Obtain size of neighbor bin
						int bin2_size = bin2.size( );

						// Calculate force on bin by neighboring bin
						for( int l = 0; l < bin1_size; l++ )
						{
							for( int m = 0; m < bin2_size; m++ )
							{
								apply_force( bin1[ l ], bin2[ m ], &dmin, &davg, &navg );
							}
						}
					}
				}			
			} 
		}			   
        }
 
	// List of particles no longer in their original bin
//	std::list< particle_t > displacedParticles;
//	std::list< particle_t >::iterator it;

        //
        //  move particles
        //
        for( int i = 0; i < binsPerRow; i++ ) 
        {
		for( int j = 0; j < binsPerRow; j++ )
		{
			// Grab bin of particles
			std::vector< particle_t >& bin = particleBins[ i * binsPerRow + j ];
			
			// Obtain the size of the bin
			int bin_size = bin.size( );
			int k = 0;
			for( ; k < bin_size; )
			{
				move( bin[ k ] );

				// Get new position; (x,y)
				int x = ( int ) bin[ k ].x / sizeOfBin;
				int y = ( int ) bin[ k ].y / sizeOfBin;
				
				// Check if the particle is still in the original bin
				// if it is not move it to its new bin otherwise continue
				if( x == i && y == j )
				{
					k++;
				}
				else
				{
					displacedParticles.push_back( bin[ k ] );
					bin[ k ] = bin[ --bin_size ];		
				}
			}

			bin.resize( k );
		}	
	}
        
	//
	//  Handle particles that moved outside of their bin
	//
	
	for( it = displacedParticles.begin(); it != displacedParticles.end(); ++it )
	//for(int i = 0; i < displacedParticles.size(); i++)
	{
		//int x = ( int ) displacedParticles[ i ].x / sizeOfBin;
		//int y = ( int ) displacedParticles[ i ].y / sizeOfBin;
		int x = ( int ) it->x / sizeOfBin;
		int y = ( int ) it->y / sizeOfBin;

		particleBins[ x * binsPerRow + y].push_back(*it);
	//	particleBins[ x * binsPerRow + y ].push_back( displacedParticles[ i ] );
	}	

	displacedParticles.clear();

	if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
