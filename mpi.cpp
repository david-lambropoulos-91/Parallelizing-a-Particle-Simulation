#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <iostream>
#include <string.h>
#include "Quadtree.cpp"

using namespace std;

typedef vector< particle_t > gridBin;
typedef vector< gridBin > particleGrid;


//
//create bins with length of cutoff
//
void create_bins( gridBin bins[ ], particle_t* particles, int n, int binsPerRow, int first ) 
{

    //put particles in bins according to their locations
    for ( int j = 0; j < n; j++ ) 
{
        int x = ( int ) floor( particles[ j ].x / binsize );
        int y = ( int ) floor( particles[ j ].y / binsize );
        bins[ x + ( y - first ) * binsPerRow ].push_back( particles[ j ] );
    }
}

void compute_forces( gridBin localBins[ ], int first_real_bin, int last_real_bin, int binsPerRow, double dmin, double davg, int navg )
{
    for ( int biter = first_real_bin; biter < last_real_bin; biter++ )
        {
            gridBin bin = localBins[ biter ];
            int particles_per_bin = bin.size( );

            for ( int j = 0; j < particles_per_bin; j++ ) 
{
                bin[ j ].ax = bin[ j ].ay = 0;
            }
            
vector< int > x_range;
            vector< int > y_range;
            
x_range.push_back( 0 );
            y_range.push_back( 0 );

            if ( biter >= binsPerRow ) 
{
                y_range.push_back( -1 );
            }
            if ( biter < binsPerRow * ( binsPerRow - 1 ) ) 
{
                y_range.push_back( 1 );
            }
            if ( biter % binsPerRow != 0 ) 
{
                x_range.push_back( -1 );
            }
            if ( biter % binsPerRow != binsPerRow - 1 ) 
{
                x_range.push_back( 1 );
            }
            
//This should manage the ranges such that the halo region is searched.
            for (int a = 0; a < x_range.size(); a++) 
{
                for (int b = 0; b < y_range.size(); b++) 
{
                    int bin_num = biter + x_range[a] + binsPerRow*y_range[b];

                    for ( int c = 0; c < bin.size( ); c++ ) 
{
                        for ( int d = 0; d < localBins[ bin_num ].size( ); d++ ) 
{
                            apply_force( bin[ c ], localBins[ bin_num ][ d ], &dmin, &davg, &navg );
                        }
                    }
                }
            }

            localBins[ biter ] = bin;
        }
}

//Partition particles into numProc bins based on their row location. 
void partition_bins( gridBin bins[ ], particle_t* particles, int* particles_per_process, int n, int binsPerRow, int numProc ) 
{
    
    memset ( particles_per_process, 0, sizeof( int ) * numProc );
    int procRows = ( binsPerRow + numProc - 1 ) / numProc;

    //put particles in bins according to their locations
    for ( int j = 0; j < n; j++ ) 
{
        int y = ( int ) floor( particles[ j ].y / binsize);
        int procdex = ( int )floor( y / procRows );
        bins[ procdex ].push_back( particles[ j ] );
        particles_per_process[ procdex ]++;

        //If this particle is in a halo bin, we must also import it twice.
        int boundcheck = y % procRows;
        if( boundcheck == 0 && y != 0 )
{
            //We now know that y is on a boundary, and must be included in procdex-1.
            bins[ procdex - 1 ].push_back( particles[ j ] );
            particles_per_process[ procdex - 1 ]++;
        }
        else if( boundcheck == procRows - 1 && y != binsPerRow - 1 )
{
            //We now know that y is on a boundary, and must be included in procdex+1.
            bins[ procdex + 1 ].push_back( particles[ j ] );
            particles_per_process[ procdex + 1 ]++;
        }
    }
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
    int* sizeOfPartitions;
    int* offsetOfPartitions;
    particle_t *sendBuf;
 
    //
    //  process command line parameters
    //
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
    
    //
    //  set up MPI
    //
    int numProc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &numProc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles;

if( ( particles = ( particle_t* ) malloc( n * sizeof( particle_t ) ) ) == NULL )
{
fprintf( stderr, "particles malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return -1;
}
    
if( ( offsetOfPartitions = ( int* ) malloc( sizeof( int ) * numProc) ) == NULL )
{
fprintf( stderr, "offsetOfPartitions malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return -1;
}
    
if( ( sizeOfPartitions = ( int* ) malloc( sizeof( int )*numProc) ) == NULL )
{
fprintf( stderr, "sizeOfPartitions malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return -1;
}

    gridBin bins[numProc];
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );


    int binsPerRow = (int) ceil( sqrt( density * n ) / binsize);
    int procRows = ( binsPerRow + numProc - 1 ) / numProc;
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );


    if( rank == 0 )
{
        init_particles( n, particles );
        
//partition_bins will sort the particles into numProc bins, based on their bin row index
        //While we do this, we want to also count the number of particles in each level;
        //We use this info to malloc a particle_t array to scatterv across all processes.
        partition_bins( bins, particles, sizeOfPartitions, n, binsPerRow, numProc );


        int totalSize = 0;
        offsetOfPartitions[ 0 ] = 0;
        
for ( int i = 0; i < (numProc - 1); i++ )
{
            totalSize += sizeOfPartitions[ i ];
            offsetOfPartitions[ i + 1 ] = offsetOfPartitions[ i ] + sizeOfPartitions[ i ]; 
        }

        totalSize += sizeOfPartitions[ numProc - 1 ];
        
//Initialize the large array of particles we will be sending
        if( ( sendBuf = ( particle_t* ) malloc( totalSize * sizeof( particle_t ) ) ) == NULL )
{
fprintf( stderr, "sendBuf malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return -1;
}
        
//This loop is meant to fill sendBuf with contiguous particles.
        for( int i = 0; i < numProc; i++)
{
            memcpy(&sendBuf[offsetOfPartitions[i]], bins[i].data(), sizeOfPartitions[i] * sizeof(particle_t));
        }
    }
        
MPI_Bcast( sizeOfPartitions, numProc, MPI_INT, 0, MPI_COMM_WORLD );
MPI_Bcast( offsetOfPartitions, numProc, MPI_INT, 0, MPI_COMM_WORLD );
    
//
    //  allocate storage for local partition
    //
    int nlocal = sizeOfPartitions[ rank ];
    size_t movementsize = n * sizeof( particle_t );
    
particle_t *local;

if( ( local = ( particle_t* ) malloc( nlocal * sizeof( particle_t ) ) ) == NULL )
{
fprintf( stderr, "local malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return -1;
}
    
particle_t *fromAbove;

if( ( fromAbove = ( particle_t* ) malloc( movementsize ) ) == NULL )
{
fprintf( stderr, "fromAbove malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return -1;
}

    particle_t *fromBelow;

if( ( fromBelow = ( particle_t* ) malloc( movementsize ) ) == NULL )
{
fprintf( stderr, "fromBelow malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return -1;
}

    particle_t* movingup;

if( ( movingup = ( particle_t* ) malloc( movementsize ) ) == NULL )
{
fprintf( stderr, "sendBuf malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return -1;
}

    particle_t* movingdown;

if( ( movingdown = ( particle_t* ) malloc( movementsize ) ) == NULL )
{
fprintf( stderr, "movingdown malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return -1;
}
    
    
    MPI_Scatterv( sendBuf, sizeOfPartitions, offsetOfPartitions, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    
//
    //  Create bins for local rows.
    //
particle_t* zippy;

    if( ( zippy = ( particle_t* ) malloc( sizeof( particle_t ) * n * numProc / 2 ) ) == NULL )
{
fprintf( stderr, "zippy malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return 0;
}

    int* zipoff;

if( ( zipoff = (int*) malloc(sizeof(int) * numProc) ) == NULL )
{
fprintf( stderr, "zipoff malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return 0;
}

    zipoff[ 0 ] = 0;

    for( int i = 1; i < numProc - 1; i++)
{
        zipoff[ i ] = zipoff[ i - 1 ] + n / 2;
    }

    int* howManyZips;

if( ( howManyZips = ( int* ) malloc( sizeof( int ) * numProc ) ) == NULL )
{
fprintf( stderr, "howManyZips malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return 0;
}

    particle_t* localzip;

if( ( localzip = ( particle_t* )malloc( sizeof( particle_t ) * n / 2 ) ) == NULL )
{
fprintf( stderr, "localzip malloc NULL at line %d in file %s\n", __LINE__, __FILE__ );
return 0;
}

    int first = min( rank * procRows, binsPerRow);
    int last = min( ( rank + 1 ) * procRows, binsPerRow );
    first--;
    int last_real_bin = last - first;
    last++;
    int first_real_bin = 1;

    if(rank == 0)
{
        first_real_bin--;
        last_real_bin--;
        first++;
    }
    else if( rank == numProc )
{
        last--;
    }
    
first_real_bin *= binsPerRow;
    last_real_bin *= binsPerRow;

    int processorBins = last - first; //General case
    int sizeOfLocalBin = processorBins * binsPerRow;

    gridBin localBins[ sizeOfLocalBin ];

    create_bins( localBins, local, nlocal, binsPerRow, first );



    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
        // 
        //  collect all global data locally (not good idea to do)
        //
        
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1)
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  compute all forces
        //
        //  Do one set of computations for each bin.
compute_forces( localBins, first_real_bin, last_real_bin, binsPerRow, dmin, davg, navg );

        MPI_Barrier( MPI_COMM_WORLD );
        if( find_option( argc, argv, "-no" ) == -1 )
        {
          
          MPI_Reduce( &davg, &rdavg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
          MPI_Reduce( &navg, &rnavg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
          MPI_Reduce( &dmin, &rdmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );

 
          if ( rank == 0 )
  {
            //
            // Computing statistical data
            //
            if ( rnavg ) 
{
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if ( rdmin < absmin ) absmin = rdmin;
          }
        }

        //
        //  move particles
        //
        //  We must first move particles as per usual, except our temp array of particles to exit our system
        //  must now be sent to the relevant folks.
        gridBin moveUp;
        gridBin moveDown;
        gridBin temp_move;
        
        int ob1 = -binsPerRow*procRows;
        int ob2 = (last_real_bin+ binsPerRow);
        int zip = 0;
        for ( int biter = first_real_bin; biter < last_real_bin; biter++ )
        {//Insert logic here
            gridBin bin = localBins[ biter ];
            
int size = bin.size( );
            for ( int p = 0; p < size; ) 
{
                move( bin[ p ] );

                int x = ( int ) floor( bin[ p ].x / binsize );
                int y = ( int ) floor( bin[ p ].y / binsize );
                int loc = ( y - first ) * binsPerRow + x;
                if ( loc != biter )
                {
                    if (loc < ob1 || loc > ob2 )
{
                        if( zip > n / 2 )
{
                            fprintf( stderr, "Error, buffer overflow localzip\n" );
                            return -1;
                        }
                        localzip[ zip ] = bin[ p ];
                        zip++;
                    }
                    else if( loc < binsPerRow )
{
                        moveUp.push_back( bin[ p ] );
}
                    else if( loc > last_real_bin )
{
                        moveDown.push_back( bin[ p ] );
}
                    else
{
                        temp_move.push_back( bin[ p ] );
                    }

                    bin.erase( bin.begin( ) + p );
                    size--;
                }
                else
{
                    p++;
                }
            }
            localBins[ biter ] = bin;
        }
        
        //Handle the zippy particles first.
        MPI_Barrier( MPI_COMM_WORLD );        
        MPI_Allgather( &zip, 1, MPI_INT, howManyZips, 1, MPI_INT, MPI_COMM_WORLD );
        MPI_Allgatherv( localzip, zip, PARTICLE, zippy, howManyZips, zipoff, PARTICLE, MPI_COMM_WORLD );

        for( int i = 0; i < numProc; i++ )
{
                if( howManyZips[ i ] > n / 2 )
{
                    fprintf( stderr, "Error, buffer overflow.\n" );
                    return -1;
                }
            for( int j = 0; j < howManyZips[ i ]; j++ )
{
                int x = ( int ) floor( zippy[ zipoff[ i ] + j ].x / binsize );
                int y = ( int ) floor( zippy[ zipoff[ i ] + j ].y / binsize );
                int biq = x + ( y - first ) * binsPerRow;
                if( biq >= first_real_bin && biq < last_real_bin )
{
                    localBins[ x + ( y - first ) * binsPerRow ].push_back( zippy[ zipoff[ i ] + j ] );
}
}
        }
        MPI_Barrier( MPI_COMM_WORLD );

        int fa, fb;
        fa = fb = 0;
        if( rank > 0 )
{
            MPI_Send( moveUp.data( ), moveUp.size( ), PARTICLE, rank - 1, rank, MPI_COMM_WORLD );
        }
        
if( rank < numProc - 1 )
{
            MPI_Status stat;
            MPI_Recv( fromBelow, n, PARTICLE, rank + 1, rank + 1, MPI_COMM_WORLD, &stat );
            MPI_Get_count( &stat, PARTICLE, &fb );
            MPI_Send( moveDown.data( ), moveDown.size( ), PARTICLE, rank + 1, rank, MPI_COMM_WORLD );
        }

        if(rank > 0){
            MPI_Status stat;
            MPI_Recv(fromAbove, n, PARTICLE, rank-1, rank-1, MPI_COMM_WORLD, &stat);
            MPI_Get_count(&stat, PARTICLE, &fa);
        }

        //Now we wish to recieve the data of all processes which have moved particles into our system.
        //We also wish to send the data of particles which have exited our system to our neighboring processes.
        MPI_Barrier( MPI_COMM_WORLD );
        for( int i = 0; i < fa; i++ )
{
            int x = ( int ) floor( fromAbove[ i ].x / binsize );
            int y = ( int ) floor( fromAbove[ i ].y / binsize );
            int ob1 = - binsPerRow*procRows;
            int ob2 = ( last_real_bin + binsPerRow );
            int loc = ( y - first ) * binsPerRow + x;
            if ( loc < ob1 || loc > ob2 )
{
                //printf("rank %d out of bounds0 \n", rank);
}
            else
{
localBins[ x + ( y - first ) * binsPerRow ].push_back( fromAbove[ i ] );
}
        }
        for( int i = 0; i < fb; i++ )
{
            int x = ( int ) floor( fromBelow[ i ].x / binsize );
            int y = ( int ) floor( fromBelow[ i ].y / binsize );
            int ob1 = - binsPerRow* procRows;
            int ob2 = ( last_real_bin + binsPerRow );
            int loc = ( y - first ) * binsPerRow + x;
            
if ( loc < ob1 || loc > ob2 )
{
                //printf("rank %d out of bounds1 \n", rank);
}
            else
{
localBins[ x + ( y - first ) * binsPerRow ].push_back( fromBelow[ i ] );
}
        }


//We accomplish this by clearing our topmost and bottommost bins, 
//With the exception of the rank 0 and rank numProc processes, handle those with care,
//And then sending the second topmost and second bottommost bin rows to our neighbors
//And then receiving the updated halo regions.


        for( int i = 0; i < first_real_bin; i++ )
{
            localBins[ i ].clear( );
        }
        for( int i = last_real_bin; i < last_real_bin+binsPerRow; i++ )
{
            localBins[ i ].clear( );
        }

        //Send the first REAL row to our neighbors.
        moveUp.clear( );
        moveDown.clear( );
        int upsize = 0;
        int downsize = 0;
        if(rank > 0)
{
            for(int eob = first_real_bin; eob < first_real_bin + binsPerRow; eob++)
{
                memcpy(&movingup[upsize], localBins[eob].data(), localBins[eob].size() * sizeof(particle_t));
                upsize += localBins[eob].size();
            }   
        }
        if(rank < numProc-1){
            for(int boe = sizeOfLocalBin; boe < last_real_bin; boe++)
{
                memcpy(&movingdown[downsize], localBins[boe].data(), localBins[boe].size() * sizeof(particle_t));
                downsize += localBins[boe].size();
            }   
        }

        MPI_Barrier(MPI_COMM_WORLD);

        fa = fb = 0;
        if(rank > 0)
{
            MPI_Send(movingup, upsize, PARTICLE, rank-1, rank, MPI_COMM_WORLD);
        }
        if(rank < numProc-1)
{
            MPI_Status stat;
            MPI_Recv(fromBelow, n, PARTICLE, rank+1, rank+1, MPI_COMM_WORLD, &stat);
            MPI_Get_count(&stat, PARTICLE, &fb);
            MPI_Send(movingdown, downsize, PARTICLE, rank+1, rank, MPI_COMM_WORLD);
        }
        if(rank > 0)
{
            MPI_Status stat;
            MPI_Recv(fromAbove, n, PARTICLE, rank-1, rank-1, MPI_COMM_WORLD, &stat);
            MPI_Get_count(&stat, PARTICLE, &fa);
        }
        
//We must now handle the data received differently though; we have to rebin it into our halo regions.
        if(rank > 0)
{
            for (int j = 0; j < fa; j++) 
{
                int x = floor(fromAbove[j].x / binsize);
                int y = floor(fromAbove[j].y / binsize);
                int ob1 = -binsPerRow*procRows;
                int ob2 = (last_real_bin+ binsPerRow);
                int loc = (y-first) * binsPerRow + x;

                if (loc < ob1 || loc > ob2 )
{
                    //printf("rank %d out of bounds2 \n", rank);
}
                else
{
localBins[x + (y-first) * binsPerRow].push_back(fromAbove[j]);
}
            }
}
        if(rank < numProc-1)
{
            for (int j = 0; j < fb; j++) 
{
                int x = floor(fromBelow[j].x / binsize);
                int y = floor(fromBelow[j].y / binsize);
                int ob1 = -binsPerRow*procRows;
                int ob2 = (last_real_bin+ binsPerRow);
                int loc = (y-first) * binsPerRow + x;

                if (loc < ob1 || loc > ob2 )
{
                    //printf("rank %d out of bounds3 \n", rank);
}
                else
{
localBins[x + (y-first) * binsPerRow].push_back(fromBelow[j]);
}
            }
}

        //We also have to rebin the particles we moved within our region.
            int tempmovesize = temp_move.size();
        for(int j = 0; j < tempmovesize; j++)
{
            int x = floor(temp_move[j].x / binsize);
            int y = floor(temp_move[j].y / binsize);
            int ob1 = -binsPerRow*procRows;
            int ob2 = (last_real_bin+ binsPerRow);
            int loc = (y-first) * binsPerRow + x;
            if (loc < ob1 || loc > ob2 )
                printf("rank %d out of bounds4 \n", rank);
            else
            localBins[x + (y-first) * binsPerRow].push_back(temp_move[j]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) 
{  
      printf( "n = %d, simulation time = %g seconds", n, simulation_time);

      if( find_option( argc, argv, "-no" ) == -1  && rank == 0)
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
        fprintf(fsum,"%d %d %g\n",n,numProc,simulation_time);
    }
  
    //
    //  release resources
    //
    free( offsetOfPartitions );
    free( sizeOfPartitions );
    free( local );
    free( particles );
    free( fromAbove );
    free( fromBelow );
    free( movingup );
    free( movingdown );
    free( howManyZips );
    free( zippy );
    free( localzip );
    free( zipoff );
    if(rank == 0)
{
        free( sendBuf );
        if( fsave )
            fclose( fsave );
        if ( fsum )
            fclose( fsum );
    }

    
    MPI_Finalize( );
    
    return 0;
}
