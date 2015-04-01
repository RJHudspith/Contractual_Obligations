/**
   @file baryons.c
   @brief baryon contraction code
 */

#include "common.h"

#include "basis_conversions.h" // nrel_rotate_slice()
#include "contractions.h"      // meson contract or whatever
#include "correlators.h"       // allocate_corrs() && free_corrs()
#include "gammas.h"            // make_gammas()
#include "GLU_timer.h"         // print_time()
#include "io.h"                // for read_prop()
#include "read_propheader.h"   // for read_propheader()


// This acts C*gamma_mu on a spinor
// Note: You need C*gamma_mu->source and (C*gamma_mu)^T->sink
struct spinor
Cgamma_mu( const struct spinor S, const int mu )
{
	int c1, c2, i, j;
	struct spinor TSNK, TSRC;

	for( c1 = 0 ; c1 < NC ; c1++ ) {
		for( c2 = 0 ; c2 < NC ; c2++ ) {

			// Sink indices first
			switch( mu ){
			case 0:
				for( i = 0 ; i < NS ; i++){
					TSNK.D[i][0].C[c1][c2] = (-1 + 0 * I ) * S.D[i][3].C[c1][c2];
					TSNK.D[i][1].C[c1][c2] = ( 1 + 0 * I ) * S.D[i][2].C[c1][c2];
					TSNK.D[i][2].C[c1][c2] = ( 1 + 0 * I ) * S.D[i][1].C[c1][c2];
					TSNK.D[i][3].C[c1][c2] = (-1 + 0 * I ) * S.D[i][0].C[c1][c2];
				}
				break;
			case 1:
				for( i = 0 ; i < NS ; i++){
					TSNK.D[i][0].C[c1][c2] = ( 0 - 1 * I ) * S.D[i][2].C[c1][c2];
					TSNK.D[i][1].C[c1][c2] = ( 0 + 1 * I ) * S.D[i][3].C[c1][c2];
					TSNK.D[i][2].C[c1][c2] = ( 0 - 1 * I ) * S.D[i][0].C[c1][c2];
					TSNK.D[i][3].C[c1][c2] = ( 0 + 1 * I ) * S.D[i][1].C[c1][c2];
				}
				break;
			case 2:
				for( i = 0 ; i < NS ; i++){
					TSNK.D[i][0].C[c1][c2] = (-1 + 0 * I ) * S.D[i][2].C[c1][c2];
					TSNK.D[i][1].C[c1][c2] = (-1 + 0 * I ) * S.D[i][3].C[c1][c2];
					TSNK.D[i][2].C[c1][c2] = (-1 + 0 * I ) * S.D[i][0].C[c1][c2];
					TSNK.D[i][3].C[c1][c2] = (-1 + 0 * I ) * S.D[i][1].C[c1][c2];
				}
				break;
			case 3:
				for( i = 0 ; i < NS ; i++){
					TSNK.D[i][0].C[c1][c2] = ( 0 - 1 * I ) * S.D[i][3].C[c1][c2];
					TSNK.D[i][1].C[c1][c2] = ( 0 - 1 * I ) * S.D[i][2].C[c1][c2];
					TSNK.D[i][2].C[c1][c2] = ( 0 - 1 * I ) * S.D[i][1].C[c1][c2];
					TSNK.D[i][3].C[c1][c2] = ( 0 - 1 * I ) * S.D[i][0].C[c1][c2];
				}
				break;
			case 4:
				for( i = 0 ; i < NS ; i++){
					TSNK.D[i][0].C[c1][c2] = ( 0 - 1 * I ) * S.D[i][1].C[c1][c2];
					TSNK.D[i][1].C[c1][c2] = ( 0 - 1 * I ) * S.D[i][0].C[c1][c2];
					TSNK.D[i][2].C[c1][c2] = ( 0 + 1 * I ) * S.D[i][3].C[c1][c2];
					TSNK.D[i][3].C[c1][c2] = ( 0 + 1 * I ) * S.D[i][2].C[c1][c2];
				}
				break;
			case 5:
				for( i = 0 ; i < NS ; i++){
					TSNK.D[i][0].C[c1][c2] = (-1 + 0 * I ) * S.D[i][1].C[c1][c2];
					TSNK.D[i][1].C[c1][c2] = ( 1 + 0 * I ) * S.D[i][0].C[c1][c2];
					TSNK.D[i][2].C[c1][c2] = (-1 + 0 * I ) * S.D[i][3].C[c1][c2];
					TSNK.D[i][3].C[c1][c2] = ( 1 + 0 * I ) * S.D[i][2].C[c1][c2];
				}
				break;
			}


			// Source indices last
			switch( mu ){
            case 0:
				for( j = 0 ; j < NS ; j++ ) {
					TSRC.D[0][j].C[c1][c2] = ( 1 + 0 * I ) * TSNK.D[3][j].C[c1][c2];
					TSRC.D[1][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[2][j].C[c1][c2];
					TSRC.D[2][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[1][j].C[c1][c2];
					TSRC.D[3][j].C[c1][c2] = ( 1 + 0 * I ) * TSNK.D[0][j].C[c1][c2];
				}
				break;
            case 1:
				for( j = 0 ; j < NS ; j++ ) {
					TSRC.D[0][j].C[c1][c2] = ( 0 + 1 * I ) * TSNK.D[2][j].C[c1][c2];
					TSRC.D[1][j].C[c1][c2] = ( 0 - 1 * I ) * TSNK.D[3][j].C[c1][c2];
					TSRC.D[2][j].C[c1][c2] = ( 0 + 1 * I ) * TSNK.D[0][j].C[c1][c2];
					TSRC.D[3][j].C[c1][c2] = ( 0 - 1 * I ) * TSNK.D[1][j].C[c1][c2];
				}
				break;
            case 2:
				for( j = 0 ; j < NS ; j++ ) {
					TSRC.D[0][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[2][j].C[c1][c2];
					TSRC.D[1][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[3][j].C[c1][c2];
					TSRC.D[2][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[0][j].C[c1][c2];
					TSRC.D[3][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[1][j].C[c1][c2];
				}
				break;
            case 3:
				for( j = 0 ; j < NS ; j++ ) {
					TSRC.D[0][j].C[c1][c2] = ( 0 - 1 * I ) * TSNK.D[3][j].C[c1][c2];
					TSRC.D[1][j].C[c1][c2] = ( 0 - 1 * I ) * TSNK.D[2][j].C[c1][c2];
					TSRC.D[2][j].C[c1][c2] = ( 0 - 1 * I ) * TSNK.D[1][j].C[c1][c2];
					TSRC.D[3][j].C[c1][c2] = ( 0 - 1 * I ) * TSNK.D[0][j].C[c1][c2];
				}
				break;
            case 4:
				for( j = 0 ; j < NS ; j++ ) {
					TSRC.D[0][j].C[c1][c2] = ( 0 - 1 * I ) * TSNK.D[1][j].C[c1][c2];
					TSRC.D[1][j].C[c1][c2] = ( 0 - 1 * I ) * TSNK.D[0][j].C[c1][c2];
					TSRC.D[2][j].C[c1][c2] = ( 0 + 1 * I ) * TSNK.D[3][j].C[c1][c2];
					TSRC.D[3][j].C[c1][c2] = ( 0 + 1 * I ) * TSNK.D[2][j].C[c1][c2];
				}
				break;
            case 5:
				for( j = 0 ; j < NS ; j++ ) {
					TSRC.D[0][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[1][j].C[c1][c2];
					TSRC.D[1][j].C[c1][c2] = ( 1 + 0 * I ) * TSNK.D[0][j].C[c1][c2];
					TSRC.D[2][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[3][j].C[c1][c2];
					TSRC.D[3][j].C[c1][c2] = ( 1 + 0 * I ) * TSNK.D[2][j].C[c1][c2];
				}
				break;
			}

		}
	}

	return TSRC;
}


// This carries out the color cross product and traces one set of Dirac indices.
// The results forms a diquark-type object
struct spinor
cross_color_trace( const struct spinor S ,
		   		   const struct spinor CgS )
{
	int c1, c2, c3, i, j, d;
	struct spinor TSRC;
	struct spinor T3SNK[3];

	// Initialize the threespinor object
	for( i = 0 ; i < NS ; i++ ) {
		for( j = 0 ; j < NS ; j++ ) {
			for( c1 = 0 ; c1 < NC ; c1++ ) {
				for( c2 = 0 ; c2 < NC ; c2++ ) {
					for( d = 0 ; d < NC ; d++ ) {
						T3SNK[d].D[i][j].C[c1][c2] = 0.0;
					}
				}
			}
		}
	}

	// Sink cross color and trace, this leaves only one set of dirac indices
	for( i = 0 ; i < NS ; i++ ) {
		for( j = 0 ; j < NS ; j++ ) {
			for( c1 = 0 ; c1 < NC ; c1++ ) {
				for( c2 = 0 ; c2 < NC ; c2++ ) {
					// Dirac sink trace
					for( d = 0 ; d < NS ; d++ ) {
						T3SNK[0].D[i][j].C[c1][c2] += S.D[i][d].C[c1][1]*CgS.D[j][d].C[c2][2] - S.D[i][d].C[c1][2]*CgS.D[j][d].C[c2][1];
						T3SNK[1].D[i][j].C[c1][c2] += S.D[i][d].C[c1][2]*CgS.D[j][d].C[c2][0] - S.D[i][d].C[c1][0]*CgS.D[j][d].C[c2][2];
						T3SNK[2].D[i][j].C[c1][c2] += S.D[i][d].C[c1][0]*CgS.D[j][d].C[c2][1] - S.D[i][d].C[c1][1]*CgS.D[j][d].C[c2][0];
					}
				}
			}
		}
	}

    // Source cross color
	for( i = 0 ; i < NS ; i++ ) {
        for( j = 0 ; j < NS ; j++ ) {
            for ( c3 = 0 ; c3 < 3 ; c3++ ){
				TSRC.D[i][j].C[0][c3] = T3SNK[c3].D[i][j].C[1][2] - T3SNK[c3].D[i][j].C[2][1]; 
				TSRC.D[i][j].C[1][c3] = T3SNK[c3].D[i][j].C[2][0] - T3SNK[c3].D[i][j].C[0][2]; 
				TSRC.D[i][j].C[2][c3] = T3SNK[c3].D[i][j].C[0][1] - T3SNK[c3].D[i][j].C[1][0]; 
			}            
		}
	}
	
	return TSRC;
}


// This contracts the diquark with the remaining propagator
complex double
baryon_contract( struct spinor DiQ ,
		   		 struct spinor S ,
				 int d0,
				 int d1,
				 int d2,
				 int d3 )
{
	int c1, c2, i, j;
	double corrr = 0.0 , corri = 0.0 ;

	for( c1 = 0 ; c1 < NC ; c1++ ) {
		for( c2 = 0 ; c2 < NC ; c2++ ) {

		corrr += creal( DiQ.D[d0][d1].C[c1][c2] ) * creal( S.D[d2][d3].C[c1][c2] ) 
			   - cimag( DiQ.D[d0][d1].C[c1][c2] ) * cimag( S.D[d2][d3].C[c1][c2] );

		corri += creal( DiQ.D[d0][d1].C[c1][c2] ) * cimag( S.D[d2][d3].C[c1][c2] ) 
			   + cimag( DiQ.D[d0][d1].C[c1][c2] ) * creal( S.D[d2][d3].C[c1][c2] );
	
		}
	}

	return corrr + I * corri;
}

// Reset a spinor, replace with zero_spinor from spinor_ops at some stage?
struct spinor
reset_spinor(){

	int i,j,c1,c2;
	struct spinor TSP;

	for( i = 0 ; i < NS ; i++ ){
		for( j = 0 ; j < NS ; j++ ){
			for( c1 = 0 ; c1 < NC ; c1++ ){
				for( c2 = 0 ; c2 < NC ; c2++ ){
		
					TSP.D[i][j].C[c1][c2] = 0.0 ;

				}
			}
		}
	}

	return TSP;
}



/* The Baryon interpolating operator is of the form: B = eps_123  q1 *( q2 Cg_23 q3 )
This gives us access to the particles:

					  pos1 pos2 pos3

			Proton  =  u  ( u Cg5 d )
   			Neutron =  d  ( d Cg5 u ) -> equivalent to proton in isospin limit
 				
			Lambda  =  u  ( d Cg5 s ) = d ( u Cg5 s )
			Sigma   =  u  ( u Cg5 s )
 			Cascade =  s  ( s Cg5 u )

			Omega   =  s  ( s Cg5 s )
			Delta   =  u  ( u Cgi u ) -> decuplett

The contraction terms are:

					    perm1    perm2

			Proton  =  (0,1,2) - (1,0,2)
			Neutron =  (0,1,2) - (1,0,2)

			Lambda  =  (0,1,2)
			Sigma   =  (0,1,2) - (1,0,2)
			Cascade =  (0,1,2) - (1,0,2)

			Omega   =  (0,1,2) - (1,0,2) + (1,2,0) - (0,2,1) + (2,0,1) - (2,1,0)

			Delta =  2 * (0,1,2) - 4 * (1,0,2)

*/

int
baryons_diagonal( struct propagator prop ,
		  const char *outfile )
{

  // Define our output correlators, with 6 channels and 16 components
  struct correlator **Buds_corr = allocate_corrs( 6 , NSNS ) ;
  struct correlator **Buud_corr = allocate_corrs( 6 , NSNS ) ;
  struct correlator **Buuu_corr = allocate_corrs( 6 , NSNS ) ;

  // and our spinors
  struct spinor *S1 ;
  if( posix_memalign( (void**)&S1 , 16 , 
		      VOL3 * sizeof( struct spinor ) ) != 0 ) {
      free_corrs( Buds_corr , 6 , NSNS ) ; 
      free_corrs( Buud_corr , 6 , NSNS ) ; 
      free_corrs( Buuu_corr , 6 , NSNS ) ; 
	  free( S1 ) ;
      printf( "[BARYONS] memalign failure \n" ) ;
      return FAILURE ;
  }

  int t ;
  // Time slice loop 
  for( t = 0 ; t < L0 ; t++ ) {

    // read in the file
    if( read_prop( prop , S1 ) == FAILURE ) {
      free_corrs( Buds_corr , 6 , NSNS ) ; 
      free_corrs( Buud_corr , 6 , NSNS ) ; 
      free_corrs( Buuu_corr , 6 , NSNS ) ; 
      free( S1 ) ;
      return FAILURE ;
    }

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    if( prop.basis == CHIRAL && prop.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    } 

    //
    // Begin contractions
    // 

    int GSRC = 0 ;
    // parallelise the furthest out loop
    #pragma omp parallel for private(GSRC)
    for( GSRC = 0 ; GSRC < 6 ; GSRC++ ) {

		 // Define some intermediate spinors
 	 	complex double term1, term2;
  		struct spinor CgS;	
  		struct spinor DiQ;	
  		register double complex Buds ;
  		register double complex Buud ;
  		register double complex Buuu ;

		int odc=0;
		int OD1;
    	for( OD1 = 0 ; OD1 < 1 ; OD1++ ) {

			int OD2;
    		for( OD2 = 0 ; OD2 < 1 ; OD2++ ) {

				Buds = 0.0 ;
				Buud = 0.0 ;
				Buuu = 0.0 ;

				// loop spatial hypercube
				int site ;
				for( site = 0 ; site < VOL3 ; site++ ) {
		
				term1 = 0.0 ;
				term2 = 0.0 ;

				CgS = reset_spinor();
				DiQ = reset_spinor();

				// Act C*gamma_mu on source and sink of spinor S1
				CgS = Cgamma_mu( S1[ site ] , GSRC );

				// Cross color product and sink Dirac trace
				DiQ = cross_color_trace( S1[ site ], CgS );

				// Contract with the final propagator and trace out the source Dirac indices
				// A polarization must still be picked for the two open Dirac indices offline
				int dirac;
				for( dirac = 0 ; dirac < NS ; dirac++ ){
		 			term1 += baryon_contract( DiQ, S1[ site ], dirac , dirac , OD1 , OD2 );
		 			term2 += baryon_contract( DiQ, S1[ site ], OD1 , dirac , dirac , OD2 );
				}

//#ifdef DEBUG
//	if( site < 2 ){
//		printf("First term (x=%d):  %le %le\n", site, creal(term1), cimag(term1) );
//		printf("Second term (x=%d): %le %le\n", site, creal(term2), cimag(term2) );
//	}
//#endif

				// Form the uds-, uud- and uuu-type baryons
				Buds += term1 ;
				Buud += term1 + term2 ;
				Buuu += 2 * term1 + 4 * term2 ;

//#ifdef DEBUG
//	if( site == 0 ){
//		printf("{Printing OD1=%d OD2=%d into index %d of corr files\n", OD1,OD2,odc);
//	}
//#endif
				}

				// Fill baryon correlator array
				Buds_corr[ GSRC ][ odc ].C[ t ] = Buds ;
				Buud_corr[ GSRC ][ odc ].C[ t ] = Buud ;
				Buuu_corr[ GSRC ][ odc ].C[ t ] = Buuu ;
			
				// Move up the open Dirac counter
				odc += 1 ;
			}
		}

    }

    // status of the computation
    printf("\r[BARYONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

#ifdef DEBUG
  debug_baryons( "Baryon: uds-type" , (const struct correlator**)Buds_corr ) ;
  debug_baryons( "Baryon: uud-type" , (const struct correlator**)Buud_corr ) ;
  debug_baryons( "Baryon: uuu-type" , (const struct correlator**)Buuu_corr ) ;
#endif


  // write out the correlator
  char outstr[ 256 ] ;
  sprintf( outstr , "%s" , outfile ) ;
  //write_correlators( outstr , (const struct correlator**)Buds_corr ,NSNS , NSNS ) ;

  free_corrs( Buds_corr , 6 , NSNS ) ;
  free_corrs( Buud_corr , 6 , NSNS ) ;
  free_corrs( Buuu_corr , 6 , NSNS ) ;

  free( S1 ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;
}





