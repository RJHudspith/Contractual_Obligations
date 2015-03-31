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

// This acts C*gamma_mu on a spinor, atm: Cg5
struct spinor
Cgamma_mu( const struct spinor S )
{
	int c1, c2, i, j;
	struct spinor TSNK, TSRC;


#ifdef DEBUG
//	printf("Input to Cg_mu: %le %le\n",creal(S.D[0][0].C[0][0]),cimag(S.D[0][0].C[0][0]));
#endif

	for( c1 = 0 ; c1 < NC ; c1++ ) {
		for( c2 = 0 ; c2 < NC ; c2++ ) {

			// Sink indices first
			for( i = 0 ; i < NS ; i++){
				TSNK.D[i][0].C[c1][c2] = (-1 + 0 * I ) * S.D[i][1].C[c1][c2];
				TSNK.D[i][1].C[c1][c2] = ( 1 + 0 * I ) * S.D[i][0].C[c1][c2];
				TSNK.D[i][2].C[c1][c2] = (-1 + 0 * I ) * S.D[i][3].C[c1][c2];
				TSNK.D[i][3].C[c1][c2] = ( 1 + 0 * I ) * S.D[i][2].C[c1][c2];
			}

#ifdef DEBUG
//	if ( c1 == 0 && c2 == 0 ){ 
//		printf("Sink done: %le %le\n",creal(TSNK.D[0][0].C[0][0]),cimag(TSNK.D[0][0].C[0][0]));
//	}
#endif
			// Source indices last
			for( j = 0 ; j < NS ; j++ ) {
				TSRC.D[0][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[1][j].C[c1][c2];
				TSRC.D[1][j].C[c1][c2] = ( 1 + 0 * I ) * TSNK.D[0][j].C[c1][c2];
				TSRC.D[2][j].C[c1][c2] = (-1 + 0 * I ) * TSNK.D[3][j].C[c1][c2];
				TSRC.D[3][j].C[c1][c2] = ( 1 + 0 * I ) * TSNK.D[2][j].C[c1][c2];
			}
		}
	}

#ifdef DEBUG
//	printf("Source done: %le %le\n",creal(TSRC.D[0][0].C[0][0]),cimag(TSRC.D[0][0].C[0][0]));
#endif

	return TSRC;
}


// This acts C*gamma_mu on a spinor, atm: Cg5
struct spinor
cross_color_trace( const struct spinor S ,
		   		   const struct spinor CgS )
{
	int c1, c2, c3, i, j, d;
	struct spinor TSRC;
	struct threespinor T3SNK;

	// Sink cross color and trace, this leaves only one set of dirac indices
	for( i = 0 ; i < NS ; i++ ) {
		for( j = 0 ; j < NS ; j++ ) {
			for( c1 = 0 ; c1 < NC ; c1++ ) {
				for( c2 = 0 ; c2 < NC ; c2++ ) {
					// Dirac sink trace
					for( d = 0 ; d < NS ; d++ ) {
						T3SNK.D[i][j].C[c1][c2][0] += S.D[i][d].C[c1][1]*CgS.D[j][d].C[c2][2] - S.D[i][d].C[c1][2]*CgS.D[j][d].C[c2][1];
						T3SNK.D[i][j].C[c1][c2][1] += S.D[i][d].C[c1][2]*CgS.D[j][d].C[c2][0] - S.D[i][d].C[c1][0]*CgS.D[j][d].C[c2][2];
						T3SNK.D[i][j].C[c1][c2][2] += S.D[i][d].C[c1][0]*CgS.D[j][d].C[c2][1] - S.D[i][d].C[c1][1]*CgS.D[j][d].C[c2][0];
					}
				}
			}
		}
	}

#ifdef DEBUG
//     printf("Sink cross color trace 0: %d %d %d %d %le %le\n",3,3,2,2,creal(T3SNK.D[3][3].C[2][2][0]),cimag(T3SNK.D[3][3].C[2][2][0]));
#endif

    // Source cross color
	for( i = 0 ; i < NS ; i++ ) {
        for( j = 0 ; j < NS ; j++ ) {
            for ( c3 = 0 ; c3 < 3 ; c3++ ){
				TSRC.D[i][j].C[0][c3] = T3SNK.D[i][j].C[1][2][c3] - T3SNK.D[i][j].C[2][1][c3]; 
				TSRC.D[i][j].C[1][c3] = T3SNK.D[i][j].C[2][0][c3] - T3SNK.D[i][j].C[0][2][c3]; 
				TSRC.D[i][j].C[2][c3] = T3SNK.D[i][j].C[0][1][c3] - T3SNK.D[i][j].C[1][0][c3]; 
			}            
		}
	}
	
	
#ifdef DEBUG
     printf("Sink cross color trace: %d %d %d %d %le %le\n",3,3,2,2,creal(TSRC.D[3][3].C[2][2]),cimag(TSRC.D[3][3].C[2][2]));
#endif

	return TSRC;
}


// This contracts the diquark with the remaining propagator
complex double
baryon_contract( const struct spinor DiQ ,
		   		 const struct spinor S ,
				 const int d0,
				 const int d1,
				 const int d2,
				 const int d3 )
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
#ifdef DEBUG
//	printf("%le %le\n",corrr,corri);
#endif
	}


	return corrr + I * corri;
}


/* The Baryon interpolating operator is of the form: B = eps_123  q1 *( q2 Cg_23 q3 )
This gives us access to the particles:

					  pos1 pos2 pos3

			Proton  =  u  ( u Cg5 d )
   			Neutron =  d  ( d Cg5 u ) -> equivalent to proton in isospin limit
 				
			Lambda  =  u  ( d Cg5 s ) = d ( u Cg5 s )
			Sigma   =  u  ( u Cg5 s )
 			Cascade =  s  ( s Cg5 u )

			Omega   =  s  ( s Cgi s )
			Delta   =  u  ( u Cgi u )

The contraction terms are:

					    perm1    perm2

			Proton  =  (0,1,2) - (1,0,2)
			Neutron =  (0,1,2) - (1,0,2)

			Lambda  =  (0,1,2)
			Sigma   =  (0,1,2) - (1,0,2)
			Cascade =  (0,1,2) - (1,0,2)

			Omega   =  (0,1,2) - (1,0,2) + (1,2,0) - (0,2,1) + (2,0,1) - (2,1,0)
			Delta   =  (0,1,2) - (1,0,2) + (1,2,0) - (0,2,1) + (2,0,1) - (2,1,0)
*/

int
baryons_diagonal( struct propagator prop ,
		  const char *outfile )
{
  // allocate the basis, maybe extern this as it is important ...
  struct gamma *GAMMAS = malloc( NSNS * sizeof( struct gamma ) ) ;

  // precompute the gamma basis
  if( make_gammas( GAMMAS , prop.basis ) == FAILURE ) {
    free( GAMMAS ) ;
    return FAILURE ;
  }

  // and our spinors
  struct spinor *S1 ;
  if( posix_memalign( (void**)&S1 , 16 , 
		      VOL3 * sizeof( struct spinor ) ) != 0 ) {
    free( S1 ) ; free( GAMMAS ) ;
    printf( "[BARYONS] memalign failure \n" ) ;
    return FAILURE ;
  }

  struct spinor CgS;	
  struct spinor DiQ;	

  struct correlator **bcorr = allocate_corrs( NSNS , NSNS ) ;

  int t ;
  // Time slice loop 
  for( t = 0 ; t < 1 ; t++ ) {

    // read in the file
    if( read_prop( prop , S1 ) == FAILURE ) {
      free_corrs( bcorr , NSNS , NSNS ) ; 
      free( GAMMAS ) ; free( S1 ) ;
      return FAILURE ;
    }

    // if we are doing nonrel-chiral mesons we switch chiral to nrel
    if( prop.basis == CHIRAL && prop.basis == NREL ) {
      nrel_rotate_slice( S1 ) ;
    } 

    //
    // Contractions go here
    // 

    int GSRC = 0 ;
    // parallelise the furthest out loop
    //#pragma omp parallel for private(GSRC)
    for( GSRC = 5 ; GSRC < 6 ; GSRC++ ) {

      int GSNK ;
      //for( GSNK = 0 ; GSNK < NSNS ; GSNK++ ) {
	  GSNK = GSRC;

		register double complex corr0 = 0.0 ;
		register double complex corr1 = 0.0 ;

		// loop spatial hypercube
		int site ;
		for( site = 0 ; site < 1 ; site++ ) {
	
			// Act C*gamma_mu on source and sink of spinor S1
			CgS = Cgamma_mu( S1[ site ] );

#ifdef DEBUG
	if( site == 0 ){
		printf("Got Cg_mu: %le %le\n",creal(CgS.D[0][0].C[0][0]),cimag(CgS.D[0][0].C[0][0]));
	}
#endif
			// Cross color product and sink Dirac trace
			DiQ = cross_color_trace( S1[ site ], CgS );

#ifdef DEBUG
	if( site == 0 ){
		printf("Got sink trace and cross color: %le %le\n",creal(DiQ.D[3][3].C[2][2]),cimag(DiQ.D[3][3].C[2][2]));
	}
#endif
			// Contract with the final propagator and trace out the source Dirac indices
			int dirac;
			for( dirac = 0 ; dirac < NS ; dirac++ ){
		 		corr0 += baryon_contract( DiQ, S1[ site ], dirac , dirac , 0 , 0 );
		 		corr1 += baryon_contract( DiQ, S1[ site ], 0 , dirac , dirac , 0 );
			}

#ifdef DEBUG
	if( site <= 10 ){
		printf("First term: %le %le\n",  creal(corr0), cimag(corr0) );
		printf("Second term: %le %le\n", creal(corr1), cimag(corr1) );
	}
#endif


		}

		// uds-type baryon correlator
		bcorr[ GSRC ][ GSNK ].C[ t ] = corr0 ;


      // GSNK}
    }

    // status of the computation
    printf("\r[BARYONS] done %.f %%", (t+1)/((L0)/100.) ) ; 
    fflush( stdout ) ;
  }
  printf( "\n" ) ;

#ifdef DEBUG
//  debug_mesons( "Baryons" , (const struct correlator**)bcorr ) ;
#endif


  // write out the correlator
  char outstr[ 256 ] ;
  sprintf( outstr , "%s" , outfile ) ;
  write_correlators( outstr , (const struct correlator**)bcorr ,
		     NSNS , NSNS ) ;

  free_corrs( bcorr , NSNS , NSNS ) ;

  free( S1 ) ;

  // rewind file and read header again
  rewind( prop.file ) ; read_propheader( &prop ) ;

  // tell us how long it all took
  print_time( ) ;

  return SUCCESS ;
}

