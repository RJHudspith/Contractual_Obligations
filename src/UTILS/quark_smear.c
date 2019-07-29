/**
   @file quark_smear.c
   @brief perform quark source and sink smearing
 */
#include "common.h"

#include "geometry.h"       // compute_spacing()
#include "grad_2.h"         // gradsq()
#include "halfspinor_ops.h" // zero_halfspinor
#include "spinor_ops.h"     // spinor_zero_site()

// performs N iterations of some smearing when the gauge
// field is fixed to Coulomb gauge we don't need to multiply
// by the link matrices, making this cheaper
#if 0
static void
CGquark_smear( struct halfspinor *S ,
	       struct halfspinor *S1 ,
	       const size_t t ,
	       const struct source_info Source )
{
  const double fac = Source.smalpha/Source.Nsmear ;
  size_t n ;
  for( n = 0 ; n < Source.Nsmear ; n++ ) {
    size_t i ;
    #pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      S1[i] = S[i] ;
      size_t mu ;
      for( mu = 0 ; mu < ND-1 ; mu++ ) {
	const size_t xpm = lat[i].neighbor[mu] ;
	const size_t xmm = lat[i].back[mu] ;
	halfspinor_Saxpy( &S1[i] , S[xpm] , fac ) ;
	halfspinor_Saxpy( &S1[i] , S[xmm] , fac ) ;
      }
      halfspinor_Saxpy( &S1[i] , S[i] , -2*(ND-1)*fac ) ;
    }
    // set S = S1, S1 = S
    #pragma omp single
    {
      struct halfspinor *t = S ;
      S  = S1 ;
      S1 = t ;
    }
  }
  return ;
}
#endif

// perform some sink smearing
void
sink_smear( struct spinor **S ,
	    struct spinor **S1 ,
	    const size_t t ,
	    const struct cut_info CUTINFO ,
	    const size_t Np )
{
  if( CUTINFO.nsink == 0 || S1 == NULL || lat == NULL ) return ;

  const double fac    = CUTINFO.sink_alpha/(double)CUTINFO.nsink ;
  const double fac_U0 = fac / CUTINFO.sink_U0 ;

  // loop number of smearing iterations
  size_t r ;
  for( r = 0 ; r < CUTINFO.nsink ; r++ ) {

    // loop combined spatial volume and number of props
    size_t iNp ;
    #pragma omp for private(iNp)
    for( iNp = 0 ; iNp < LCU*Np ; iNp++ ) {

      const size_t i = iNp%LCU ;
      const size_t n = iNp/LCU ;
      
      struct spinor *PS  = S[n] ;
      struct spinor *PS1 = S1[n] ;
      
      PS1[i] = PS[i] ;

      struct spinor tmp ;
      size_t mu ;
      const size_t Uidx = i + LCU*t ;
      for( mu = 0 ; mu < ND-1 ; mu++ ) {
	const size_t xpmu = lat[i].neighbor[mu] ;
	const size_t xmmu = lat[i].back[mu] ;
	const size_t Ubck = lat[ Uidx ].back[mu] ;

	gauge_spinor( &tmp , lat[ Uidx ].O[mu] , PS[xpmu] ) ;
	spinor_Saxpy( &PS1[i] , fac_U0 , tmp ) ;
	  
	gaugedag_spinor( &tmp , lat[ Ubck ].O[mu] , PS[xmmu] ) ;
	spinor_Saxpy( &PS1[i] , fac_U0 , tmp ) ;
      }
      spinor_Saxpy( &PS1[i] , -2*(ND-1)*fac , PS[i] ) ;
    }
    // shallow pointer swap
    #pragma omp single
    {
      struct spinor *temp ;
      size_t n ;
      for( n = 0 ; n < Np ; n++ ) {
	temp  = S[n] ;
	S[n]  = S1[n] ;
	S1[n] = temp ;
      }
    }
  }

  return ;
}

// performs N iterations of some smearing
// works like
// S = exp( asmear\grad^2 ) S
// by approximating the exp as the iteration
// S = ( 1 + asmear grad^2/Nsmear )^Nsmear S
// This is very like the C_0 term of NRQCD
void
source_smear( struct halfspinor *S ,
	      struct halfspinor *S1 ,
	      const size_t t ,
	      const struct source_info Source )
{
  const double fac = Source.smalpha/Source.Nsmear ;
  size_t n ;
  for( n = 0 ; n < Source.Nsmear ; n++ ) {
    size_t i ;
    #pragma omp for private(i)
    for( i = 0 ; i < LCU ; i++ ) {
      struct halfspinor der ;
      S1[i] = S[i] ;
      gradsq( &der , S , i , t ) ;
      halfspinor_Saxpy( &S1[i] , der , fac ) ;
    }
    // set S = S1, S1 = S
    #pragma omp single
    {
      struct halfspinor *t = S ;
      S  = S1 ;
      S1 = t ;
    }
  }
  return ;
}

// sum over spatial indices a spinor
void
sum_spatial_sep( struct spinor *SUM_r2 ,
		 const struct measurements M ,
		 const size_t site1 )
{
  size_t n , r ;
#if (defined __AVX__) && (ND==4)
  __m256d sum[M.Nprops][ 8*NCNC ] , *pt ; // spinor
  double *pB ;
  size_t k ;
  // set the sum to zero
  for( n = 0 ; n < M.Nprops ; n++ ) {
    pt = sum[n] ;
    for( k = 0 ; k < NCNC ; k++ ) {
      *pt = _mm256_setzero_pd() ; pt++ ;
      *pt = _mm256_setzero_pd() ; pt++ ;
      *pt = _mm256_setzero_pd() ; pt++ ;
      *pt = _mm256_setzero_pd() ; pt++ ;
      *pt = _mm256_setzero_pd() ; pt++ ;
      *pt = _mm256_setzero_pd() ; pt++ ;
      *pt = _mm256_setzero_pd() ; pt++ ;
      *pt = _mm256_setzero_pd() ; pt++ ;
    }
  }
  // sum over each spatial separation
  for( r = 0 ; r < (size_t)M.NR ; r++ ) {
    const size_t site2 = compute_spacing( M.rlist[r].MOM ,
					  site1 , ND-1 ) ;

    for( n = 0 ; n < M.Nprops ; n++ ) {
      pB = (double*)M.S[n][site2].D ;
      register __m256d B ;
      pt = sum[n] ;
      // inline avx sum
      for( k = 0 ; k < NCNC ; k++ ) {
	B = _mm256_loadu_pd( pB ) ; pB+=4 ;
	*pt = _mm256_add_pd( *pt , B ) ; pt++ ;
	B = _mm256_loadu_pd( pB ) ; pB+=4 ;
	*pt = _mm256_add_pd( *pt , B ) ; pt++ ;
	B = _mm256_loadu_pd( pB ) ; pB+=4 ;
	*pt = _mm256_add_pd( *pt , B ) ; pt++ ;
	B = _mm256_loadu_pd( pB ) ; pB+=4 ;
	*pt = _mm256_add_pd( *pt , B ) ; pt++ ;
	B = _mm256_loadu_pd( pB ) ; pB+=4 ;
	*pt = _mm256_add_pd( *pt , B ) ; pt++ ;
	B = _mm256_loadu_pd( pB ) ; pB+=4 ;
	*pt = _mm256_add_pd( *pt , B ) ; pt++ ;
	B = _mm256_loadu_pd( pB ) ; pB+=4 ;
	*pt = _mm256_add_pd( *pt , B ) ; pt++ ;
	B = _mm256_loadu_pd( pB ) ; pB+=4 ;
	*pt = _mm256_add_pd( *pt , B ) ; pt++ ;
      }
    }
  }
  // copy sums back unaliged copy
  for( n = 0 ; n < M.Nprops ; n++ ) {
    pt = sum[n] ;
    pB = (double*)SUM_r2[n].D ;
    for( k = 0 ; k < NCNC ; k++ ) {
      _mm256_storeu_pd( pB , *pt ) ; pB+=4 ; pt++ ;
      _mm256_storeu_pd( pB , *pt ) ; pB+=4 ; pt++ ;
      _mm256_storeu_pd( pB , *pt ) ; pB+=4 ; pt++ ;
      _mm256_storeu_pd( pB , *pt ) ; pB+=4 ; pt++ ;
      _mm256_storeu_pd( pB , *pt ) ; pB+=4 ; pt++ ;
      _mm256_storeu_pd( pB , *pt ) ; pB+=4 ; pt++ ;
      _mm256_storeu_pd( pB , *pt ) ; pB+=4 ; pt++ ;
      _mm256_storeu_pd( pB , *pt ) ; pB+=4 ; pt++ ;
    }
  }
#else
  // set the sum to zero
  for( n = 0 ; n < M.Nprops ; n++ ) {
    spinor_zero_site( &SUM_r2[ n ] ) ;
  }
  // sum over each spatial separation
  for( r = 0 ; r < (size_t)M.NR ; r++ ) {
    const size_t site2 = compute_spacing( M.rlist[r].MOM , site1 ,
					  ND-1 ) ;
    for( n = 0 ; n < M.Nprops ; n++ ) {
      add_spinors( &SUM_r2[n] , M.S[n][site2] ) ;
      // we could really be creative here and put all sorts of functions in.
      // This one below weights further away points exponentially, although the
      // gauge field kinda does that already. 
      // spinor_Saxpy( &SUM_r2[n] , exp( -M.rlist[r].nsq*0.1 ) , M.S[n][site2] ) ;
    }
  }
#endif

  return ;
}
