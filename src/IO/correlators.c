/**
   @file correlators.c
   @brief now holds the correlator definitions and IO
 */
#include "common.h"

#include "crc32c.h"       // crc 
#include "correlators.h"  // so that I can alphabetise
#include "cut_output.h"   // for write_mom_veclist

#ifdef HAVE_IMMINTRIN_H
  #include "SSE2_OPS.h"
#endif

// does a DFT with +/- M -> sum_mom
static int
DFT_correlator( struct measurements *M ,
		const size_t stride1 ,
		const size_t stride2 ,
		const size_t tshifted )

{
  size_t idx ;
#pragma omp for private(idx) schedule(dynamic)
  for( idx = 0 ; idx < stride1*stride2*M->nmom[0] ; idx++ ) {
    
    const size_t gidx = idx%(stride1*stride2) ;
    const size_t p = idx/(stride1*stride2) ;
    const size_t i = gidx/stride2 ;
    const size_t j = gidx%stride2 ;
    size_t site ;
    
    #ifdef HAVE_IMMINTRIN_H
    const __m128d *in = (const __m128d*)M -> in[ gidx ] ;
    const __m128d *eipx = (const __m128d*)M -> dft_mom[ p ] ;
    register __m128d s = _mm_setzero_pd() ;
    for( site = 0 ; site < LCU ; site++ ) {
      s = _mm_add_pd( s , SSE2_MUL( *in , *eipx ) ) ;
      in++ ; eipx++ ;
    }
    // cast into void
    double complex sum = 0.0 ;
    _mm_store_pd( (void*)&sum , s ) ;
    #else
    register double complex sum = 0.0 ;
    for( site = 0 ; site < LCU ; site++ ) {
      sum +=
	M -> in[ gidx ][ site ] * M -> dft_mom[ p ][ site ] ;
    }
    #endif
    M -> corr[ i ][ j ].mom[ p ].C[ tshifted ] = sum ;
  }
  return SUCCESS ;
}

// FFT using FFTW OR we just do the zero momentum sum
static int
FFT_correlator( struct measurements *M ,
		const size_t stride1 ,
		const size_t stride2 ,
		const size_t tshifted )
{
#ifdef HAVE_FFTW3_H
  fftw_plan *fwd = (fftw_plan*)M -> forward ;
#endif
  // momentum projection
  size_t idx ;
#pragma omp for private(idx) schedule(dynamic)
  for( idx = 0 ; idx < stride1*stride2 ; idx++ ) {
    const size_t i = idx/stride2 ;
    const size_t j = idx%stride2 ;
    size_t p ;
    #ifdef HAVE_FFTW3_H
    fftw_execute( fwd[ idx ] ) ;
    for( p = 0 ; p < (size_t)M -> nmom[ 0 ] ; p++ ) {
      M -> corr[ i ][ j ].mom[ p ].C[ tshifted ] =
	M -> out[ idx ][ M -> list[ p ].idx ] ;
    }
    #else
    register double complex sum = 0.0 ;
    for( p = 0 ; p < LCU ; p++ ) {
      sum += M -> in[ idx ][ p ] ;
    }
    M -> corr[ i ][ j ].mom[ 0 ].C[ tshifted ] = sum ;
    #endif
  }
  return SUCCESS ;
}

// allocation of dispersion relation correlation function
struct mcorr **
allocate_momcorrs( const size_t length1 , 
		   const size_t length2 ,
		   const size_t nmom )
{
  // check if stuff is non-zero
  if( length1 == 0 || length2 == 0 || nmom == 0 ) {
    fprintf( stderr , "[IO] corr lengths are 0 :: ( %zu , %zu , %zu ) \n" ,
	     length1 , length2 , nmom ) ;
    return NULL ;
  }
  struct mcorr **mcorr = malloc( length1 * sizeof( struct mcorr* ) ) ;
  size_t i , j , p ;
  for( i = 0 ; i < length1 ; i++ ) {
    mcorr[ i ] = malloc( length2 * sizeof( struct mcorr ) ) ;
    for( j = 0 ; j < length2 ; j++ ) {
      mcorr[ i ][ j ].mom = malloc( nmom * sizeof( struct correlator ) ) ;
      for( p = 0 ; p < nmom ; p++ ) {
	mcorr[ i ][ j ].mom[ p ].C = malloc( LT * sizeof( double complex ) ) ;
      }
    }
  }
  return mcorr ;
}

// compute the momentum-projected correlation function
int
compute_correlator( struct measurements *M , 
		    const size_t stride1 , 
		    const size_t stride2 ,
		    const size_t tshifted )
{
  if( M -> configspace == GLU_TRUE ) {    
    size_t idx ;
    #pragma omp for private(idx) schedule(dynamic)
    for( idx = 0 ; idx < stride1*stride2 ; idx++ ) {
      const size_t i = idx/stride2 ;
      const size_t j = idx%stride2 ;
      size_t p ;
      for( p = 0 ; p < (size_t)M -> nmom[ 0 ] ; p++ ) {
	M -> corr[ i ][ j ].mom[ p ].C[ tshifted ] =
	  M -> in[ idx ][ M -> list[ p ].idx ] ;
      }
    }
  } else if( M -> is_wall_mom == GLU_TRUE ||
	     M -> is_dft == GLU_TRUE ) {
    DFT_correlator( M , stride1 , stride2 , tshifted ) ;
  } else {
    FFT_correlator( M , stride1 , stride2 , tshifted ) ;
  }
  return SUCCESS ;
}

// momcorr freer
void
free_momcorrs( struct mcorr **mcorr , 
	       const size_t length1 ,
	       const size_t length2 ,
	       const size_t nmom ) 
{
  size_t i , j , p ;
  for( i = 0 ; i < length1 ; i++ ) {
    for( j = 0 ; j < length2 ; j++ ) {
      for( p = 0 ; p < nmom ; p++ ) {
	free( mcorr[ i ][ j ].mom[ p ].C ) ;
      }
      free( mcorr[ i ][ j ].mom ) ;
    }
    free( mcorr[ i ] ) ;
  }
  free( mcorr ) ;
  return ;
}

// write the full correlator matrix
void
write_momcorr( const char *outfile ,
	       const struct mcorr **corr ,
	       const struct veclist *list ,
	       const double twist[ ND ] ,
	       const size_t NSRC ,
	       const size_t NSNK ,
	       const int *nmom , 
	       const char *type )
{
  // write out the correlator
  char outstr[ strlen(outfile)+strlen(type)+2 ] ;
  if( !strcmp( type , "" ) ) { 
    sprintf( outstr , "%s" , outfile ) ;
  } else {
    sprintf( outstr , "%s.%s" , outfile , type ) ;
  }

  fprintf( stdout , "[IO] writing correlation matrix to %s \n" , outstr ) ;

  FILE *output_file = fopen( outstr , "wb" ) ;

  uint32_t magic[ 1 ] = { CORR_MAGIC } ; // THIS SPELLS CORR in ascii

  fwrite( magic , sizeof( uint32_t ) , 1 , output_file ) ;

  write_mom_veclist( output_file , twist , nmom , list , ND-1 ) ;

  uint32_t NMOM[ 1 ] = { nmom[0] } ;
  
  fwrite( NMOM , sizeof( uint32_t ) , 1 , output_file ) ;

  uint32_t L0[ 1 ] = { LT } , cksuma = 0 , cksumb = 0 ;

  size_t p ;
  for( p = 0 ; p < (size_t)nmom[0] ; p++ ) {
    
    uint32_t NGSRC[ 1 ] = { (uint32_t)NSRC } ;
    uint32_t NGSNK[ 1 ] = { (uint32_t)NSNK } ;
    
    fwrite( NGSRC , sizeof( uint32_t ) , 1 , output_file ) ;
    fwrite( NGSNK , sizeof( uint32_t ) , 1 , output_file ) ;

    size_t GSRC , GSNK ;
    for( GSRC = 0 ; GSRC < NSRC ; GSRC++ ) {
      for( GSNK = 0 ; GSNK < NSNK ; GSNK++ ) {
	fwrite( L0 , sizeof( uint32_t ) , 1 , output_file ) ;
	fwrite( corr[GSRC][GSNK].mom[p].C , sizeof( double complex ) , LT , output_file ) ; 
	// accumulate the newer, fancier checksum
	DML_checksum_accum_crc32c( &cksuma , &cksumb ,
				   p + nmom[0] * ( GSNK + NSNK * GSRC ) , 
				   corr[GSRC][GSNK].mom[p].C , 
				   sizeof( double complex ) * LT ) ;
      }
    }
  }
  // write out both checksums
  uint32_t csum[ 2 ] = { cksuma , cksumb } ;
  fwrite( csum , sizeof( uint32_t ) , 2 , output_file ) ;

  fclose( output_file ) ;

  return ;
}

// little wrapper function
void
write_momcorr_WW( const struct measurements M ,
		  const char *outfile ,
		  const size_t NSRC ,
		  const size_t NSNK )
{
  // write out the ND-1 momentum-injected correlator and maybe the wall
  write_momcorr( outfile , (const struct mcorr**)M.corr , M.list ,
		 M.sum_twist , NSRC , NSNK , M.nmom , "" ) ;
  write_momcorr( outfile , (const struct mcorr**)M.wwcorr , M.wwlist ,
		 M.sum_twist , NSRC , NSNK , M.wwnmom , "ww" ) ;
  return ;
}
