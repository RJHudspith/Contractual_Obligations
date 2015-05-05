/**
   @file io.c
   @file readers and alike
 */

#include "common.h"

#include "crc32.h"        // checksum calc
#include "GLU_bswap.h"    // byteswaps
#include "matrix_ops.h"   // matrix equivs
#include "spinor_ops.h"   // zero the spinor

// fill our spinor
static void
fill_spinor( struct spinor *__restrict S ,
	     void *tmp ,
	     const int ND1 , 
	     const int d1shift ,
	     const int tmpsize )
{
  if( tmpsize == sizeof( float complex ) ) {
    const float complex *ftmp = (const float complex*)tmp ;
    int d1d2 ;
    #pragma omp parallel for private(d1d2)
    for( d1d2 = 0 ; d1d2 < ( ND1 * ND1 ) ; d1d2++ ) {
      const int d1 = d1d2 / ND1 + d1shift ;
      const int d2 = d1d2 % ND1 + d1shift ;
      // unroll the matching
      colormatrix_equiv_f2d( (double complex*)S -> D[ d1 ][ d2 ].C ,
			     ftmp + d1d2 * NCNC ) ;
    }
  } else {
    const double complex *ftmp = (const double complex*)tmp ;
    int d1d2 ;
    #pragma omp parallel for private(d1d2)
    for( d1d2 = 0 ; d1d2 < ( ND1 * ND1 ) ; d1d2++ ) {
      const int d1 = d1d2 / ND1 + d1shift ;
      const int d2 = d1d2 % ND1 + d1shift ;
      // unroll the matching
      colormatrix_equiv( (double complex*)S -> D[ d1 ][ d2 ].C ,
			 ftmp + d1d2 * NCNC ) ;
    }
  }
  return ;
}

// the question is ... Who checks the checksum?
int
check_checksum( FILE *fprop )
{
  // spinsize
  const int spinsize = NSNS * NCNC ;

  double *prop_buf = malloc( spinsize*2 * sizeof( double ) ) ;

  // crc accumulators
  uint32_t CRCsum29 = 0 , CRCsum31 = 0 ;

  // ok so this must be a geometry/dirac thing
  int site;
  for( site = 0 ; site < LVOLUME ; site++ ) {
    if( fread( prop_buf , sizeof(double) , spinsize*2 , fprop ) != spinsize*2 ) {
      printf( "[IO] fread failure ... exiting \n" ) ;
      return FAILURE ;
    }
    // we need to know what end is up, this is actually kinda tricky
    DML_checksum_accum( &CRCsum29 , &CRCsum31 , site , 
			(char*)prop_buf , 
			2 * spinsize * sizeof(double) ) ;
  }
  
  // we are at the end of the file now so we can fscanf for the value I will have to grok the output
  uint32_t rCRCsum29 , rCRCsum31 ;
  if( fscanf( fprop , "%x %x" , &rCRCsum29 , &rCRCsum31 ) != 2 ) {
    printf( "[IO] file read failure at the checksums \n" ) ;
    free( prop_buf ) ;
    return FAILURE ;
  }

  if( CRCsum29 != rCRCsum29 || CRCsum31 != rCRCsum31 ) {
    printf( "[IO] mismatched checksums \n" ) ;
    printf( "[IO] Computed Checksums %x %x\n" , CRCsum29 , CRCsum31 ) ;
    printf( "[IO] File Read Checksums %x %x\n" , rCRCsum29 , rCRCsum31 ) ;
    free( prop_buf ) ;
    return FAILURE ;
  }

  return SUCCESS ;
}

// Read light propagator on a time slice 
// should we accumulate the checksum? Probably
static int 
read_chiralprop( struct propagator prop ,
		 struct spinor *S )
{
  const int spinsize = NCNC * NSNS ;

  // do we need to byte swap?
  const GLU_bool must_swap = prop.endian != WORDS_BIGENDIAN ? \
    GLU_TRUE : GLU_FALSE ;

  // single precision storage
  float complex *ftmp = NULL ;
  if( prop.precision == SINGLE ) {
    ftmp = ( float complex* )malloc( spinsize * sizeof( float complex ) ) ;
  }

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // depending on the precision, read in the data
    if( prop.precision == SINGLE ) {
      if( fread( ftmp , sizeof( float complex ) , spinsize , prop.file ) != 
	  spinsize ) {
	printf( "[IO] chiral propagator failure single prec (%d) \n" , i ) ;
	free( ftmp ) ;
	return FAILURE ;
      }
      if( must_swap ) bswap_32( 2 * spinsize , ftmp ) ;
      // cast our ftmp to double complex spinor
      fill_spinor( &S[i] , ftmp , NS , 0 , sizeof( float complex ) ) ;
    } else {
      // Read in propagator on a timeslice elements of our struct should be byte-compatible
      if( fread( S[i].D , sizeof( double complex ) , spinsize , prop.file ) != 
	  spinsize ) {
	printf( "[IO] chiral propagator failure double prec (%d) \n" , i ) ;
	return FAILURE ;
      }
      if( must_swap ) bswap_64( 2 * spinsize , S[i].D ) ;
    }
  }

  // free the single prec memory if used
  if( prop.precision == SINGLE ) {
    free( ftmp ) ;
  }

  return SUCCESS ;
}

// Read NRQCD propagator time slice 
static int 
read_nrprop( struct propagator prop , 
	     struct spinor *S )
{
  const int NR_NS = NS >> 1 ;
  const int spinsize = NCNC * NR_NS * NR_NS ;

  // read in site-by-site
  const GLU_bool must_swap = prop.endian != WORDS_BIGENDIAN ? \
    GLU_TRUE : GLU_FALSE ;

  // temporaries depending on prop precision
  double complex *tmp = NULL ;
  float complex *ftmp = NULL ;
  if( prop.precision == SINGLE ) {
    ftmp = malloc( spinsize * sizeof( float complex ) ) ;
  } else {
    tmp = malloc( spinsize * sizeof( double complex ) ) ;
  }

  // zeros our spinor over LCU
  spinor_zero( S ) ;

  int i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // single precision read 
    if( prop.precision == SINGLE ) {
      // Read in tslice 
      if( fread( ftmp , sizeof( float complex ) , spinsize , prop.file ) != 
	  spinsize ) {
	printf( "[IO] nrel propagator read failure \n" ) ;
	free( ftmp ) ;
	return FAILURE ;
      }
      if( must_swap ) bswap_32( 2 * spinsize , ftmp ) ;

      // fill the lower indices ( for example 4D :: 2,2 2,3 3,2 3,3 ) of propagator
      fill_spinor( &S[i] , ftmp , NR_NS , NR_NS , sizeof(float complex) ) ;
    } else {
      // Read in tslice 
      if( fread( tmp , sizeof(double complex) , spinsize , prop.file ) != 
	  spinsize ) {
	printf( "[IO] nrel propagator read failure \n" ) ;
	free( tmp ) ;
	return FAILURE ;
      }
      if( must_swap ) bswap_64( 2 * spinsize , tmp ) ;

      // fill the lower indices 2,2 2,3 3,2 3,3 of propagator
      fill_spinor( &S[i] , tmp , NR_NS , NR_NS , sizeof(double complex) ) ;
    }
  }

  // free the temp storage
  free( tmp ) ;
  free( ftmp ) ;

  return SUCCESS ;
}

// thin wrapper for propagator reading
int
read_prop( struct propagator prop ,
	   struct spinor *S )
{
  switch( prop.basis ) {
  case CHIRAL :
    return read_chiralprop( prop , S ) ;
  case NREL :
    return read_nrprop( prop , S ) ;
  }
  return FAILURE ;
}

