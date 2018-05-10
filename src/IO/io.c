/**
   @file io.c
   @file readers and alike
 */
#include "common.h"

#include "crc32.h"        // checksum calc
#include "gammas.h"       // gamma matrix technology
#include "GLU_bswap.h"    // byteswaps
#include "io.h"           // alphabetising
#include "matrix_ops.h"   // matrix equivs
#include "polyakov.h"     // static quark computation
#include "spinor_ops.h"   // zero the spinor

// fill our spinor
static void
fill_spinor( struct spinor *__restrict S ,
	     void *tmp ,
	     const size_t ND1 , 
	     const size_t d1shift ,
	     const size_t d2shift ,
	     const size_t tmpsize )
{
  if( tmpsize == sizeof( float complex ) ) {
    const float complex *ftmp = (const float complex*)tmp ;
    size_t d1d2 ;
    for( d1d2 = 0 ; d1d2 < ( ND1 * ND1 ) ; d1d2++ ) {
      const size_t d1 = d1d2 / ND1 + d1shift ;
      const size_t d2 = d1d2 % ND1 + d2shift ;
      // unroll the matching
      colormatrix_equiv_f2d( (double complex*)S -> D[ d1 ][ d2 ].C ,
			     ftmp + d1d2 * NCNC ) ;
    }
  } else {
    const double complex *ftmp = (const double complex*)tmp ;
    size_t d1d2 ;
    for( d1d2 = 0 ; d1d2 < ( ND1 * ND1 ) ; d1d2++ ) {
      const size_t d1 = d1d2 / ND1 + d1shift ;
      const size_t d2 = d1d2 % ND1 + d2shift ;
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
  const size_t spinsize = NSNS * NCNC ;

  double *prop_buf = malloc( spinsize*2 * sizeof( double ) ) ;

  // crc accumulators
  uint32_t CRCsum29 = 0 , CRCsum31 = 0 ;

  // ok so this must be a geometry/dirac thing
  size_t site;
  for( site = 0 ; site < LVOLUME ; site++ ) {
    if( fread( prop_buf , sizeof(double) , spinsize*2 , fprop ) != spinsize*2 ) {
      fprintf( stderr , "[IO] fread failure ... exiting \n" ) ;
      free( prop_buf ) ;
      return FAILURE ;
    }
    // we need to know what end is up, this is actually kinda tricky
    DML_checksum_accum( &CRCsum29 , &CRCsum31 , site , 
			(char*)prop_buf , 
			2 * spinsize * sizeof(double) ) ;
  }

  free( prop_buf ) ;
  
  // we are at the end of the file now so we can fscanf for the value I will have to grok the output
  uint32_t rCRCsum29 , rCRCsum31 ;
  if( fscanf( fprop , "%x %x" , &rCRCsum29 , &rCRCsum31 ) != 2 ) {
    fprintf( stderr , "[IO] file read failure at the checksums \n" ) ;
    return FAILURE ;
  }

  if( CRCsum29 != rCRCsum29 || CRCsum31 != rCRCsum31 ) {
    fprintf( stderr , "[IO] mismatched checksums \n" ) ;
    fprintf( stderr , "[IO] Computed Checksums %x %x\n" , 
	     CRCsum29 , CRCsum31 ) ;
    fprintf( stderr , "[IO] File Read Checksums %x %x\n" , 
	     rCRCsum29 , rCRCsum31 ) ;
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
  if( LCU%IO_NBLOCK != 0 ) {
    fprintf( stderr , "[IO] Local volume not a multiple of IO blocking"
	     " (IO_BLOCK %d) vs. (LCU %zu)\n" , IO_NBLOCK , LCU ) ;
    return FAILURE ;
  }
  const size_t spinsize = NCNC * NSNS ;

  // do we need to byte swap?
  const GLU_bool must_swap = prop.endian != WORDS_BIGENDIAN ? \
    GLU_TRUE : GLU_FALSE ;

  // single precision storage
  float complex *ftmp = NULL ;
  if( prop.precision == SINGLE ) {
    ftmp = ( float complex* )malloc( IO_NBLOCK*spinsize * sizeof( float complex ) ) ;
  }

  size_t i , mu ;
  for( i = 0 ; i < LCU ; i += IO_NBLOCK ) {
    // depending on the precision, read in the data
    if( prop.precision == SINGLE ) {
      if( fread( ftmp , sizeof( float complex ) ,
		 IO_NBLOCK*spinsize , prop.file ) != IO_NBLOCK*spinsize ) {
	fprintf( stderr , "[IO] chiral propagator failure single prec (%zu)\n" ,
		 i ) ;
	free( ftmp ) ;
	return FAILURE ;
      }
      if( must_swap ) bswap_32( 2 * IO_NBLOCK * spinsize , ftmp ) ;
      // cast our ftmp to double complex spinor
      for( mu = 0 ; mu < IO_NBLOCK ; mu++ ) {
	fill_spinor( &S[i+mu] , ftmp+mu*spinsize ,
		     NS , 0 , 0 , sizeof( float complex ) ) ;
      }
    } else {
      // Read in propagator on a timeslice elements of our struct should be byte-compatible
      for( mu = 0 ; mu < IO_NBLOCK ; mu++ ) {
	if( fread( S[i+mu].D , sizeof( double complex ) , spinsize , prop.file ) != 
	    spinsize ) {
	  fprintf( stderr , "[IO] chiral propagator failure double prec (%zu)\n" ,
		   i ) ;
	  return FAILURE ;
	}
	if( must_swap ) bswap_64( 2 * spinsize , S[i+mu].D ) ;
      }
    }
  }
  // free the possibly allocated floating-point storage
  free( ftmp ) ;

  return SUCCESS ;
}

// Read NRQCD propagator time slice 
static int 
read_nrprop( struct propagator prop , 
	     struct spinor *S , 
	     const proptype basis )
{
  // non rel prop is only the top corner
  const size_t NR_NS = NS >> 1 ;
  const size_t spinsize = NCNC * NR_NS * NR_NS ;

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

  size_t i ;
  for( i = 0 ; i < LCU ; i++ ) {
    // single precision read 
    if( prop.precision == SINGLE ) {
      // Read in tslice 
      if( fread( ftmp , sizeof( float complex ) , spinsize , prop.file ) != 
	  spinsize ) {
	fprintf( stderr , "[IO] nrel propagator read failure \n" ) ;
	free( ftmp ) ;
	return FAILURE ;
      }
      if( must_swap ) bswap_32( 2 * spinsize , ftmp ) ;

      // convention dictates that forward is the bottom right and backward is top left
      if( basis == NREL_FWD ) {
	fill_spinor( &S[i] , ftmp , NR_NS , NR_NS , NR_NS , sizeof(float complex) ) ;
      } else {
	fill_spinor( &S[i] , ftmp , NR_NS , 0 , 0 , sizeof(float complex) ) ;
      }
      //
    } else {
      // Read in tslice 
      if( fread( tmp , sizeof(double complex) , spinsize , prop.file ) != 
	  spinsize ) {
	fprintf( stderr , "[IO] nrel propagator read failure \n" ) ;
	free( tmp ) ;
	return FAILURE ;
      }
      if( must_swap ) bswap_64( 2 * spinsize , tmp ) ;

      // convention dictates that forward is the bottom right and backward is top left
      if( basis == NREL_FWD ) {
	fill_spinor( &S[i] , tmp , NR_NS , NR_NS , NR_NS , sizeof(double complex) ) ;
      } else {
	fill_spinor( &S[i] , tmp , NR_NS , 0 , 0 , sizeof(double complex) ) ;
      }
      //
    }
  }

  // free the temp storage
  free( tmp ) ;
  free( ftmp ) ;

  return SUCCESS ;
}

int
compute_nrprop( struct propagator prop , 
		struct spinor *S ,
		const size_t t )
{
  // copy from the prop
  size_t i ;
  spinor_zero( S ) ;
  for( i = 0 ; i < LCU ; i++ ) {
    if( prop.NRQCD.backward == GLU_TRUE ) {
      colormatrix_equiv( (void*)S[i].D[0][0].C , (void*)prop.H[i+LCU*t].D[0] ) ;
      colormatrix_equiv( (void*)S[i].D[0][1].C , (void*)prop.H[i+LCU*t].D[1] ) ;
      colormatrix_equiv( (void*)S[i].D[1][0].C , (void*)prop.H[i+LCU*t].D[2] ) ;
      colormatrix_equiv( (void*)S[i].D[1][1].C , (void*)prop.H[i+LCU*t].D[3] ) ;
    } else {
      colormatrix_equiv( (void*)S[i].D[2][2].C , (void*)prop.H[i+LCU*t].D[0] ) ;
      colormatrix_equiv( (void*)S[i].D[2][3].C , (void*)prop.H[i+LCU*t].D[1] ) ;
      colormatrix_equiv( (void*)S[i].D[3][2].C , (void*)prop.H[i+LCU*t].D[2] ) ;
      colormatrix_equiv( (void*)S[i].D[3][3].C , (void*)prop.H[i+LCU*t].D[3] ) ;
    }
  }
  return SUCCESS ;
}

// read timeslice above of Nprops
int
read_ahead( struct propagator *prop ,
	    struct spinor **S , 
	    int *error_code ,
	    const size_t Nprops ,
	    const size_t t )
{ 
  // loops for IO
#pragma omp master
  {
    if( read_prop( prop[0] , S[0] , t ) == FAILURE ) {
      *error_code = FAILURE ;
    }
  }
  size_t mu ;
  for( mu = 1 ; mu < Nprops ; mu++ ) {
#pragma omp single nowait
    {
      if( read_prop( prop[mu] , S[mu] , t ) == FAILURE ) {
	*error_code = FAILURE ;
      }
    }
  }
  return 0 ;
}

// thin wrapper for propagator reading
int
read_prop( struct propagator prop ,
	   struct spinor *S ,
	   const size_t t )
{
  switch( prop.basis ) {
  case CHIRAL :
    return read_chiralprop( prop , S ) ;
  case NREL_CORR :
    return compute_nrprop( prop , S , t ) ;
  case NREL_FWD :
  case NREL_BWD :
    return read_nrprop( prop , S , prop.basis ) ;
  }
  return FAILURE ;
}

