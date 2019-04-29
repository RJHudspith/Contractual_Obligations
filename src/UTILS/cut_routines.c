/*
    Copyright 2013 Renwick James Hudspith

    This file (cut_routines.c) is part of GLU.

    GLU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GLU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GLU.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   @file cut_routines.c
   @brief simple routines for creating our cut Lie-fields in momentum space

  The main routines for performing the cuts in momenta
  i.e. limiting the lattice artefacts by
  a ) Setting a maximum and minumum momentum range
  b ) Performing a cylinder cut so no hard, on-axis momenta contribute
  c ) Performing a conical cut which is like a weaker cylinder cut.

  I have three different cylinder cutting routines, they all agree. AS is the
  strictest. CJD and DF are the same.
 */

#include "common.h"

#include "geometry.h"  // for the get mom functions

/**
   @var rats
   @brief asymmetry ratios of the directions
 */
double rats[ ND ] ; // asymmetry ratios

/**
   @var small
   @brief smallest lattice side length
 */
static double norm = 0.0 ;
static int small = 1 ; // smallest lattice size

/**
   @enum list_creation
   @brief either add a momentum to the list or do not
 */
enum{ DO_NOT_ADD , ADD_TO_LIST } list_creation ;

/**
   @enum momenta_saved
   @brief have we saved a file with these momenta?
 */
enum{ MOMENTUM_CONFIG , NO_MOMENTUM_CONFIG } momenta_saved ;

// typical distinct permutation vectors for the DFT
static const int perms_zerz[ 4 ][ 3 ] = { { +1 , +1 , +1 } ,
					  { +1 , +1 , -1 } ,
					  { +1 , -1 , +1 } ,
					  { +1 , -1 , -1 } } ;
// permutations where there is a single zero in prototype
static const int perms_onez[ 6 ][ 3 ] = { { +1 ,  0 , +1 } ,
					  { +1 , +1 ,  0 } ,
					  {  0 , +1 , +1 } ,
					  { +1 ,  0 , -1 } ,
					  { +1 , -1 ,  0 } ,
					  {  0 , +1 , -1 } } ;
// permutations where there are two zeros in prototype
static const int perms_twoz[ 3 ][ 3 ] = { { +1 ,  0 ,  0 } ,
					  {  0 , +1 ,  0 } ,
					  {  0 ,  0 , +1 } } ;
// permutations where it is all zero is clearly just zero
static const int perms_thrz[ 1 ][ 3 ] = { {  0 ,  0 ,  0 } } ;


// calculate similar orbits i.e measure anisotropy
void
simorb_ratios( const int DIMS ,
	       const GLU_bool configspace )
{
  size_t mu ;
  if( configspace == GLU_TRUE ) {
    for( mu = 0 ; mu < (size_t)DIMS ; mu++ ) {
      rats[mu] = 1 ;
    }  
    return ;
  }
  // similar orbits in terms of the smallest dimension ...
  // find the smallest
  small = Latt.dims[0] ;
  for( mu = 1 ; mu < (size_t)DIMS ; mu++ ) {
    if( Latt.dims[ mu ] < (size_t)small ) {
      small = Latt.dims[ mu ] ;
    }
  }
  norm = 0. ;
  for( mu = 0 ; mu < (size_t)ND ; mu++ ) {
    rats[ mu ] = (double)small / (double)Latt.dims[ mu ] ; 
    fprintf( stdout , "[CUTS] rats :: %f %f %f\n" , 
	     rats[mu] , (double)small , (double)Latt.dims[ mu ] ) ;
    norm += ( rats[mu] * rats[mu] ) ;
  }
  norm = sqrt( norm ) ;
  return ;
}

////////// Cylinder cutting procedurals //////////
// gets the body diagonal vectors for our lattice
static inline void
get_diagonal( double n[ ND ] ,
	      const int i ,
	      const int DIMS )
{
  int mu , subvol = 1 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( mu < DIMS ) {
      n[ mu ] = ( ( i - i % subvol ) / subvol ) % 2 ;
      if( n[ mu ] == 0 ) {
	n[ mu ] -- ;
      }
      subvol *= 2 ;
    } else {// set it to 0?
      n[ mu ] = 0 ;
    }
  }
  return ;
}

/**
   @fn static int cylinder_DF( const double q[ ND ] , const int DIMS , const double cyl_width )
   @brief this one comes from de-forcrand <a href="http://arxiv.org/abs/hep-lat/0112043"> paper </a>
   @param q :: physical momentum
   @param DIMS :: dimension of the problem i.e. ND or ND-1
   @param cyl_width :: width of the cylinder

  and uses the formula 

  \f[
  \delta q = | q_\mu - ( q_\mu n_\mu ) * n_\mu |
  \f]

  and chooses momenta in a width :: \f$ width * 2 * \pi / ( small )\f$ 

  @return #ADD_TO_LIST or #DO_NOT_ADD
 */
// generic cylinder calc
static int 
cylinder_DF( const double q[ ND ] ,
	     const int DIMS ,
	     const double cyl_width )
{
  // test that this satisfies the correct cylinder cutting procedure
  const double norm = 1. / ( 2. * sqrt( DIMS ) ) ;
  const int diagonals = 2 << ( ND - 1 ) ;
 
  // generix for loop over the diagonals
  int mu ;
  double x[ ND ] ;
  for( mu = 0 ; mu < diagonals ; mu ++ ) {
    // inline for getting a diagonal for lexi order
    get_diagonal( x , mu , DIMS ) ;

    double scalar_prod = 0. ; 
    int nu ;
    for( nu = 0 ; nu < ND ; nu ++ ) {
      scalar_prod += q[ nu ] * x[ nu ] ;
    }
    scalar_prod *= norm ;

    double mod = 0.0 ; 
    for( nu = 0 ; nu < ND ; nu ++ ) {
      register const double temp = q[ nu ] - scalar_prod * x[ nu ] ;
      mod += temp * temp ;
    }

    if( sqrt( mod ) <= cyl_width ) {
      return ADD_TO_LIST ;
    }
  }
  return DO_NOT_ADD ;
}

/////////// HYPERCUBIC and PSQ cuts /////////////

// classic psq cut
static double 
gen_calc_psq( const int k[ ND ] , 
	      const int DIMS )
{
  int mu ;
  double res = 0. ;
  for( mu = 0 ; mu < DIMS ; mu++ ) {
    register double temp = rats[mu] * k[ mu ] ;
    res += temp * temp ;
  }
  return res ;
}

// standard hypercube cut 
static GLU_bool
gen_calc_hyp( const int k[ ND ] , 
	      const double root_mxmom , 
	      const int DIMS )
{
  int mu , flag = 0 ;
  for( mu = 0 ; mu < DIMS ; mu++ ) {
    if( fabs( (double)k[ mu ] * rats[ mu ] ) <= root_mxmom ) {
      flag ++ ;
    }
  }
  if( flag != DIMS ) {
    return GLU_FALSE ;
  } else {
    return GLU_TRUE ;
  }
}

// does it pass the test?
static GLU_bool
safe_momenta( const int n[ ND ] ,
	      struct cut_info CUTINFO ,
	      const int DIMS )
{
  const double latt_width = TWOPI * CUTINFO.cyl_width / small ;
  switch( CUTINFO.type ) {
  case HYPERCUBIC_CUT :
    return gen_calc_hyp( n , sqrt( CUTINFO.max_mom ) , DIMS ) ;
  case PSQ_CUT :
    if( gen_calc_psq( n , DIMS ) <= CUTINFO.max_mom ) {
      return GLU_TRUE ;
    }
    break ;
  case CYLINDER_CUT :
    if( gen_calc_psq( n , DIMS ) <= CUTINFO.max_mom ) {
      double k[ ND ] ; 
      compute_p( k , n , DIMS ) ;
      // CJD and DF are equivalent
      #ifdef CYLINDER_AS
      if( cylinder_AS( k , DIMS , latt_width ) == 1 ) {
	return GLU_TRUE ;
      }
      #else
      if( cylinder_DF( k , DIMS , latt_width ) == 1 ) {
	return GLU_TRUE ;
      }
      #endif
    }
    break ;
  default :
    if( gen_calc_psq( n , DIMS ) <= CUTINFO.max_mom ) {
      return GLU_TRUE ;
    }
  }
  return GLU_FALSE ;
}

// same thing with our new veclist
// gets the momentum list ..
static int 
get_mom_veclist( struct veclist_int *kept , 
		 const struct cut_info CUTINFO ,
		 const int LOOP ,
		 const int DIMS )
{ 
  simorb_ratios( DIMS , GLU_FALSE ) ;

  int mu ;
  fprintf( stdout , "[CUTS] " ) ;
  for( mu = 0 ; mu < DIMS ; mu++ ) {
    fprintf( stdout , "simorb[ %f ] " , rats[mu] ) ;
  }
  fprintf( stdout , "\n" ) ;

  // loop the correct length
  int i , in = 0 ; 
  for( i = 0 ; i < LOOP ; i++ ) {
    int n[ ND ] ;
    get_mom_pipi( n , i , DIMS ) ; 
    // if it is not allowed we do not add it in
    if( safe_momenta( n , CUTINFO , 
		      DIMS ) == GLU_TRUE ) {
      kept[in].idx = i ; 
      in++ ; 
    }
  }
  
  fprintf( stdout , "[CUTS] Kept vs reject %f vs %f \n" , 
	   100 * in/( double )LOOP, 
	   100 * ( LOOP - in )/( double )LOOP ) ; 

  // test the list in the shifted bz .. passes unless too many momenta are used
#pragma omp parallel for private(i)
  for( i = 0 ; i < in ; i++ ) {

    int k[ ND ] , sum[ ND ] ; 
    get_mom_pipi( k , kept[ i ].idx , DIMS ) ; 
    get_mom_pipi( sum , kept[ in - i - 1 ].idx , DIMS ) ; 

    int mu ;
    for( mu = 0 ; mu < DIMS ; mu++ ) {
      if( ( k[mu] + sum[mu] ) != 0 ) {
	fprintf( stderr , "NON +/- Symmetric Momenta @ %d \n" , i ) ;
	    
	fprintf( stderr , "(" ) ;
	for( mu = 0 ; mu < ND ; mu++ ){
	  fprintf( stderr , " %d " , k[mu] ) ;
	}
	fprintf( stderr , ") != (" ) ;
	for( mu = 0 ; mu < ND ; mu++ ) {
	  fprintf( stderr , " %d " , sum[mu] );
	}   
	fprintf( stderr , ")" ) ;
      }
    }
  }

//  This loop takes the momenta defined in the -Pi -> Pi BZ , that have been
//  accepted by the cut routine.
//  It calls get_mom_pipi which gets the value of the momenta in the 
// -pi , pi BZ
//  It then shifts this to the 0 -> 2Pi BZ that is the FFTW output. There are
//  other ways to do this, but I chose to use the periodicity of the FFT. 
//  Meaning that the -Pi -> 0 portion is the same as the Pi -> 2Pi bit.
//  We rewrite "kept" with the FFTW BZ position of the momenta hence we can just
//  call A[kept[i]].O[mu] will be a cut momenta.
#pragma omp parallel for private(i)
  for( i = 0 ; i < in ; i++ ) {
    //int k[ ND ] ;
    get_mom_pipi( kept[i].MOM , kept[i].idx , DIMS ) ; 
    kept[i].idx = get_site_pipiBZ( kept[i].MOM , DIMS ) ;
    // compute the fourier-mode squared
    int mu ;
    kept[i].nsq = 0 ;
    for( mu = 0 ; mu < DIMS ; mu++ ) {
      kept[i].nsq += ( kept[i].MOM[ mu ] * kept[i].MOM[ mu ] ) ;
    }
  }

  return in ; 
}

// does it pass the test?
static GLU_bool
passes_cuts( const int i , 
	     const struct cut_info CUTINFO ,
	     const int DIMS )
{
  // vector from the origin
  int n[ ND ] ;
  get_vec_from_origin( n , i , DIMS ) ;
  return safe_momenta( n , CUTINFO , DIMS ) ;
}

// gets the vector list
static size_t
get_veclist( struct veclist_int *kept , 
	     const struct cut_info CUTINFO , 
	     const size_t LOOP ,
	     const size_t DIMS )
{
  simorb_ratios( DIMS , GLU_TRUE ) ;
  size_t i , size = 0 ;
  for( i = 0 ; i < LOOP ; i++ ) {
    if( passes_cuts( i , CUTINFO , DIMS ) == GLU_TRUE ) {
      kept[ size ].idx = i ;
      get_vec_from_origin( kept[ size ].MOM , i , DIMS ) ;
      size_t mu ;
      kept[ size ].nsq = 0 ;
      for( mu = 0 ; mu < DIMS ; mu++ ) {
	kept[ size ].nsq += ( kept[ size ].MOM[ mu ] * kept[ size ].MOM[ mu ] ) ;
      }
      size ++ ;
    }
  }
  return size ;
}

// Computes the momentum list
struct veclist_int*
compute_veclist_int( int *list_size , 
		     const struct cut_info CUTINFO ,
		     const size_t DIMS ,
		     const GLU_bool CONFIGSPACE )
{
  int in[1] = { 1 } ;

  // loop up to DIMS
  int LOOP = 1 , i ;
  for( i = 0 ; i < (int)DIMS ; i++ ) {
    LOOP *= Latt.dims[ i ] ;
  }

  struct veclist_int *kept = malloc( LOOP * sizeof( struct veclist_int ) ) ;
  if( CONFIGSPACE == GLU_TRUE ) {
    in[0] = get_veclist( kept , CUTINFO , LOOP , DIMS ) ;
  } else {
    in[0] = get_mom_veclist( kept , CUTINFO , LOOP , DIMS ) ;
  }
  
  // stop 0-byte malloc
  if( in[0] < 1 ) {
    return NULL ;
  }
  
  *list_size = in[0] ;

  return kept ;
}

// Computes the momentum list
struct veclist*
compute_veclist( int *list_size , 
		 const struct cut_info CUTINFO ,
		 const size_t DIMS ,
		 const GLU_bool CONFIGSPACE )
{
  struct veclist_int *kept = compute_veclist_int( list_size , 
						  CUTINFO ,
						  DIMS ,
						  CONFIGSPACE ) ;
  struct veclist *list = malloc( (*list_size) * sizeof( struct veclist ) ) ;
  size_t i ;
  for( i = 0 ; i < (size_t)(*list_size) ; i ++  ) {
    list[i].nsq = (double)kept[i].nsq ;
    size_t mu ;
    for( mu = 0 ; mu < DIMS ; mu++ ) {
      list[i].MOM[mu] = (double)kept[i].MOM[mu] ;
    }
    list[i].idx = kept[i].idx ;
  }
  free( kept ) ;
  return list ;
}

// 
struct veclist*
DFT_mom_veclist( int *list_size ,
		 const struct cut_info CUTINFO ,
		 const int DIMS )
{
  // count the number of zeros in proto_mom
  size_t Nzeros = 0 ;
  size_t p , mu ;
  for( mu = 0 ; mu < ND-1 ; mu++ ) {
    if( fabs( CUTINFO.proto_mom[mu] ) < NRQCD_TOL ) {
      Nzeros++ ;
    }
  }

  size_t Nperm = 1 ;
  switch( Nzeros ) {
    case 0 : Nperm = 4 ; break ;
    case 1 : Nperm = 6 ; break ;
    case 2 : Nperm = 3 ; break ;
    case 3 : Nperm = 1 ; break ;
  default : Nperm = 1 ; break ;
  }
  
  const size_t Ntot = Nperm * CUTINFO.Nalphas ;
  // this is the list of independent momenta: note that it is
  // half what you would expect due to explicit p,-p symmetry

  struct veclist *list = calloc( Ntot , sizeof( struct veclist ) ) ;
  size_t theta ;
  for( theta = 0 ; theta < CUTINFO.Nalphas ; theta++ ) {
    for( p = 0 ; p < Nperm ; p++ ) {
      const size_t idx = p + Nperm*theta ;
      list[idx].idx = idx ;
      list[idx].nsq = 0.0 ;
      int loc_perms[ ND-1 ] ;
      switch( Nzeros ) {
      case 0 : memcpy( loc_perms , perms_zerz[p] , (ND-1)*sizeof(int) ) ; break ;
      case 1 : memcpy( loc_perms , perms_onez[p] , (ND-1)*sizeof(int) ) ; break ;
      case 2 : memcpy( loc_perms , perms_twoz[p] , (ND-1)*sizeof(int) ) ; break ;
      case 3 : memcpy( loc_perms , perms_thrz[p] , (ND-1)*sizeof(int) ) ; break ;
      default :
	memcpy( loc_perms , perms_thrz[p] , (ND-1)*sizeof(int) ) ; break ;
      }
      for( mu = 0 ; mu < (size_t)DIMS ; mu++ ) {
	list[idx].MOM[ mu ] = ( loc_perms[mu] * CUTINFO.thetas[theta] ) ;
	list[idx].nsq += list[idx].MOM[mu] * list[idx].MOM[mu] ;
      }
    }
  }
  list_size[ 0 ] = Ntot ;
  return list ;
}

// single wall-momentum
struct veclist*
wall_mom_veclist( int *list_size ,
		  const double sum_mom[ ND-1 ] ,
		  const int DIMS )
{
  struct veclist *list = calloc( 1 , sizeof( struct veclist ) ) ;
  size_t mu ;
  list[ 0 ].idx = 0 ;
  list[ 0 ].nsq = 0.0 ;
  for( mu = 0 ; mu < (size_t)DIMS ; mu++ ) {
    list[ 0 ].MOM[ mu ] = ( sum_mom[mu] ) ;
    list[ 0 ].nsq += list[0].MOM[mu] * list[0].MOM[mu] ;
  }
  list_size[ 0 ] = 1 ;
  return list ;
}

// passes a zero'd veclist 
struct veclist*
zero_veclist( int *list_size ,
	      const size_t DIMS ,
	      const GLU_bool CONFIGSPACE )
{
  struct veclist *list = calloc( 1 , sizeof( struct veclist ) ) ;
  list[ 0 ].idx = 0 ;
  size_t mu ; 
  for( mu = 0 ; mu < DIMS ; mu++ ) {
    list[ 0 ].MOM[ mu ] = 0 ;
  }
  list[ 0 ].nsq = 0 ;
  list_size[ 0 ] = 1 ;
  return list ;
}
