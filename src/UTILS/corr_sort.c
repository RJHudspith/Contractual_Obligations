/**
   @file corr_sort.c
   @brief type-agnostic merge and heap sorts
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "corr_sort.h" // for alphabetisation

// tuning for when to use insertion sort
#define INSERTION_TUNE ( 2 << 9 )

// insertion sort
static int
insertion_sort( void *list1 , void *list2 , 
		const size_t base1 , const size_t base2 ,
		const size_t N ,
		int(*compare)( const void *a , const void *b ) )
{
  char insert1[ base1 ] , insert2[ base2 ] ; // temporary storage
  int i , hole ;
  // standard insertion sort
  for( i = 1 ; i < N ; i++ ) {
    memcpy( insert1 , (char*)list1 + i*base1 , base1 ) ;
    memcpy( insert2 , (char*)list2 + i*base2 , base2 ) ;
    hole = i ;
    while( ( hole > 0 ) && ( compare( insert1 , (char*)list1 + base1*(hole - 1) ) ) ) {
      memcpy( (char*)list1 + base1*hole , (char*)list1 + base1*(hole - 1) , base1 ) ;
      memcpy( (char*)list2 + base2*hole , (char*)list2 + base2*(hole - 1) , base2 ) ;
      hole-- ;
    }
    memcpy( (char*)list1 + base1*hole , insert1 , base1 ) ;
    memcpy( (char*)list2 + base2*hole , insert2 , base2 ) ;
  }
  return 0 ;
}

// merge two lists
static int
merge( void *list1 , void *list2 ,
       const void *lower1 , const void *lower2 , const size_t low_size ,
       const void *upper1 , const void *upper2 , const size_t upp_size ,
       const size_t base1 , const size_t base2 , 
       int(*compare)( const void *a , const void *b ) )
{
  int a = 0 , b = 0 , counter = 0 ;
  while( a < low_size || b < upp_size ) {
    if( a < low_size ) {
      if( b < upp_size ) {
	if( compare( (char*)lower1 + base1*a , (char*)upper1 + base1*b ) ) {
	  memcpy( (char*)list1 + base1*counter , (char*)lower1 + base1*a , base1 ) ;
	  memcpy( (char*)list2 + base2*counter , (char*)lower2 + base2*a , base2 ) ;
	  a++ ;
	} else {
	  memcpy( (char*)list1 + base1*counter , (char*)upper1 + base1*b , base1 ) ;
	  memcpy( (char*)list2 + base2*counter , (char*)upper2 + base2*b , base2 ) ;
	  b++ ;
	}
      } else {
	memcpy( (char*)list1 + base1*counter , (char*)lower1 + base1*a , base1 ) ;
	memcpy( (char*)list2 + base2*counter , (char*)lower2 + base2*a , base2 ) ;
	a++ ;
      }
    } else {
      memcpy( (char*)list1 + base1*counter , (char*)upper1 + base1*b , base1 ) ;
      memcpy( (char*)list2 + base2*counter , (char*)upper2 + base2*b , base2 ) ;
      b++ ;
    }
    counter++ ;
  }
  free( (void*)lower1 ) ; free( (void*)upper1 ) ; 
  free( (void*)lower2 ) ; free( (void*)upper2 ) ;
  return 0 ;
}

/// heapsort
static void 
siftDown( void *list1 , void *list2 ,
	  size_t root , 
	  const size_t bottom ,
	  const size_t base1 , const size_t base2 ,
	  int(*compare)( const void *a , const void *b ) )
{
  size_t done = 0 , maxChild ; 
  while ( ( root*2 <= bottom ) && ( !done ) ) {
    if( root*2 == bottom ) {
      maxChild = root * 2 ; 
    } else if( compare( (char*)list1 + ( root * 2 + 1 ) *base1 , 
			(char*)list1 + root * 2 * base1 ) ) {
      maxChild = root * 2 ; 
    } else {
      maxChild = root * 2 + 1 ; 
    }
    if( compare( (char*)list1 + root*base1 , (char*)list1 + maxChild*base1 ) ) {
      swap( list1 , base1 , root , maxChild ) ;
      swap( list2 , base2 , root , maxChild ) ;
      root = maxChild ; 
    } else {
      done = 1 ;
    } 
  }
  return ;
}

// heapsorts the data wrt to x index
void 
heap_sort( void *list1 , void *list2 ,
	   const size_t base1 , const size_t base2 , 
	   const size_t N , 
	   int(*compare)( const void *a , const void *b ) )
{
  int i ; 
  for ( i = ( N >>1 ) ;  i >= 0 ;  i-- ) {
    siftDown( list1 , list2 , i , N-1 , base1 , base2 , compare ) ; 
  }
  
  // swapsies
  for ( i = N-1 ;  i >= 1 ;  i-- ) {
    // swap the x and the y
    swap( list1 , base1 , 0 , i ) ;
    swap( list2 , base2 , 0 , i ) ;
    // and go down the rabbit hole
    siftDown( list1 , list2 , 0 , i-1 , base1 , base2 , compare ) ; 
  }
  return ;
}

int
lt_dbl( const void *a , 
	const void *b )
{
  return (*(const double*)a < *(const double*)b) ; 
}

int
lt_int( const void *a , 
	const void *b )
{
  return (*(const int*)a < *(const int *)b ) ;
}

// merge_sort
int
merge_sort( void *list1 , void *list2 ,
	    const size_t base1 , const size_t base2 ,
	    const size_t size ,
	    int(*compare)( const void *a , const void *b ) )
{
  // recursively break the list up into even buckets
  if( size <= INSERTION_TUNE ) { 
    return insertion_sort( list1 , list2 , 
			   base1 , base2 , 
			   size , compare ) ; 
  }

  // break list into two sublists
  const size_t low_size = size >> 1 ;
  const size_t upp_size = size - low_size ;

  // this really shouldn't happen as the insertion
  // sort should capture such behaviour, is mostly
  // to remove the static anlayzer warning
  if( low_size < 1 || upp_size < 1 ) return 0 ;
  
  // temporary lists
  char *lower1 = malloc( low_size * base1 ) ;
  char *lower2 = malloc( low_size * base2 ) ;
  char *upper1 = malloc( upp_size * base1 ) ;
  char *upper2 = malloc( upp_size * base2 ) ;

  memcpy( lower1 , list1 , low_size * base1 ) ;
  memcpy( lower2 , list2 , low_size * base2 ) ;
  memcpy( upper1 , (char*)list1 + low_size * base1 , upp_size * base1 ) ;
  memcpy( upper2 , (char*)list2 + low_size * base2 , upp_size * base2 ) ;

  // recursively break up lower and upper into smaller lists
  merge_sort( lower1 , lower2 , low_size , base1 , base2 , compare ) ;
  merge_sort( upper1 , upper2 , upp_size , base1 , base2 , compare ) ;

  return merge( list1 , list2 , 
		lower1 , lower2 , low_size , 
		upper1 , upper2 , upp_size , 
		base1 , base2 , compare ) ;
}

// type-agnostic swap
void
swap( void *list , 
      const size_t base , 
      const size_t idx1 , 
      const size_t idx2 )
{
  char t[ base ] ;
  memcpy( t , (char*)list + idx1 * base , base ) ;
  memcpy( (char*)list + idx1 * base , (char*)list + idx2 * base , base ) ;
  memcpy( (char*)list + idx2 * base , t , base ) ;
}

#undef INSERTION_TUNE
