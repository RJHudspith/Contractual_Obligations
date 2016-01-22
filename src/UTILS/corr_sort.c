/**
   @file corr_sort.c
   @brief type-agnostic merge and heap sorts
 */
#include <stdlib.h> // malloc
#include <string.h> // memcpy

#include "corr_sort.h" // for alphabetisation

// tuning for when to use insertion sort
static size_t INSERTION_TUNE = 2 << 9 ;

// insertion sort
static int
insertion_sort( void *list1 , void *list2 , 
		const size_t base1 , const size_t base2 ,
		const size_t lower , const size_t upper ,
		int(*compare)( const void *a , const void *b ) )
{
  char insert1[ base1 ] , insert2[ base2 ] ; // temporary storage
  int i , hole ;
  // standard insertion sort
  for( i = lower ; i < upper ; i++ ) {
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
       void *tlist1 , void *tlist2 ,
       const size_t lower , const size_t middle , const size_t upper ,
       const size_t base1 , const size_t base2 , 
       int(*compare)( const void *a , const void *b ) )
{
  size_t i , a = lower , b = middle ;
  for( i = lower ; i < upper ; i++ ) {
    // potentially have an a candidate to put in list
    if( a < middle ) {
      // if it satisfies the compare poke it in
      if( b >= upper || compare( (char*)tlist1 + base1*a , 
				 (char*)tlist1 + base1*b ) ) {
	goto swapa ;
      // otherwise put a b in the list
      } else {
	goto swapb ;
      }
      // fill the list up with bees
    } else {
      goto swapb ;
    }
  swapa :
    memcpy( (char*)list1 + base1*i , (char*)tlist1 + base1*a , base1 ) ;
    memcpy( (char*)list2 + base2*i , (char*)tlist2 + base2*a , base2 ) ;
    a++ ;
    continue ;
  swapb :
    memcpy( (char*)list1 + base1*i , (char*)tlist1 + base1*b , base1 ) ;
    memcpy( (char*)list2 + base2*i , (char*)tlist2 + base2*b , base2 ) ;
    b++ ;
    continue ;
  }
  // copy merged list back into workspace lists
  memcpy( (char*)tlist1 + base1*lower , (char*)list1 + base1*lower , 
	  (upper-lower) * base1 ) ;
  memcpy( (char*)tlist2 + base2*lower , (char*)list2 + base2*lower , 
	  (upper-lower ) * base2 ) ;
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

// double comparison example
int
lt_dbl( const void *a , 
	const void *b )
{
  return (*(const double*)a < *(const double*)b) ; 
}

// int comparison exampled
int
lt_int( const void *a , 
	const void *b )
{
  return (*(const int*)a < *(const int *)b ) ;
}

// merge_sort
int
mergesort( void *list1 , void *list2 ,
	   void *tlist1 , void *tlist2 ,
	   const size_t lower , const size_t upper ,
	   const size_t base1 , const size_t base2 ,
	   int(*compare)( const void *a , const void *b ) )
{
  // recursively break the list up into even buckets
  if( (upper-lower) <= INSERTION_TUNE ) { 
    // sort temporary list as this is the one we merge back into 
    // the full lists, this is why pivot_tune < size
    return insertion_sort( tlist1 , tlist2 , base1 , base2 , 
			   lower , upper , compare ) ; 
  }

  // return when we only have one element in the list or if 
  // we are using the insertion sort to speed things up
  if( ( upper - lower ) == 1 ) return 0 ;

  // the middle of this list
  const size_t middle = ( lower + upper ) >> 1 ;

  // recursively break up lower and upper into smaller lists
  mergesort( list1 , list2 , tlist1 , tlist2 , 
	     lower , middle , base1 , base2 , compare ) ;

  mergesort( list1 , list2 , tlist1 , tlist2 , 
	     middle , upper , base1 , base2 , compare ) ;

  // merge the two lists
  return merge( list1 , list2 , tlist1 , tlist2 ,
		lower , middle , upper , base1 , base2 , compare ) ;
}

// new merge sort allocates temporary space and uses that
// for the sorted lists
int
merge_sort( void *list1 , void *list2 ,
	    const size_t base1 , const size_t base2 ,
	    const size_t size ,
	    int(*compare)( const void *a , const void *b ) )
{
  // allocate two temporary lists
  char *tlist1 = malloc( size * base1 ) ;
  char *tlist2 = malloc( size * base2 ) ;

  // copy into the two workspace lists
  memcpy( tlist1 , list1 , size * base1 ) ;
  memcpy( tlist2 , list2 , size * base2 ) ;

  // make sure we do at least one split
  if( size < INSERTION_TUNE ) {
    INSERTION_TUNE = size-1 ;
  }

  // do the merging
  mergesort( list1 , list2 , tlist1 , tlist2 ,
	     0 , size , base1 , base2 , compare ) ;

  // free the lists
  free( tlist1 ) ;
  free( tlist2 ) ;

  return 0 ;
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
