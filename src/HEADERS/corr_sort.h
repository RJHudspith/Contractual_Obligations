/**
   @file corr_sort.h
   @brief prototype declarations for generic sorting routines
 */
#ifndef CORR_SORT_H
#define CORR_SORT_H

/**
   @fn void heap_sort( void *list1 , void *list2 , const size_t base1 , const size_t base2 , const size_t N , int(*compare)( const void *a , const void *b ) )
   @brief perform a heap sort of list 1, swapping list2 accordingly
 */
void 
heap_sort( void *list1 , void *list2 ,
	   const size_t base1 , const size_t base2 , 
	   const size_t N , 
	   int(*compare)( const void *a , const void *b ) ) ;

/**
   @fn int lt_dbl( const void *a , const void *b )
   @brief less than double function
 */
int
lt_dbl( const void *a , 
	const void *b ) ;

/**
   @fn int lt_int( const void *a , const void *b )
   @brief less than integer function
 */
int
lt_int( const void *a , 
	const void *b ) ;

/**
   @fn int merge_sort( void *list1 , void *list2 , const size_t base1 , const size_t base2 , const size_t size , int(*compare)( const void *a , const void *b ) )
   @brief merge sort @list1 swapping @list2 accordingly
 */
int
merge_sort( void *list1 , void *list2 ,
	    const size_t base1 , const size_t base2 ,
	    const size_t size ,
	    int(*compare)( const void *a , const void *b ) ) ;

/**
   @fn void swap( void *list , const size_t base , const size_t idx1 , const size_t idx2 )
   @brief swap elements of list in idx1 and idx2
 */
void
swap( void *list , 
      const size_t base , 
      const size_t idx1 , 
      const size_t idx2 ) ;

#endif
