/**
   @file progress_bar.c
   @brief little progress bar to tell us how we are doing
 */
#include <stdio.h>

static size_t prev_nequals = 0 ;

void
progress_bar( const size_t state ,
	      const size_t goal )
{
  // idea is to have a buffer of some chars
  if( state == 0 ) {
    fprintf( stdout , "\n____________________________"
	     "__________________________\n\n" ) ;
    fprintf( stdout , "0%%         25%%          50%%"
	     "          75%%          100%%\n" ) ;
  } else if( state == ( goal - 1 ) ) {
    fprintf( stdout , "\n____________________________"
	     "__________________________\n" ) ;
    fprintf( stdout , "\n\n" ) ;
  } else {
    // tell us how many equals to put down
    const size_t nequals = (size_t)( state*( 62./(double)goal ) ) ;
    size_t i ;
    for( i = 0 ; i < ( nequals - prev_nequals ) ; i++ ) {
      fprintf( stdout , "=" ) ; fflush( stdout ) ;
    }
    prev_nequals = nequals ;
  }
  return ;
}
