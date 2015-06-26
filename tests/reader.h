/**
   @file reader.h
   @brief declarations for file reading
 */
#ifndef TESTS_READER_H
#define TESTS_READER_H

/**
   @fn int find_desired_mom( const int **momentum , const int *moms , const int NMOM )
   @brief finds the user-specified momentum in our momentum list returnin the index
 */
int
find_desired_mom( const int **momentum , 
		  const int *moms , 
		  const int NMOM ) ;

/**
   @fn void write_momlist( const int **momentum , const int NMOM )
   @brief write out the full momentum list
 */
void
write_momlist( const int **momentum ,
	       const int NMOM ) ;

/**
   @fn struct mcorr** process_file( int ***momentum , FILE *infile , uint32_t NGSRC[1] , uint32_t NGSNK[1] , uint32_t NMOM[1] )
   @brief read file, allocate mcorr and pack it
 */
struct mcorr**
process_file( int ***momentum ,
	      FILE *infile ,
	      uint32_t NGSRC[1] ,
	      uint32_t NGSNK[1] ,
	      uint32_t NMOM[1] ) ;

#endif
