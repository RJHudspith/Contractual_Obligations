/**
   @file crc32.h
   @brief prototype function for computing the checksum
 */
#ifndef CRC32_H
#define CRC32_H

/**
   @fn void CKSUM_GET( uint32_t *total_crc , uint32_t *total_bytes )
   @brief get the accumulated checksum in total_crc
 */
void 
CKSUM_GET( uint32_t *total_crc , 
	   uint32_t *total_bytes ) ;

/**
   @fn void CKSUM_ADD( void *memptr , const uint32_t nbytes )
   @brief add up the checksum starting at memptr for nbytes
   @warning increments static memory address in the file
 */
void 
CKSUM_ADD( void *memptr , 
	   const uint32_t nbytes ) ;

/**
   @fn void DML_checksum_accum( uint32_t *checksuma , uint32_t *checksumb , const uint32_t rank, char *buf, size_t size )
   @brief compute the CRC32 checksum
 */
void 
DML_checksum_accum( uint32_t *checksuma , 
		    uint32_t *checksumb , 
		    const uint32_t rank, 
		    char *buf, 
		    size_t size ) ;

#endif
