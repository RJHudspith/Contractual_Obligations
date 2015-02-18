/**
   @file crc32.h
   @brief prototype function for computing the checksum
 */
#ifndef CRC32_H
#define CRC32_H

/**
   @fn void DML_checksum_accum( uint32_t *checksuma , uint32_t *checksumb , const uint32_t rank, unsigned char *buf, size_t size )
   @brief compute the CRC32 checksum
 */
void 
DML_checksum_accum( uint32_t *checksuma , 
		    uint32_t *checksumb , 
		    const uint32_t rank, 
		    unsigned char *buf, 
		    size_t size ) ;

/**
   @fn void CKSUM_ADD( void *memptr , const uint32_t nbytes )
   @brief BQCD's crc accumulator
   @param memptr :: start of the memory that we crc
   @param nbytes :: total number of bytes allocated in memptr
 */
void 
CKSUM_ADD( void *memptr , 
	   const uint32_t nbytes ) ;

/**
   @fn void CKSUM_GET( uint32_t *total_crc, uint32_t *total_bytes ) ;
   @brief BQCD's crc
 */
void 
CKSUM_GET( uint32_t *total_crc, 
	   uint32_t *total_bytes ) ;

#endif
