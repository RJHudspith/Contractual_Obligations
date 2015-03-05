/**
   @file crc32.h
   @brief prototype function for computing the checksum
 */
#ifndef CRC32_H
#define CRC32_H

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
