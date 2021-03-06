/**
   @file currents.h
   @brief conserved-local and local-local current prototype functions
 */
#ifndef CURRENTS_H
#define CURRENTS_H

/**
   @fn void contract_conserved_local_site( struct PIdata *DATA_AA , struct PIdata *DATA_VV ,  const struct site *lat , const struct spinor *S1 , const struct spinor *S1UP , const struct spinor *S2 , const struct spinor *S2UP , const struct gamma *GAMMAS , const size_t AGMAP[ ND ] , const size_t VGMAP[ ND ] , const size_t x , const size_t t ) 
   @brief conserved-local Wilson current at a single site x + LCU * t
 */
void
contract_conserved_local_site( struct PIdata *DATA_AA ,
			       struct PIdata *DATA_VV ,
			       const struct site *lat ,
			       const struct spinor *S1 ,
			       const struct spinor *S1UP ,
			       const struct spinor *S2 ,
			       const struct spinor *S2UP ,
			       const struct gamma *GAMMAS ,
			       const size_t AGMAP[ ND ] ,
			       const size_t VGMAP[ ND ] ,
			       const size_t x ,
			       const size_t t ) ;


/**
   @fn void contract_local_local( struct PIdata *DATA_AA , struct PIdata *DATA_VV , const struct spinor *S1 , const struct spinor *S2 , const struct gamma *GAMMAS , const size_t AGMAP[ ND ] , const size_t VGMAP[ ND ] , const size_t x , const size_t t ) 
   @brief local-local vector and axial currents
*/
void
contract_local_local_site( struct PIdata *DATA_AA ,
			   struct PIdata *DATA_VV ,
			   const struct spinor *S1 ,
			   const struct spinor *S2 ,
			   const struct gamma *GAMMAS ,
			   const size_t AGMAP[ ND ] ,
			   const size_t VGMAP[ ND ] ,
			   const size_t x ,
			   const size_t t ) ;

#endif
