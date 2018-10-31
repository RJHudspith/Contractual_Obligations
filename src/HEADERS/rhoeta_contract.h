/**
   @file rhoeta_contract.h
   @brief prototype declarations for rho-eta contractions 
 */
#ifndef RHOETA_CONTRACT_H
#define RHOETA_CONTRACT_H

/**
   @fn double complex rhoeta_contract( const struct spinor S , const struct spinor Sadj , const size_t GSRC , const size_t GSNK , const struct gamma *GAMMAS )
   @brief contract the (connected only!) rho-eta dimeson for 1 flavour
   @warning computes only the connected pieces!
 */
double complex
rhoeta_contract( const struct spinor S ,
		 const struct spinor Sadj ,
		 const size_t GSRC ,
		 const size_t GSNK ,
		 const struct gamma *GAMMAS ) ;

#endif
