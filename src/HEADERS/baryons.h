/**
   @file baryons.h
   @brief baryon contraction code prototype functions
 */

#ifndef BARYONS_H
#define BARYONS_H

/**
   @fn int baryons_diagonal( struct propagator prop , const char *outfile )
   @brief flavour diagonal baryon contractions

   The Baryon interpolating operator is of the form: B = eps_123  q1 *( q2 Cg_23 q3 )
   This gives us access to the particles:

   pos1 pos2 pos3

   Proton  =  u  ( u Cg5 d )
   Neutron =  d  ( d Cg5 u ) -> equivalent to proton in isospin limit
 				
   Lambda  =  u  ( d Cg5 s ) = d ( u Cg5 s )
   Sigma   =  u  ( u Cg5 s )
   Cascade =  s  ( s Cg5 u )

   Omega   =  s  ( s Cg5 s )
   Delta   =  u  ( u Cgi u ) -> decuplett

   The contraction terms are:

   perm1    perm2

   Proton  =  (0,1,2) - (1,0,2)
   Neutron =  (0,1,2) - (1,0,2)

   Lambda  =  (0,1,2)
   Sigma   =  (0,1,2) - (1,0,2)
   Cascade =  (0,1,2) - (1,0,2)

   Omega   =  (0,1,2) - (1,0,2) + (1,2,0) - (0,2,1) + (2,0,1) - (2,1,0)

   Delta =  2 * (0,1,2) - 4 * (1,0,2)

   @return #SUCCESS or #FAILURE
 */
int
baryons_diagonal( struct propagator prop ,
		  const char *outfile ) ;

#endif
