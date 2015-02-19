/**
   @file gammas.c
   @brief reshuffler for the gamma matrices

   int gamma_matrix: 
   	gamma matrices for mesons, returns reshuffled index and sign
 
       0 -> gamma_0 (time-component)
       1 -> gamma_1
       2 -> gamma_2
       3 -> gamma_3
       4 -> unity
       5 -> gamma_5
       6 -> gamma_0 gamma_5
       7 -> gamma_1 gamma_5
       8 -> gamma_2 gamma_5
       9 -> gamma_3 gamma_5
       10-> gamma_0 gamma_1
       11-> gamma_0 gamma_2
       12-> gamma_0 gamma_3
       13-> gamma_1 gamma_2
       14-> gamma_1 gamma_3
       15-> gamma_2 gamma_3    
 */
#include "common.h" // needed for definition of NS

// map is e^{in\Pi/4}
// 1  -> 0
// I  -> 1
// -1 -> 2
// -I -> 3
int 
gamma_matrix( int g[ NS ] , 
	      int ig[ NS ] , 
	      const int mu )
{
  switch (mu) {
  case 0: /* gamma_0 */
    ig[0] = 2; g[0] = 2 ;
    ig[1] = 3; g[1] = 2 ;
    ig[2] = 0; g[2] = 2 ;
    ig[3] = 1; g[3] = 2 ;
    return 0;
  case 1: /* gamma_1 */
    ig[0] = 3; g[0] =  3 ;
    ig[1] = 2; g[1] =  3 ;
    ig[2] = 1; g[2] =  1 ;
    ig[3] = 0; g[3] =  1 ;
    return 0;
  case 2: /* gamma_2 */
    ig[0] = 3; g[0] = 2 ;
    ig[1] = 2; g[1] = 0 ;
    ig[2] = 1; g[2] = 0 ;
    ig[3] = 0; g[3] = 2 ;
    return 0;
  case 3: /* gamma_3 */
    ig[0] = 2; g[0] = 3 ;
    ig[1] = 3; g[1] = 1 ;
    ig[2] = 0; g[2] = 1 ;
    ig[3] = 1; g[3] = 3 ;
    return 0;
  case 4: /* unity */
    ig[0] = 0; g[0] = 0 ;
    ig[1] = 1; g[1] = 0 ;
    ig[2] = 2; g[2] = 0 ;
    ig[3] = 3; g[3] = 0 ;
    return 0;
  case 5: /* gamma_5 */
    ig[0] = 0; g[0] = 0 ;
    ig[1] = 1; g[1] = 0 ;
    ig[2] = 2; g[2] = 2 ;
    ig[3] = 3; g[3] = 2 ;
    return 0;
  case 6: /* gamma_0 gamma_5 */
    ig[0] = 2; g[0] = 2 ;
    ig[1] = 3; g[1] = 2 ;
    ig[2] = 0; g[2] = 1 ;
    ig[3] = 1; g[3] = 1 ;
    return 0;
  case 7: /* gamma_1 gamma_5 */
    ig[0] = 3; g[0] = 3 ;
    ig[1] = 2; g[1] = 3 ;
    ig[2] = 1; g[2] = 3 ;
    ig[3] = 0; g[3] = 3 ;
    return 0;
  case 8: /* gamma_2 gamma_5 */
    ig[0] = 3; g[0] = 2 ;
    ig[1] = 2; g[1] = 0 ;
    ig[2] = 1; g[2] = 2 ;
    ig[3] = 0; g[3] = 0 ;
    return 0;
  case 9: /* gamma_3 gamma_5 */
    ig[0] = 2; g[0] = 3 ;
    ig[1] = 3; g[1] = 1 ;
    ig[2] = 0; g[2] = 3 ;
    ig[3] = 1; g[3] = 1 ;
    return 0;
  case 10: /* gamma_0 gamma_1 */
    ig[0] = 1; g[0] = 3 ;
    ig[1] = 0; g[1] = 3;
    ig[2] = 3; g[2] = 1;
    ig[3] = 2; g[3] = 1;
    return 0;
  case 11: /* gamma_0 gamma_2 */
    ig[0] = 1; g[0] = 2 ;
    ig[1] = 0; g[1] = 0 ;
    ig[2] = 3; g[2] = 0 ;
    ig[3] = 2; g[3] = 2 ;
    return 0;
  case 12: /* gamma_0 gamma_3 */
    ig[0] = 0; g[0] = 3 ;
    ig[1] = 1; g[1] = 1 ;
    ig[2] = 2; g[2] = 1 ;
    ig[3] = 3; g[3] = 3 ;
    return 0;
  case 13: /* gamma_1 gamma_2 */
    ig[0] = 0; g[0] = 1 ;
    ig[1] = 1; g[1] = 3 ;
    ig[2] = 2; g[2] = 1 ;
    ig[3] = 3; g[3] = 3 ;
    return 0;
  case 14: /* gamma_1 gamma_3 */
    ig[0] = 1; g[0] = 2 ;
    ig[1] = 0; g[1] = 0 ;
    ig[2] = 3; g[2] = 2 ;
    ig[3] = 2; g[3] = 0 ;
    return 0;
  case 15: /* gamma_2 gamma_3 */
    ig[0] = 1; g[0] = 1 ;
    ig[1] = 0; g[1] = 1 ;
    ig[2] = 3; g[2] = 1 ;
    ig[3] = 2; g[3] = 1 ;
    return 0;
  default:
    return -1;
  }
  // should never get here
  return 0;
}
