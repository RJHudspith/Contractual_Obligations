/**
   @file structs.h
   @brief storage for all the structs
 */

#ifndef STRUCTS_H
#define STRUCTS_H

/**
   @struct colormatrix
   @brief color matrix
*/
struct colormatrix {
  double complex C[ NC ][ NC ] ;
} ;

/**
   @struct spinmatrix
   @brief dirac components
*/
struct spinmatrix {
  double complex D[ NS ][ NS ] ;
} ;

/**
   @struct spinor
   @brief spinor is a spinmatrix of colormatrices
*/
struct spinor{
  struct colormatrix D[ NS ][ NS ] ;
} ;

/**
   @struct correlator
   @brief correlator data storage
 */
struct correlator{
  double complex *C ;
} ;

/**
   @struct lilgamma
   @brief tiny structure for the gammas
 */
struct lilgamma {
  unsigned int g : 2 ;
  unsigned int ig : 2 ;
} ;

/**
   @struct gamma
   @brief gamma matrix type
   uint8_t so it behaves better in cache
*/
struct gamma{
  uint8_t g[ NS ] ;
  uint8_t ig[ NS ] ;
} ;

/**
   @struct PIdata
   @brief VPF storage
 */
struct PIdata {
  double complex PI[ ND ][ ND ] ;
} ;

/**
   @struct site
   @brief the gauge field format
 */
struct site
{
  double complex O[ ND ][ NCNC ] ;
  int neighbor[ ND ] ;
  int back[ ND ] ;
} ;

/**
   @enum config_size
   @brief for writing out files
   loops number of spaces for this format uses NC and NCNC
 */
enum config_size
  { LOOP_SMALL = NCNC - 1 ,
    LOOP_GAUGE = 2 * ( NC - 1 ) * NC ,
    LOOP_NCxNC = 2 * NCNC } ;

/**
   @struct QCDheader
   @brief contains the NERSC header information
   This is only used in chklat_stuff.c and should probably be moved there
 */
struct QCDheader{
  int ntoken ; 
  char **token ; 
  char **value ; 
};

/**
   @struct head_data
   @brief information taken from the header
   @param endianess :: tells us which end is up
   @param config_type :: NERSC, HIREP, SCIDAC, LIME ... etc
   @param precision :: what precision the configs were saved in
   @param plaquette :: the read value of the average plaquette
   @param trace :: the read value of the average link trace
   @param checksum :: the read value of the checksum
 */
struct head_data{
  GLU_endian endianess ;
  GLU_output config_type ;
  file_prec precision ;
  double plaquette ;
  double trace ;
  uint32_t checksum ;
  uint32_t checksumb ;
} ;

/**
   @struct latt_info
   @brief (useful?) lattice information
   @param dims[mu] :: lattice dimensions in c-order, x moves quickest
   @param Volume :: lattice volume
   @param Lcu :: volume with one (slowest moving) dimension removed
   @param Lsq :: volume with two (slowest moving) dimensions removed
   @param flow :: the configuration number
   @param gf_alpha :: the tuning parameter of our gauge fixer. 1.0/12.0 is good.
   @param sm_alpha :: the smearing parameters from the input file
   @param head :: what header type we use
   @param Seed :: the seed we use for our RNG 
 */
struct latt_info{
  int dims[ ND ] ; // dimensions in x,y,z,t order opposite to FFTW
  int Lsq ; // dims[0] * dims[1] x,y plane
  int Lcu ; // dims[2] * Lsq x,y,z cubic volume
  int Volume ; // lattice Volume
  int flow ; // config number , gets passed around a bit
  header_mode head ;// Which header type are we using
  double twiddles[ ND ] ; // fourier transform twiddles
} ;

/**
   @struct meson_info
   @brief meson contraction information
   @param map :: maps contractions indices
   @param source :: of type #sourcetype
   @param proptype1 :: type of prop1
   @param proptype2 :: type of prop2
   @param outfile :: output file name
 */
struct meson_info {
  int map[2] ;
  sourcetype source ;
  proptype proptype1 ;
  proptype proptype2 ;
  char outfile[ 256 ] ;
} ;

#endif
