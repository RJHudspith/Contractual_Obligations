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
   @struct cut_info
   @brief cutting information storage
   @param dir :: either spatial or temporal cuts allowed for now
   @param type :: psq,hypercubic,cylinder or conical
   @param max_mom :: maximum allowed p^2 for the vector of ints definition
   @param max_t :: maximum T allowed in static potential
   @param where :: where is our file outputted to?
   @param definition :: are our gauge fields logarithmic or AntiHermitian_projly defined?
   @param angle :: conical angle from the p=0.
   @param cyl_width :: width of the cylinder in lattice units
 */
struct cut_info{
  momentum_cut_def type ; // enumerated cutting type
  int max_mom ; // maximum momentum allowed for the cut
  double cyl_width ; // cylinder with
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

/**
   @struct VPF_info
   @brief VPF contractions information
   @param map :: maps contractions indices
   @param current :: of type #current_type
   @param proptype1 :: type of prop1
   @param proptype2 :: type of prop2
   @param outfile :: output file name
 */
struct VPF_info {
  int map[2] ;
  current_type current ;
  proptype proptype1 ;
  proptype proptype2 ;
  char outfile[ 256 ] ;
} ;

/**
   @struct WME_info
   @brief WME information
   @param map :: maps contractions indices
   @param proptype1 :: type of prop0
   @param proptype2 :: type of prop1
   @param proptype3 :: type of prop2
   @param proptype4 :: type of prop3
   @param outfile :: output file name
 */
struct WME_info {
  int map[4] ;
  proptype proptype1 ;
  proptype proptype2 ;
  proptype proptype3 ;
  proptype proptype4 ;
  char outfile[ 256 ] ;
} ;

/**
   @struct veclist
   @brief storage for the momenta
 */
struct veclist {
  int idx ;
  int MOM[ ND ] ;
} ;

/**
   @struct propfile
   @brief propagator file storage
 */
struct propfile {
  char filename[ GLU_STR_LENGTH ] ;
} ;

/**
   @struct input_info
   @brief input data struct 
 */
struct input_info {
  struct propfile *prop_files ;
  int nprops ;
  struct meson_info *mesons ;
  int nmesons ;
  struct VPF_info *VPF ;
  int nVPF ;
  struct WME_info *wme ;
  int nWME ;
  struct cut_info CUTINFO ;
  int dims[ ND ] ;
} ;

#endif
