/**
   @file structs.h
   @brief storage for all the structs
 */
#ifndef STRUCTS_H
#define STRUCTS_H

#define ALIGNMENT (16)

/**
   @struct baryon_info
   @brief baryon contraction info
   @param map :: contraction map indices
   @param outfile :: output file name
 */
struct baryon_info {
  size_t map[ 3 ] ;
  char outfile[ 256 ] ;
} ;

/**
   @struct block
   @brief wrapper for flattened spinmatrix
 */
struct block {
  double complex M[ NSNS ] __attribute__((aligned(ALIGNMENT))) ;
} ;

/**
   @struct colormatrix
   @brief color matrix
*/
struct colormatrix {
  double complex C[ NC ][ NC ] __attribute__((aligned(ALIGNMENT))) ;
} ;

/**
   @struct halfspinor
   @brief heavy prop
 */
struct halfspinor {
  double complex D[ ND ][ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
} ;

/**
   @struct halfspinorf
   @brief heavy prop
 */
struct halfspinor_f {
  float complex D[ ND ][ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
} ;

/**
   @struct correlator
   @brief correlator data storage
 */
struct correlator{
  double complex *C ;
} ;

/**
   @param type :: cut type we are doing
   @param max_mom :: maximum n^2
   @param cyl_width :: cylinder width
   @param configspace :: configuration space hook
 */
struct cut_info{
  momentum_cut_def type ; // enumerated cutting type
  size_t max_mom ; // maximum momentum allowed for the cut
  double cyl_width ; // cylinder with
  GLU_bool configspace ; // are we doing configuration space analysis
  size_t max_r2 ; // maximum r^2 we will sum to
  double proto_mom[ ND-1 ] ;
  double thetas[ 16 ] ;
  size_t Nalphas ;
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
   @struct input_info
   @brief input data struct 
 */
struct input_info {
  size_t nprops ;
  struct baryon_info *baryons ;
  size_t nbaryons ;
  struct meson_info *diquarks ;
  size_t ndiquarks ;
  struct meson_info *mesons ;
  size_t nmesons ;
  struct penta_info *pentas ;
  size_t npentas ;
  struct tetra_info *tetras ;
  size_t ntetras ;
  struct VPF_info *VPF ;
  size_t nVPF ;
  struct WME_info *wme ;
  size_t nWME ;
  struct cut_info CUTINFO ;
  size_t dims[ ND ] ;
} ;

/**
   @struct inputs
   @brief input file information is packed in this struct
 */
struct inputs {
  char TOKEN[ GLU_STR_LENGTH ] ;
  char VALUE[ GLU_STR_LENGTH ] ;
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
  size_t dims[ ND ] ; // dimensions in x,y,z,t order opposite to FFTW
  size_t Lsq ; // dims[0] * dims[1] x,y plane
  size_t Lcu ; // dims[2] * Lsq x,y,z cubic volume
  size_t Volume ; // lattice Volume
  size_t flow ; // config number , gets passed around a bit
  header_mode head ;// Which header type are we using
  double twiddles[ ND ] ; // fourier transform twiddles
  uint32_t Nthreads ; // number of threads
  uint32_t Seed ;
} ;

// structure of measurements
struct measurements {
  struct spinor **S ;
  struct spinor **Sf ;
  struct spinor *SUM ;
  struct gamma *GAMMAS ;
  struct veclist *list ;
  struct veclist *wwlist ;
  struct veclist_int *rlist ;
  int NR ;
  int *nmom ;
  int *wwnmom ;
  double complex **in ;
  double complex **out ;
#ifdef HAVE_FFTW3_H
  fftw_plan *forward , *backward ;
#else
  int *forward , *backward ;
#endif
  struct mcorr **corr ;
  struct mcorr **wwcorr ;
  GLU_bool is_wall ;
  GLU_bool configspace ;
  double sum_mom[ ND ] ;
  double sum_twist[ ND ] ;
  GLU_bool is_wall_mom ;
  double complex **dft_mom ;
  GLU_bool is_dft ;
  size_t Nprops ;
} ;

/**
   @struct meson_info
   @brief meson contraction information
   @param map :: maps contractions indices
   @param outfile :: output file name
 */
struct meson_info {
  size_t map[2] ;
  char outfile[ 256 ] ;
} ;

/**
   @struct mcorr
   @brief storage for ( p_{ND-1} , t ) correlation functions
 */
struct mcorr {
  struct correlator *mom ;
} ;

/**
   @struct NRQCD_params
   @brief various NRQCD params storage
 */
struct NRQCD_params {
  double U0 ;  // tadpole
  // correction terms
  double C0 , C1 , C2 , C3 , C4 , C5 , C6 , C7 , C8 , C9EB , C10EB , C11 ;
  double M_0 ; // bare heavy quark mass
  size_t N ;   // number of hamiltonian applications
  GLU_bool backward ; // direction of propagator
} ;

/**
   @struct PIdata
   @brief VPF storage
 */
struct PIdata {
  double complex PI[ ND ][ ND ] ;
} ;

/**
   @struct propagator
   @brief container for the propagator
 */
struct propagator {
  FILE *file ;
  proptype basis ;
  sourcetype source ;
  size_t origin[ ND ] ;
  boundaries bound[ ND ] ;
  double twist[ ND ] ;
  double mom_source[ ND ] ;
  struct halfspinor_f *H ;
  struct NRQCD_params NRQCD ;
  fp_precision precision ;
  endianness endian ;
  smearing smear ;
  size_t Nsmear ;
  size_t t ;
} ;

/**
   @struct QCDheader
   @brief contains the NERSC header information
   This is only used in chklat_stuff.c and should probably be moved there
 */
struct QCDheader {
  int ntoken ; 
  char **token ; 
  char **value ; 
} ;

/**
   @struct site
   @brief the gauge field format
 */
struct site {
  double complex O[ ND ][ NCNC ] __attribute__((aligned(ALIGNMENT))) ;
  int neighbor[ ND ] ;
  int back[ ND ] ;
} ;

/**
   @struct spinmatrix
   @brief dirac components
*/
struct spinmatrix {
  double complex D[ NS ][ NS ] __attribute__((aligned(ALIGNMENT))) ;
} ;

/**
   @struct spinor
   @brief spinor is a spinmatrix of colormatrices
*/
struct spinor{
  struct colormatrix D[ NS ][ NS ] __attribute__((aligned(ALIGNMENT))) ;
} ;

/**
   @struct Ospinor
   @brief opposite-ordering spinor
 */
struct Ospinor {
  struct spinmatrix C[ NC ][ NC ] __attribute__((aligned(ALIGNMENT))) ;
} ;

/**
   @struct penta_info
   @brief pentaquark contraction info
   @param map :: contraction map indices
   @param outfile :: output file name
 */
struct penta_info {
  size_t map[3] ;
  char outfile[ 256 ] ;
} ;

/**
   @struct tetra_info
   @brief tetraquark contraction info
   @param map :: contraction map indices
   @param outfile :: output file name
 */
struct tetra_info {
  size_t map[4] ;
  char outfile[ 256 ] ;
} ;

/**
   @struct veclist
   @brief storage for the momenta
 */
struct veclist_int {
  size_t idx ;
  int MOM[ ND ] ;
  int nsq ;
} ;

/**
   @struct veclist
   @brief storage for the momenta
 */
struct veclist {
  size_t idx ;
  double MOM[ ND ] ;
  double nsq ;
} ;

/**
   @struct VPF_info
   @brief VPF contractions information
   @param map :: maps contractions indices
   @param current :: of type #current_type
   @param outfile :: output file name
 */
struct VPF_info {
  size_t map[2] ;
  current_type current ;
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
  size_t map[4] ;
  proptype proptype1 ;
  proptype proptype2 ;
  proptype proptype3 ;
  proptype proptype4 ;
  char outfile[ 256 ] ;
} ;

// allocation of a field
struct field {
  double complex **O ;
} ;

// little struct for the NRQCD temporaries
struct NRQCD_fields {
  struct halfspinor *S ;
  struct halfspinor *S1 ;
  struct halfspinor *S2 ;
  struct halfspinor *H ;
  struct field *Fmunu ;
} ;

#endif
