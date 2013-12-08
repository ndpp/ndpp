module constants
  implicit none

  ! ============================================================================
  ! VERSIONING NUMBERS

  ! OpenMC major, minor, and release numbers
  integer, parameter :: VERSION_MAJOR   = 1
  integer, parameter :: VERSION_MINOR   = 0
  integer, parameter :: VERSION_RELEASE = 0

  ! Revision numbers for binary files
  !!! Will need something like this for my output files.
!   integer, parameter :: REVISION_SOURCE     = 0

  ! ============================================================================
  ! ADJUSTABLE PARAMETERS

  ! NOTE: This is the only section of the constants module that should ever be
  ! adjusted. Modifying constants in other sections may cause the code to fail.

  ! Monoatomic ideal-gas scattering treatment threshold
  real(8), parameter :: FREEGAS_THRESHOLD_DEFAULT = 400.0

  ! Significance level for confidence intervals
  real(8), parameter :: CONFIDENCE_LEVEL = 0.95_8

  ! Used for surface current tallies
  real(8), parameter :: TINY_BIT = 1e-8_8

  ! User for precision in geometry
  real(8), parameter :: FP_PRECISION = 1e-14_8
  real(8), parameter :: FP_REL_PRECISION = 1e-5_8
  real(8), parameter :: FP_COINCIDENT = 1e-12_8

  ! Maximum number of words in a single line, length of line, and length of
  ! single word
  integer, parameter :: MAX_WORDS    = 500
  integer, parameter :: MAX_LINE_LEN = 250
  integer, parameter :: MAX_WORD_LEN = 150
  integer, parameter :: MAX_FILE_LEN = 255

  ! ============================================================================
  ! PHYSICAL CONSTANTS

  real(8), parameter ::            &
       PI             = 3.1415926535898_8, & ! pi
       MASS_NEUTRON   = 1.0086649156,      & ! mass of a neutron
       MASS_PROTON    = 1.00727646677,     & ! mass of a proton
       AMU            = 1.66053873e-27,    & ! 1 amu in kg
       JOULES_PER_MEV = 1.60217657E-13_8,  & ! Joules from MeV
       N_AVOGADRO     = 0.602214179,       & ! Avogadro's number in 10^24/mol
       K_BOLTZMANN    = 8.617342e-11,      & ! Boltzmann constant in MeV/K
       INFINITY       = huge(0.0_8),       & ! positive infinity
       ZERO           = 0.0_8,             &
       ONE            = 1.0_8,             &
       TWO            = 2.0_8

  ! ============================================================================
  ! NDPP OBJECT CONSTANTS
  ! Type of output to produce for the scattering portion of NDPP
  integer, parameter ::           &
       SCATT_TYPE_LEGENDRE  =  0, & ! Produce Legendre moments of the distro
       SCATT_TYPE_TABULAR   =  1    ! Produce Tabular Representation of distro

  ! The Index in scattdata % mu_bounds(:,g) which contains the low and high
  ! boundaries of integration of the angular variable.
  integer, parameter :: &
       MU_LO =  1,      &
       MU_HI =  2

  ! HDF5 output type (ascii and binary defined below)
  integer, parameter :: &
       H5     = 3,      &
       NO_OUT = 4,      &
       HUMAN  = 5

  ! Value for freegas_cutoff in cross_sections.xml which means 'use the global'
  real(8), parameter :: &
       GLOBAL_FREEGAS_CUTOFF = -TWO

  ! Value for freegas_cutoff (everywhere) which means 'apply everywhere'
  real(8), parameter :: &
       INFINITE_FREEGAS_CUTOFF = -ONE

  ! Various Default values for input parameters, ints, reals, bool, chars:
  integer, parameter :: &
       MU_BINS_DEFAULT     = 2001, & ! # of angles for numerical integration
       SCATT_ORDER_DEFAULT = 5,    & ! Legendre order default, if omitted
       THREADS_DEFAULT     = -1,   & ! Default number of threads (-1 means use envvar)
       SCATT_TYPE_DEFAULT  = SCATT_TYPE_LEGENDRE ! Default scatt results type

  real(8), parameter :: &
       PRINT_TOL_DEFAULT   = 1.0E-8_8, & ! Value at which scattering matrices
                                         ! values are not printed
       THIN_TOL_DEFAULT    = 1.0E-2_8    ! Value for grid thinning

  logical, parameter :: &
       INTEGRATE_CHI_DEFAULT = .false., & ! Integrate Chi?
       USE_FREEGAS_DEFAULT   = .true.     ! Use free-gas treatment?


  ! ============================================================================
  ! SCATTDATA OBJECT CONSTANTS
  integer, parameter :: MAX_LEGENDRE_ORDER   = 10

  ! Number of equiprobable bins
  integer, parameter :: NUM_EP = 32
  ! Reciprocal of NUM_EP
  real(8), parameter :: R_NUM_EP = ONE / 32.0_8

  ! Flag to interpolate on angular distributions with nearest neighbor or
  ! linear interpolation
  logical, parameter :: INTERP_NEAREST = .false.
  ! logical, parameter :: INTERP_NEAREST = .true.

  ! Fraction of maximum s(a,b) value to use as cutoff for determining
  ! the range of integration
  real(8), parameter :: SAB_THRESHOLD = 1.0E-6_8

  ! Brent root finding algorithm threshold for the mu variable
  real(8), parameter :: BRENT_MU_THRESH = 1.0E-6_8

  ! Adaptive Simpsons integration (of mu) tolerance
  real(8), parameter :: ADAPTIVE_MU_TOL = 1.0E-7_8

  ! Adaptive Simpsons integration (of mu) maximum recursion depth/iterations
  integer, parameter :: ADAPTIVE_MU_ITS = 15

  ! Adaptive Simpsons integration (of Eout) tolerance
  real(8), parameter :: ADAPTIVE_EOUT_TOL = 1.0E-8_8

  ! Adaptive Simpsons integration (of Eout) maximum recursion depth/iterations
  integer, parameter :: ADAPTIVE_EOUT_ITS = 15

  ! Number of Eout points per bin for thermal scatter collisions
  integer, parameter :: SAB_EPTS_PER_BIN = 3  ! 0 Implies no expansion

  ! ============================================================================
  ! CROSS SECTION RELATED CONSTANTS

  ! Interpolation flag
  integer, parameter ::   &
       HISTOGRAM     = 1, & ! y is constant in x
       LINEAR_LINEAR = 2, & ! y is linear in x
       LINEAR_LOG    = 3, & ! y is linear in ln(x)
       LOG_LINEAR    = 4, & ! ln(y) is linear in x
       LOG_LOG       = 5    ! ln(y) is linear in ln(x)

  ! Particle type
  integer, parameter :: &
       NEUTRON  = 1, &
       PHOTON   = 2, &
       ELECTRON = 3

  ! Angular distribution type
  integer, parameter :: &
       ANGLE_ISOTROPIC = 1, & ! Isotropic angular distribution
       ANGLE_32_EQUI   = 2, & ! 32 equiprobable bins
       ANGLE_TABULAR   = 3    ! Tabular angular distribution

  ! Secondary energy mode for S(a,b) inelastic scattering
  integer, parameter :: &
       SAB_SECONDARY_EQUAL  = 0, & ! Equally-likely outgoing energy bins
       SAB_SECONDARY_SKEWED = 1    ! Skewed outgoing energy bins

  ! Elastic mode for S(a,b) elastic scattering
  integer, parameter :: &
       SAB_ELASTIC_DISCRETE = 3, & ! Sample from discrete cosines
       SAB_ELASTIC_EXACT    = 4    ! Exact treatment for coherent elastic

  ! Reaction types
  integer, parameter :: &
       TOTAL_XS    = 1, &
       ELASTIC     = 2, &
       N_LEVEL     = 4, &
       MISC        = 5, &
       N_2ND       = 11, &
       N_2N        = 16, &
       N_3N        = 17, &
       N_FISSION   = 18, &
       N_F         = 19, &
       N_NF        = 20, &
       N_2NF       = 21, &
       N_NA        = 22, &
       N_N3A       = 23, &
       N_2NA       = 24, &
       N_3NA       = 25, &
       N_NP        = 28, &
       N_N2A       = 29, &
       N_2N2A      = 30, &
       N_ND        = 32, &
       N_NT        = 33, &
       N_N3HE      = 34, &
       N_ND2A      = 35, &
       N_NT2A      = 36, &
       N_4N        = 37, &
       N_3NF       = 38, &
       N_2NP       = 41, &
       N_3NP       = 42, &
       N_N2P       = 44, &
       N_NPA       = 45, &
       N_N1        = 51, &
       N_N40       = 90, &
       N_NC        = 91, &
       N_DISAPPEAR = 101, &
       N_GAMMA     = 102, &
       N_P         = 103, &
       N_D         = 104, &
       N_T         = 105, &
       N_3HE       = 106, &
       N_A         = 107, &
       N_2A        = 108, &
       N_3A        = 109, &
       N_2P        = 111, &
       N_PA        = 112, &
       N_T2A       = 113, &
       N_D2A       = 114, &
       N_PD        = 115, &
       N_PT        = 116, &
       N_DA        = 117

  ! ACE table types
  integer, parameter :: &
       ACE_NEUTRON   = 1, & ! continuous-energy neutron
       ACE_THERMAL   = 2, & ! thermal S(a,b) scattering data
       ACE_DOSIMETRY = 3    ! dosimetry cross sections

  ! Fission neutron emission (nu) type
  integer, parameter ::   &
       NU_NONE       = 0, & ! No nu values (non-fissionable)
       NU_POLYNOMIAL = 1, & ! Nu values given by polynomial
       NU_TABULAR    = 2    ! Nu values given by tabular distribution

  ! Cross section filetypes
  integer, parameter :: &
       ASCII  = 1, & ! ASCII cross section file
       BINARY = 2    ! Binary cross section file

  ! Probability table parameters
  integer, parameter :: &
       URR_CUM_PROB = 1, &
       URR_TOTAL    = 2, &
       URR_ELASTIC  = 3, &
       URR_FISSION  = 4, &
       URR_N_GAMMA  = 5, &
       URR_HEATING  = 6

  ! Maximum number of partial fission reactions
  integer, parameter :: PARTIAL_FISSION_MAX = 4

  ! ============================================================================
  ! MISCELLANEOUS CONSTANTS

  ! indicates that an array index hasn't been set
  integer, parameter :: NONE = 0

  ! Codes for read errors -- better hope these numbers are never used in an
  ! input file!
  integer, parameter :: ERROR_INT  = -huge(0)
  real(8), parameter :: ERROR_REAL = -huge(0.0_8) * 0.917826354_8

  ! Unit numbers
  integer, parameter :: UNIT_NDPP  = 20 ! unit # for writing ndpp_lib.xml file
  integer, parameter :: UNIT_NUC   = 21 ! unit # for writing nuclide library

end module constants
