MODULE MyConstants
!///////////////////////////////////////////////////////////////////////////////
! $Id: MyConstants.f90 155 2017-05-23 17:11:03Z jyli $
!-------------------------------------------------------------------------------
! D E S C R I P T I O N :
!    Define commonly used constants.  
! A U TH O R: 
!    Jason Li (SSAI)
!///////////////////////////////////////////////////////////////////////////////


USE DataTypeDef

IMPLICIT NONE

!...............................................................................
! return status codes: more or less aligned with HDF5 status codes
!...............................................................................
INTEGER, PARAMETER :: FAILURE_STATE = -1
INTEGER, PARAMETER :: SUCCESS_STATE =  0
INTEGER, PARAMETER :: WARNING_STATE =  1

!...............................................................................
! size constants: borrowed from sdptk and CONFIG_READER
!...............................................................................
! see PGS_PC.f of sdptk toolkit
!  PGSd_PC_FILE_NAME_MAX:
INTEGER(I4B), PUBLIC, PARAMETER :: MAX_FILE_NAME_LENGTH = 256 

! PGSd_PC_PATH_LENGTH_MAX:
INTEGER(I4B), PUBLIC, PARAMETER :: MAX_PATH_NAME_LENGTH = 768

! PGSd_PC_FILE_PATH_MAX = PGSd_PC_PATH_LENGTH_MAX + PGSd_PC_FILE_NAME_MAX 
INTEGER(I4B), PUBLIC, PARAMETER :: MAX_FILE_PATH_NAME_LENGTH = 1024


! see GetConfig.inc of the CONFIG_READER:
! CFG_KEY_LEN and CFG_VAL_LEN: as defined in CONFIG_READER GetConfig.inc
INTEGER(I4B), PUBLIC, PARAMETER :: MAX_STRING_LENGTH = 260


!...............................................................................
! unsigned INTEGERs (as defined in OMI):
!...............................................................................
INTEGER(I4B), PUBLIC, PARAMETER :: FILLVALUE_UI1B = 255    ! = 2^8  - 1
INTEGER(I4B), PUBLIC, PARAMETER :: FILLVALUE_UI2B = 65535  ! = 2^16 - 1

!...............................................................................
! signed INTEGERs:
!...............................................................................
INTEGER(I1B), PUBLIC, PARAMETER :: FILLVALUE_I1B = - HUGE(1_I1B)
INTEGER(I2B), PUBLIC, PARAMETER :: FILLVALUE_I2B = - HUGE(1_I2B)
INTEGER(I4B), PUBLIC, PARAMETER :: FILLVALUE_I4B = - HUGE(1_I4B) 

!...............................................................................
! float type (as defined in OMI):
!...............................................................................
REAL(SP), PUBLIC, PARAMETER :: FILLVALUE_SP = -2.0**100    
REAL(DP), PUBLIC, PARAMETER :: FILLVALUE_DP = -2.0d0**100  

!...............................................................................
! Math constants:
!...............................................................................
REAL(SP), PARAMETER :: PI      = 3.141592653589793238462643383279502884197_sp
REAL(SP), PARAMETER :: PIO2    = 1.57079632679489661923132169163975144209858_sp
REAL(SP), PARAMETER :: TWOPI   = 6.283185307179586476925286766559005768394_sp
REAL(SP), PARAMETER :: SQRT2   = 1.41421356237309504880168872420969807856967_sp
REAL(SP), PARAMETER :: EULER   = 0.5772156649015328606065120900824024310422_sp
REAL(DP), PARAMETER :: PI_D    = 3.141592653589793238462643383279502884197_dp
REAL(DP), PARAMETER :: PIO2_D  = 1.57079632679489661923132169163975144209858_dp
REAL(DP), PARAMETER :: TWOPI_D = 6.283185307179586476925286766559005768394_dp

END MODULE MyConstants

