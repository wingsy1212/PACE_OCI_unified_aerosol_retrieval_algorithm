! 
! These are common variables from both proxy and synthetic 
! data at the original resolution. 
! Vinay Kayetha (SSAI)
!
Integer(Kind=4) :: l1b_nXtrack, l1b_nLines
Character(Len=255) :: l1b_date_time_str
Real(Kind=4), Dimension(:,:), Allocatable :: Latitude
Real(Kind=4), Dimension(:,:), Allocatable :: Longitude
Real(Kind=4), Dimension(:,:), Allocatable :: SolarAzimuthAngle
Real(Kind=4), Dimension(:,:), Allocatable :: SolarZenithAngle
Real(Kind=4), Dimension(:,:), Allocatable :: ViewingAzimuthAngle
Real(Kind=4), Dimension(:,:), Allocatable :: ViewingZenithAngle
Real(Kind=4), Dimension(:,:), Allocatable :: RelativeAzimuthAngle
Real(Kind=4), Dimension(:,:), Allocatable :: TerrainHeight
Real(Kind=4), Dimension(:,:), Allocatable :: UVAI
Real(Kind=4), Dimension(:,:), Allocatable :: Residue_1km
Real(Kind=4), Dimension(:,:,:), Allocatable :: Reflectivity_1km
!
!
! Based on comparison of OCI Ver2 L1B's and Proxy (TROPOMI measurements 
! convoluted with OCI RSR). It is determined that OCI 354 and 388 channels 
! needs an increase in reflectances by 2.5 and 0.9 percent (July 2024).
! This soft-calibration is applied in the L1B reader (OCIUAAER_L1BModule.f90). 
!
! Real(Kind=4), PARAMETER :: softcalib354 = 1.025
! Real(Kind=4), PARAMETER :: softcalib388 = 1.009
! no soft calibration.
  Real(Kind=4), PARAMETER :: softcalib354 = 1.0
  Real(Kind=4), PARAMETER :: softcalib388 = 1.0
