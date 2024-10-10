MODULE OCI_UV1km_DataModule

USE H5Util_class
USE HDF5
USE DataTypeDef
USE MyConstants

 IMPLICIT NONE

!==============================================================================
! 1km UV data
!==============================================================================
  REAL(KIND=8),        DIMENSION(:), ALLOCATABLE :: Time
!  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: Latitude, Longitude
!  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: SolarZenithAngle
!  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: SolarAzimuthAngle
!  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: ViewingZenithAngle
!  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: ViewingAzimuthAngle
!  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: RelativeAzimuthAngle
  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: TerrainPressure
  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: SnowIce_fraction
!  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: TerrainHeight
  REAL(KIND=4),    DIMENSION(:,:,:), ALLOCATABLE :: SurfaceAlbedo_Oceancorr
!  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: UVAerosolIndex
  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: CloudOpticalDepth
  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: CloudFraction
!  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: Residue
!  REAL(KIND=4),    DIMENSION(:,:,:), ALLOCATABLE :: Reflectivity
  REAL(KIND=4),    DIMENSION(:,:,:), ALLOCATABLE :: NormRadiances
  INTEGER(KIND=2),   DIMENSION(:,:), ALLOCATABLE :: SurfaceTypeNativeRes

! Functions
  PUBLIC  :: allocate_UV1km_Data
  PUBLIC  :: deallocate_UV1km_Data

  CONTAINS   
!==============================================================================

  FUNCTION allocate_UV1km_Data(nX, nT) RESULT(STATUS)
!==============================================================================
! TITLE: Allocate memory for output data
! NAME: allocate_UV1km_Data
! INPUTS:
!     nT        Number of scans (XDim)
!     nX        Number of Cross-track pixel
!==============================================================================
  IMPLICIT NONE

  INTEGER(KIND=4), INTENT(IN) :: nT, nX
  INTEGER(KIND=4)        :: STATUS

!==============================================================================
! Allocate space for Time
!==============================================================================
  ALLOCATE(Time(nT), stat=STATUS)                                           
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate Time"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Allocate space for NormRadiances
!==============================================================================
  ALLOCATE(NormRadiances(2,nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate NormRadiances"
     CALL EXIT(1)
  ENDIF 
                                
!==============================================================================
! Allocate space for SurfaceTypeNativeRes
!==============================================================================
  ALLOCATE(SurfaceTypeNativeRes(nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate SurfaceTypeNativeRes"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Allocate space for TerrainPressure
!==============================================================================
  ALLOCATE(TerrainPressure(nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate TerrainPressure"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Allocate space for SnowIce_fraction
!==============================================================================
  ALLOCATE(SnowIce_fraction(nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate SnowIce_Fraction"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Allocate space for SurfaceAlbedo_Oceancorr
!==============================================================================
  ALLOCATE(SurfaceAlbedo_Oceancorr(2,nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate SurfaceAlbedo_Oceancorr"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Allocate space for CloudOpticalDepth
!==============================================================================
  ALLOCATE(CloudOpticalDepth(nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate CloudOpticalDepth"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Allocate space for CloudFraction
!==============================================================================
  ALLOCATE(CloudFraction(nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate CloudFraction"
     CALL EXIT(1)
  ENDIF 


  RETURN

 END FUNCTION allocate_UV1km_Data
!!
!!
!!
 FUNCTION deallocate_UV1km_Data() RESULT(STATUS)
!==============================================================================
! TITLE: Deallocate memory for Input data
! NAME: deallocateInput
!==============================================================================
  IMPLICIT NONE

  INTEGER(KIND=4)        :: STATUS

!==============================================================================
! Deallocate space for CloudOpticalDepth
!==============================================================================
  IF(ALLOCATED(CloudOpticalDepth))DEALLOCATE(CloudOpticalDepth,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate CloudOpticalDepth"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Deallocate space for CloudFraction
!==============================================================================
  IF(ALLOCATED(CloudFraction))DEALLOCATE(CloudFraction,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate CloudFraction"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Deallocate space for SurfaceAlbedo_Oceancorr
!==============================================================================
  IF(ALLOCATED(SurfaceAlbedo_Oceancorr))DEALLOCATE(SurfaceAlbedo_Oceancorr,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate SurfaceAlbedo_Oceancorr"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Deallocate space for SnowIce_fraction
!==============================================================================
  IF(ALLOCATED(SnowIce_fraction))DEALLOCATE(SnowIce_fraction,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate SnowIce_fraction"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Deallocate space for TerrainPressure
!==============================================================================
  IF(ALLOCATED(TerrainPressure))DEALLOCATE(TerrainPressure,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate TerrainPressure"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Deallocate space for SurfaceTypeNativeRes
!==============================================================================
  IF(ALLOCATED(SurfaceTypeNativeRes))DEALLOCATE(SurfaceTypeNativeRes,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate SurfaceTypeNativeRes"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Deallocate space for NormRadiances
!==============================================================================
  IF(ALLOCATED(NormRadiances))DEALLOCATE(NormRadiances,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate NormRadiances"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Deallocate space for Time
!==============================================================================
  IF(ALLOCATED(Time))DEALLOCATE(Time,stat=STATUS)                           
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate Time "
     CALL EXIT(1)
  ENDIF 

  RETURN

 END FUNCTION Deallocate_UV1km_Data
!!
!!
!==============================================================================
!==============================================================================
!
!!
END MODULE OCI_UV1km_DataModule
