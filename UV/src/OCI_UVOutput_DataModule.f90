MODULE OCI_UVOutput_DataModule

USE H5Util_class
USE HDF5
USE DataTypeDef
USE MyConstants

 IMPLICIT NONE

!==============================================================================
! Output UV data
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
  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: UVAerosolIndex
  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: CloudOpticalDepth
  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: CloudFraction
  REAL(KIND=4),      DIMENSION(:,:), ALLOCATABLE :: Residue
  REAL(KIND=4),    DIMENSION(:,:,:), ALLOCATABLE :: Reflectivity
  REAL(KIND=4),    DIMENSION(:,:,:), ALLOCATABLE :: NormRadiances

! Functions
  PUBLIC  :: allocate_UVOutput_Data
  PUBLIC  :: deallocate_UVOutput_Data
  PUBLIC  :: inti_UVOutput_Data

  CONTAINS   
!==============================================================================

  FUNCTION allocate_UVOutput_Data(nX, nT) RESULT(STATUS)
!==============================================================================
! TITLE: Allocate memory for output data
! NAME: allocate_UVOutput_Data
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
! Allocate space for UVAerosolIndex
!==============================================================================
  ALLOCATE(UVAerosolIndex(nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate UVAerosolIndex"
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

!==============================================================================
! Allocate space for Residue
!==============================================================================
  ALLOCATE(Residue(nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate Residue"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Allocate space for Reflectivity
!==============================================================================
  ALLOCATE(Reflectivity(2,nX,nT), stat=STATUS)                               
  IF (STATUS < 0) THEN
     PRINT *,"Error : Unable to allocate Reflectivity "
     CALL EXIT(1)
  ENDIF 

!  STATUS = AER_S_SUCCESS

  RETURN

 END FUNCTION allocate_UVOutput_Data
!!
!!
!!
 FUNCTION deallocate_UVOutput_Data() RESULT(STATUS)
!==============================================================================
! TITLE: Deallocate memory for Input data
! NAME: deallocateInput
!==============================================================================
  IMPLICIT NONE

  INTEGER(KIND=4)        :: STATUS

!==============================================================================
! Deallocate space for Reflectivity
!==============================================================================
  IF(ALLOCATED(Reflectivity))DEALLOCATE(Reflectivity,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate Reflectivity"
     CALL EXIT(1)
  ENDIF 

!==============================================================================
! Deallocate space for Residue
!==============================================================================
  IF(ALLOCATED(Residue))DEALLOCATE(Residue,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate Residue"
     CALL EXIT(1)
  ENDIF 

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
! Deallocate space for UVAerosolIndex
!==============================================================================
  IF(ALLOCATED(UVAerosolIndex))DEALLOCATE(UVAerosolIndex,stat=STATUS)         
  IF (STATUS < 0) THEN
     PRINT *,"Error :  Unable to deallocate UVAerosolIndex"
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

 END FUNCTION Deallocate_UVOutput_Data
!!
!!
!==============================================================================
!==============================================================================
!!
SUBROUTINE inti_UVOutput_Data 

USE MyConstants, ONLY: FILLVALUE_SP

IMPLICIT NONE

INTEGER :: STATUS

STATUS = 1

  IF(ALLOCATED(Reflectivity))               Reflectivity = FILLVALUE_SP         
  IF(ALLOCATED(Residue))                         Residue = FILLVALUE_SP        
  IF(ALLOCATED(CloudOpticalDepth))     CloudOpticalDepth = FILLVALUE_SP         
  IF(ALLOCATED(CloudFraction))             CloudFraction = FILLVALUE_SP         
  IF(ALLOCATED(UVAerosolIndex))           UVAerosolIndex = FILLVALUE_SP        
  IF(ALLOCATED(SurfaceAlbedo_Oceancorr))SurfaceAlbedo_Oceancorr = FILLVALUE_SP         
  IF(ALLOCATED(SnowIce_fraction))       SnowIce_fraction = FILLVALUE_SP         
  IF(ALLOCATED(TerrainPressure))         TerrainPressure = FILLVALUE_SP        
  IF(ALLOCATED(Time))                               Time = FILLVALUE_SP                          

END SUBROUTINE inti_UVOutput_Data
!!
!==============================================================================
!==============================================================================
!!
END MODULE OCI_UVOutput_DataModule
