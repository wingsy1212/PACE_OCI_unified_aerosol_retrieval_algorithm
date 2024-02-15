 MODULE Get_OceanLUT_H5module 
!==============================================================================
!
! FILENAME:
!     Get_OceanLUT_module.f90
!
! DESCRIPTION:
!     This module read the look-up table in he4
!
! AUTHORS:
!     Changwoo Ahn / Science Systems and Applications, Inc.
!
! HISTORY: Mar 28, 2014
!==============================================================================

 IMPLICIT NONE

 INTEGER(KIND = 4), PARAMETER :: nwave_ocean = 2 ! 380, 340 nm
 INTEGER(KIND = 4), PARAMETER :: nsza_ocean = 16
 INTEGER(KIND = 4), PARAMETER :: nvza_ocean = 16
 INTEGER(KIND = 4), PARAMETER :: nraa_ocean = 16
 
 REAL(KIND=4), DIMENSION(nwave_ocean, nsza_ocean, nvza_ocean, nraa_ocean) :: oceanler
 REAL(KIND=4), DIMENSION(nsza_ocean) :: sza_oceanset
 REAL(KIND=4), DIMENSION(nvza_ocean) :: vza_oceanset
 REAL(KIND=4), DIMENSION(nraa_ocean) :: raa_oceanset
 
 CONTAINS

!==============================================================================
!==============================================================================

SUBROUTINE Read_OceanLUTparams(lut_fn, oceanler, nwave_ocean, &
                           nsza_ocean, nvza_ocean, nraa_ocean,&
                       sza_oceanset, vza_oceanset, raa_oceanset)
 USE HDF5
 USE OCIUAAER_Config_Module
   
 IMPLICIT NONE

!
  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
  INTEGER(KIND=4) :: STATUS, hdf_err
!
  REAL(KIND=4), DIMENSION(nwave_ocean, nsza_ocean, nvza_ocean, nraa_ocean), &
                                       INTENT(OUT) :: oceanler
  REAL(KIND=4), DIMENSION(nsza_ocean), INTENT(OUT) :: sza_oceanset
  REAL(KIND=4), DIMENSION(nvza_ocean), INTENT(OUT) :: vza_oceanset
  REAL(KIND=4), DIMENSION(nraa_ocean), INTENT(OUT) :: raa_oceanset
  INTEGER(KIND=4):: nwave_ocean,nsza_ocean,nvza_ocean,nraa_ocean
!
! File, group, dataset and attribute names.
  CHARACTER(LEN=256) :: group_name, dataset_name, attribute_name
 
! HID_T type integers.
  INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id
 
! HSIZE_T type integer.
  INTEGER(HSIZE_T), DIMENSION(1) :: dims1
  INTEGER(HSIZE_T), DIMENSION(4) :: dims4

! Open SnowIce climatology file.
  CALL H5Fopen_f(cfg%uv_nc4, H5F_ACC_RDONLY_F, file_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!  PRINT *, 'Now reading... ', lut_fn
 
! Open geolocation group.
  group_name = "/ocncorr"
  CALL H5Gopen_f(file_id, group_name, group_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read Main dataset.
  dataset_name = "OceanLER_ai"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims4 = SHAPE(oceanler)
  CALL H5Dread_f(dataset_id, datatype_id, oceanler, dims4, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  
! Read SZA dataset.
  dataset_name = "NumberOfSolarZenithAngle"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims1 = SHAPE(sza_oceanset)
  CALL H5Dread_f(dataset_id, datatype_id, sza_oceanset, dims1, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
 

! Read VZA dataset.
  dataset_name = "NumberOfViewingZenithAngle"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims1 = SHAPE(vza_oceanset)
  CALL H5Dread_f(dataset_id, datatype_id, vza_oceanset, dims1, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90


! Read RAA dataset.
  dataset_name = "NumberOfRelativeAzimuthAngle"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims1 = SHAPE(raa_oceanset)
  CALL H5Dread_f(dataset_id, datatype_id, raa_oceanset, dims1, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  CALL H5fclose_f(file_id, hdf_err)
! 
  RETURN
 
! Very simple error handling (way too simple).
90 WRITE(6,99)
99 FORMAT("Error in subroutine snowice_Reader!")
  STOP
!
 END SUBROUTINE Read_OceanLUTparams

END MODULE Get_OceanLUT_H5module 
