 MODULE GetLUT_H5module 
!==============================================================================
!
! FILENAME:
!     GetLUT_module.f90
!
! DESCRIPTION:
!     This module read the look-up table in he4
!
! AUTHORS:
!     Changwoo Ahn / Science Systems and Applications, Inc.
!
! HISTORY: April 20, 2012
!==============================================================================

 IMPLICIT NONE

 INTEGER(KIND = 4), PARAMETER :: nwaves = 2 ! 354, 380 nm
 INTEGER(KIND = 4), PARAMETER :: nzhgt = 5 
 INTEGER(KIND = 4), PARAMETER :: nw0 = 21
 INTEGER(KIND = 4), PARAMETER :: naod = 7
 INTEGER(KIND = 4), PARAMETER :: nsza1 = 7
 INTEGER(KIND = 4), PARAMETER :: nraa = 11
 INTEGER(KIND = 4), PARAMETER :: nvza = 14
 
 REAL(KIND=4), DIMENSION(nwaves,nzhgt,nw0,naod,nsza1,nraa,nvza):: radp10, radp6
 REAL(KIND=4), DIMENSION(nwaves,nzhgt,nw0,naod,nsza1,nvza) :: trp10, trp6
 REAL(KIND=4), DIMENSION(nwaves,nzhgt,nw0,naod) :: sbp10, sbp6

 CONTAINS

 SUBROUTINE ReadLUTparams(lut_fn)

 USE HDF5
   
 IMPLICIT NONE
!
  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
  INTEGER(KIND=4) :: STATUS, hdf_err
       
! File, group, dataset and attribute names.
  CHARACTER(LEN=256) :: group_name, dataset_name, attribute_name
 
! HID_T type integers.
  INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id
 
! HSIZE_T type integer.
  INTEGER(HSIZE_T), DIMENSION(4) :: dims4
  INTEGER(HSIZE_T), DIMENSION(6) :: dims6
  INTEGER(HSIZE_T), DIMENSION(7) :: dims7
   
! Open AerosolLUT climatology file.
  CALL H5Fopen_f(lut_fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!  PRINT *, 'Now reading... ', lut_fn
 
! Open geolocation group.
  group_name = "Look-up Tables for OMAERUV"
  CALL H5Gopen_f(file_id, group_name, group_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

  dataset_name = "RadianceP10"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims7 = SHAPE(radp10)
  CALL H5Dread_f(dataset_id, datatype_id, radp10, dims7, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  
! Read RadianceP1000 dataset.
  dataset_name = "RadianceP10"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims7 = SHAPE(radp10)
  CALL H5Dread_f(dataset_id, datatype_id, radp10, dims7, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read RadianceP0600 dataset.
  dataset_name = "RadianceP6"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims7 = SHAPE(radp10)
  CALL H5Dread_f(dataset_id, datatype_id, radp6, dims7, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read RadianceP1000 dataset.
  dataset_name = "TransmittanceP10"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims6 = SHAPE(trp10)
  CALL H5Dread_f(dataset_id, datatype_id, trp10, dims6, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read TransmittanceP0600 dataset.
  dataset_name = "TransmittanceP6"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims6 = SHAPE(trp6)
  CALL H5Dread_f(dataset_id, datatype_id, trp6, dims6, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read SBARP1000 dataset.
  dataset_name = "SphericalAlbedoP10"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims4 = SHAPE(sbp10)
  CALL H5Dread_f(dataset_id, datatype_id, sbp10, dims4, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read SBARP0600 dataset.
  dataset_name = "SphericalAlbedoP6"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims4 = SHAPE(sbp6)
  CALL H5Dread_f(dataset_id, datatype_id, sbp6, dims4, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
 
!
  CALL H5fclose_f(file_id, hdf_err)
! 
  RETURN
 
! Very simple error handling (way too simple).
90 WRITE(6,99)
99 FORMAT("Error in subroutine Aerosol LUT Reader!")
  STOP
  

 END SUBROUTINE ReadLUTparams

END MODULE GetLUT_H5module 
