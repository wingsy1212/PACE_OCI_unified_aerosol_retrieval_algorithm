 MODULE GetLUT_LER_AI_H5module 
!==============================================================================
! DESCRIPTION:
!     This module reads the AI_LER LUT H5 file.
!==============================================================================

 IMPLICIT NONE

 REAL(KIND=4), DIMENSION(:,:,:,:,:), ALLOCATABLE :: rad354_ler, rad388_ler
 REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: rad_lin354_ai_ler, rad_lin388_ai_ler
!
 REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: SURFALBSET_ler, PRESSURESET_ler
 INTEGER(KIND = 4)                       :: nsalb_ler, nplev_ler, nsza_ler

 CONTAINS

!
!==================================================================
!==================================================================
!
 SUBROUTINE ReadLUT_LER_AIparams(lut_fn)

 USE HDF5

 IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
  INTEGER(KIND=4) :: nplev, nsza, nvza, nraa
  INTEGER(KIND=4) :: hdf_err, STATUS 

! HID_T type integers.
  INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id

! HSIZE_T type integer.
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  INTEGER(HSIZE_T), DIMENSION(5) :: dimsN
  CHARACTER(LEN=256) :: group_name, dataset_name, attribute_name
    
! Open the file. Error check.    
  lut_fn = cfg%uv_nc4
  CALL H5Fopen_f(lut_fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!  PRINT *, 'Now reading... ', lut_fn

! Open DATA group.
  group_name = "/ai_ler"
  CALL H5Gopen_f(file_id, group_name, group_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

  nplev_ler = 4
  nsalb_ler = 10
  nsza_ler = 9
  nvza = 14
  nraa = 11

  !! Allocate memory
  ALLOCATE( rad354_ler( nplev_ler,nsalb_ler, nsza_ler,nvza,nraa ), &
            rad388_ler( nplev_ler,nsalb_ler, nsza_ler,nvza,nraa ), &
            SURFALBSET_ler(nsalb_ler), PRESSURESET_ler(nplev_ler),  STAT = STATUS)
  IF (STATUS < 0) THEN 
      PRINT *,'Error : Allocation of variables for HE4 LUT AI failed.'
      CALL EXIT(1)
  ENDIF     
!
  dataset_name = "MolecularAtmosphere354"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dimsN = SHAPE(rad354_ler)
  CALL H5Dread_f(dataset_id, datatype_id, rad354_ler, dimsN, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "MolecularAtmosphere388"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dimsN = SHAPE(rad388_ler)
  CALL H5Dread_f(dataset_id, datatype_id, rad388_ler, dimsN, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "NumberOfSceneReflectivity"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims = SHAPE(surfalbset_ler)
  CALL H5Dread_f(dataset_id, datatype_id, surfalbset_ler, dims, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "NumberOfPressure"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims = SHAPE(pressureset_ler)
  CALL H5Dread_f(dataset_id, datatype_id, pressureset_ler, dims, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  ALLOCATE(rad_lin354_ai_ler(nplev_ler*nsalb_ler*nsza_ler*nvza*nraa),&
           rad_lin388_ai_ler(nplev_ler*nsalb_ler*nsza_ler*nvza*nraa) )
!
  rad_lin354_ai_ler(:) = RESHAPE(rad354_ler, (/nplev_ler*nsalb_ler*nsza_ler*nvza*nraa/) )
  rad_lin388_ai_ler(:) = RESHAPE(rad388_ler, (/nplev_ler*nsalb_ler*nsza_ler*nvza*nraa/) )
!
  CALL H5fclose_f(file_id, hdf_err)

  RETURN

! Very simple error handling (way too simple).
90 WRITE(6,99)
99 FORMAT("Error in subroutine LUT_AI_LER_Reader!")
  STOP

 END SUBROUTINE ReadLUT_LER_AIparams
!
!==================================================================
!==================================================================
!
END MODULE GetLUT_LER_AI_H5module 
