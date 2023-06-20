 MODULE GetLUT_Miecloud_AI_H5module 
!==============================================================================
! DESCRIPTION:
!     This module reads AI_Miecloud LUT H5 file. 
!==============================================================================

 IMPLICIT NONE
 
 REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: rad_lin354_ai, rad_lin388_ai
 REAL(KIND=4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: rad354, rad388
 REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: SURFALBSET, CODSET
 INTEGER(KIND = 4)                       :: nplev_mie, nsalb_mie, ncod_mie, nsza_mie, nvza_mie, nraa_mie
!
 REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: MieRad_lin2p3_ai, MieRad_lin388_ai
 REAL(KIND=4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: MieRad2p3, MieRad388
 REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: plev2p3_nodes, salb2p3_nodes, cod2p3_nodes
 REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: sza2p3_nodes, vza2p3_nodes, raa2p3_nodes
 INTEGER(KIND = 4)                       :: nplev_swir, nsalb_swir, ncod_swir, nsza_swir, nvza_swir, nraa_swir
!

 CONTAINS
!
!==================================================================
!==================================================================
!
 SUBROUTINE ReadLUTAIparams(lut_fn) 

 USE HDF5

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
  INTEGER(KIND=4) :: nsza, hdf_err, STATUS
!
! HID_T type integers.
  INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id

! HSIZE_T type integer.
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  INTEGER(HSIZE_T), DIMENSION(6) :: dimsN
  CHARACTER(LEN=256) :: group_name, dataset_name, attribute_name
    
! Open the file. Error check.    
  CALL H5Fopen_f(lut_fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!  PRINT *, 'Now reading... ', lut_fn

! Open DATA group.
  group_name = "/Look-up Tables for TROPOMAER_AI"
  CALL H5Gopen_f(file_id, group_name, group_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

  nplev_mie = 3
  nsalb_mie = 5
  ncod_mie = 9
  nsza_mie = 9
  nvza_mie = 14
  nraa_mie = 11

  !! Allocate memory
  ALLOCATE( rad354( nplev_mie,nsalb_mie, ncod_mie, nsza_mie,nvza_mie,nraa_mie ), &
            rad388( nplev_mie,nsalb_mie, ncod_mie, nsza_mie,nvza_mie,nraa_mie ), &
            SURFALBSET( nsalb_mie ), CODSET( ncod_mie),  STAT = STATUS)
  IF (STATUS < 0) THEN 
      PRINT *,'Error : Allocation of variables for AI MieCloud LUT failed.'
      CALL EXIT(1)
  ENDIF     
!
  dataset_name = "MieCloud354"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dimsN = SHAPE(rad354)
  CALL H5Dread_f(dataset_id, datatype_id, rad354, dimsN, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "MieCloud388"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dimsN = SHAPE(rad388)
  CALL H5Dread_f(dataset_id, datatype_id, rad388, dimsN, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "NumberOfSurfaceReflectivity"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims = SHAPE(surfalbset)
  CALL H5Dread_f(dataset_id, datatype_id, surfalbset, dims, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "NumberOfCloudOpticalDepth"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims = SHAPE(codset)
  CALL H5Dread_f(dataset_id, datatype_id, codset, dims, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  ALLOCATE(rad_lin354_ai(nplev_mie*nsalb_mie*ncod_mie*nsza_mie*nvza_mie*nraa_mie),&
           rad_lin388_ai(nplev_mie*nsalb_mie*ncod_mie*nsza_mie*nvza_mie*nraa_mie) )
!
  rad_lin354_ai(:) = RESHAPE(rad354, (/nplev_mie*nsalb_mie*ncod_mie*nsza_mie*nvza_mie*nraa_mie/) )
  rad_lin388_ai(:) = RESHAPE(rad388, (/nplev_mie*nsalb_mie*ncod_mie*nsza_mie*nvza_mie*nraa_mie/) )
!
  CALL H5fclose_f(file_id, hdf_err)

  RETURN

! Very simple error handling (way too simple).
90 WRITE(6,99)
99 FORMAT("Error in subroutine LUT_AI_LER_Reader!")
  STOP

 END SUBROUTINE ReadLUTAIparams



!
!==================================================================
!==================================================================
!
 SUBROUTINE Read_SWIRLUTAIparams(lut_fn) 

 USE HDF5

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
  INTEGER(KIND=4) :: nsza, hdf_err, STATUS
!
! HID_T type integers.
  INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id

! HSIZE_T type integer.
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  INTEGER(HSIZE_T), DIMENSION(6) :: dimsN
  CHARACTER(LEN=256) :: group_name, dataset_name, attribute_name
    
! Open the file. Error check.    
  CALL H5Fopen_f(lut_fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!  PRINT *, 'Now reading... ', lut_fn

! Open DATA group.
  group_name = "/Look-up Tables for UV-SWIR AI"
  CALL H5Gopen_f(file_id, group_name, group_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

  nplev_swir = 2
  nsalb_swir = 5
  ncod_swir = 9
  nsza_swir = 9
  nvza_swir = 14
  nraa_swir = 11

  !! Allocate memory
  ALLOCATE( MieRad2p3( nplev_swir,nsalb_swir, ncod_swir, nsza_swir,nvza_swir,nraa_swir ), &
            MieRad388( nplev_swir,nsalb_swir, ncod_swir, nsza_swir,nvza_swir,nraa_swir ), &
            plev2p3_nodes(nplev_swir), salb2p3_nodes(nsalb_swir), cod2p3_nodes(ncod_swir),  &
	    sza2p3_nodes(nsza_swir), vza2p3_nodes(nvza_swir), raa2p3_nodes(nraa_swir), STAT = STATUS)
	    
  IF (STATUS < 0) THEN 
      PRINT *,'Error : Allocation of variables for AI swirCloud LUT failed.'
      CALL EXIT(1)
  ENDIF     
!
  plev2p3_nodes(1:nplev_swir) = [1013.25, 600.0]
  salb2p3_nodes(1:nsalb_swir) = [0.00, 0.10, 0.20, 0.30, 0.40]
   cod2p3_nodes(1: ncod_swir) = [0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 100.0]
   sza2p3_nodes(1: nsza_swir) = [0.0, 20.0, 40.0, 60.0, 66.0, 72.0, 80.0, 84.0, 88.0]
   vza2p3_nodes(1: nvza_swir) = [0, 12, 18, 26, 32, 36, 40, 46, 50, 54, 56, 60, 66, 72]
   raa2p3_nodes(1: nraa_swir) = [0, 30, 60, 90, 120, 150, 160, 165, 170, 175, 180]
!
  dataset_name = "MieCloud2p3"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dimsN = SHAPE(MieRad2p3)
  CALL H5Dread_f(dataset_id, datatype_id, MieRad2p3, dimsN, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "MieCloud388"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dimsN = SHAPE(MieRad388)
  CALL H5Dread_f(dataset_id, datatype_id, MieRad388, dimsN, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  ALLOCATE(MieRad_lin2p3_ai(nplev_swir*nsalb_swir*ncod_swir*nsza_swir*nvza_swir*nraa_swir),&
           MieRad_lin388_ai(nplev_swir*nsalb_swir*ncod_swir*nsza_swir*nvza_swir*nraa_swir) )
!
  MieRad_lin2p3_ai(:) = RESHAPE(MieRad2p3, (/nplev_swir*nsalb_swir*ncod_swir*nsza_swir*nvza_swir*nraa_swir/) )
  MieRad_lin388_ai(:) = RESHAPE(MieRad388, (/nplev_swir*nsalb_swir*ncod_swir*nsza_swir*nvza_swir*nraa_swir/) )
!
  CALL H5fclose_f(file_id, hdf_err)

  RETURN

! Very simple error handling (way too simple).
90 WRITE(6,99)
99 FORMAT("Error in subroutine LUT_AI_LER_Reader!")
  STOP

 END SUBROUTINE Read_SWIRLUTAIparams


!
!==================================================================
!==================================================================
!
END MODULE GetLUT_Miecloud_AI_H5module 
