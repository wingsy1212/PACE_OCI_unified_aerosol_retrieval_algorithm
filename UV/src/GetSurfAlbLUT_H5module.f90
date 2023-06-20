 MODULE GetSurfAlbLUT_H5module 
!==============================================================================
! DESCRIPTION:
!     This module reads the UV surface albedo LUT H5 file.
!==============================================================================

 IMPLICIT NONE

 REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: SRFLER354, SRFLER388 
 REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: GRDLAT, GRDLON

 CONTAINS
!!
!!
!==============================================================================
!==============================================================================
!
 SUBROUTINE ReadSurfAlbLUTparams(lut_fn, SRFLER354, SRFLER388, GRDLAT, GRDLON)

 USE HDF5

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
  INTEGER(KIND=4) :: STATUS, hdf_err
  INTEGER(KIND=4) :: nmonth, nlats, nlons
!
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: SRFLER354, SRFLER388
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE, INTENT(OUT) :: GRDLAT, GRDLON

! HID_T type integers.
  INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id

! HSIZE_T type integer.
  INTEGER(HSIZE_T), DIMENSION(2) :: dims2
  INTEGER(HSIZE_T), DIMENSION(3) :: dims3
  CHARACTER(LEN=256) :: group_name, dataset_name, attribute_name

! Open the file. Error check.    
  CALL H5Fopen_f(lut_fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!  PRINT *, 'Now reading... ', lut_fn

! Open DATA group.
  group_name = "/OMI Surface Albedo Look-up Tables for OMAERUV"
  CALL H5Gopen_f(file_id, group_name, group_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

  nmonth = 12
  nlats = 720
  nlons = 1440

  !! Allocate memory
  ALLOCATE( SRFLER354( nmonth, nlons, nlats ), &
            SRFLER388( nmonth, nlons, nlats ), &
            GRDLAT( nlons, nlats ), &
            GRDLON( nlons, nlats ), STAT = STATUS)
  IF (STATUS < 0) THEN 
      PRINT *,'Error : Allocation of variables for Surface Albedo LUT failed.'
      CALL EXIT(1)
  ENDIF     
!
  dataset_name = "SRFLER354"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims3 = SHAPE(SRFLER354)
  CALL H5Dread_f(dataset_id, datatype_id, SRFLER354, dims3, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "SRFLER388"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims3 = SHAPE(SRFLER388)
  CALL H5Dread_f(dataset_id, datatype_id, SRFLER388, dims3, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "GRDLAT"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims2 = SHAPE(GRDLAT)
  CALL H5Dread_f(dataset_id, datatype_id, GRDLAT, dims2, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  dataset_name = "GRDLON"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims2 = SHAPE(GRDLON)
  CALL H5Dread_f(dataset_id, datatype_id, GRDLON, dims2, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  CALL H5fclose_f(file_id, hdf_err)
!
  RETURN

! Very simple error handling (way too simple).
90 WRITE(6,99)
99 FORMAT("Error in subroutine SurfAlb_LUT_Reader!")
  STOP

 END SUBROUTINE ReadSurfAlbLUTparams
!!
!!
!!
!==============================================================================
!==============================================================================
!
FUNCTION Get_UVSurfaceAlbedo(mon,lat,lon,surfAlb) RESULT(STATUS)
!==============================================================================
! TITLE:   Read the OCI Surface Albedo Lookup Table
! NAME:    GetOCI_SurfaceAlbedo
! INPUTS:  mon, lat, lon
! OUTPUTS: surfAlb(:)
!==============================================================================
  IMPLICIT NONE

  INTEGER(KIND=4),            INTENT(IN)  :: mon
  REAL(KIND=4),               INTENT(IN)  :: lat, lon
  REAL(KIND=4), DIMENSION(:), INTENT(OUT) :: surfAlb
  INTEGER(KIND=4)                         :: tmpmon
  INTEGER(KIND=4) :: ilon, ilat
!
  INTEGER(KIND=4)    :: STATUS

! -- Extract the surfae albedo at 354 and 388 nm from Pawan's surface LER data set ---
   IF (lat .lt. 0.0) then 
       ilat = 360 + int(4.0*lat + 0.000001)
   ENDIF
   IF (lat .ge. 0.0) then 
       ilat = 361 + int(4.0*lat - 0.000001)
   ENDIF
   IF (lon .lt. 0.0) then 
       ilon = 720 + int(4.0*lon + 0.000001)
   ENDIF
   IF (lon .ge. 0.0) then 
       ilon = 721 + int(4.0*lon - 0.000001)
   ENDIF
!
   if (ilat .gt. 720)  ilat = 720
   if (ilat .lt. 1) ilat = 1
   if (ilon .gt. 1440) ilon = 1440
   if (ilon .lt. 1) ilon = 1
!
   surfAlb(1) = SRFLER354(mon, ilon,ilat)
   if ((surfAlb(1).LT.0.0) .OR. (surfAlb(1).GT.1.0)) surfAlb(1) = -9999.
!   
   surfAlb(2) = SRFLER388(mon, ilon,ilat)
   if ((surfAlb(2).LT.0.0) .OR. (surfAlb(2).GT.1.0)) surfAlb(2) = -9999.
!
  STATUS = 1

  RETURN

 END FUNCTION Get_UVSurfaceAlbedo
!!
!==============================================================================
!==============================================================================
!
END MODULE GetSurfAlbLUT_H5module 
