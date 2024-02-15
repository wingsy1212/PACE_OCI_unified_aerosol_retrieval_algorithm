MODULE Get_SnowIce_module

 IMPLICIT NONE

 INTEGER(KIND = 4), PARAMETER :: nlon = 360
 INTEGER(KIND = 4), PARAMETER :: nlat = 180
 INTEGER(KIND = 4), PARAMETER :: nmonth = 12

 REAL(KIND=4), DIMENSION(nlon, nlat, nmonth) :: snowice
 REAL(KIND=4), DIMENSION(nlon) :: swice_Longitude         ! Value of longitude for EPIC pixel of interest.
 REAL(KIND=4), DIMENSION(nlat) :: swice_Latitude          ! Value of latitude for EPIC pixel of interest.
!
CONTAINS

!==============================================================================
! Snow_Ice climtology data set  Reader
!==============================================================================
! 
SUBROUTINE snowice_Reader(lut_fn, snowice, swice_Longitude, swice_Latitude)

  USE HDF5
 
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
  INTEGER(KIND=4) :: STATUS, hdf_err
       
  REAL(KIND=4), DIMENSION(nlon, nlat, nmonth), INTENT(OUT) :: snowice
  REAL(KIND=4), DIMENSION(nlon), INTENT(OUT) :: swice_Longitude         
  REAL(KIND=4), DIMENSION(nlat), INTENT(OUT) :: swice_Latitude         
 
! File, group, dataset and attribute names.
  CHARACTER(LEN=256) :: group_name, dataset_name, attribute_name
 
! HID_T type integers.
  INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id
 
! HSIZE_T type integer.
  INTEGER(HSIZE_T), DIMENSION(1) :: dims1
  INTEGER(HSIZE_T), DIMENSION(2) :: dims2
  INTEGER(HSIZE_T), DIMENSION(3) :: dims3
   
! Open SnowIce climatology file.
  lut_fn = cfg%uv_nc4
  CALL H5Fopen_f(lut_fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!  PRINT *, 'Now reading... ', lut_fn
 
! Open geolocation group.
  group_name = "/snowice"
  CALL H5Gopen_f(file_id, group_name, group_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read Snow_Ice dataset.
  dataset_name = "SNOW_ICE"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims3 = SHAPE(snowice)
  CALL H5Dread_f(dataset_id, datatype_id, snowice, dims3, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! --- In SNOW_ICE.h5 file, GRID_LAT and GRID_LON names are reversed.
!     Therefore, one should read latitude for longitude or vice versa !!!!

! Read longitude dataset.
  dataset_name = "GRID_LON"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims1 = SHAPE(swice_Latitude)
  CALL H5Dread_f(dataset_id, datatype_id, swice_Latitude, dims1, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
 
! Read latitude dataset.
  dataset_name = "GRID_LAT"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims1 = SHAPE(swice_Longitude)
  CALL H5Dread_f(dataset_id, datatype_id, swice_Longitude, dims1, hdf_err)
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
 
 END SUBROUTINE snowice_Reader
!
!==================================================================
!==================================================================

  FUNCTION Get_snowice_fraction(mon,lat,lon, swicefrac) RESULT(STATUS)

! TITLE: Get the snow/ice fraction value for a given month
! NAME:  Get_snowice_fraction
! INPUTS: mon, lat, lon
! OUTPUTS: swicefrac

  IMPLICIT NONE

  INTEGER(KIND=4),            INTENT(IN)  :: mon
  REAL(KIND=4),               INTENT(IN)  :: lat, lon
  REAL(KIND=4),               INTENT(OUT) :: swicefrac
  INTEGER(KIND=4)                         :: tmpmon
  INTEGER(KIND=4) :: ilon, ilat
  INTEGER(KIND=4)    :: STATUS

! -- Extract the closest grid indices from Snow_Ice data set ---
   ilon = int(lon + 179.5) + 1
   ilat = int(90. +lat) + 1
!
   IF (ilat .GT. 180)  ilat = 180
   IF (ilon .GT. 360)  ilon = 360
!
   swicefrac = snowice(ilon,ilat, mon)
   if ((swicefrac.LT.0.0) .OR. (swicefrac.GT.100.0)) swicefrac = -9999.0
!PRINT *, mon, lat,lon,ilat,ilon,swicefrac
  STATUS = 1

  RETURN

 END FUNCTION Get_snowice_fraction
!
!==================================================================
!==================================================================
!
END MODULE Get_SnowIce_module
