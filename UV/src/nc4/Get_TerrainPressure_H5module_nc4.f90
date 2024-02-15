 MODULE Get_TerrainPressure_H5module 
!==============================================================================
!
! FILENAME:
!     Get_TerrainPressure_H5module.f90
!
! DESCRIPTION:
!     This module read the OCI-derived surface albedo look-up table in he4
!
!==============================================================================

 IMPLICIT NONE

 INTEGER(KIND = 4), PARAMETER :: nlon = 3600
 INTEGER(KIND = 4), PARAMETER :: nlat = 1800

 REAL(KIND=4), DIMENSION(nlon, nlat) :: tpres
 REAL(KIND=4), DIMENSION(nlon) :: terrain_Longitude
 REAL(KIND=4), DIMENSION(nlat) :: terrain_Latitude

 CONTAINS
!
!==================================================================
!==================================================================
!

SUBROUTINE ReadTerrainPressparams(lut_fn, tpres, terrain_Longitude, terrain_Latitude)

 USE HDF5
 USE OCIUAAER_Config_Module

 IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
  INTEGER(KIND=4) :: fid, swid, STATUS
  INTEGER(KIND=4) :: err
!
  REAL(KIND=4), DIMENSION(nlon, nlat), INTENT(OUT) :: tpres
  REAL(KIND=4), DIMENSION(nlon), INTENT(OUT) :: terrain_Longitude
  REAL(KIND=4), DIMENSION(nlat), INTENT(OUT) :: terrain_Latitude

! File, group, dataset and attribute names.
  CHARACTER(LEN=256) :: group_name, dataset_name, attribute_name

! HID_T type integers.
  INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id

! HSIZE_T type integer.
  INTEGER(HSIZE_T), DIMENSION(1) :: dims1
  INTEGER(HSIZE_T), DIMENSION(2) :: dims2

! Regular four-byte integer.
  INTEGER(KIND=4) :: hdf_err

! Open SnowIce climatology file.
  lut_fn = cfg%uv_nc4
  CALL H5Fopen_f(lut_fn, H5F_ACC_RDONLY_F, file_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!  PRINT *, 'Now reading... ', lut_fn

! Open DATA group.
  group_name = "/surfprs"
  CALL H5Gopen_f(file_id, group_name, group_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read TERR_PRES dataset.
  dataset_name = "TERR_PRES"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims2 = SHAPE(tpres)
  CALL H5Dread_f(dataset_id, datatype_id, tpres, dims2, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read longitude dataset.
  dataset_name = "LON"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims1 = SHAPE(terrain_Longitude)
  CALL H5Dread_f(dataset_id, datatype_id, terrain_Longitude, dims1, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90

! Read latitude dataset. 
  dataset_name = "LAT"
  CALL H5Dopen_f(group_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dget_type_f(dataset_id, datatype_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  dims1 = SHAPE(terrain_Latitude)
  CALL H5Dread_f(dataset_id, datatype_id, terrain_Latitude, dims1, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
  CALL H5Dclose_f(dataset_id, hdf_err)
  IF (hdf_err .NE. 0) GO TO 90
!
  CALL H5fclose_f(file_id, hdf_err)
!
  RETURN

! Very simple error handling (way too simple).
90 WRITE(6,99)
99 FORMAT("Error in subroutine TerrainPressure_Reader!")
  STOP

 END SUBROUTINE ReadTerrainPressparams
!
!==================================================================
!==================================================================
!

  FUNCTION Get_TerrainPressure(lat,lon, thispres) RESULT(STATUS)
!==============================================================================
!
! TITLE:
!     Get a terrain pressure with a given geolocation 
! NAME:
!     Get_TerrainPressure
! INPUTS:
!      lat, lon
! OUTPUTS:
!     thispres
!==============================================================================
  IMPLICIT NONE

! Number of wavelengths and Wavelength arrays
  REAL(KIND=4), INTENT(IN)  :: lat, lon
  REAL(KIND=4), INTENT(OUT) :: thispres

  INTEGER(KIND=4) :: tmpmon
  INTEGER(KIND=4) :: ilon, ilat, STATUS

! -- Extract a terrain pressure value from the TERRAIN PRRESSURE data set ---
   if (lon .lt. 0.0) then 
      ilon = 1800 + int(10*lon + 0.00001)
   endif
   if (lon .ge. 0.0) then
      ilon = 1801 + int(10*lon - 0.00001)
   endif
!
   if (lat .lt. 0.0) then 
      ilat = 900 + int(10*lat + 0.00001)
   endif
   if (lat .ge. 0.0) then 
      ilat = 901 + int(10*lat - 0.00001)
   endif
!   
   if (ilat .gt. 1800) ilat = 1800
   if (ilat .lt. 1)  ilat = 1
   if (ilon .gt. 3600) ilon = 3600
   if (ilon .lt. 1) ilon = 1
!
   thispres = tpres(ilon,ilat)
   if ( (thispres.LT.0.0) .OR. (thispres.GT.1014.0) ) thispres = -1.2676506E+30
   
  STATUS = 1

  RETURN

 END FUNCTION Get_TerrainPressure
!
!==================================================================
!==================================================================
!
END MODULE Get_TerrainPressure_H5module 
