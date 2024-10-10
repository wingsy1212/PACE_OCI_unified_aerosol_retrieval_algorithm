MODULE Get_TerrainPressure_H5module_nc4
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

      use netcdf
      USE OCIUAAER_Config_Module

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: lut_fn
      !
      REAL(KIND=4), DIMENSION(nlon, nlat), INTENT(OUT) :: tpres
      REAL(KIND=4), DIMENSION(nlon), INTENT(OUT) :: terrain_Longitude
      REAL(KIND=4), DIMENSION(nlat), INTENT(OUT) :: terrain_Latitude

      integer, dimension (1) :: start1, edge1, stride1
      integer, dimension (2) :: start2, edge2, stride2

      integer               ::  status
      character(len=255)    ::  sds_name
      character(len=255)    ::  dset_name
      character(len=255)    ::  attr_name
      character(len=255)    ::  group_name

      integer               ::  nc_id
      integer               ::  dim_id
      integer               ::  dset_id
      integer               ::  grp_id

      status = nf90_open(cfg%uv_nc4, nf90_nowrite, nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to open UV lut_nc4 file: ", status
         return
      end if

      group_name = 'surfprs'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if
      dset_name = 'TERR_PRES'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start2  = (/ 1,1 /)
      edge2   = SHAPE(tpres)
      stride2 = (/ 1,1 /)
      status = nf90_get_var(grp_id, dset_id, tpres, start=start2, &
         stride=stride2, count=edge2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'LON'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start1  = (/ 1 /)
      edge1   = SHAPE(terrain_Longitude)
      stride1 = (/ 1 /)
      status = nf90_get_var(grp_id, dset_id, terrain_Longitude, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'LAT'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = SHAPE(terrain_Latitude)
      status = nf90_get_var(grp_id, dset_id, terrain_Latitude, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      status = nf90_close(nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to close lut_nc4 file: ", status
         return
      end if


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
END MODULE Get_TerrainPressure_H5module_nc4 
