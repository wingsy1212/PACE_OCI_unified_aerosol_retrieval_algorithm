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

      use netcdf
      USE OCIUAAER_Config_Module
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: lut_fn
       
      REAL(KIND=4), DIMENSION(nlon, nlat, nmonth), INTENT(OUT) :: snowice
      REAL(KIND=4), DIMENSION(nlon), INTENT(OUT) :: swice_Longitude
      REAL(KIND=4), DIMENSION(nlat), INTENT(OUT) :: swice_Latitude
 
      integer, dimension (1) :: start1, edge1, stride1
      integer, dimension (3) :: start3, edge3, stride3

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

      group_name = 'snowice'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if
      dset_name = 'SNOW_ICE'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start3  = (/ 1,1,1 /)
      edge3   = SHAPE(snowice)
      stride3 = (/ 1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, snowice, start=start3, &
         stride=stride3, count=edge3)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'GRID_LON'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start1  = (/ 1 /)
      edge1   = SHAPE(swice_Latitude)
      stride1 = (/ 1 /)
      status = nf90_get_var(grp_id, dset_id, swice_Latitude, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'GRID_LAT'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = SHAPE(swice_Longitude)
      status = nf90_get_var(grp_id, dset_id, swice_Longitude, start=start1, &
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
