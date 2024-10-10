MODULE GetSurfAlbLUT_H5module_nc4
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

      use netcdf
      USE OCIUAAER_Config_Module

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: lut_fn
      INTEGER(KIND=4) :: nmonth, nlats, nlons
      !
      REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: SRFLER354, SRFLER388
      REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE, INTENT(OUT) :: GRDLAT, GRDLON


      integer, dimension (2) :: start2, edge2, stride2
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
      group_name = 'surfalb'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'SRFLER354'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start3  = (/ 1,1,1 /)
      edge3   = SHAPE(SRFLER354)
      stride3 = (/ 1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, SRFLER354, start=start3, &
         stride=stride3, count=edge3)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'SRFLER388'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge3   = SHAPE(SRFLER388)
      status = nf90_get_var(grp_id, dset_id, SRFLER388, start=start3, &
         stride=stride3, count=edge3)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'GRDLAT'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start2  = (/ 1,1 /)
      edge2   = SHAPE(GRDLAT)
      stride2 = (/ 1,1 /)
      status = nf90_get_var(grp_id, dset_id, GRDLAT, start=start2, &
         stride=stride2, count=edge2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'GRDLON'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge2   = SHAPE(GRDLON)
      status = nf90_get_var(grp_id, dset_id, GRDLON, start=start2, &
         stride=stride2, count=edge2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      status = nf90_close(nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to close lut_nc4 file: ", status
         return
      end if

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
END MODULE GetSurfAlbLUT_H5module_nc4
