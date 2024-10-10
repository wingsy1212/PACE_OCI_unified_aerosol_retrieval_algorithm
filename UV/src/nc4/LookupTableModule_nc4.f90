MODULE LookupTableModule_nc4
   !==============================================================================
   ! FILENAME:
   !     LookupTableModule.f90
   !
   ! DESCRIPTION:
   !     Read look up tables and ancillary data ( binary in the mean time)
   !==============================================================================

   IMPLICIT NONE

   ! Look up table  dimensions
   INTEGER(KIND=4), PARAMETER :: nwave =  2, & !Number of Wavelengths in model
      nzae =  5, & !Number of Aerosol Height Level
      nw0model = 21, & !Number of Profiles in Aerosol mod.
      ntau =  7, & !Number of Optical Depth Values
      nsza =  7, & !Number of Solar Zenith Angles
      nphi = 11, & !Number of Relative Az Angles
      ntheta = 14, & !Number of Satellite Zenith Angles
      nw0sel =  7    !Number of Profiles selected for
                                               !  current calculation

   ! Wavelength table (in nm) *** wave ***
   REAL(KIND=4), DIMENSION(3) :: wave_table
   DATA wave_table/354.0,388.0,500.0/

   ! Aerosol Height table (in km) *** zae ***
   REAL(KIND=4), DIMENSION(5) :: zae_table
   DATA zae_table/0.0,1.5,3.0,6.0,10.0/

   ! Solar Zenith Angle table *** sza ***
   REAL(KIND=4), DIMENSION(7) :: sza_table
   DATA sza_table/0.,20.,40.,60.,66., 72.,80./

   REAL(KIND=4), DIMENSION(9) :: sza_table_ep
   DATA sza_table_ep/0.,20.,40.,60.,66., 72.,80., 84., 88./

   ! Relative Azimuth Angle table *** phi ***
   REAL(KIND=4), DIMENSION(11) :: phi_table
   DATA phi_table/0.,30.,60.,90.,120.,150.,160., 165.,170.,175.,180./

   REAL(KIND=4), DIMENSION(11) :: phi_table_ep
   DATA phi_table_ep /0.,30.,60.,90.,120.,150.,160.,165.,170.,175.,180./

   ! Satellite Zenith Angle table *** theta ***
   REAL(KIND=4), DIMENSION(14):: theta_table
   DATA theta_table/0.,12.,18.,26.,32.,36.,40.,46.,50.,54.,56.,60., 66., 72./

   REAL(KIND=4), DIMENSION(14):: theta_table_ep
   DATA theta_table_ep /0.,12.,18.,26.,32.,36.,40.,46.,50.,54.,56.,60., 66., 72./

   ! Some Month/Day/Year conversion data
   INTEGER(KIND=4), PARAMETER         :: nday_lim=24
   REAL (KIND=4), DIMENSION(nday_lim) :: month_day_lim
   DATA month_day_lim/1,31,32,59,60,90,91,120,121,151,152,181,182,&
      212,213,243,244,273,274,304,305,334,335,366/

   INTEGER(KIND=4), PARAMETER         :: nmonth=12
   REAL (KIND=4), DIMENSION(nmonth)   :: half_month_day
   DATA half_month_day/15,45,74,105,135,166,196,227,258,288,319,349/

   ! Radiances at given Latitude, Longitude and Day of Year
   INTEGER(KIND=4), PARAMETER :: nlat=180,nlon=360,nlon_uv=360

   !  Radiance matrices used by the new interpolation scheme
   REAL (KIND=4), DIMENSION(:), ALLOCATABLE :: radlin_p10, radlin_p06

   !  Transmittance matrices used by the new interpolation scheme
   REAL (KIND=4), DIMENSION(:), ALLOCATABLE :: tralin_p10, tralin_p06
 
   !  Coeficients and indices needed for the interpolation (3D & 2D)
   REAL (KIND=4) :: cofs(64),cofs_2d(16)
   INTEGER (KIND=4) :: indsol,indscn,indphi, indsol_2d,indscn_2d

   ! Memory used to store results of look up table interpolations
   REAL(KIND=4),DIMENSION(:,:,:,:),ALLOCATABLE :: radiance_geo,          &
      radiance_geo_p10,      &
      radiance_geo_p06,      &
      transmittance_geo,     &
      transmittance_geo_p10, &
      transmittance_geo_p06, &
      sphalb,                &
      radiance_toa
   REAL(KIND=4),DIMENSION(:,:,:,:),ALLOCATABLE :: radiance_geo_fmode,          &
      radiance_geo_p10_fmode,      &
      radiance_geo_p06_fmode,      &
      transmittance_geo_fmode,     &
      transmittance_geo_p10_fmode, &
      transmittance_geo_p06_fmode, &
      sphalb_fmode,                &
      radiance_toa_fmode
   REAL(KIND=4),DIMENSION(:,:,:,:),ALLOCATABLE :: radiance_geo_cmode,          &
      radiance_geo_p10_cmode,      &
      radiance_geo_p06_cmode,      &
      transmittance_geo_cmode,     &
      transmittance_geo_p10_cmode, &
      transmittance_geo_p06_cmode, &
      sphalb_cmode,                &
      radiance_toa_cmode

   ! Surface catagory
   INTEGER(KIND=4), PARAMETER :: nlat_highRes=1080, nlon_highRes=2160
   INTEGER(KIND=4), DIMENSION(nlat_highRes,nlon_highRes) :: ref_sfc

   ! Aerosol level heights at given Latitude, Longitude and Day of Year
   REAL (KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: zaer_cm

   ! CALIOP Aerosol level heights at given Latitude, Longitude and Day of Year
   REAL (KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: zaer_cm2

   ! AIRS CO at given Latitude, Longitude and Day of Year
   REAL (KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: airsco_cm, airsco_cm2


   ! Functions
   PUBLIC :: Read_ancillaryLUTs
   PUBLIC :: allocateLUT
   PUBLIC :: deallocateLUT

CONTAINS
   !
   !==============================================================================
   !==============================================================================

   SUBROUTINE Read_ancillaryLUTs()
      use netcdf
      USE OCIUAAER_Config_Module
      implicit none

      !character(*),           intent (in)      :: sfc_file, AIRSCO_clm_file
      integer                                  :: status
      integer, dimension (2)                   :: start2, edge2, stride2
      integer, dimension (3)                   :: start3, edge3, stride3

      character(len=255)     ::  lut_fn
      character(len=255)     ::  sds_name
      character(len=255)    ::  dset_name
      character(len=255)    ::  attr_name
      character(len=255)    ::  group_name

      integer               ::  nc_id
      integer               ::  dim_id
      integer               ::  dset_id
      integer               ::  grp_id
      integer               ::  unitnum, unitnum1
      
      Print*, 'Now reading...', cfg%uv_nc4
      status = nf90_open(cfg%uv_nc4, nf90_nowrite, nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
         return
      end if

      group_name = 'sfc'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'data'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start2  = (/ 1,1 /)
      edge2   = (/ nlat_highRes, nlon_highRes /)
      stride2 = (/ 1,1 /)
      status = nf90_get_var(grp_id, dset_id, ref_sfc, start=start2, &
         stride=stride2, count=edge2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!	Print*, 'Success sfc'

      group_name = 'airsco_clm'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'data'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start3  = (/ 1,1,1 /)
      edge3   = (/ nlon, nlat, nmonth /)
      stride3 = (/ 1,1,1 /) 
      status = nf90_get_var(grp_id, dset_id, airsco_cm, start=start3, &
         stride=stride3, count=edge3)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!      airsco_cm = 0
!    Print*, 'Success with airsco_clm', airsco_cm(10, 10, 10)
    
      group_name = 'zaer'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'data'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start3  = (/ 1,1,1 /)
      edge3   = (/ nmonth, nlon, nlat  /)
      stride3 = (/ 1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, zaer_cm, start=start3, &
         stride=stride3, count=edge3)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!      zaer_cm = 0
!	Print*, 'Success with zaer'
	
      group_name = 'zaer2'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'data'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start3  = (/ 1,1,1 /)
      edge3   = (/ nmonth, nlon, nlat  /)
      stride3 = (/ 1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, zaer_cm2, start=start3, &
         stride=stride3, count=edge3)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!      zaer_cm2 = 0
!	Print*, 'Success with zaer_cm2'

      status = nf90_close(nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to close lut_nc4 file: ", status
         return
      end if

      RETURN

   END SUBROUTINE Read_ancillaryLUTs
 
   !
   !==============================================================================
   !==============================================================================
   !
   FUNCTION GetAerosolLayerHeight(mon,lat,lon,aeroHeight) RESULT(status)

      IMPLICIT NONE

      ! Number of wavelengths and Wavelength arrays
      INTEGER(KIND=4), INTENT(IN)  :: mon
      REAL(KIND=4),    INTENT(IN)  :: lat, lon
      REAL(KIND=4),    INTENT(OUT) :: aeroHeight

      ! For Climatology Lookup table retrieval
      INTEGER(KIND=4) :: latin_visclm, lonin_visclm, status

      STATUS = -1
      latin_visclm = 0.0d00
      lonin_visclm = 0.0d00

      ! Extract Aerosol Layer Height of pixel from climatological database
      IF(lat .GE.  0.0)THEN
         latin_visclm = 89.5 - (int(lat) + 0.5) + 1.0
      ENDIF

      IF(lat .LT. 0.0)THEN
         latin_visclm = 89.5 - (int(lat) - 0.5) + 1.0
      ENDIF

      IF(lon .GE. 0.0)THEN
         lonin_visclm = 179.5 + (int(lon) + 0.5) + 1.0
      ENDIF

      IF(lon .LT. 0.0)THEN
         lonin_visclm = 179.5 + (int(lon) - 0.5) + 1.0
      ENDIF
      if (lonin_visclm .lt. 1) lonin_visclm = 1
      if (lonin_visclm .gt. 360) lonin_visclm = 360
      if (latin_visclm .lt. 1) latin_visclm = 1
      if (latin_visclm .gt. 180) latin_visclm = 180

      aeroHeight= zaer_cm(mon,lonin_visclm, 180-latin_visclm+1)

      STATUS = 1

      RETURN

   END FUNCTION GetAerosolLayerHeight

   !
   !==============================================================================
   !==============================================================================
   !
   FUNCTION GetAerosolLayer_CALIOPHeight(mon,lat,lon,aeroHeight) RESULT(status)

      IMPLICIT NONE

      ! Number of wavelengths and Wavelength arrays
      INTEGER(KIND=4), INTENT(IN)  :: mon
      REAL(KIND=4),    INTENT(IN)  :: lat, lon
      REAL(KIND=4),    INTENT(OUT) :: aeroHeight

      ! For Climatology Lookup table retrieval
      INTEGER(KIND=4) :: latin_visclm, lonin_visclm, STATUS

      STATUS = -1
      latin_visclm = 0.0d00
      lonin_visclm = 0.0d00

      ! Extract Aerosol Layer Height of pixel from climatological database
      IF(lat .GE.  0.0)THEN
         latin_visclm = 89.5 - (int(lat) + 0.5) + 1.0
      ENDIF

      IF(lat .LT. 0.0)THEN
         latin_visclm = 89.5 - (int(lat) - 0.5) + 1.0
      ENDIF

      IF(lon .GE. 0.0)THEN
         lonin_visclm = 179.5 + (int(lon) + 0.5) + 1.0
      ENDIF

      IF(lon .LT. 0.0)THEN
         lonin_visclm = 179.5 + (int(lon) - 0.5) + 1.0
      ENDIF
      if (lonin_visclm .lt. 1) lonin_visclm = 1
      if (lonin_visclm .gt. 360) lonin_visclm = 360
      if (latin_visclm .lt. 1) latin_visclm = 1
      if (latin_visclm .gt. 180) latin_visclm = 180

      aeroHeight= zaer_cm2(mon,lonin_visclm, 180-latin_visclm+1)

      STATUS = 1

      RETURN

   END FUNCTION GetAerosolLayer_CALIOPHeight
 
   !
   !==============================================================================
   !==============================================================================
   !
   FUNCTION GetAIRSCO_Clm(mon,lat,lon,airsco_thismonth) RESULT(status)

      IMPLICIT NONE

      ! Number of wavelengths and Wavelength arrays
      INTEGER(KIND=4), INTENT(IN)  :: mon
      REAL(KIND=4),    INTENT(IN)  :: lat, lon
      REAL(KIND=4),    INTENT(OUT) :: airsco_thismonth

      ! For Climatology Lookup table retrieval
      INTEGER(KIND=4) :: latin_visclm, lonin_visclm, status

      STATUS = -1
      latin_visclm = 0.0d00
      lonin_visclm = 0.0d00

      ! Extract AIRS CO of pixel from climatological database
      IF(lat .GE.  0.0)THEN
         latin_visclm = 89.5 - (int(lat) + 0.5) + 1.0
      ENDIF

      IF(lat .LT. 0.0)THEN
         latin_visclm = 89.5 - (int(lat) - 0.5) + 1.0
      ENDIF

      IF(lon .GE. 0.0)THEN
         lonin_visclm = 179.5 + (int(lon) + 0.5) + 1.0
      ENDIF

      IF(lon .LT. 0.0)THEN
         lonin_visclm = 179.5 + (int(lon) - 0.5) + 1.0
      ENDIF
      if (lonin_visclm .lt. 1) lonin_visclm = 1
      if (lonin_visclm .gt. 360) lonin_visclm = 360
      if (latin_visclm .lt. 1) latin_visclm = 1
      if (latin_visclm .gt. 180) latin_visclm = 180

      airsco_thismonth = airsco_cm(lonin_visclm, 180-latin_visclm+1,mon)

      STATUS = 1

      RETURN

   END FUNCTION GetAIRSCO_Clm
 
   !
   !==============================================================================
   !==============================================================================
   !
   FUNCTION GetSurfaceType(lat,lon,surfType) RESULT(STATUS)
      !
      ! TITLE:    Read the Surface Type database
      ! NAME:     GetSurfaceType
      ! INPUTS:   lat, lon
      ! OUTPUTS:  surfType
      !
      IMPLICIT NONE

      REAL(KIND=4),               INTENT(IN)  :: lat, lon
      INTEGER(KIND=2),            INTENT(OUT) :: surfType

      ! For Climatology Lookup table retrieval
      INTEGER(KIND=4) :: ilon, ilat, STATUS

      STATUS = 1

      !  Use the high resolution surface category reference catalog (1/6 degree bin, 1080 x 2160 )
      IF (lon .lt. 0.0) ilon = 2160 + INT(6.0*lon + 0.000001)
      IF (lon .ge. 0.0) ilon = 1 + INT(6.0*lon - 0.000001)
      IF (lat .lt. 0.0) ilat = 541 - INT(6.0*lat + 0.000001)
      IF (lat .ge. 0.0) ilat = 540 - INT(6.0*lat - 0.000001)
      !
      if (ilat .gt. 1080) ilat = 1080
      if (ilat .lt. 1) ilat = 1
      if (ilon .gt. 2160) ilon = 2160
      if (ilon .lt. 1) ilon = 1
      !
      surfType = ref_sfc(ilat,ilon)

      STATUS = 1

      RETURN

   END FUNCTION GetSurfaceType
 
   !
   !==============================================================================
   !==============================================================================
   !
   FUNCTION allocateLUT() RESULT(STATUS)
      !
      ! TITLE: Allocate memory for look up tables and arrays
      ! NAME:  allocateLUT
      ! INPUTS:  none
      ! OUTPUTS: none
      !

      IMPLICIT NONE

      ! for error messaging
      INTEGER(KIND=4)        :: STATUS

      STATUS = -1

      ! Allocate space for radiance_geo_p10
      ALLOCATE(radiance_geo_p10(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(radiance_geo_p10_fmode(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(radiance_geo_p10_cmode(nwave,nzae,nw0sel,ntau),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) ' Unable to allocate for radiance_geo_p10'
         CALL EXIT(1)
      ENDIF

      ! Allocate space for transmittance_geo_p10
      ALLOCATE(transmittance_geo_p10(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(transmittance_geo_p10_fmode(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(transmittance_geo_p10_cmode(nwave,nzae,nw0sel,ntau),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) ' Unable to allocate for transmittance_geo_p10'
         CALL EXIT(1)
      ENDIF

      ! Allocate space for radiance_geo_p06
      ALLOCATE(radiance_geo_p06(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(radiance_geo_p06_fmode(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(radiance_geo_p06_cmode(nwave,nzae,nw0sel,ntau),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) ' Unable to allocate for radiance_geo_p06'
         CALL EXIT(1)
      ENDIF

      ! Allocate space for transmittance_geo_p06
      ALLOCATE(transmittance_geo_p06(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(transmittance_geo_p06_fmode(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(transmittance_geo_p06_cmode(nwave,nzae,nw0sel,ntau),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) ' Unable to allocate for transmittance_geo_p06'
         CALL EXIT(1)
      ENDIF

      ! Allocate space for radiance_geo
      ALLOCATE(radiance_geo(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(radiance_geo_fmode(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(radiance_geo_cmode(nwave,nzae,nw0sel,ntau),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) ' Unable to allocate for radiance_geo'
         CALL EXIT(1)
      ENDIF

      ! Allocate space for transmittance_geo
      ALLOCATE(transmittance_geo(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(transmittance_geo_fmode(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(transmittance_geo_cmode(nwave,nzae,nw0sel,ntau),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) ' Unable to allocate for transmittance_geo'
         CALL EXIT(1)
      ENDIF

      ! Allocate space for radiance_toa
      ALLOCATE(radiance_toa(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(radiance_toa_fmode(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(radiance_toa_cmode(nwave,nzae,nw0sel,ntau),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) ' Unable to allocate for radiance_toa'
         CALL EXIT(1)
      ENDIF

      ! Allocate space for sphalb
      ALLOCATE(sphalb(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(sphalb_fmode(nwave,nzae,nw0sel,ntau),stat=status)
      ALLOCATE(sphalb_cmode(nwave,nzae,nw0sel,ntau),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) 'Unable to allocate for sphalb '
         CALL EXIT(1)
      ENDIF

      ! Allocate space for zaer_cm (REANALYSIS HGT)
      ALLOCATE(zaer_cm(nmonth,nlon,nlat),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) 'Unable to allocate for zaer_cm '
         CALL EXIT(1)
      ENDIF

      ! Allocate space for zaer_cm2 (CALIOP HGT)
      ALLOCATE(zaer_cm2(nmonth,nlon,nlat),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) 'Unable to allocate for zaer_cm2 '
         CALL EXIT(1)
      ENDIF

      ! Allocate space for airsco_cm (AIRS CO monthly climatology)
      ALLOCATE(airsco_cm(nlon,nlat,nmonth),stat=status)

      IF(status /= 0) THEN
         WRITE(*,*) 'Unable to allocate for airsco_cm '
         CALL EXIT(1)
      ENDIF


      STATUS = 1

      RETURN

   END FUNCTION allocateLUT
   !
   !==============================================================================
   !==============================================================================
   !
   FUNCTION deallocateLUT() RESULT(STATUS)
      !
      ! TITLE:   Deallocate memory for look up tables and arrays
      ! NAME:    deallocateLUT
      ! INPUTS:  none
      ! OUTPUTS: none
      !

      IMPLICIT NONE

      ! for error messaging
      INTEGER(KIND=4)    :: STATUS

      STATUS = -1

      !==============================================================================
      ! Deallocate space for airsco_cm
      !==============================================================================
      IF(ALLOCATED(airsco_cm))DEALLOCATE(airsco_cm,stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate airsco_cm'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for zaer_cm2
      !==============================================================================
      IF(ALLOCATED(zaer_cm2))DEALLOCATE(zaer_cm2,stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate zaer_cm2'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for zaer_cm
      !==============================================================================
      IF(ALLOCATED(zaer_cm))DEALLOCATE(zaer_cm,stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate zaer_cm'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for sphalb
      !==============================================================================
      IF(ALLOCATED(sphalb))DEALLOCATE(sphalb,stat=status)
      IF(ALLOCATED(sphalb_fmode))DEALLOCATE(sphalb_fmode,stat=status)
      IF(ALLOCATED(sphalb_cmode))DEALLOCATE(sphalb_fmode,stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate sphalb'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for radiance_toa
      !==============================================================================
      IF(ALLOCATED(radiance_toa))DEALLOCATE(radiance_toa,stat=status)
      IF(ALLOCATED(radiance_toa_fmode))DEALLOCATE(radiance_toa_fmode,stat=status)
      IF(ALLOCATED(radiance_toa_cmode))DEALLOCATE(radiance_toa_fmode,stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate radiance_toa'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for transmittance_geo
      !==============================================================================
      IF(ALLOCATED(transmittance_geo))DEALLOCATE(transmittance_geo,stat=status)
      IF(ALLOCATED(transmittance_geo_fmode))DEALLOCATE(transmittance_geo_fmode,stat=status)
      IF(ALLOCATED(transmittance_geo_cmode))DEALLOCATE(transmittance_geo_cmode,stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate transmittance_geo'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for radiance_geo
      !==============================================================================
      IF(ALLOCATED(radiance_geo))DEALLOCATE(radiance_geo,stat=status)
      IF(ALLOCATED(radiance_geo_fmode))DEALLOCATE(radiance_geo_fmode,stat=status)
      IF(ALLOCATED(radiance_geo_cmode))DEALLOCATE(radiance_geo_cmode,stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate radiance_geo'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for transmittance_geo_p06
      !==============================================================================
      IF(ALLOCATED(transmittance_geo_p06))DEALLOCATE(transmittance_geo_p06, &
         stat=status)
      IF(ALLOCATED(transmittance_geo_p06_fmode))DEALLOCATE(transmittance_geo_p06_fmode, &
         stat=status)
      IF(ALLOCATED(transmittance_geo_p06_cmode))DEALLOCATE(transmittance_geo_p06_cmode, &
         stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate transmittance_geo_p06'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for radiance_geo_p06
      !==============================================================================
      IF(ALLOCATED(radiance_geo_p06))DEALLOCATE(radiance_geo_p06,stat=status)
      IF(ALLOCATED(radiance_geo_p06_fmode))DEALLOCATE(radiance_geo_p06_fmode,stat=status)
      IF(ALLOCATED(radiance_geo_p06_cmode))DEALLOCATE(radiance_geo_p06_cmode,stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate radiance_geo_p06'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for transmittance_geo_p10
      !==============================================================================
      IF(ALLOCATED(transmittance_geo_p10))DEALLOCATE(transmittance_geo_p10, &
         stat=status)
      IF(ALLOCATED(transmittance_geo_p10_fmode))DEALLOCATE(transmittance_geo_p10_fmode, &
         stat=status)
      IF(ALLOCATED(transmittance_geo_p10_cmode))DEALLOCATE(transmittance_geo_p10_cmode, &
         stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate transmittance_geo_p10'
         RETURN
      ENDIF

      !==============================================================================
      ! Deallocate space for radiance_geo_p10
      !==============================================================================
      IF(ALLOCATED(radiance_geo_p10))DEALLOCATE(radiance_geo_p10,stat=status)
      IF(ALLOCATED(radiance_geo_p10_fmode))DEALLOCATE(radiance_geo_p10_fmode,stat=status)
      IF(ALLOCATED(radiance_geo_p10_cmode))DEALLOCATE(radiance_geo_p10_cmode,stat=status)
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate radiance_geo_p10'
         RETURN
      ENDIF

      !******************************************************************************
      !  Deallocate new variables needed for the new interpolation scheme

      !------------------------------------------------------------------------------
      !  Deallocate space for radlin_p10
      !------------------------------------------------------------------------------
      IF(ALLOCATED(radlin_p10))DEALLOCATE(radlin_p10,stat=status)          !  ?
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate radlin_p10'
         RETURN
      ENDIF

      !------------------------------------------------------------------------------
      !  Deallocate space for radlin_p06
      !------------------------------------------------------------------------------
      IF(ALLOCATED(radlin_p06))DEALLOCATE(radlin_p06,stat=status)          !  ?
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate radlin_p06'
         RETURN
      ENDIF

      !------------------------------------------------------------------------------
      !  Deallocate space for tralin_p10
      !------------------------------------------------------------------------------
      IF(ALLOCATED(tralin_p10))DEALLOCATE(tralin_p10,stat=status)          !  ?
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate tralin_p10'
         RETURN
      ENDIF

      !------------------------------------------------------------------------------
      !  Deallocate space for tralin_p06
      !------------------------------------------------------------------------------
      IF(ALLOCATED(tralin_p06))DEALLOCATE(tralin_p06,stat=status)          !  ?
      IF(status /= 0) THEN
         WRITE (*, '(A)') ' Unable to deallocate tralin_p06'
         RETURN
      ENDIF
      !******************************************************************************


      STATUS = 1

      RETURN

   END FUNCTION deallocateLUT
!
!==============================================================================
!==============================================================================
!
END MODULE LookupTableModule_nc4
