MODULE Get_ssaclimLUT2_H5module_nc4
   !==============================================================================
   !
   ! FILENAME:
   !     Get_ssaval_modlue_Hiren.f90
   !
   ! DESCRIPTION:
   !     This module read the OMI ssa at 388 nm table (omaca_ssa388_v111.he4) in he4
   !
   ! AUTHORS:
   !     Hiren Jethva / USRA-GESTAR
   !==============================================================================

   USE InterpolationModule
   use netcdf
   USE OCIUAAER_Config_Module
 
   IMPLICIT NONE

   REAL(KIND=4), DIMENSION(:,:,:,:,:), ALLOCATABLE        :: ssa388_smk, ssa388_dst
   REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE            :: ssa388_clim_smk, ssa388_clim_dst
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE                :: ssalon_table, ssalat_table

CONTAINS

SUBROUTINE Read_ssaclimLUTparams(lut_fn, ssa388_smk, ssa388_dst, &
                               ssa388_clim_smk, ssa388_clim_dst, &
			       ssalon_table, ssalat_table)

 IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
  INTEGER(KIND=4) :: fid, swid, STATUS, err
  INTEGER(KIND=4) :: nxgrids,nygrids,nyears, nmonths, ndays
!
  REAL(KIND=4), DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(OUT) :: ssa388_smk, ssa388_dst
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: ssa388_clim_smk, ssa388_clim_dst
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ssalon_table, ssalat_table

      integer, dimension (1) :: start1, edge1, stride1
      integer, dimension (3) :: start3, edge3, stride3
      integer, dimension (5) :: start5, edge5, stride5

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

      group_name = 'ssa388_clim'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if
      nxgrids = 360
      nygrids = 180
      nyears  = 19
      nmonths = 12
      ndays   = 31

      !! Allocate memory
      ALLOCATE( ssa388_smk( nxgrids, nygrids, nyears, nmonths, ndays), &
                ssa388_dst( nxgrids, nygrids, nyears, nmonths, ndays), &
           ssa388_clim_smk( nxgrids, nygrids, nmonths), &
           ssa388_clim_dst( nxgrids, nygrids, nmonths), &
              ssalon_table( nxgrids), &
	      ssalat_table( nygrids), STAT = status)

      IF (STATUS < 0) THEN
         PRINT *,'Error : Allocation of variables for SSA388 Clim LUT failed.'
         CALL EXIT(1)
      ENDIF

      dset_name = 'NumberOfXGrids'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start1  = (/ 1 /)
      edge1   = SHAPE(ssalon_table)
      stride1 = (/ 1 /)
      status = nf90_get_var(grp_id, dset_id, ssalon_table, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!      
      dset_name = 'NumberOfYGrids'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = SHAPE(ssalat_table)
      status = nf90_get_var(grp_id, dset_id, ssalat_table, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!
      dset_name = 'OMI_Gridded_SSA388_Smoke_For_OMACA'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start5  = (/ 1,1,1,1,1 /)
      edge5   = SHAPE(ssa388_smk)
      stride5 = (/ 1,1,1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, ssa388_smk, start=start5, &
         stride=stride5, count=edge5)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!      
      dset_name = 'OMI_Gridded_SSA388_Dust_For_OMACA'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge5   = SHAPE(ssa388_dst)
      status = nf90_get_var(grp_id, dset_id, ssa388_dst, start=start5, &
         stride=stride5, count=edge5)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!
      dset_name = 'OMI_Gridded_SSA388_Smoke_Clim_For_OMACA'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start3  = (/ 1,1,1 /)
      edge3   = SHAPE(ssa388_clim_smk)
      stride3 = (/ 1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, ssa388_clim_smk, start=start3, &
         stride=stride3, count=edge3)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!
      dset_name = 'OMI_Gridded_SSA388_Dust_Clim_For_OMACA'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge3   = SHAPE(ssa388_clim_dst)
      status = nf90_get_var(grp_id, dset_id, ssa388_clim_dst, start=start3, &
         stride=stride3, count=edge3)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
!
      status = nf90_close(nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to close lut_nc4 file: ", status
         return
      end if

   END SUBROUTINE Read_ssaClimLUTparams


  FUNCTION GetOMI_SSAlbedoClimValue(yr, mon, date, inlat, inlon, inatype, ssa354val, ssa388val, ssa500val) RESULT(status)
!==============================================================================
! TITLE:
!     Find the OMI Single Scattering Albedo Value (388 nm) from the Lookup Table
!
! NAME:
!     GetOMI_SSAlbedoClimValue
!
! INPUTS:
!     year, mon, date, lat, lon, atype
!
! OUTPUTS:
!     ssa388val_smk(:)
!     ssa388val_dst(:)
!==============================================================================

  IMPLICIT NONE

!==============================================================================
!  Bounding values in table for input parameters
!==============================================================================
  REAL(KIND=4)    :: theta1,theta2,sza1,sza2,phi1,phi2,w1,w2
  REAL(KIND=4)    :: fractionalLat, fractionalLon 

!==============================================================================
!  Indices for bounding values in table for input parameters
!==============================================================================
  INTEGER(KIND=4) :: jtheta, jtheta2, jsza, jsza2, jphi, jphi2, jw1


!==============================================================================
! Number of wavelengths and Wavelength arrays
!==============================================================================
  INTEGER(KIND=4),            INTENT(IN)  :: yr, mon, date
  REAL(KIND=4),               INTENT(IN)  :: inlat, inlon
  INTEGER(KIND=2),            INTENT(IN)  :: inatype
  REAL(KIND=4),               INTENT(OUT) :: ssa354val, ssa388val, ssa500val
!
  INTEGER(KIND=4)                         :: yridx, mnthidx, dtidx
  INTEGER(KIND=4)                         :: iXGrd, jYGrd, status

!==============================================================================
! For Climatology Lookup table retrieval
!==============================================================================
  INTEGER(KIND=4) :: ilon, ilat

!   Added wavelength conversion coefficients for SSA on May 13, 2016---H. Jethva (USRA)
    REAL(KIND=4), DIMENSION(5)              :: coef_ssa354_smk, coef_ssa354_dst
    REAL(KIND=4), DIMENSION(5)              :: coef_ssa500_smk, coef_ssa500_dst
    DATA coef_ssa354_smk/10.067429,-41.106924,65.361537,-44.802955,11.480581/ 
    DATA coef_ssa500_smk/-25.404624,123.45677,-219.03895,172.84429,-50.858004/

    DATA coef_ssa354_dst/3.0061727,-12.641012,23.569893,-18.655039,5.7199183/
    DATA coef_ssa500_dst/10.715972,-49.210211,88.437284,-68.556941,19.613582/


  status = -1

!Daily/Monthly/Clim SSA388 search for data withing 2004-2022

!*************************************!
IF (yr .GE. 2004 .AND. yr .LE. 2022) THEN

yridx = yr-2003
mnthidx = mon
dtidx = date

! Initialize ssa388val
  ssa354val=0.
  ssa388val=0.
  ssa500val=0.

  DO iXGrd = 1, 360
     DO jYGrd = 1, 180

     IF(inlon .GE. ssalon_table(iXGrd)-0.5 .AND. inlon .LE. ssalon_table(iXGrd)+0.5 .AND. &
        inlat .GE. ssalat_table(jYGrd)-0.5 .AND. inlat .LE. ssalat_table(jYGrd)+0.5)THEN

        IF(inatype .EQ. 1 .AND. ssa388_smk(iXGrd,jYGrd,yridx,mnthidx,dtidx) .GT. 0)THEN
           ssa388val = ssa388_smk(iXGrd,jYGrd,yridx,mnthidx,dtidx)
        ENDIF
        IF(inatype .EQ. 2 .AND. ssa388_dst(iXGrd,jYGrd,yridx,mnthidx,dtidx) .GT. 0)THEN
           ssa388val = ssa388_dst(iXGrd,jYGrd,yridx,mnthidx,dtidx)
        ENDIF


        IF(inatype .EQ. 1 .AND. ssa388_smk(iXGrd,jYGrd,yridx,mnthidx,dtidx) .LE. 0)THEN
           ssa388val = ssa388_clim_smk(iXGrd,jYGrd,mnthidx)
        ENDIF

        IF(inatype .EQ. 2 .AND. ssa388_dst(iXGrd,jYGrd,yridx,mnthidx,dtidx) .LE. 0)THEN
           ssa388val = ssa388_clim_dst(iXGrd,jYGrd,mnthidx)
        ENDIF

     ENDIF

    ENDDO
  ENDDO


ENDIF   !IF(yr .GE. 2004 .AND. yr .LE. 2022)THEN
!*************************************!



!Use of climatology value of 'ssa388val' beyond the temporal coverage of daily regional dataset used in the previous step 

!*******************!
IF(yr .GE. 2023)THEN

mnthidx = mon

! Initialize ssa388val
  ssa388val=0.
  ssa354val=0.
  ssa500val=0.

DO iXGrd = 1, 360
     DO jYGrd = 1, 180

     IF(inlon .GE. ssalon_table(iXGrd)-0.5 .AND. inlon .LE. ssalon_table(iXGrd)+0.5 .AND. &
        inlat .GE. ssalat_table(jYGrd)-0.5 .AND. inlat .LE. ssalat_table(jYGrd)+0.5)THEN

        IF(inatype .EQ. 1.) THEN
           ssa388val = ssa388_clim_smk(iXGrd,jYGrd,mnthidx)
        ENDIF

        IF(inatype .EQ. 2.)THEN
           ssa388val = ssa388_clim_dst(iXGrd,jYGrd,mnthidx)
        ENDIF

     ENDIF

    ENDDO
  ENDDO

ENDIF
!*******************!


!****************************************!
! Conversion of 'ssa388val' to 'ssa354val'
!****************************************!
  IF(inatype .EQ. 1)THEN
    ssa354val = coef_ssa354_smk(1) + &
                coef_ssa354_smk(2)*ssa388val + &
	        coef_ssa354_smk(3)*ssa388val**2.0 + &
	        coef_ssa354_smk(4)*ssa388val**3.0 + &
	        coef_ssa354_smk(5)*ssa388val**4.0
  ENDIF

  IF(inatype .EQ. 2)THEN
    ssa354val = coef_ssa354_dst(1) + &
                coef_ssa354_dst(2)*ssa388val + &
		coef_ssa354_dst(3)*ssa388val**2.0 + &
		coef_ssa354_dst(4)*ssa388val**3.0 + &
		coef_ssa354_dst(5)*ssa388val**4.0
  ENDIF


!****************************************!
! Conversion of 'ssa388val' to 'ssa500val'
!****************************************!
  IF(inatype .EQ. 1)THEN
    ssa500val = coef_ssa500_smk(1) + &
                coef_ssa500_smk(2)*ssa388val + &
		coef_ssa500_smk(3)*ssa388val**2.0 + &
		coef_ssa500_smk(4)*ssa388val**3.0 + &
		coef_ssa500_smk(5)*ssa388val**4.0
  ENDIF

  IF(inatype .EQ. 2)THEN
    ssa500val = coef_ssa500_dst(1) + &
                coef_ssa500_dst(2)*ssa388val + &
		coef_ssa500_dst(3)*ssa388val**2.0 + &
		coef_ssa500_dst(4)*ssa388val**3.0 + &
		coef_ssa500_dst(5)*ssa388val**4.0
  ENDIF


  status = 1

  RETURN

 END FUNCTION GetOMI_SSAlbedoClimValue

END MODULE Get_ssaclimLUT2_H5module_nc4
