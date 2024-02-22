MODULE Get_OceanLUT_H5module
   !==============================================================================
   !
   ! FILENAME:
   !     Get_OceanLUT_module.f90
   !
   ! DESCRIPTION:
   !     This module read the look-up table in he4
   !
   ! AUTHORS:
   !     Changwoo Ahn / Science Systems and Applications, Inc.
   !
   ! HISTORY: Mar 28, 2014
   !==============================================================================

   IMPLICIT NONE

   INTEGER(KIND = 4), PARAMETER :: nwave_ocean = 2 ! 380, 340 nm
   INTEGER(KIND = 4), PARAMETER :: nsza_ocean = 16
   INTEGER(KIND = 4), PARAMETER :: nvza_ocean = 16
   INTEGER(KIND = 4), PARAMETER :: nraa_ocean = 16
 
   REAL(KIND=4), DIMENSION(nwave_ocean, nsza_ocean, nvza_ocean, nraa_ocean) :: oceanler
   REAL(KIND=4), DIMENSION(nsza_ocean) :: sza_oceanset
   REAL(KIND=4), DIMENSION(nvza_ocean) :: vza_oceanset
   REAL(KIND=4), DIMENSION(nraa_ocean) :: raa_oceanset
 
CONTAINS

   !==============================================================================
   !==============================================================================

   SUBROUTINE Read_OceanLUTparams(lut_fn, oceanler, nwave_ocean, &
      nsza_ocean, nvza_ocean, nraa_ocean,&
      sza_oceanset, vza_oceanset, raa_oceanset)
      use netcdf
      USE OCIUAAER_Config_Module
   
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: lut_fn
      !
      REAL(KIND=4), DIMENSION(nwave_ocean, nsza_ocean, nvza_ocean, nraa_ocean), &
         INTENT(OUT) :: oceanler
      REAL(KIND=4), DIMENSION(nsza_ocean), INTENT(OUT) :: sza_oceanset
      REAL(KIND=4), DIMENSION(nvza_ocean), INTENT(OUT) :: vza_oceanset
      REAL(KIND=4), DIMENSION(nraa_ocean), INTENT(OUT) :: raa_oceanset
      INTEGER(KIND=4):: nwave_ocean,nsza_ocean,nvza_ocean,nraa_ocean
         !
      integer, dimension (1) :: start1, edge1, stride1
      integer, dimension (4) :: start4, edge4, stride4

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

      group_name = 'ocncorr'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'OceanLER_ai'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start4  = (/ 1,1,1,1 /)
      edge4   = SHAPE(oceanler)
      stride4 = (/ 1,1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, oceanler, start=start4, &
         stride=stride4, count=edge4)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'NumberOfSolarZenithAngle'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start1  = (/ 1 /)
      edge1   = SHAPE(sza_oceanset)
      stride1 = (/ 1 /)
      status = nf90_get_var(grp_id, dset_id, sza_oceanset, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'NumberOfViewingZenithAngle'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = SHAPE(vza_oceanset)
      status = nf90_get_var(grp_id, dset_id, vza_oceanset, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'NumberOfRelativeAzimuthAngle'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = SHAPE(raa_oceanset)
      status = nf90_get_var(grp_id, dset_id, raa_oceanset, start=start1, &
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
   !
   END SUBROUTINE Read_OceanLUTparams

END MODULE Get_OceanLUT_H5module 
