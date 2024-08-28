MODULE GetLUT_H5module_nc4
   !==============================================================================
   !
   ! FILENAME:
   !     GetLUT_module.f90
   !
   ! DESCRIPTION:
   !     This module read the look-up table in he4
   !
   ! AUTHORS:
   !     Changwoo Ahn / Science Systems and Applications, Inc.
   !
   ! HISTORY: April 20, 2012
   !==============================================================================

   IMPLICIT NONE

   INTEGER(KIND = 4), PARAMETER :: nwaves = 2 ! 354, 380 nm
   INTEGER(KIND = 4), PARAMETER :: nzhgt = 5
   INTEGER(KIND = 4), PARAMETER :: nw0 = 21
   INTEGER(KIND = 4), PARAMETER :: naod = 7
   INTEGER(KIND = 4), PARAMETER :: nsza1 = 7
   INTEGER(KIND = 4), PARAMETER :: nraa = 11
   INTEGER(KIND = 4), PARAMETER :: nvza = 14
 
   REAL(KIND=4), DIMENSION(nwaves,nzhgt,nw0,naod,nsza1,nraa,nvza):: radp10, radp6
   REAL(KIND=4), DIMENSION(nwaves,nzhgt,nw0,naod,nsza1,nvza) :: trp10, trp6
   REAL(KIND=4), DIMENSION(nwaves,nzhgt,nw0,naod) :: sbp10, sbp6

CONTAINS

   SUBROUTINE ReadLUTparams(lut_fn)

      use netcdf
      USE OCIUAAER_Config_Module
   
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: lut_fn

      integer, dimension (4) :: start4, edge4, stride4
      integer, dimension (6) :: start6, edge6, stride6
      integer, dimension (7) :: start7, edge7, stride7

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

      group_name = 'aer'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if
      dset_name = 'RadianceP10'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start7  = (/ 1,1,1,1,1,1,1 /)
      edge7   = SHAPE(radp10)
      stride7 = (/ 1,1,1,1,1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, radp10, start=start7, &
         stride=stride7, count=edge7)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'RadianceP6'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge7   = SHAPE(radp6)
      status = nf90_get_var(grp_id, dset_id, radp6, start=start7, &
         stride=stride7, count=edge7)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'TransmittanceP10'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start6  = (/ 1,1,1,1,1,1 /)
      edge6   = SHAPE(trp10)
      stride6 = (/ 1,1,1,1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, trp10, start=start6, &
         stride=stride6, count=edge6)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'TransmittanceP6'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge6   = SHAPE(trp6)
      status = nf90_get_var(grp_id, dset_id, trp6, start=start6, &
         stride=stride6, count=edge6)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'SphericalAlbedoP10'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start4  = (/ 1,1,1,1 /)
      edge4   = SHAPE(sbp10)
      stride4 = (/ 1,1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, sbp10, start=start4, &
         stride=stride4, count=edge4)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'SphericalAlbedoP6'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge4   = SHAPE(sbp6)
      status = nf90_get_var(grp_id, dset_id, sbp6, start=start4, &
         stride=stride4, count=edge4)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      status = nf90_close(nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to close lut_nc4 file: ", status
         return
      end if

   END SUBROUTINE ReadLUTparams

END MODULE GetLUT_H5module_nc4
