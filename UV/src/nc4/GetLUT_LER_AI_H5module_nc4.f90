MODULE GetLUT_LER_AI_H5module
   !==============================================================================
   ! DESCRIPTION:
   !     This module reads the AI_LER LUT H5 file.
   !==============================================================================

   IMPLICIT NONE

   REAL(KIND=4), DIMENSION(:,:,:,:,:), ALLOCATABLE :: rad354_ler, rad388_ler
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: rad_lin354_ai_ler, rad_lin388_ai_ler
   !
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: SURFALBSET_ler, PRESSURESET_ler
   INTEGER(KIND = 4)                       :: nsalb_ler, nplev_ler, nsza_ler

CONTAINS

   !
   !==================================================================
   !==================================================================
   !
   SUBROUTINE ReadLUT_LER_AIparams(lut_fn)

      use netcdf
      USE OCIUAAER_Config_Module

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: lut_fn
      INTEGER(KIND=4) :: nplev, nsza, nvza, nraa
      !
      integer, dimension (1) :: start1, edge1, stride1
      integer, dimension (5) :: start5, edge5, stride5

      integer               ::  status
      character(len=255)    ::  sds_name
      character(len=255)    ::  dset_name
      character(len=255)    ::  attr_name
      character(len=255)    ::  group_name

      integer               ::  nc_id
      integer               ::  dim_id
      integer               ::  dset_id
      integer               ::  grp_id

      nplev_ler = 4
      nsalb_ler = 10
      nsza_ler = 9
      nvza = 14
      nraa = 11

      !! Allocate memory
      ALLOCATE( rad354_ler( nplev_ler,nsalb_ler, nsza_ler,nvza,nraa ), &
         rad388_ler( nplev_ler,nsalb_ler, nsza_ler,nvza,nraa ), &
         SURFALBSET_ler(nsalb_ler), PRESSURESET_ler(nplev_ler),  STAT = STATUS)
      IF (STATUS < 0) THEN
         PRINT *,'Error : Allocation of variables for HE4 LUT AI failed.'
         CALL EXIT(1)
      ENDIF

      status = nf90_open(cfg%uv_nc4, nf90_nowrite, nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to open UV lut_nc4 file: ", status
         return
      end if

      group_name = 'ai_ler'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'MolecularAtmosphere354'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start5  = (/ 1,1,1,1,1 /)
      edge5   = SHAPE(rad354_ler)
      stride5 = (/ 1,1,1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, rad354_ler, start=start5, &
         stride=stride5, count=edge5)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'MolecularAtmosphere388'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge5   = SHAPE(rad388_ler)
      status = nf90_get_var(grp_id, dset_id, rad388_ler, start=start5, &
         stride=stride5, count=edge5)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'NumberOfSceneReflectivity'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start1  = (/ 1 /)
      edge1   = SHAPE(surfalbset_ler)
      stride1 = (/ 1 /)
      status = nf90_get_var(grp_id, dset_id, surfalbset_ler, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'NumberOfPressure'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = SHAPE(pressureset_ler)
      status = nf90_get_var(grp_id, dset_id, pressureset_ler, start=start1, &
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
      ALLOCATE(rad_lin354_ai_ler(nplev_ler*nsalb_ler*nsza_ler*nvza*nraa),&
         rad_lin388_ai_ler(nplev_ler*nsalb_ler*nsza_ler*nvza*nraa) )
      !
      rad_lin354_ai_ler(:) = RESHAPE(rad354_ler, (/nplev_ler*nsalb_ler*nsza_ler*nvza*nraa/) )
      rad_lin388_ai_ler(:) = RESHAPE(rad388_ler, (/nplev_ler*nsalb_ler*nsza_ler*nvza*nraa/) )
      !
   END SUBROUTINE ReadLUT_LER_AIparams
!
!==================================================================
!==================================================================
!
END MODULE GetLUT_LER_AI_H5module 
