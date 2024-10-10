MODULE GetLUT_Miecloud_AI_H5module_nc4
   !==============================================================================
   ! DESCRIPTION:
   !     This module reads AI_Miecloud LUT H5 file.
   !==============================================================================

   IMPLICIT NONE
 
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: rad_lin354_ai, rad_lin388_ai
   REAL(KIND=4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: rad354, rad388
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: SURFALBSET, CODSET
   INTEGER(KIND = 4)                       :: nplev_mie, nsalb_mie, ncod_mie, nsza_mie, nvza_mie, nraa_mie
   !
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: MieRad_lin2p3_ai, MieRad_lin388_ai
   REAL(KIND=4), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: MieRad2p3, MieRad388
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: plev2p3_nodes, salb2p3_nodes, cod2p3_nodes
   REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: sza2p3_nodes, vza2p3_nodes, raa2p3_nodes
   INTEGER(KIND = 4)                       :: nplev_swir, nsalb_swir, ncod_swir, nsza_swir, nvza_swir, nraa_swir
!

CONTAINS
   !
   !==================================================================
   !==================================================================
   !
   SUBROUTINE ReadLUTAIparams(lut_fn)

      use netcdf
      USE OCIUAAER_Config_Module

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: lut_fn
      !
      integer, dimension (1) :: start1, edge1, stride1
      integer, dimension (6) :: start6, edge6, stride6

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

      nplev_mie = 3
      nsalb_mie = 5
      ncod_mie = 9
      nsza_mie = 9
      nvza_mie = 14
      nraa_mie = 11

      !! Allocate memory
      ALLOCATE( rad354( nplev_mie,nsalb_mie, ncod_mie, nsza_mie,nvza_mie,nraa_mie ), &
         rad388( nplev_mie,nsalb_mie, ncod_mie, nsza_mie,nvza_mie,nraa_mie ), &
         SURFALBSET( nsalb_mie ), CODSET( ncod_mie),  STAT = STATUS)
      IF (STATUS < 0) THEN
         PRINT *,'Error : Allocation of variables for AI MieCloud LUT failed.'
         CALL EXIT(1)
      ENDIF

      group_name = 'ai_mie'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if
      dset_name = 'MieCloud354'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start6  = (/ 1,1,1,1,1,1 /)
      edge6   = SHAPE(rad354)
      stride6 = (/ 1,1,1,1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, rad354, start=start6, &
         stride=stride6, count=edge6)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'MieCloud388'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge6   = SHAPE(rad388)
      status = nf90_get_var(grp_id, dset_id, rad388, start=start6, &
         stride=stride6, count=edge6)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'NumberOfSurfaceReflectivity'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start1  = (/ 1 /)
      edge1   = SHAPE(surfalbset)
      stride1 = (/ 1 /)
      status = nf90_get_var(grp_id, dset_id, surfalbset, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'NumberOfCloudOpticalDepth'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = SHAPE(codset)
      status = nf90_get_var(grp_id, dset_id, codset, start=start1, &
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
      ALLOCATE(rad_lin354_ai(nplev_mie*nsalb_mie*ncod_mie*nsza_mie*nvza_mie*nraa_mie),&
         rad_lin388_ai(nplev_mie*nsalb_mie*ncod_mie*nsza_mie*nvza_mie*nraa_mie) )
      !
      rad_lin354_ai(:) = RESHAPE(rad354, (/nplev_mie*nsalb_mie*ncod_mie*nsza_mie*nvza_mie*nraa_mie/) )
      rad_lin388_ai(:) = RESHAPE(rad388, (/nplev_mie*nsalb_mie*ncod_mie*nsza_mie*nvza_mie*nraa_mie/) )

   END SUBROUTINE ReadLUTAIparams



   !
   !==================================================================
   !==================================================================
   !
   SUBROUTINE Read_SWIRLUTAIparams(lut_fn)

      use netcdf
      USE OCIUAAER_Config_Module

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: lut_fn
      !

      integer, dimension (6) :: start6, edge6, stride6

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

      nplev_swir = 2
      nsalb_swir = 5
      ncod_swir = 9
      nsza_swir = 9
      nvza_swir = 14
      nraa_swir = 11

      !! Allocate memory
      ALLOCATE( MieRad2p3( nplev_swir,nsalb_swir, ncod_swir, nsza_swir,nvza_swir,nraa_swir ), &
         MieRad388( nplev_swir,nsalb_swir, ncod_swir, nsza_swir,nvza_swir,nraa_swir ), &
         plev2p3_nodes(nplev_swir), salb2p3_nodes(nsalb_swir), cod2p3_nodes(ncod_swir),  &
         sza2p3_nodes(nsza_swir), vza2p3_nodes(nvza_swir), raa2p3_nodes(nraa_swir), STAT = STATUS)

      IF (STATUS < 0) THEN
         PRINT *,'Error : Allocation of variables for AI swirCloud LUT failed.'
         CALL EXIT(1)
      ENDIF
      !
      plev2p3_nodes(1:nplev_swir) = [1013.25, 600.0]
      salb2p3_nodes(1:nsalb_swir) = [0.00, 0.10, 0.20, 0.30, 0.40]
      cod2p3_nodes(1: ncod_swir) = [0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 100.0]
      sza2p3_nodes(1: nsza_swir) = [0.0, 20.0, 40.0, 60.0, 66.0, 72.0, 80.0, 84.0, 88.0]
      vza2p3_nodes(1: nvza_swir) = [0, 12, 18, 26, 32, 36, 40, 46, 50, 54, 56, 60, 66, 72]
      raa2p3_nodes(1: nraa_swir) = [0, 30, 60, 90, 120, 150, 160, 165, 170, 175, 180]

      group_name = 'ai_mie'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'MieCloud354'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start6  = (/ 1,1,1,1,1,1 /)
      edge6   = SHAPE(MieRad2p3)
      stride6 = (/ 1,1,1,1,1,1 /)
      status = nf90_get_var(grp_id, dset_id, MieRad2p3, start=start6, &
         stride=stride6, count=edge6)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'MieCloud388'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge6   = SHAPE(MieRad388)
      status = nf90_get_var(grp_id, dset_id, MieRad388, start=start6, &
         stride=stride6, count=edge6)
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
      ALLOCATE(MieRad_lin2p3_ai(nplev_swir*nsalb_swir*ncod_swir*nsza_swir*nvza_swir*nraa_swir),&
         MieRad_lin388_ai(nplev_swir*nsalb_swir*ncod_swir*nsza_swir*nvza_swir*nraa_swir) )
      !
      MieRad_lin2p3_ai(:) = RESHAPE(MieRad2p3, (/nplev_swir*nsalb_swir*ncod_swir*nsza_swir*nvza_swir*nraa_swir/) )
      MieRad_lin388_ai(:) = RESHAPE(MieRad388, (/nplev_swir*nsalb_swir*ncod_swir*nsza_swir*nvza_swir*nraa_swir/) )
      !

   END SUBROUTINE Read_SWIRLUTAIparams


!
!==================================================================
!==================================================================
!
END MODULE GetLUT_Miecloud_AI_H5module_nc4
