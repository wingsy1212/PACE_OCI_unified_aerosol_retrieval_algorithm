 MODULE Get_omacaLUT7dim_H5module 
!==============================================================================
!
! FILENAME:
!     Get_omacaLUT7dim_module.f90
!
! DESCRIPTION:
!     This module read the Above Cloud Aerosol look-up table in he4
!
! AUTHORS:
!     Changwoo Ahn / Science Systems and Applications, Inc.
!
! HISTORY: Jul 10, 2014
!==============================================================================
 USE InterpolationModule
 use netcdf
 USE OCIUAAER_Config_Module
 
 IMPLICIT NONE

 REAL(KIND=4) :: sflux_lutaac
 DATA sflux_lutaac /3.1415926/

 INTEGER(KIND = 4), PARAMETER :: nwav_aac = 2
 INTEGER(KIND = 4), PARAMETER :: naod_aac = 7
 INTEGER(KIND = 4), PARAMETER :: ncod_aac = 8
 INTEGER(KIND = 4), PARAMETER :: nsza_aac = 7
 INTEGER(KIND = 4), PARAMETER :: nraa_aac = 11
 INTEGER(KIND = 4), PARAMETER :: nvza_aac = 14
 INTEGER(KIND = 4), PARAMETER :: nssa_aac = 7
 INTEGER(KIND = 4), PARAMETER :: nsalb_aac = 5

    REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: wavetbl_aac, aodtbl_aac, &
                     codtbl_aac, szatbl_aac, &
                     vzatbl_aac, raatbl_aac, &
		     ssatbl_aac, salbtbl_aac

!  Smoke over land & ocean variables...
    REAL (KIND=4), DIMENSION(:), ALLOCATABLE :: rad388_smokelin1013zhgt3_aac,rad388_smokelin1013zhgt4_aac,&
                                                rad388_smokelin1013zhgt5_aac,rad388_smokelin1013zhgt6_aac,&
                                                rad388_smokelin800zhgt5_aac,rad388_smokelin800zhgt6_aac,&
                                                rad388_smokelin800zhgt7_aac,rad388_smokelin800zhgt8_aac,&
                                                uvaimie_smokelin1013zhgt3_aac,uvaimie_smokelin1013zhgt4_aac,&
                                                uvaimie_smokelin1013zhgt5_aac,uvaimie_smokelin1013zhgt6_aac,&
                                                uvaimie_smokelin800zhgt5_aac,uvaimie_smokelin800zhgt6_aac,&
                                                uvaimie_smokelin800zhgt7_aac,uvaimie_smokelin800zhgt8_aac

!  Dust over land variables...
   REAL (KIND=4), DIMENSION(:), ALLOCATABLE ::  rad388_dustlin1013zhgt3_aac,rad388_dustlin1013zhgt4_aac,&
                                                rad388_dustlin1013zhgt5_aac,rad388_dustlin1013zhgt6_aac,&
                                                rad388_dustlin800zhgt5_aac,rad388_dustlin800zhgt6_aac,&
                                                rad388_dustlin800zhgt7_aac,rad388_dustlin800zhgt8_aac,&
                                                uvaimie_dustlin1013zhgt3_aac,uvaimie_dustlin1013zhgt4_aac,&
                                                uvaimie_dustlin1013zhgt5_aac,uvaimie_dustlin1013zhgt6_aac,&
                                                uvaimie_dustlin800zhgt5_aac,uvaimie_dustlin800zhgt6_aac,&
                                                uvaimie_dustlin800zhgt7_aac,uvaimie_dustlin800zhgt8_aac

!  Dust over ocean variables...
   REAL (KIND=4), DIMENSION(:), ALLOCATABLE ::  rad388_dustolin1013zhgt3_aac,rad388_dustolin1013zhgt4_aac,&
                                                rad388_dustolin1013zhgt5_aac,rad388_dustolin1013zhgt6_aac,&
                                                rad388_dustolin800zhgt5_aac,rad388_dustolin800zhgt6_aac,&
                                                rad388_dustolin800zhgt7_aac,rad388_dustolin800zhgt8_aac,&
                                                uvaimie_dustolin1013zhgt3_aac,uvaimie_dustolin1013zhgt4_aac,&
                                                uvaimie_dustolin1013zhgt5_aac,uvaimie_dustolin1013zhgt6_aac,&
                                                uvaimie_dustolin800zhgt5_aac,uvaimie_dustolin800zhgt6_aac,&
                                                uvaimie_dustolin800zhgt7_aac,uvaimie_dustolin800zhgt8_aac

!
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: fint_rad388_aac, fint_uvaimie_aac


 CONTAINS

!
!

 SUBROUTINE Read_aac_LUTparams(lut_fn, &
 			       rad388_smokelin1013zhgt3_aac,rad388_smokelin1013zhgt4_aac,&
                               rad388_smokelin1013zhgt5_aac,rad388_smokelin1013zhgt6_aac,&
                               rad388_smokelin800zhgt5_aac,rad388_smokelin800zhgt6_aac,&
                               rad388_smokelin800zhgt7_aac,rad388_smokelin800zhgt8_aac,&
                               uvaimie_smokelin1013zhgt3_aac,uvaimie_smokelin1013zhgt4_aac,&
                               uvaimie_smokelin1013zhgt5_aac,uvaimie_smokelin1013zhgt6_aac,&
                               uvaimie_smokelin800zhgt5_aac,uvaimie_smokelin800zhgt6_aac,&
                               uvaimie_smokelin800zhgt7_aac,uvaimie_smokelin800zhgt8_aac,&

!	                       ;DUST Over-land LUT Parameters...
                               rad388_dustlin1013zhgt3_aac,rad388_dustlin1013zhgt4_aac,&
                               rad388_dustlin1013zhgt5_aac,rad388_dustlin1013zhgt6_aac,&
                               rad388_dustlin800zhgt5_aac,rad388_dustlin800zhgt6_aac,&
                               rad388_dustlin800zhgt7_aac,rad388_dustlin800zhgt8_aac,&
                               uvaimie_dustlin1013zhgt3_aac,uvaimie_dustlin1013zhgt4_aac,&
                               uvaimie_dustlin1013zhgt5_aac,uvaimie_dustlin1013zhgt6_aac,&
                               uvaimie_dustlin800zhgt5_aac,uvaimie_dustlin800zhgt6_aac,&
                               uvaimie_dustlin800zhgt7_aac,uvaimie_dustlin800zhgt8_aac,&

!	                       ;DUST Over-ocean LUT Parameters...
                               rad388_dustolin1013zhgt3_aac,rad388_dustolin1013zhgt4_aac,&
                               rad388_dustolin1013zhgt5_aac,rad388_dustolin1013zhgt6_aac,&
                               rad388_dustolin800zhgt5_aac,rad388_dustolin800zhgt6_aac,&
                               rad388_dustolin800zhgt7_aac,rad388_dustolin800zhgt8_aac,&
                               uvaimie_dustolin1013zhgt3_aac,uvaimie_dustolin1013zhgt4_aac,&
                               uvaimie_dustolin1013zhgt5_aac,uvaimie_dustolin1013zhgt6_aac,&
                               uvaimie_dustolin800zhgt5_aac,uvaimie_dustolin800zhgt6_aac,&
                               uvaimie_dustolin800zhgt7_aac,uvaimie_dustolin800zhgt8_aac,&
                               wavetbl_aac, aodtbl_aac, codtbl_aac, szatbl_aac,vzatbl_aac, &
                               raatbl_aac,ssatbl_aac, salbtbl_aac)
      

 IMPLICIT NONE

!  INTEGER(KIND=4) :: nwav_aac,naod_aac,ncod_aac,nsza_aac,nraa_aac,nvza_aac,nssa_aac, nsalb_aac
  INTEGER (KIND=4) :: linIndex_aac(2), ialo, iz, id

! File, group, dataset and attribute names.
  CHARACTER(LEN=*), INTENT(IN) :: lut_fn
! HSIZE_T type integer.
  INTEGER :: wavdims(nwav_aac), aoddims(naod_aac), coddims(ncod_aac)
  INTEGER :: szadims(nsza_aac), vzadims(nvza_aac), raadims(nraa_aac)
  INTEGER :: ssadims(nssa_aac), salbdims(nsalb_aac)
  INTEGER, DIMENSION(7) :: dims7

  
! Regular four-byte integer.
  INTEGER(KIND=4) :: hdf_err, status
  
!
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: rad388_smokelin1013zhgt3_aac,rad388_smokelin1013zhgt4_aac,&
                                                             rad388_smokelin1013zhgt5_aac,rad388_smokelin1013zhgt6_aac,&
                                                             rad388_smokelin800zhgt5_aac,rad388_smokelin800zhgt6_aac,&
                                                             rad388_smokelin800zhgt7_aac,rad388_smokelin800zhgt8_aac,&
                                                             uvaimie_smokelin1013zhgt3_aac,uvaimie_smokelin1013zhgt4_aac,&
                                                             uvaimie_smokelin1013zhgt5_aac,uvaimie_smokelin1013zhgt6_aac,&
                                                             uvaimie_smokelin800zhgt5_aac,uvaimie_smokelin800zhgt6_aac,&
                                                             uvaimie_smokelin800zhgt7_aac,uvaimie_smokelin800zhgt8_aac

!
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: rad388_dustlin1013zhgt3_aac,rad388_dustlin1013zhgt4_aac,&
                                                             rad388_dustlin1013zhgt5_aac,rad388_dustlin1013zhgt6_aac,&
                                                             rad388_dustlin800zhgt5_aac,rad388_dustlin800zhgt6_aac,&
                                                             rad388_dustlin800zhgt7_aac,rad388_dustlin800zhgt8_aac,&
                                                             uvaimie_dustlin1013zhgt3_aac,uvaimie_dustlin1013zhgt4_aac,&
                                                             uvaimie_dustlin1013zhgt5_aac,uvaimie_dustlin1013zhgt6_aac,&
                                                             uvaimie_dustlin800zhgt5_aac,uvaimie_dustlin800zhgt6_aac,&
                                                             uvaimie_dustlin800zhgt7_aac,uvaimie_dustlin800zhgt8_aac

!
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: rad388_dustolin1013zhgt3_aac,rad388_dustolin1013zhgt4_aac,&
                                                             rad388_dustolin1013zhgt5_aac,rad388_dustolin1013zhgt6_aac,&
                                                             rad388_dustolin800zhgt5_aac,rad388_dustolin800zhgt6_aac,&
                                                             rad388_dustolin800zhgt7_aac,rad388_dustolin800zhgt8_aac,&
                                                             uvaimie_dustolin1013zhgt3_aac,uvaimie_dustolin1013zhgt4_aac,&
                                                             uvaimie_dustolin1013zhgt5_aac,uvaimie_dustolin1013zhgt6_aac,&
                                                             uvaimie_dustolin800zhgt5_aac,uvaimie_dustolin800zhgt6_aac,&
                                                             uvaimie_dustolin800zhgt7_aac,uvaimie_dustolin800zhgt8_aac

    REAL(KIND=4), DIMENSION(:), ALLOCATABLE,  INTENT(OUT) :: wavetbl_aac, aodtbl_aac,&
                                                             codtbl_aac, szatbl_aac, &
                                                             vzatbl_aac, raatbl_aac, ssatbl_aac, salbtbl_aac
!  For SMOKE Over Land & Ocean
    REAL(KIND=4), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE  :: rad388_1013zhgt3_aac_smk, rad388_1013zhgt4_aac_smk, &
    							    rad388_1013zhgt5_aac_smk, rad388_1013zhgt6_aac_smk, &
                                                            rad388_800zhgt5_aac_smk, rad388_800zhgt6_aac_smk, &
							    rad388_800zhgt7_aac_smk, rad388_800zhgt8_aac_smk, &
                                                            uvaimie_1013zhgt3_aac_smk, uvaimie_1013zhgt4_aac_smk, &
							    uvaimie_1013zhgt5_aac_smk, uvaimie_1013zhgt6_aac_smk, &
                                                            uvaimie_800zhgt5_aac_smk, uvaimie_800zhgt6_aac_smk, &
							    uvaimie_800zhgt7_aac_smk, uvaimie_800zhgt8_aac_smk

              
!  For DUST Over Land
   REAL(KIND=4), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE  ::  rad388_1013zhgt3_aac_dst, rad388_1013zhgt4_aac_dst, &
   							    rad388_1013zhgt5_aac_dst, rad388_1013zhgt6_aac_dst, &
                                                            rad388_800zhgt5_aac_dst, rad388_800zhgt6_aac_dst, &
							    rad388_800zhgt7_aac_dst, rad388_800zhgt8_aac_dst, &
                                                            uvaimie_1013zhgt3_aac_dst, uvaimie_1013zhgt4_aac_dst, &
							    uvaimie_1013zhgt5_aac_dst, uvaimie_1013zhgt6_aac_dst, &
                                                            uvaimie_800zhgt5_aac_dst, uvaimie_800zhgt6_aac_dst, &
							    uvaimie_800zhgt7_aac_dst, uvaimie_800zhgt8_aac_dst

  
! For DUST Over Ocean
   REAL(KIND=4), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE  ::  rad388_1013zhgt3_aac_dsto, rad388_1013zhgt4_aac_dsto, &
   							    rad388_1013zhgt5_aac_dsto, rad388_1013zhgt6_aac_dsto, &
                                                            rad388_800zhgt5_aac_dsto, rad388_800zhgt6_aac_dsto, &
							    rad388_800zhgt7_aac_dsto, rad388_800zhgt8_aac_dsto, &
                                                            uvaimie_1013zhgt3_aac_dsto, uvaimie_1013zhgt4_aac_dsto, &
							    uvaimie_1013zhgt5_aac_dsto, uvaimie_1013zhgt6_aac_dsto, &
                                                            uvaimie_800zhgt5_aac_dsto, uvaimie_800zhgt6_aac_dsto, &
							    uvaimie_800zhgt7_aac_dsto, uvaimie_800zhgt8_aac_dsto
 
!--------------------------------------------------------------------------------------------------------------------------------------------------
      integer, dimension (1) :: start1, edge1, stride1
      integer, dimension (7) :: start7, edge7, stride7

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

      group_name = 'omacalut'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

     ALLOCATE( wavetbl_aac(nwav_aac), aodtbl_aac(naod_aac), &
                codtbl_aac(ncod_aac), szatbl_aac(nsza_aac), &
                vzatbl_aac(nvza_aac), raatbl_aac(nraa_aac), &
                ssatbl_aac(nssa_aac),salbtbl_aac(nsalb_aac) )

      dset_name = 'NumberOfWavelength'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start1  = (/ 1 /)
      edge1   = (/ nwav_aac /)
      stride1 = (/ 1 /)
      status = nf90_get_var(grp_id, dset_id, wavetbl_aac, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'NumberOfAerosolOpticalDepth'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = (/ naod_aac /)
      status = nf90_get_var(grp_id, dset_id, aodtbl_aac, start=start1, &
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
      edge1   = (/ ncod_aac /)
      status = nf90_get_var(grp_id, dset_id, codtbl_aac, start=start1, &
         stride=stride1, count=edge1)
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
      edge1   = (/ nsza_aac /)
      status = nf90_get_var(grp_id, dset_id, szatbl_aac, start=start1, &
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
      edge1   = (/ nvza_aac /)
      status = nf90_get_var(grp_id, dset_id, vzatbl_aac, start=start1, &
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
      edge1   = (/ nraa_aac /)
      status = nf90_get_var(grp_id, dset_id, raatbl_aac, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'NumberOfSingleScatteringAlbedoDust'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = (/ nssa_aac /)
      status = nf90_get_var(grp_id, dset_id, ssatbl_aac, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'NumberOfSurfaceAlbedo'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge1   = (/ nsalb_aac /)
      status = nf90_get_var(grp_id, dset_id, salbtbl_aac, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

    !! Allocate memory

!   FOR SMOKE
    ALLOCATE(rad388_1013zhgt3_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt4_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt5_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt6_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt3_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt4_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt5_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt6_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)

    ALLOCATE(rad388_800zhgt5_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt6_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt7_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt8_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt5_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt6_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt7_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt8_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)

! FOR DUST OVER LAND

    ALLOCATE(rad388_1013zhgt3_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt4_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt5_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt6_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt3_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt4_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt5_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt6_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)

    ALLOCATE(rad388_800zhgt5_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt6_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt7_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt8_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt5_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt6_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt7_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt8_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)


! FOR DUST OVER OCEAN
    ALLOCATE(rad388_1013zhgt3_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt4_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt5_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt6_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt3_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt4_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt5_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt6_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)


    ALLOCATE(rad388_800zhgt5_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt6_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt7_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt8_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt5_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt6_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt7_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt8_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)

    ALLOCATE(wavetbl_aac(nwav_aac), aodtbl_aac(naod_aac), codtbl_aac(ncod_aac),  szatbl_aac(nsza_aac), &
              vzatbl_aac(nvza_aac), raatbl_aac(nraa_aac), ssatbl_aac(nssa_aac),salbtbl_aac(nsalb_aac), &
	      STAT = status)

  
! -------------------------------------------------------------------------------------------------------------------------
! FOR SMOKE
!
!
   dset_name = 'Radiance388_SMOKE_PRESLEV1013_AERHGT3p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start7  = (/ 1,1,1,1,1,1,1 /)
   edge7   = SHAPE(rad388_1013zhgt3_aac_smk)
   stride7 = (/ 1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt3_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt3_aac_smk = rad388_1013zhgt3_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV1013_AERHGT3p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt3_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt3_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt3_aac_smk = uvaimie_1013zhgt3_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV1013_AERHGT4p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt4_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt4_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt4_aac_smk = rad388_1013zhgt4_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV1013_AERHGT4p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt4_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt4_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt4_aac_smk = uvaimie_1013zhgt4_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV1013_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt5_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt5_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt5_aac_smk = rad388_1013zhgt5_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV1013_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt5_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt5_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt5_aac_smk = uvaimie_1013zhgt5_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV1013_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt6_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt6_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt6_aac_smk = rad388_1013zhgt6_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV1013_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt6_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt6_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt6_aac_smk = uvaimie_1013zhgt6_aac_smk

!
! --- Read 800 hPa files ----------------------------------------------------------------------
!
   dset_name = 'Radiance388_SMOKE_PRESLEV0800_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt5_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt5_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt5_aac_smk = rad388_800zhgt5_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV0800_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt5_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt5_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt5_aac_smk = uvaimie_800zhgt5_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV0800_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt6_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt6_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt6_aac_smk = rad388_800zhgt6_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV0800_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt6_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt6_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt6_aac_smk = uvaimie_800zhgt6_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV0800_AERHGT7p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt7_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt7_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt7_aac_smk = rad388_800zhgt7_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV0800_AERHGT7p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt7_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt7_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt7_aac_smk = uvaimie_800zhgt7_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV0800_AERHGT8p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt8_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt8_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt8_aac_smk = rad388_800zhgt8_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV0800_AERHGT8p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt8_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt8_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt8_aac_smk = uvaimie_800zhgt8_aac_smk

!

! -------------------------------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------------------------------
! FOR DUST
!

   dset_name = 'Radiance388_DUST_PRESLEV1013_AERHGT3p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt3_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt3_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt3_aac_dst = rad388_1013zhgt3_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV1013_AERHGT3p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt3_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt3_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt3_aac_dst = uvaimie_1013zhgt3_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV1013_AERHGT4p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt4_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt4_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt4_aac_dst = rad388_1013zhgt4_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV1013_AERHGT4p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt4_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt4_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt4_aac_dst = uvaimie_1013zhgt4_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV1013_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt5_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt5_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt5_aac_dst = rad388_1013zhgt5_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV1013_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt5_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt5_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt5_aac_dst = uvaimie_1013zhgt5_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV1013_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt6_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt6_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt6_aac_dst = rad388_1013zhgt6_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV1013_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt6_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt6_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt6_aac_dst = uvaimie_1013zhgt6_aac_dst

!
! --- Read 800 hPa files ----------------------------------------------------------------------
!
   dset_name = 'Radiance388_DUST_PRESLEV0800_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt5_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt5_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt5_aac_dst = rad388_800zhgt5_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV0800_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt5_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt5_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt5_aac_dst = uvaimie_800zhgt5_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV0800_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt6_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt6_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt6_aac_dst = rad388_800zhgt6_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV0800_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt6_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt6_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt6_aac_dst = uvaimie_800zhgt6_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV0800_AERHGT7p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt7_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt7_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt7_aac_dst = rad388_800zhgt7_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV0800_AERHGT7p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt7_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt7_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt7_aac_dst = uvaimie_800zhgt7_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV0800_AERHGT8p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt8_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt8_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt8_aac_dst = rad388_800zhgt8_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV0800_AERHGT8p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt8_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt8_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt8_aac_dst = uvaimie_800zhgt8_aac_dst
 

! -------------------------------------------------------------------------------------------------------------------------
! FOR DUSTO OVER OCEAN
!

   dset_name = 'Radiance388_DUSTO_PRESLEV1013_AERHGT3p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt3_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt3_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt3_aac_dsto = rad388_1013zhgt3_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV1013_AERHGT3p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt3_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt3_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt3_aac_dsto = uvaimie_1013zhgt3_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV1013_AERHGT4p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt4_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt4_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt4_aac_dsto = rad388_1013zhgt4_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV1013_AERHGT4p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt4_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt4_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt4_aac_dsto = uvaimie_1013zhgt4_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV1013_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt5_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt5_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt5_aac_dsto = rad388_1013zhgt5_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV1013_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt5_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt5_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt5_aac_dsto = uvaimie_1013zhgt5_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV1013_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt6_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt6_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt6_aac_dsto = rad388_1013zhgt6_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV1013_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt6_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt6_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt6_aac_dsto = uvaimie_1013zhgt6_aac_dsto

!
! --- Read 800 hPa files ----------------------------------------------------------------------
!
   dset_name = 'Radiance388_DUSTO_PRESLEV0800_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt5_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt5_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt5_aac_dsto = rad388_800zhgt5_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV0800_AERHGT5p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt5_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt5_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt5_aac_dsto = uvaimie_800zhgt5_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV0800_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt6_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt6_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt6_aac_dsto = rad388_800zhgt6_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV0800_AERHGT6p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt6_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt6_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt6_aac_dsto = uvaimie_800zhgt6_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV0800_AERHGT7p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt7_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt7_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt7_aac_dsto = rad388_800zhgt7_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV0800_AERHGT7p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt7_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt7_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt7_aac_dsto = uvaimie_800zhgt7_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV0800_AERHGT8p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt8_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt8_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt8_aac_dsto = rad388_800zhgt8_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV0800_AERHGT8p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt8_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt8_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt8_aac_dsto = uvaimie_800zhgt8_aac_dsto
 
    
!
! -- Do Lagrange interpolation for Radiance fields ----------------------
  ALLOCATE(rad388_smokelin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_smokelin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_smokelin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_smokelin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)
   !
  ALLOCATE(rad388_smokelin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_smokelin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_smokelin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_smokelin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)

! DUST Over-land Parameters...
  ALLOCATE(rad388_dustlin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustlin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustlin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustlin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)
   !
  ALLOCATE(rad388_dustlin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustlin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustlin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustlin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)

! DUST Over-ocean Parameters...
  ALLOCATE(rad388_dustolin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustolin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustolin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustolin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)
!
  ALLOCATE(rad388_dustolin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustolin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustolin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           rad388_dustolin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)


  rad388_smokelin1013zhgt3_aac(:)  = RESHAPE(rad388_1013zhgt3_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin1013zhgt4_aac(:)  = RESHAPE(rad388_1013zhgt4_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin1013zhgt5_aac(:)  = RESHAPE(rad388_1013zhgt5_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin1013zhgt6_aac(:)  = RESHAPE(rad388_1013zhgt6_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt3_aac(:)  = RESHAPE(uvaimie_1013zhgt3_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt4_aac(:)  = RESHAPE(uvaimie_1013zhgt4_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt5_aac(:)  = RESHAPE(uvaimie_1013zhgt5_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt6_aac(:)  = RESHAPE(uvaimie_1013zhgt6_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
!
  rad388_smokelin800zhgt5_aac(:)  = RESHAPE(rad388_800zhgt5_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin800zhgt6_aac(:)  = RESHAPE(rad388_800zhgt6_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin800zhgt7_aac(:)  = RESHAPE(rad388_800zhgt7_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin800zhgt8_aac(:)  = RESHAPE(rad388_800zhgt8_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt5_aac(:)  = RESHAPE(uvaimie_800zhgt5_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt6_aac(:)  = RESHAPE(uvaimie_800zhgt6_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt7_aac(:)  = RESHAPE(uvaimie_800zhgt7_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt8_aac(:)  = RESHAPE(uvaimie_800zhgt8_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )

! DUST Over-land Parameters...
  rad388_dustlin1013zhgt3_aac(:)  = RESHAPE(rad388_1013zhgt3_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin1013zhgt4_aac(:)  = RESHAPE(rad388_1013zhgt4_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin1013zhgt5_aac(:)  = RESHAPE(rad388_1013zhgt5_aac_dst, &
 	 				(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin1013zhgt6_aac(:)  = RESHAPE(rad388_1013zhgt6_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt3_aac(:)  = RESHAPE(uvaimie_1013zhgt3_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt4_aac(:)  = RESHAPE(uvaimie_1013zhgt4_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt5_aac(:)  = RESHAPE(uvaimie_1013zhgt5_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt6_aac(:)  = RESHAPE(uvaimie_1013zhgt6_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
!
  rad388_dustlin800zhgt5_aac(:)  = RESHAPE(rad388_800zhgt5_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin800zhgt6_aac(:)  = RESHAPE(rad388_800zhgt6_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin800zhgt7_aac(:)  = RESHAPE(rad388_800zhgt7_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin800zhgt8_aac(:)  = RESHAPE(rad388_800zhgt8_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt5_aac(:)  = RESHAPE(uvaimie_800zhgt5_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt6_aac(:)  = RESHAPE(uvaimie_800zhgt6_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt7_aac(:)  = RESHAPE(uvaimie_800zhgt7_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt8_aac(:)  = RESHAPE(uvaimie_800zhgt8_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )


! DUST Over-ocean Parameters...
  rad388_dustolin1013zhgt3_aac(:)  = RESHAPE(rad388_1013zhgt3_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin1013zhgt4_aac(:)  = RESHAPE(rad388_1013zhgt4_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin1013zhgt5_aac(:)  = RESHAPE(rad388_1013zhgt5_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin1013zhgt6_aac(:)  = RESHAPE(rad388_1013zhgt6_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt3_aac(:)  = RESHAPE(uvaimie_1013zhgt3_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt4_aac(:)  = RESHAPE(uvaimie_1013zhgt4_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt5_aac(:)  = RESHAPE(uvaimie_1013zhgt5_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt6_aac(:)  = RESHAPE(uvaimie_1013zhgt6_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
!
  rad388_dustolin800zhgt5_aac(:)  = RESHAPE(rad388_800zhgt5_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin800zhgt6_aac(:)  = RESHAPE(rad388_800zhgt6_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin800zhgt7_aac(:)  = RESHAPE(rad388_800zhgt7_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin800zhgt8_aac(:)  = RESHAPE(rad388_800zhgt8_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt5_aac(:)  = RESHAPE(uvaimie_800zhgt5_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt6_aac(:)  = RESHAPE(uvaimie_800zhgt6_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt7_aac(:)  = RESHAPE(uvaimie_800zhgt7_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt8_aac(:)  = RESHAPE(uvaimie_800zhgt8_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )

 DEALLOCATE(rad388_1013zhgt3_aac_smk, rad388_1013zhgt4_aac_smk, rad388_1013zhgt5_aac_smk, rad388_1013zhgt6_aac_smk,&
            rad388_800zhgt5_aac_smk, rad388_800zhgt6_aac_smk, rad388_800zhgt7_aac_smk, rad388_800zhgt8_aac_smk,&
            uvaimie_1013zhgt3_aac_smk, uvaimie_1013zhgt4_aac_smk, uvaimie_1013zhgt5_aac_smk, uvaimie_1013zhgt6_aac_smk,&
            uvaimie_800zhgt5_aac_smk, uvaimie_800zhgt6_aac_smk, uvaimie_800zhgt7_aac_smk, uvaimie_800zhgt8_aac_smk, STAT=status)

 DEALLOCATE(rad388_1013zhgt3_aac_dst, rad388_1013zhgt4_aac_dst, rad388_1013zhgt5_aac_dst, rad388_1013zhgt6_aac_dst,&
            rad388_800zhgt5_aac_dst, rad388_800zhgt6_aac_dst, rad388_800zhgt7_aac_dst, rad388_800zhgt8_aac_dst,&
            uvaimie_1013zhgt3_aac_dst, uvaimie_1013zhgt4_aac_dst, uvaimie_1013zhgt5_aac_dst, uvaimie_1013zhgt6_aac_dst,&
            uvaimie_800zhgt5_aac_dst, uvaimie_800zhgt6_aac_dst, uvaimie_800zhgt7_aac_dst, uvaimie_800zhgt8_aac_dst, STAT=status)

 DEALLOCATE(rad388_1013zhgt3_aac_dsto, rad388_1013zhgt4_aac_dsto, rad388_1013zhgt5_aac_dsto, rad388_1013zhgt6_aac_dsto,&
            rad388_800zhgt5_aac_dsto, rad388_800zhgt6_aac_dsto, rad388_800zhgt7_aac_dsto, rad388_800zhgt8_aac_dsto,&
            uvaimie_1013zhgt3_aac_dsto, uvaimie_1013zhgt4_aac_dsto, uvaimie_1013zhgt5_aac_dsto, uvaimie_1013zhgt6_aac_dsto,&
            uvaimie_800zhgt5_aac_dsto, uvaimie_800zhgt6_aac_dsto, uvaimie_800zhgt7_aac_dsto, uvaimie_800zhgt8_aac_dsto, STAT=status)

!

90   status = nf90_close(nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to close lut_nc4 file: ", status
      return
   end if
!
 END SUBROUTINE Read_aac_LUTparams


 SUBROUTINE Interpol_aac_LUTparams(ocean,atype,inpterr, inzhgt, inssa, insalb354, insalb388, &
                                   fint_rad388_aac, fint_uvaimie_aac)
!
   USE LookupTableModule

   IMPLICIT NONE

   INTEGER(KIND=4) :: status, version, ocean
   INTEGER(KIND=2), INTENT(IN)  :: atype
   REAL(KIND=4),    INTENT(IN)  ::  inpterr, inzhgt, inssa, insalb354, insalb388


!
   REAL(KIND=4)     ::  inpterr2
!
   REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: fint_rad388_aac, fint_uvaimie_aac
!
   REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE        :: rad_smokeout_aac, rad_dustout_aac, rad_dustoout_aac 
   REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE        :: uvaimie_smokeout_aac, uvaimie_dustout_aac, uvaimie_dustoout_aac
!
   REAL(KIND=4), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: intrad388_1013_aac, intrad388_800_aac,&
                                                           intuvaimie_1013_aac, intuvaimie_800_aac
!
  REAL(KIND=4), DIMENSION(:,:,:,:), ALLOCATABLE      :: int2_rad388_1013_aac, int2_rad388_800_aac, &
                                                        int2_uvaimie_1013_aac, int2_uvaimie_800_aac
!
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE      :: int3_rad388_1013_aac, int3_rad388_800_aac, &
                                                      int3_uvaimie_1013_aac, int3_uvaimie_800_aac

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE      :: int4_rad388_1013_aac, int4_rad388_800_aac, &
                                                    int4_uvaimie_1013_aac, int4_uvaimie_800_aac


!  Bounding values in table for input parameters
!==============================================================================
  REAL(KIND=4)    :: theta1,theta2,sza1,sza2,phi1,phi2,w1,w2
  REAL(KIND=4)    :: fractionalZhgt, fractionalSSA, meanfrac, wt
  REAL(KIND=4)    :: fractionalSALB354, fractionalSALB388

!==============================================================================
!  Indices for bounding values in table for input parameters
!==============================================================================
  INTEGER(KIND=4) :: jtheta, jtheta2, jsza, jsza2, jphi, jphi2, jw1

  REAL(KIND=4), DIMENSION(7) :: ssa388val_smk = (/0.780184, 0.807809, 0.845545, 0.887624, 0.935943, 0.964274, 1.000/)
  REAL(KIND=4), DIMENSION(7) :: ssa388val_dst = (/0.765130, 0.828882, 0.870567, 0.902707, 0.948519, 0.972289, 1.000/)
  REAL(KIND=4), DIMENSION(4) :: zhgt1013val = (/3.0, 4.0, 5.0, 6.0/)
  REAL(KIND=4), DIMENSION(5) :: salbval = (/0.0, 0.05, 0.10, 0.15, 0.20/)

  INTEGER(KIND=4) :: zhgtlw, zhgtup, ssalw, ssaup, salblw354, salbup354, salblw388, salbup388
  INTEGER(KIND=4) ::  ialo, iz

!
   ALLOCATE(intrad388_1013_aac(nsalb_aac,4,nssa_aac,naod_aac,ncod_aac), &
             intrad388_800_aac(nsalb_aac,4,nssa_aac,naod_aac,ncod_aac) )

   ALLOCATE(intuvaimie_1013_aac(nsalb_aac,4,nssa_aac,naod_aac,ncod_aac), &
             intuvaimie_800_aac(nsalb_aac,4,nssa_aac,naod_aac,ncod_aac) )
!
   ALLOCATE(    rad_smokeout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE(     rad_dustout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE(    rad_dustoout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE(uvaimie_smokeout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE( uvaimie_dustout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE(uvaimie_dustoout_aac(nssa_aac,naod_aac,ncod_aac) )



!=============================================================================
! Find the two index values for input zhgt vlaues
!=============================================================================
  status=FindTableEntry(inzhgt, zhgt1013val, 4, theta1,theta2, zhgtlw, zhgtup,&
                                                                 fractionalZhgt)
  if (zhgtlw .ge. 4) zhgtup = zhgtlw
  if (zhgtlw .lt. 1) then
      zhgtlw = 1
      zhgtup = zhgtlw
  endif

!=============================================================================
! Find the two index values for input SSA vlaues
!=============================================================================
  IF(atype.eq.1) THEN
  status=FindTableEntry(inssa, ssa388val_smk, 7, theta1,theta2, ssalw, ssaup,&
                                                           fractionalSSA)
  ENDIF

  IF(atype.eq.2) THEN
  status=FindTableEntry(inssa, ssa388val_dst, 7, theta1,theta2, ssalw, ssaup,&
                                                           fractionalSSA)
  ENDIF

  IF (ssalw .ge. 7) ssaup = ssalw
  IF (ssalw .lt. 1) THEN
      ssalw = 1
      ssaup = ssalw
  ENDIF

!=============================================================================
! Find the two index values for input SALB354 vlaues
!=============================================================================
  status=FindTableEntry(insalb354, salbval, 5, theta1,theta2, salblw354, salbup354,&
                                                                  fractionalSALB354)

!=============================================================================
! Find the two index values for input SALB388 vlaues
!=============================================================================
  status=FindTableEntry(insalb388, salbval, 5, theta1,theta2, salblw388, salbup388,&
                                                                  fractionalSALB388)


  IF(atype.eq.1) THEN
   !
   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin1013zhgt3_aac, rad_smokeout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin1013zhgt4_aac, rad_smokeout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin1013zhgt5_aac, rad_smokeout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin1013zhgt6_aac, rad_smokeout_aac)
     intrad388_1013_aac(ialo,zhgtlw,:,:,:) = rad_smokeout_aac
  
     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin1013zhgt3_aac, rad_smokeout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin1013zhgt4_aac, rad_smokeout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin1013zhgt5_aac, rad_smokeout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin1013zhgt6_aac, rad_smokeout_aac)
     intrad388_1013_aac(ialo,zhgtup,:,:,:) = rad_smokeout_aac
   ENDDO   !DO ialo = salblw388, salbup388
   

   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin1013zhgt3_aac, uvaimie_smokeout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin1013zhgt4_aac, uvaimie_smokeout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin1013zhgt5_aac, uvaimie_smokeout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin1013zhgt6_aac, uvaimie_smokeout_aac)
     intuvaimie_1013_aac(ialo,zhgtlw,:,:,:) = uvaimie_smokeout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin1013zhgt3_aac, uvaimie_smokeout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin1013zhgt4_aac, uvaimie_smokeout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin1013zhgt5_aac, uvaimie_smokeout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin1013zhgt6_aac, uvaimie_smokeout_aac)
     intuvaimie_1013_aac(ialo,zhgtup,:,:,:) = uvaimie_smokeout_aac
   ENDDO   !DO ialo = salblw388, salbup388


 ! -- 800 hPa --------------------------------------------------
   !
   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin800zhgt5_aac, rad_smokeout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin800zhgt6_aac, rad_smokeout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin800zhgt7_aac, rad_smokeout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin800zhgt8_aac, rad_smokeout_aac)
     intrad388_800_aac(ialo,zhgtlw,:,:,:) = rad_smokeout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin800zhgt5_aac, rad_smokeout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin800zhgt6_aac, rad_smokeout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin800zhgt7_aac, rad_smokeout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_smokelin800zhgt8_aac, rad_smokeout_aac)
     intrad388_800_aac(ialo,zhgtup,:,:,:) = rad_smokeout_aac
   ENDDO   !DO ialo = salblw388, salbup388


   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin800zhgt5_aac, uvaimie_smokeout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin800zhgt6_aac, uvaimie_smokeout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin800zhgt7_aac, uvaimie_smokeout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin800zhgt8_aac, uvaimie_smokeout_aac)
     intuvaimie_800_aac(ialo,zhgtlw,:,:,:) = uvaimie_smokeout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin800zhgt5_aac, uvaimie_smokeout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin800zhgt6_aac, uvaimie_smokeout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin800zhgt7_aac, uvaimie_smokeout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_smokelin800zhgt8_aac, uvaimie_smokeout_aac)
     intuvaimie_800_aac(ialo,zhgtup,:,:,:) = uvaimie_smokeout_aac
   ENDDO   !DO ialo = salblw388, salbup388

 ENDIF


 IF(atype.eq.2.and.ocean.eq.0) THEN     !over-land pixel
   !
   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin1013zhgt3_aac, rad_dustout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin1013zhgt4_aac, rad_dustout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin1013zhgt5_aac, rad_dustout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin1013zhgt6_aac, rad_dustout_aac)
     intrad388_1013_aac(ialo,zhgtlw,:,:,:) = rad_dustout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin1013zhgt3_aac, rad_dustout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin1013zhgt4_aac, rad_dustout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin1013zhgt5_aac, rad_dustout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin1013zhgt6_aac, rad_dustout_aac)
     intrad388_1013_aac(ialo,zhgtup,:,:,:) = rad_dustout_aac
   ENDDO   !DO ialo = salblw388, salbup388

   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1)  &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin1013zhgt3_aac, uvaimie_dustout_aac)
     IF(zhgtlw.EQ.2)  &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin1013zhgt4_aac, uvaimie_dustout_aac)
     IF(zhgtlw.EQ.3)  &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin1013zhgt5_aac, uvaimie_dustout_aac)
     IF(zhgtlw.EQ.4)  &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin1013zhgt6_aac, uvaimie_dustout_aac)
     intuvaimie_1013_aac(ialo,zhgtlw,:,:,:) = uvaimie_dustout_aac

     IF(zhgtup.EQ.1)  &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin1013zhgt3_aac, uvaimie_dustout_aac)
     IF(zhgtup.EQ.2)  &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin1013zhgt4_aac, uvaimie_dustout_aac)
     IF(zhgtup.EQ.3)  &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin1013zhgt5_aac, uvaimie_dustout_aac)
     IF(zhgtup.EQ.4)  &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin1013zhgt6_aac, uvaimie_dustout_aac)
     intuvaimie_1013_aac(ialo,zhgtup,:,:,:) = uvaimie_dustout_aac
   ENDDO   !DO ialo = salblw388, salbup388

 ! -- 800 hPa --------------------------------------------------
   !
   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin800zhgt5_aac, rad_dustout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin800zhgt6_aac, rad_dustout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin800zhgt7_aac, rad_dustout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin800zhgt8_aac, rad_dustout_aac)
     intrad388_800_aac(ialo,zhgtlw,:,:,:) = rad_dustout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin800zhgt5_aac, rad_dustout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin800zhgt6_aac, rad_dustout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin800zhgt7_aac, rad_dustout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustlin800zhgt8_aac, rad_dustout_aac)
     intrad388_800_aac(ialo,zhgtup,:,:,:) = rad_dustout_aac
   ENDDO   !DO ialo = salblw388, salbup388

   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin800zhgt5_aac, uvaimie_dustout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin800zhgt6_aac, uvaimie_dustout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin800zhgt7_aac, uvaimie_dustout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin800zhgt8_aac, uvaimie_dustout_aac)
     intuvaimie_800_aac(ialo,zhgtlw,:,:,:) = uvaimie_dustout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin800zhgt5_aac, uvaimie_dustout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin800zhgt6_aac, uvaimie_dustout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin800zhgt7_aac, uvaimie_dustout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustlin800zhgt8_aac, uvaimie_dustout_aac)
     intuvaimie_800_aac(ialo,zhgtup,:,:,:) = uvaimie_dustout_aac
   ENDDO   !DO ialo = salblw388, salbup388
 ENDIF


IF(atype.eq.2.and.ocean.eq.1) THEN     !over-ocean pixel
   !
   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin1013zhgt3_aac, rad_dustoout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin1013zhgt4_aac, rad_dustoout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin1013zhgt5_aac, rad_dustoout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin1013zhgt6_aac, rad_dustoout_aac)
     intrad388_1013_aac(ialo,zhgtlw,:,:,:) = rad_dustoout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin1013zhgt3_aac, rad_dustoout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin1013zhgt4_aac, rad_dustoout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin1013zhgt5_aac, rad_dustoout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin1013zhgt6_aac, rad_dustoout_aac)
     intrad388_1013_aac(ialo,zhgtup,:,:,:) = rad_dustoout_aac
   ENDDO   !DO ialo = salblw388, salbup388

   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin1013zhgt3_aac, uvaimie_dustoout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin1013zhgt4_aac, uvaimie_dustoout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin1013zhgt5_aac, uvaimie_dustoout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin1013zhgt6_aac, uvaimie_dustoout_aac)
     intuvaimie_1013_aac(ialo,zhgtlw,:,:,:) = uvaimie_dustoout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin1013zhgt3_aac, uvaimie_dustoout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin1013zhgt4_aac, uvaimie_dustoout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin1013zhgt5_aac, uvaimie_dustoout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin1013zhgt6_aac, uvaimie_dustoout_aac)
     intuvaimie_1013_aac(ialo,zhgtup,:,:,:) = uvaimie_dustoout_aac
   ENDDO   !DO ialo = salblw388, salbup388


! -- 800 hPa --------------------------------------------------
   !
   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin800zhgt5_aac, rad_dustoout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin800zhgt6_aac, rad_dustoout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin800zhgt7_aac, rad_dustoout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin800zhgt8_aac, rad_dustoout_aac)
     intrad388_800_aac(ialo,zhgtlw,:,:,:) = rad_dustoout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin800zhgt5_aac, rad_dustoout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin800zhgt6_aac, rad_dustoout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin800zhgt7_aac, rad_dustoout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, rad388_dustolin800zhgt8_aac, rad_dustoout_aac)
     intrad388_800_aac(ialo,zhgtup,:,:,:) = rad_dustoout_aac
   ENDDO   !DO ialo = salblw388, salbup388

   DO ialo = salblw388, salbup388
     IF(zhgtlw.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin800zhgt5_aac, uvaimie_dustoout_aac)
     IF(zhgtlw.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin800zhgt6_aac, uvaimie_dustoout_aac)
     IF(zhgtlw.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin800zhgt7_aac, uvaimie_dustoout_aac)
     IF(zhgtlw.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin800zhgt8_aac, uvaimie_dustoout_aac)
     intuvaimie_800_aac(ialo,zhgtlw,:,:,:) = uvaimie_dustoout_aac

     IF(zhgtup.EQ.1) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin800zhgt5_aac, uvaimie_dustoout_aac)
     IF(zhgtup.EQ.2) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin800zhgt6_aac, uvaimie_dustoout_aac)
     IF(zhgtup.EQ.3) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin800zhgt7_aac, uvaimie_dustoout_aac)
     IF(zhgtup.EQ.4) &
     status = InterpRadiance_aac(ialo,1,naod_aac,1,ncod_aac,1,nssa_aac, uvaimie_dustolin800zhgt8_aac, uvaimie_dustoout_aac)
     intuvaimie_800_aac(ialo,zhgtup,:,:,:) = uvaimie_dustoout_aac
   ENDDO   !DO ialo = salblw388, salbup388
 ENDIF

! --------------------------------------------------------------------------------

!!  -- Do linear interpolation on SSA --
 ALLOCATE(int2_rad388_1013_aac(nsalb_aac,4,naod_aac,ncod_aac), &
           int2_rad388_800_aac(nsalb_aac,4,naod_aac,ncod_aac), &
         int2_uvaimie_1013_aac(nsalb_aac,4,naod_aac,ncod_aac), &
	  int2_uvaimie_800_aac(nsalb_aac,4,naod_aac,ncod_aac) )

!
 DO ialo = salblw388, salbup388
    int2_rad388_1013_aac(ialo,zhgtlw,:,:) = intrad388_1013_aac(ialo,zhgtlw,ssalw,:,:)*(1.0 - fractionalSSA) + &
                                            intrad388_1013_aac(ialo,zhgtlw,ssaup,:,:)*fractionalSSA
    int2_rad388_1013_aac(ialo,zhgtup,:,:) = intrad388_1013_aac(ialo,zhgtup,ssalw,:,:)*(1.0 - fractionalSSA) + &
                                            intrad388_1013_aac(ialo,zhgtup,ssaup,:,:)*fractionalSSA
 ENDDO   ! DO ialo = salblw354, salbup354
!
 DO ialo = salblw388, salbup388
    int2_uvaimie_1013_aac(ialo,zhgtlw,:,:) = intuvaimie_1013_aac(ialo,zhgtlw,ssalw,:,:)*(1.0 - fractionalSSA) + &
                                             intuvaimie_1013_aac(ialo,zhgtlw,ssaup,:,:)*fractionalSSA
    int2_uvaimie_1013_aac(ialo,zhgtup,:,:) = intuvaimie_1013_aac(ialo,zhgtup,ssalw,:,:)*(1.0 - fractionalSSA) + &
                                             intuvaimie_1013_aac(ialo,zhgtup,ssaup,:,:)*fractionalSSA
 ENDDO   ! DO ialo = salblw388, salbup388
!
 DO ialo = salblw388, salbup388
    int2_rad388_800_aac(ialo,zhgtlw,:,:) = intrad388_800_aac(ialo,zhgtlw,ssalw,:,:)*(1.0 - fractionalSSA) + &
                                           intrad388_800_aac(ialo,zhgtlw,ssaup,:,:)*fractionalSSA
    int2_rad388_800_aac(ialo,zhgtup,:,:) = intrad388_800_aac(ialo,zhgtup,ssalw,:,:)*(1.0 - fractionalSSA) + &
                                           intrad388_800_aac(ialo,zhgtup,ssaup,:,:)*fractionalSSA
 ENDDO   ! DO ialo = salblw388, salbup388
!
 DO ialo = salblw388, salbup388
    int2_uvaimie_800_aac(ialo,zhgtlw,:,:) = intuvaimie_800_aac(ialo,zhgtlw,ssalw,:,:)*(1.0 - fractionalSSA) + &
                                            intuvaimie_800_aac(ialo,zhgtlw,ssaup,:,:)*fractionalSSA
    int2_uvaimie_800_aac(ialo,zhgtup,:,:) = intuvaimie_800_aac(ialo,zhgtup,ssalw,:,:)*(1.0 - fractionalSSA) + &
                                            intuvaimie_800_aac(ialo,zhgtup,ssaup,:,:)*fractionalSSA
 ENDDO   ! DO ialo = salblw388, salbup388


!!  -- Do linear interpolation on ZHGT --
 ALLOCATE(int3_rad388_1013_aac(nsalb_aac,naod_aac,ncod_aac), &
           int3_rad388_800_aac(nsalb_aac,naod_aac,ncod_aac), &
         int3_uvaimie_1013_aac(nsalb_aac,naod_aac,ncod_aac), &
	  int3_uvaimie_800_aac(nsalb_aac,naod_aac,ncod_aac) )

 DO ialo = salblw388, salbup388
     int3_rad388_1013_aac(ialo,:,:) =  int2_rad388_1013_aac(ialo,zhgtlw,:,:)*(1.0 - fractionalZhgt) + &
                                       int2_rad388_1013_aac(ialo,zhgtup,:,:)*fractionalZhgt
      int3_rad388_800_aac(ialo,:,:) =  int2_rad388_800_aac(ialo,zhgtlw,:,:)*(1.0 - fractionalZhgt) + &
                                       int2_rad388_800_aac(ialo,zhgtup,:,:)*fractionalZhgt
    int3_uvaimie_1013_aac(ialo,:,:) =  int2_uvaimie_1013_aac(ialo,zhgtlw,:,:)*(1.0 - fractionalZhgt) + &
                                       int2_uvaimie_1013_aac(ialo,zhgtup,:,:)*fractionalZhgt
     int3_uvaimie_800_aac(ialo,:,:) =  int2_uvaimie_800_aac(ialo,zhgtlw,:,:)*(1.0 - fractionalZhgt) + &
                                       int2_uvaimie_800_aac(ialo,zhgtup,:,:)*fractionalZhgt
 ENDDO

!
 ALLOCATE(int4_rad388_1013_aac(naod_aac,ncod_aac), &
           int4_rad388_800_aac(naod_aac,ncod_aac), &
         int4_uvaimie_1013_aac(naod_aac,ncod_aac), &
	  int4_uvaimie_800_aac(naod_aac,ncod_aac) )
!
  int4_rad388_1013_aac(:,:) =  int3_rad388_1013_aac(salblw388,:,:)*(1.0 - fractionalSALB388) + &
                               int3_rad388_1013_aac(salbup388,:,:)*fractionalSALB388
  int4_uvaimie_1013_aac(:,:) = int3_uvaimie_1013_aac(salblw388,:,:)*(1.0 - fractionalSALB388) + &
                               int3_uvaimie_1013_aac(salbup388,:,:)*fractionalSALB388
  !
  int4_rad388_800_aac(:,:) =  int3_rad388_800_aac(salblw388,:,:)*(1.0 - fractionalSALB388) + &
                              int3_rad388_800_aac(salbup388,:,:)*fractionalSALB388
  int4_uvaimie_800_aac(:,:) = int3_uvaimie_800_aac(salblw388,:,:)*(1.0 - fractionalSALB388) + &
                              int3_uvaimie_800_aac(salbup388,:,:)*fractionalSALB388

! -- Final interpolated Radiance on pressure --------------------
 ALLOCATE(fint_rad388_aac(naod_aac,ncod_aac), fint_uvaimie_aac(naod_aac,ncod_aac) )

!==============================================================================
! Determine the Terrain Pressure Fraction
!==============================================================================
  if (inpterr .LT. 800.0) then 
     inpterr2 = 800.0
  else
     inpterr2 = inpterr
  endif
  wt = (LOG(1013.0) - LOG(inpterr2)) / (LOG(1013.0) - LOG(800.0))
!write(*,*) 'inpterr2 , wt : ', inpterr2 , wt

  fint_rad388_aac(:,:) = int4_rad388_1013_aac(:,:)*(1.0 -wt) + int4_rad388_800_aac(:,:)*wt
  fint_uvaimie_aac(:,:) = int4_uvaimie_1013_aac(:,:)*(1.0 -wt) + int4_uvaimie_800_aac(:,:)*wt

!
!write(*,*) 'fint_rad388_aac(:,6) : ', fint_rad388_aac(:,6)

  DEALLOCATE(rad_smokeout_aac, rad_dustout_aac, rad_dustoout_aac, &
             uvaimie_smokeout_aac, uvaimie_dustout_aac, uvaimie_dustoout_aac, & 
             intrad388_1013_aac, intrad388_800_aac,  &
             intuvaimie_1013_aac, intuvaimie_800_aac, STAT=status)
!
  DEALLOCATE(int2_rad388_1013_aac, int2_rad388_800_aac, &
             int2_uvaimie_1013_aac, int2_uvaimie_800_aac, STAT=status)
!
  DEALLOCATE(int3_rad388_1013_aac, int3_rad388_800_aac, &
             int3_uvaimie_1013_aac, int3_uvaimie_800_aac, STAT=status)
!
  DEALLOCATE(int4_rad388_1013_aac, int4_rad388_800_aac, &
            int4_uvaimie_1013_aac, int4_uvaimie_800_aac, STAT=status)

 END SUBROUTINE Interpol_aac_LUTparams


 FUNCTION InterpRadiance_aac(iwave,nw01,nw02,ntau1,ntau2,nssa1,nssa2, inRad_aac, outRad_aac) RESULT(status)

!==============================================================================
!
!       This subroutine interpolates radiance values once the Lookup tables
!       have been reduced in size and converted to vectorial form in order to
!       speed up the program by using more efficiently the indices
!
!       Written by Marcos Andrade based on f77 TOMS code (2007)
!
!==============================================================================
  USE LookupTableModule

  IMPLICIT NONE

  INTEGER(KIND=4),                INTENT(IN)    :: iwave
  INTEGER(KIND=4),                INTENT(IN)    :: nw01,nw02
  INTEGER(KIND=4),                INTENT(IN)    :: ntau1,ntau2
  INTEGER(KIND=4),                INTENT(IN)    :: nssa1,nssa2
!
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(IN)   :: inRad_aac
  REAL(KIND=4), DIMENSION(:,:,:),         INTENT(OUT)   :: outRad_aac

  INTEGER(KIND=4)    :: iw0, jw0, itau, issa
  INTEGER(KIND=4)    :: ierr, polyInterp
  REAL(KIND=4)       :: dy

  INTEGER :: i,k,m,n
  INTEGER (KIND=4) :: isza,ithe,iphi, status
  INTEGER (KIND=4) :: idxgral(64), idxaux, model

  status = -1

 n = 0
 do k=1,4
    ithe = indscn+k-1
    do m=1,4
       iphi = indphi+m-1
       do i=1,4
          isza = indsol+i-1
          n = n+1
!               5 -> the number of actual layer heights used on the calculations

          idxgral(n) = iwave+ (isza-1)*5*nssa_aac*naod_aac*ncod_aac + &
                       (iphi-1)*5*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac + &
                       (ithe-1)*5*nssa_aac*naod_aac*ncod_aac*nsza_aac
       enddo
    enddo
 enddo

DO issa = nssa1,nssa2
   DO iw0 = nw01,nw02
      jw0 = iw0 - nw01 + 1
      DO itau = ntau1,ntau2

         idxaux = (issa-1)*5 + (iw0-1)*5*nssa_aac + (itau-1)*5*nssa_aac*naod_aac

         outRad_aac(issa, jw0,itau) = SUM( inRad_aac(idxgral + idxaux)*cofs )
      ENDDO  ! DO itau = 1,ncod_aac
   ENDDO  ! DO iw0  = 1,naod_aac
ENDDO  ! DO issa = 1, nssa_aac

 status = 1
 
 RETURN

 END FUNCTION InterpRadiance_aac


END MODULE Get_omacaLUT7dim_H5module 
