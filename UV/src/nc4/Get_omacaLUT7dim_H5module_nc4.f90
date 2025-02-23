 MODULE Get_omacaLUT7dim_H5module_nc4 
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
 INTEGER(KIND = 4), PARAMETER :: ncod_aac = 10
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
  REAL (KIND=4), DIMENSION(:), ALLOCATABLE ::  rad388_smokelin1013zhgt3_aac,rad388_smokelin1013zhgt4_aac,&
                                                rad388_smokelin1013zhgt5_aac,rad388_smokelin1013zhgt6_aac,&
                                                rad388_smokelin1013zhgt9_aac,rad388_smokelin1013zhgt12_aac,&
						rad388_smokelin1013zhgt15_aac,&
                                                rad388_smokelin800zhgt5_aac,rad388_smokelin800zhgt6_aac,&
                                                rad388_smokelin800zhgt7_aac,rad388_smokelin800zhgt8_aac,&
                                                rad388_smokelin800zhgt11_aac,rad388_smokelin800zhgt14_aac,&
						rad388_smokelin800zhgt17_aac,&
                                                uvaimie_smokelin1013zhgt3_aac,uvaimie_smokelin1013zhgt4_aac,&
                                                uvaimie_smokelin1013zhgt5_aac,uvaimie_smokelin1013zhgt6_aac,&
                                                uvaimie_smokelin1013zhgt9_aac,uvaimie_smokelin1013zhgt12_aac,&
						uvaimie_smokelin1013zhgt15_aac,&
                                                uvaimie_smokelin800zhgt5_aac,uvaimie_smokelin800zhgt6_aac,&
                                                uvaimie_smokelin800zhgt7_aac,uvaimie_smokelin800zhgt8_aac,&
                                                uvaimie_smokelin800zhgt11_aac,uvaimie_smokelin800zhgt14_aac,&
						uvaimie_smokelin800zhgt17_aac

!  Dust over land variables...
 REAL (KIND=4), DIMENSION(:), ALLOCATABLE ::  rad388_dustlin1013zhgt3_aac,rad388_dustlin1013zhgt4_aac,&
                                                rad388_dustlin1013zhgt5_aac,rad388_dustlin1013zhgt6_aac,&
                                                rad388_dustlin1013zhgt9_aac,rad388_dustlin1013zhgt12_aac,&
						rad388_dustlin1013zhgt15_aac,&
                                                rad388_dustlin800zhgt5_aac,rad388_dustlin800zhgt6_aac,&
                                                rad388_dustlin800zhgt7_aac,rad388_dustlin800zhgt8_aac,&
                                                rad388_dustlin800zhgt11_aac,rad388_dustlin800zhgt14_aac,&
						rad388_dustlin800zhgt17_aac,&
                                                uvaimie_dustlin1013zhgt3_aac,uvaimie_dustlin1013zhgt4_aac,&
                                                uvaimie_dustlin1013zhgt5_aac,uvaimie_dustlin1013zhgt6_aac,&
                                                uvaimie_dustlin1013zhgt9_aac,uvaimie_dustlin1013zhgt12_aac,&
						uvaimie_dustlin1013zhgt15_aac,&
                                                uvaimie_dustlin800zhgt5_aac,uvaimie_dustlin800zhgt6_aac,&
                                                uvaimie_dustlin800zhgt7_aac,uvaimie_dustlin800zhgt8_aac,&
                                                uvaimie_dustlin800zhgt11_aac,uvaimie_dustlin800zhgt14_aac,&
						uvaimie_dustlin800zhgt17_aac

!  Dust over ocean variables...
   REAL (KIND=4), DIMENSION(:), ALLOCATABLE ::  rad388_dustolin1013zhgt3_aac,rad388_dustolin1013zhgt4_aac,&
                                                rad388_dustolin1013zhgt5_aac,rad388_dustolin1013zhgt6_aac,&
                                                rad388_dustolin1013zhgt9_aac,rad388_dustolin1013zhgt12_aac,&
						rad388_dustolin1013zhgt15_aac,&
                                                rad388_dustolin800zhgt5_aac,rad388_dustolin800zhgt6_aac,&
                                                rad388_dustolin800zhgt7_aac,rad388_dustolin800zhgt8_aac,&
                                                rad388_dustolin800zhgt11_aac,rad388_dustolin800zhgt14_aac,&
						rad388_dustolin800zhgt17_aac,&
                                                uvaimie_dustolin1013zhgt3_aac,uvaimie_dustolin1013zhgt4_aac,&
                                                uvaimie_dustolin1013zhgt5_aac,uvaimie_dustolin1013zhgt6_aac,&
                                                uvaimie_dustolin1013zhgt9_aac,uvaimie_dustolin1013zhgt12_aac,&
						uvaimie_dustolin1013zhgt15_aac,&
                                                uvaimie_dustolin800zhgt5_aac,uvaimie_dustolin800zhgt6_aac,&
                                                uvaimie_dustolin800zhgt7_aac,uvaimie_dustolin800zhgt8_aac,&
                                                uvaimie_dustolin800zhgt11_aac,uvaimie_dustolin800zhgt14_aac,&
						uvaimie_dustolin800zhgt17_aac

!
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: fint_rad388_aac, fint_uvaimie_aac


 CONTAINS
!
!

 SUBROUTINE Read_aac_LUTparams(lut_fn, rad388_smokelin1013zhgt3_aac,rad388_smokelin1013zhgt4_aac,&
                               rad388_smokelin1013zhgt5_aac,rad388_smokelin1013zhgt6_aac,&
                               rad388_smokelin1013zhgt9_aac,rad388_smokelin1013zhgt12_aac,&
			       rad388_smokelin1013zhgt15_aac,&
                               rad388_smokelin800zhgt5_aac,rad388_smokelin800zhgt6_aac,&
                               rad388_smokelin800zhgt7_aac,rad388_smokelin800zhgt8_aac,&
                               rad388_smokelin800zhgt11_aac,rad388_smokelin800zhgt14_aac,&
			       rad388_smokelin800zhgt17_aac,&
                               uvaimie_smokelin1013zhgt3_aac,uvaimie_smokelin1013zhgt4_aac,&
                               uvaimie_smokelin1013zhgt5_aac,uvaimie_smokelin1013zhgt6_aac,&
                               uvaimie_smokelin1013zhgt9_aac,uvaimie_smokelin1013zhgt12_aac,&
			       uvaimie_smokelin1013zhgt15_aac,&
                               uvaimie_smokelin800zhgt5_aac,uvaimie_smokelin800zhgt6_aac,&
                               uvaimie_smokelin800zhgt7_aac,uvaimie_smokelin800zhgt8_aac,&
                               uvaimie_smokelin800zhgt11_aac,uvaimie_smokelin800zhgt14_aac,&
			       uvaimie_smokelin800zhgt17_aac,&

!                              ;DUST Over-land LUT Parameters...
                               rad388_dustlin1013zhgt3_aac,rad388_dustlin1013zhgt4_aac,&
                               rad388_dustlin1013zhgt5_aac,rad388_dustlin1013zhgt6_aac,&
                               rad388_dustlin1013zhgt9_aac,rad388_dustlin1013zhgt12_aac,&
			       rad388_dustlin1013zhgt15_aac,&
                               rad388_dustlin800zhgt5_aac,rad388_dustlin800zhgt6_aac,&
                               rad388_dustlin800zhgt7_aac,rad388_dustlin800zhgt8_aac,&
                               rad388_dustlin800zhgt11_aac,rad388_dustlin800zhgt14_aac,&
			       rad388_dustlin800zhgt17_aac,&
                               uvaimie_dustlin1013zhgt3_aac,uvaimie_dustlin1013zhgt4_aac,&
                               uvaimie_dustlin1013zhgt5_aac,uvaimie_dustlin1013zhgt6_aac,&
                               uvaimie_dustlin1013zhgt9_aac,uvaimie_dustlin1013zhgt12_aac,&
			       uvaimie_dustlin1013zhgt15_aac,&
                               uvaimie_dustlin800zhgt5_aac,uvaimie_dustlin800zhgt6_aac,&
                               uvaimie_dustlin800zhgt7_aac,uvaimie_dustlin800zhgt8_aac,&
                               uvaimie_dustlin800zhgt11_aac,uvaimie_dustlin800zhgt14_aac,&
			       uvaimie_dustlin800zhgt17_aac,&

!                              ;DUST Over-ocean LUT Parameters...
                               rad388_dustolin1013zhgt3_aac,rad388_dustolin1013zhgt4_aac,&
                               rad388_dustolin1013zhgt5_aac,rad388_dustolin1013zhgt6_aac,&
                               rad388_dustolin1013zhgt9_aac,rad388_dustolin1013zhgt12_aac,&
			       rad388_dustolin1013zhgt15_aac,&
                               rad388_dustolin800zhgt5_aac,rad388_dustolin800zhgt6_aac,&
                               rad388_dustolin800zhgt7_aac,rad388_dustolin800zhgt8_aac,&
                               rad388_dustolin800zhgt11_aac,rad388_dustolin800zhgt14_aac,&
			       rad388_dustolin800zhgt17_aac,&
                               uvaimie_dustolin1013zhgt3_aac,uvaimie_dustolin1013zhgt4_aac,&
                               uvaimie_dustolin1013zhgt5_aac,uvaimie_dustolin1013zhgt6_aac,&
                               uvaimie_dustolin1013zhgt9_aac,uvaimie_dustolin1013zhgt12_aac,&
			       uvaimie_dustolin1013zhgt15_aac,&
                               uvaimie_dustolin800zhgt5_aac,uvaimie_dustolin800zhgt6_aac,&
                               uvaimie_dustolin800zhgt7_aac,uvaimie_dustolin800zhgt8_aac,&
                               uvaimie_dustolin800zhgt11_aac,uvaimie_dustolin800zhgt14_aac,&
			       uvaimie_dustolin800zhgt17_aac,&
                               wavetbl_aac, aodtbl_aac, codtbl_aac, &
			       szatbl_aac,vzatbl_aac,raatbl_aac,ssatbl_aac, salbtbl_aac)
      

 IMPLICIT NONE

!  INTEGER(KIND=4) :: nwav_aac,naod_aac,ncod_aac,nsza_aac,nraa_aac,nvza_aac,nssa_aac, nsalb_aac
  INTEGER (KIND=4) :: linIndex_aac(2), ialo, iz, id

! File, group, dataset and attribute names.
  CHARACTER(LEN=*), INTENT(IN) :: lut_fn

  INTEGER :: wavdims(nwav_aac), aoddims(naod_aac), coddims(ncod_aac)
  INTEGER :: szadims(nsza_aac), vzadims(nvza_aac), raadims(nraa_aac)
  INTEGER :: ssadims(nssa_aac), salbdims(nsalb_aac)
  INTEGER, DIMENSION(7) :: dims7

  
! Regular four-byte integer.
  INTEGER(KIND=4) :: hdf_err, status
  
!
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: rad388_smokelin1013zhgt3_aac,rad388_smokelin1013zhgt4_aac,&
                                                             rad388_smokelin1013zhgt5_aac,rad388_smokelin1013zhgt6_aac,&
                                                             rad388_smokelin1013zhgt9_aac,rad388_smokelin1013zhgt12_aac,&
							     rad388_smokelin1013zhgt15_aac,&
                                                             rad388_smokelin800zhgt5_aac,rad388_smokelin800zhgt6_aac,&
                                                             rad388_smokelin800zhgt7_aac,rad388_smokelin800zhgt8_aac,&
                                                             rad388_smokelin800zhgt11_aac,rad388_smokelin800zhgt14_aac,&
							     rad388_smokelin800zhgt17_aac,&
                                                             uvaimie_smokelin1013zhgt3_aac,uvaimie_smokelin1013zhgt4_aac,&
                                                             uvaimie_smokelin1013zhgt5_aac,uvaimie_smokelin1013zhgt6_aac,&
                                                             uvaimie_smokelin1013zhgt9_aac,uvaimie_smokelin1013zhgt12_aac,&
							     uvaimie_smokelin1013zhgt15_aac,&
                                                             uvaimie_smokelin800zhgt5_aac,uvaimie_smokelin800zhgt6_aac,&
                                                             uvaimie_smokelin800zhgt7_aac,uvaimie_smokelin800zhgt8_aac,&
                                                             uvaimie_smokelin800zhgt11_aac,uvaimie_smokelin800zhgt14_aac,&
							     uvaimie_smokelin800zhgt17_aac
!
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: rad388_dustlin1013zhgt3_aac,rad388_dustlin1013zhgt4_aac,&
                                                             rad388_dustlin1013zhgt5_aac,rad388_dustlin1013zhgt6_aac,&
                                                             rad388_dustlin1013zhgt9_aac,rad388_dustlin1013zhgt12_aac,&
							     rad388_dustlin1013zhgt15_aac,&
                                                             rad388_dustlin800zhgt5_aac,rad388_dustlin800zhgt6_aac,&
                                                             rad388_dustlin800zhgt7_aac,rad388_dustlin800zhgt8_aac,&
                                                             rad388_dustlin800zhgt11_aac,rad388_dustlin800zhgt14_aac,&
							     rad388_dustlin800zhgt17_aac,&
                                                             uvaimie_dustlin1013zhgt3_aac,uvaimie_dustlin1013zhgt4_aac,&
                                                             uvaimie_dustlin1013zhgt5_aac,uvaimie_dustlin1013zhgt6_aac,&
                                                             uvaimie_dustlin1013zhgt9_aac,uvaimie_dustlin1013zhgt12_aac,&
							     uvaimie_dustlin1013zhgt15_aac,&
                                                             uvaimie_dustlin800zhgt5_aac,uvaimie_dustlin800zhgt6_aac,&
                                                             uvaimie_dustlin800zhgt7_aac,uvaimie_dustlin800zhgt8_aac,&
                                                             uvaimie_dustlin800zhgt11_aac,uvaimie_dustlin800zhgt14_aac,&
							     uvaimie_dustlin800zhgt17_aac
!
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: rad388_dustolin1013zhgt3_aac,rad388_dustolin1013zhgt4_aac,&
                                                             rad388_dustolin1013zhgt5_aac,rad388_dustolin1013zhgt6_aac,&
                                                             rad388_dustolin1013zhgt9_aac,rad388_dustolin1013zhgt12_aac,&
							     rad388_dustolin1013zhgt15_aac,&
                                                             rad388_dustolin800zhgt5_aac,rad388_dustolin800zhgt6_aac,&
                                                             rad388_dustolin800zhgt7_aac,rad388_dustolin800zhgt8_aac,&
                                                             rad388_dustolin800zhgt11_aac,rad388_dustolin800zhgt14_aac,&
							     rad388_dustolin800zhgt17_aac,&
                                                             uvaimie_dustolin1013zhgt3_aac,uvaimie_dustolin1013zhgt4_aac,&
                                                             uvaimie_dustolin1013zhgt5_aac,uvaimie_dustolin1013zhgt6_aac,&
                                                             uvaimie_dustolin1013zhgt9_aac,uvaimie_dustolin1013zhgt12_aac,&
							     uvaimie_dustolin1013zhgt15_aac,&
                                                             uvaimie_dustolin800zhgt5_aac,uvaimie_dustolin800zhgt6_aac,&
                                                             uvaimie_dustolin800zhgt7_aac,uvaimie_dustolin800zhgt8_aac,&
                                                             uvaimie_dustolin800zhgt11_aac,uvaimie_dustolin800zhgt14_aac,&
							     uvaimie_dustolin800zhgt17_aac
!
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE,  INTENT(OUT) :: wavetbl_aac, aodtbl_aac,&
                                                             codtbl_aac, szatbl_aac, &
                                                             vzatbl_aac, raatbl_aac, ssatbl_aac, salbtbl_aac
!  For SMOKE Over Land & Ocean
    REAL(KIND=4), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE  :: rad388_1013zhgt3_aac_smk, rad388_1013zhgt4_aac_smk, &
    							    rad388_1013zhgt5_aac_smk, rad388_1013zhgt6_aac_smk,&
                                                            rad388_1013zhgt9_aac_smk, rad388_1013zhgt12_aac_smk, &
							    rad388_1013zhgt15_aac_smk,&
                                                            rad388_800zhgt5_aac_smk, rad388_800zhgt6_aac_smk, &
							    rad388_800zhgt7_aac_smk, rad388_800zhgt8_aac_smk,&
                                                            rad388_800zhgt11_aac_smk, rad388_800zhgt14_aac_smk, &
							    rad388_800zhgt17_aac_smk,&
                                                            uvaimie_1013zhgt3_aac_smk, uvaimie_1013zhgt4_aac_smk, &
							    uvaimie_1013zhgt5_aac_smk, uvaimie_1013zhgt6_aac_smk,&
                                                            uvaimie_1013zhgt9_aac_smk, uvaimie_1013zhgt12_aac_smk, &
							    uvaimie_1013zhgt15_aac_smk,&
                                                            uvaimie_800zhgt5_aac_smk, uvaimie_800zhgt6_aac_smk,&
							    uvaimie_800zhgt7_aac_smk, uvaimie_800zhgt8_aac_smk,&
                                                            uvaimie_800zhgt11_aac_smk, uvaimie_800zhgt14_aac_smk, &
							    uvaimie_800zhgt17_aac_smk
              
!  For DUST Over Land
   REAL(KIND=4), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE  ::  rad388_1013zhgt3_aac_dst, rad388_1013zhgt4_aac_dst, &
   							    rad388_1013zhgt5_aac_dst, rad388_1013zhgt6_aac_dst,&
                                                            rad388_1013zhgt9_aac_dst, rad388_1013zhgt12_aac_dst, &
							    rad388_1013zhgt15_aac_dst,&
                                                            rad388_800zhgt5_aac_dst, rad388_800zhgt6_aac_dst,&
							    rad388_800zhgt7_aac_dst, rad388_800zhgt8_aac_dst,&
                                                            rad388_800zhgt11_aac_dst, rad388_800zhgt14_aac_dst, &
							    rad388_800zhgt17_aac_dst,&
                                                            uvaimie_1013zhgt3_aac_dst, uvaimie_1013zhgt4_aac_dst,&
							    uvaimie_1013zhgt5_aac_dst, uvaimie_1013zhgt6_aac_dst,&
                                                            uvaimie_1013zhgt9_aac_dst, uvaimie_1013zhgt12_aac_dst, &
							    uvaimie_1013zhgt15_aac_dst,&
                                                            uvaimie_800zhgt5_aac_dst, uvaimie_800zhgt6_aac_dst,&
							    uvaimie_800zhgt7_aac_dst, uvaimie_800zhgt8_aac_dst,&
                                                            uvaimie_800zhgt11_aac_dst, uvaimie_800zhgt14_aac_dst, &
							    uvaimie_800zhgt17_aac_dst
  
! For DUST Over Ocean
   REAL(KIND=4), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE  ::  rad388_1013zhgt3_aac_dsto, rad388_1013zhgt4_aac_dsto, &
   							    rad388_1013zhgt5_aac_dsto, rad388_1013zhgt6_aac_dsto,&
                                                            rad388_1013zhgt9_aac_dsto, rad388_1013zhgt12_aac_dsto, &
							    rad388_1013zhgt15_aac_dsto,&
                                                            rad388_800zhgt5_aac_dsto, rad388_800zhgt6_aac_dsto, &
							    rad388_800zhgt7_aac_dsto, rad388_800zhgt8_aac_dsto,&
                                                            rad388_800zhgt11_aac_dsto, rad388_800zhgt14_aac_dsto, &
							    rad388_800zhgt17_aac_dsto,&
                                                            uvaimie_1013zhgt3_aac_dsto, uvaimie_1013zhgt4_aac_dsto, &
							    uvaimie_1013zhgt5_aac_dsto, uvaimie_1013zhgt6_aac_dsto,&
                                                            uvaimie_1013zhgt9_aac_dsto, uvaimie_1013zhgt12_aac_dsto, &
							    uvaimie_1013zhgt15_aac_dsto,&
                                                            uvaimie_800zhgt5_aac_dsto, uvaimie_800zhgt6_aac_dsto,&
							    uvaimie_800zhgt7_aac_dsto, uvaimie_800zhgt8_aac_dsto,&
                                                            uvaimie_800zhgt11_aac_dsto, uvaimie_800zhgt14_aac_dsto, &
							    uvaimie_800zhgt17_aac_dsto
 
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
!
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
!
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
!
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
!
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
!
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
!
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
!
      dset_name = 'NumberOfSingleScatteringAlbedoSmoke'
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
!
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
    ALLOCATE( rad388_1013zhgt3_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt4_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt5_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt6_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt9_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt12_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt15_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt3_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt4_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt5_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt6_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt9_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt12_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt15_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)

    ALLOCATE( rad388_800zhgt5_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_800zhgt6_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_800zhgt7_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_800zhgt8_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt11_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt14_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt17_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt5_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt6_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt7_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt8_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt11_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt14_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt17_aac_smk( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)

! FOR DUST OVER LAND

    ALLOCATE( rad388_1013zhgt3_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt4_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt5_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt6_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt9_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt12_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt15_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt3_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt4_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt5_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt6_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt9_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt12_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt15_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)

    ALLOCATE( rad388_800zhgt5_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_800zhgt6_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_800zhgt7_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_800zhgt8_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt11_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt14_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt17_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt5_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt6_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt7_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt8_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt11_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt14_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt17_aac_dst( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)


! FOR DUST OVER OCEAN
    ALLOCATE( rad388_1013zhgt3_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt4_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt5_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt6_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_1013zhgt9_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt12_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_1013zhgt15_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt3_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt4_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt5_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt6_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt9_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt12_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_1013zhgt15_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             STAT = status)


    ALLOCATE( rad388_800zhgt5_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_800zhgt6_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_800zhgt7_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
              rad388_800zhgt8_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt11_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt14_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             rad388_800zhgt17_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt5_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt6_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt7_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt8_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt11_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt14_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
             uvaimie_800zhgt17_aac_dsto( nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac ), &
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
   dset_name = 'Radiance388_SMOKE_PRESLEV1013_AERHGT9p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt9_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt9_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt9_aac_smk = rad388_1013zhgt9_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV1013_AERHGT9p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt9_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt9_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt9_aac_smk = uvaimie_1013zhgt9_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV1013_AERHGT12p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt12_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt12_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt12_aac_smk = rad388_1013zhgt12_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV1013_AERHGT12p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt12_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt12_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt12_aac_smk = uvaimie_1013zhgt12_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV1013_AERHGT15p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt15_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt15_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt15_aac_smk = rad388_1013zhgt15_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV1013_AERHGT15p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt15_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt15_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt15_aac_smk = uvaimie_1013zhgt15_aac_smk




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
   dset_name = 'Radiance388_SMOKE_PRESLEV0800_AERHGT11p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt11_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt11_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt11_aac_smk = rad388_800zhgt11_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV0800_AERHGT11p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt11_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt11_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt11_aac_smk = uvaimie_800zhgt11_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV0800_AERHGT14p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt14_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt14_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt14_aac_smk = rad388_800zhgt14_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV0800_AERHGT14p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt14_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt14_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt14_aac_smk = uvaimie_800zhgt14_aac_smk
!
   dset_name = 'Radiance388_SMOKE_PRESLEV0800_AERHGT17p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt17_aac_smk)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt17_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt17_aac_smk = rad388_800zhgt17_aac_smk / sflux_lutaac
!
   dset_name = 'UVAIMie_SMOKE_PRESLEV0800_AERHGT17p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt17_aac_smk)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt17_aac_smk, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt17_aac_smk = uvaimie_800zhgt17_aac_smk




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
   dset_name = 'Radiance388_DUST_PRESLEV1013_AERHGT9p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt9_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt9_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt9_aac_dst = rad388_1013zhgt9_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV1013_AERHGT9p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt9_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt9_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt9_aac_dst = uvaimie_1013zhgt9_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV1013_AERHGT12p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt12_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt12_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt12_aac_dst = rad388_1013zhgt12_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV1013_AERHGT12p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt12_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt12_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt12_aac_dst = uvaimie_1013zhgt12_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV1013_AERHGT15p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt15_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt15_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt15_aac_dst = rad388_1013zhgt15_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV1013_AERHGT15p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt15_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt15_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt15_aac_dst = uvaimie_1013zhgt15_aac_dst



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
!
   dset_name = 'Radiance388_DUST_PRESLEV0800_AERHGT11p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt11_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt11_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt11_aac_dst = rad388_800zhgt11_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV0800_AERHGT11p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt11_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt11_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt11_aac_dst = uvaimie_800zhgt11_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV0800_AERHGT14p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt14_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt14_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt14_aac_dst = rad388_800zhgt14_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV0800_AERHGT14p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt14_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt14_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt14_aac_dst = uvaimie_800zhgt14_aac_dst
!
   dset_name = 'Radiance388_DUST_PRESLEV0800_AERHGT17p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt17_aac_dst)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt17_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt17_aac_dst = rad388_800zhgt17_aac_dst / sflux_lutaac
!
   dset_name = 'UVAIMie_DUST_PRESLEV0800_AERHGT17p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt17_aac_dst)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt17_aac_dst, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt17_aac_dst = uvaimie_800zhgt17_aac_dst


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
   dset_name = 'Radiance388_DUSTO_PRESLEV1013_AERHGT9p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt9_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt9_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt9_aac_dsto = rad388_1013zhgt9_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV1013_AERHGT9p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt9_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt9_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt9_aac_dsto = uvaimie_1013zhgt9_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV1013_AERHGT12p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt12_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt12_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt12_aac_dsto = rad388_1013zhgt12_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV1013_AERHGT12p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt12_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt12_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt12_aac_dsto = uvaimie_1013zhgt12_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV1013_AERHGT15p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_1013zhgt15_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_1013zhgt15_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_1013zhgt15_aac_dsto = rad388_1013zhgt15_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV1013_AERHGT15p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_1013zhgt15_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_1013zhgt15_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_1013zhgt15_aac_dsto = uvaimie_1013zhgt15_aac_dsto



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
   dset_name = 'Radiance388_DUSTO_PRESLEV0800_AERHGT11p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt11_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt11_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt11_aac_dsto = rad388_800zhgt11_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV0800_AERHGT11p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt11_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt11_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt11_aac_dsto = uvaimie_800zhgt11_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV0800_AERHGT14p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt14_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt14_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt14_aac_dsto = rad388_800zhgt14_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV0800_AERHGT14p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt14_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt14_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt14_aac_dsto = uvaimie_800zhgt14_aac_dsto
!
   dset_name = 'Radiance388_DUSTO_PRESLEV0800_AERHGT17p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(rad388_800zhgt17_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, rad388_800zhgt17_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  rad388_800zhgt17_aac_dsto = rad388_800zhgt17_aac_dsto / sflux_lutaac
!
   dset_name = 'UVAIMie_DUSTO_PRESLEV0800_AERHGT17p00_aac'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   edge7   = SHAPE(uvaimie_800zhgt17_aac_dsto)
   status = nf90_get_var(grp_id, dset_id, uvaimie_800zhgt17_aac_dsto, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
  uvaimie_800zhgt17_aac_dsto = uvaimie_800zhgt17_aac_dsto

 
    
!
! -- Do Lagrange interpolation for Radiance fields ----------------------
 ALLOCATE(   rad388_smokelin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_smokelin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_smokelin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_smokelin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_smokelin1013zhgt9_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_smokelin1013zhgt12_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_smokelin1013zhgt15_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_smokelin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_smokelin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_smokelin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_smokelin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_smokelin1013zhgt9_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin1013zhgt12_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin1013zhgt15_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)
!
  ALLOCATE(  rad388_smokelin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_smokelin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_smokelin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_smokelin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_smokelin800zhgt11_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_smokelin800zhgt14_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_smokelin800zhgt17_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_smokelin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_smokelin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_smokelin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_smokelin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin800zhgt11_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin800zhgt14_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_smokelin800zhgt17_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)

! DUST Over-land Parameters...
 ALLOCATE(   rad388_dustlin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustlin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustlin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustlin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustlin1013zhgt9_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustlin1013zhgt12_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustlin1013zhgt15_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustlin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustlin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustlin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustlin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustlin1013zhgt9_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin1013zhgt12_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin1013zhgt15_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)
!
  ALLOCATE(  rad388_dustlin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustlin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustlin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustlin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustlin800zhgt11_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustlin800zhgt14_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustlin800zhgt17_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustlin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustlin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustlin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustlin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin800zhgt11_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin800zhgt14_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustlin800zhgt17_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)

! DUST Over-ocean Parameters...

 ALLOCATE(   rad388_dustolin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustolin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustolin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustolin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustolin1013zhgt9_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustolin1013zhgt12_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustolin1013zhgt15_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustolin1013zhgt3_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustolin1013zhgt4_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustolin1013zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustolin1013zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustolin1013zhgt9_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin1013zhgt12_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin1013zhgt15_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)
!
  ALLOCATE(  rad388_dustolin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustolin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustolin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
             rad388_dustolin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustolin800zhgt11_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustolin800zhgt14_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            rad388_dustolin800zhgt17_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustolin800zhgt5_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustolin800zhgt6_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustolin800zhgt7_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
            uvaimie_dustolin800zhgt8_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin800zhgt11_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin800zhgt14_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           uvaimie_dustolin800zhgt17_aac(nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac), &
           STAT = status)


  rad388_smokelin1013zhgt3_aac(:)  = RESHAPE(rad388_1013zhgt3_aac_smk,  & 
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin1013zhgt4_aac(:)  = RESHAPE(rad388_1013zhgt4_aac_smk,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin1013zhgt5_aac(:)  = RESHAPE(rad388_1013zhgt5_aac_smk,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin1013zhgt6_aac(:)  = RESHAPE(rad388_1013zhgt6_aac_smk,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin1013zhgt9_aac(:)  = RESHAPE(rad388_1013zhgt9_aac_smk,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin1013zhgt12_aac(:) = RESHAPE(rad388_1013zhgt12_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin1013zhgt15_aac(:) = RESHAPE(rad388_1013zhgt15_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt3_aac(:) = RESHAPE(uvaimie_1013zhgt3_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt4_aac(:) = RESHAPE(uvaimie_1013zhgt4_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt5_aac(:) = RESHAPE(uvaimie_1013zhgt5_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt6_aac(:) = RESHAPE(uvaimie_1013zhgt6_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt9_aac(:) = RESHAPE(uvaimie_1013zhgt9_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt12_aac(:)= RESHAPE(uvaimie_1013zhgt12_aac_smk,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin1013zhgt15_aac(:)= RESHAPE(uvaimie_1013zhgt15_aac_smk,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
!
  rad388_smokelin800zhgt5_aac(:)  = RESHAPE(rad388_800zhgt5_aac_smk,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin800zhgt6_aac(:)  = RESHAPE(rad388_800zhgt6_aac_smk,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin800zhgt7_aac(:)  = RESHAPE(rad388_800zhgt7_aac_smk,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin800zhgt8_aac(:)  = RESHAPE(rad388_800zhgt8_aac_smk,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin800zhgt11_aac(:) = RESHAPE(rad388_800zhgt11_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin800zhgt14_aac(:) = RESHAPE(rad388_800zhgt14_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_smokelin800zhgt17_aac(:) = RESHAPE(rad388_800zhgt17_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt5_aac(:) = RESHAPE(uvaimie_800zhgt5_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt6_aac(:) = RESHAPE(uvaimie_800zhgt6_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt7_aac(:) = RESHAPE(uvaimie_800zhgt7_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt8_aac(:) = RESHAPE(uvaimie_800zhgt8_aac_smk, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt11_aac(:)= RESHAPE(uvaimie_800zhgt11_aac_smk,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt14_aac(:)= RESHAPE(uvaimie_800zhgt14_aac_smk,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_smokelin800zhgt17_aac(:)= RESHAPE(uvaimie_800zhgt17_aac_smk,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )


! DUST Over-land Parameters...

  rad388_dustlin1013zhgt3_aac(:)  = RESHAPE(rad388_1013zhgt3_aac_dst,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin1013zhgt4_aac(:)  = RESHAPE(rad388_1013zhgt4_aac_dst,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin1013zhgt5_aac(:)  = RESHAPE(rad388_1013zhgt5_aac_dst,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin1013zhgt6_aac(:)  = RESHAPE(rad388_1013zhgt6_aac_dst,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin1013zhgt9_aac(:)  = RESHAPE(rad388_1013zhgt9_aac_dst,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin1013zhgt12_aac(:) = RESHAPE(rad388_1013zhgt12_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin1013zhgt15_aac(:) = RESHAPE(rad388_1013zhgt15_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt3_aac(:) = RESHAPE(uvaimie_1013zhgt3_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt4_aac(:) = RESHAPE(uvaimie_1013zhgt4_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt5_aac(:) = RESHAPE(uvaimie_1013zhgt5_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt6_aac(:) = RESHAPE(uvaimie_1013zhgt6_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt9_aac(:) = RESHAPE(uvaimie_1013zhgt9_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt12_aac(:)= RESHAPE(uvaimie_1013zhgt12_aac_dst,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin1013zhgt15_aac(:)= RESHAPE(uvaimie_1013zhgt15_aac_dst,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
!
  rad388_dustlin800zhgt5_aac(:)  = RESHAPE(rad388_800zhgt5_aac_dst,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin800zhgt6_aac(:)  = RESHAPE(rad388_800zhgt6_aac_dst,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin800zhgt7_aac(:)  = RESHAPE(rad388_800zhgt7_aac_dst,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin800zhgt8_aac(:)  = RESHAPE(rad388_800zhgt8_aac_dst,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin800zhgt11_aac(:) = RESHAPE(rad388_800zhgt11_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin800zhgt14_aac(:) = RESHAPE(rad388_800zhgt14_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustlin800zhgt17_aac(:) = RESHAPE(rad388_800zhgt17_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt5_aac(:) = RESHAPE(uvaimie_800zhgt5_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt6_aac(:) = RESHAPE(uvaimie_800zhgt6_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt7_aac(:) = RESHAPE(uvaimie_800zhgt7_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt8_aac(:) = RESHAPE(uvaimie_800zhgt8_aac_dst, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt11_aac(:)= RESHAPE(uvaimie_800zhgt11_aac_dst,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt14_aac(:)= RESHAPE(uvaimie_800zhgt14_aac_dst,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustlin800zhgt17_aac(:)= RESHAPE(uvaimie_800zhgt17_aac_dst,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )


! DUST Over-ocean Parameters...

  rad388_dustolin1013zhgt3_aac(:)  = RESHAPE(rad388_1013zhgt3_aac_dsto,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin1013zhgt4_aac(:)  = RESHAPE(rad388_1013zhgt4_aac_dsto,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin1013zhgt5_aac(:)  = RESHAPE(rad388_1013zhgt5_aac_dsto,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin1013zhgt6_aac(:)  = RESHAPE(rad388_1013zhgt6_aac_dsto,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin1013zhgt9_aac(:)  = RESHAPE(rad388_1013zhgt9_aac_dsto,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin1013zhgt12_aac(:) = RESHAPE(rad388_1013zhgt12_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin1013zhgt15_aac(:) = RESHAPE(rad388_1013zhgt15_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt3_aac(:) = RESHAPE(uvaimie_1013zhgt3_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt4_aac(:) = RESHAPE(uvaimie_1013zhgt4_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt5_aac(:) = RESHAPE(uvaimie_1013zhgt5_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt6_aac(:) = RESHAPE(uvaimie_1013zhgt6_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt9_aac(:) = RESHAPE(uvaimie_1013zhgt9_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt12_aac(:)= RESHAPE(uvaimie_1013zhgt12_aac_dsto,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin1013zhgt15_aac(:)= RESHAPE(uvaimie_1013zhgt15_aac_dsto,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
!
  rad388_dustolin800zhgt5_aac(:)  = RESHAPE(rad388_800zhgt5_aac_dsto,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin800zhgt6_aac(:)  = RESHAPE(rad388_800zhgt6_aac_dsto,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin800zhgt7_aac(:)  = RESHAPE(rad388_800zhgt7_aac_dsto,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin800zhgt8_aac(:)  = RESHAPE(rad388_800zhgt8_aac_dsto,  &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin800zhgt11_aac(:) = RESHAPE(rad388_800zhgt11_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin800zhgt14_aac(:) = RESHAPE(rad388_800zhgt14_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  rad388_dustolin800zhgt17_aac(:) = RESHAPE(rad388_800zhgt17_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt5_aac(:) = RESHAPE(uvaimie_800zhgt5_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt6_aac(:) = RESHAPE(uvaimie_800zhgt6_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt7_aac(:) = RESHAPE(uvaimie_800zhgt7_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt8_aac(:) = RESHAPE(uvaimie_800zhgt8_aac_dsto, &
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt11_aac(:)= RESHAPE(uvaimie_800zhgt11_aac_dsto,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt14_aac(:)= RESHAPE(uvaimie_800zhgt14_aac_dsto,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )
  uvaimie_dustolin800zhgt17_aac(:)= RESHAPE(uvaimie_800zhgt17_aac_dsto,&
  					(/nsalb_aac*nssa_aac*naod_aac*ncod_aac*nsza_aac*nvza_aac*nraa_aac/) )

 !DEALLOCATE(rad388_1013zhgt3_aac_smk, rad388_1013zhgt4_aac_smk, rad388_1013zhgt5_aac_smk, rad388_1013zhgt6_aac_smk,&
 !           rad388_1013zhgt9_aac_smk, rad388_1013zhgt12_aac_smk, rad388_1013zhgt15_aac_smk,&
 !           rad388_800zhgt5_aac_smk, rad388_800zhgt6_aac_smk, rad388_800zhgt7_aac_smk, rad388_800zhgt8_aac_smk,&
 !           rad388_800zhgt11_aac_smk, rad388_800zhgt14_aac_smk, rad388_800zhgt17_aac_smk,&
 !           uvaimie_1013zhgt3_aac_smk, uvaimie_1013zhgt4_aac_smk, uvaimie_1013zhgt5_aac_smk, uvaimie_1013zhgt6_aac_smk,&
 !           uvaimie_1013zhgt9_aac_smk, uvaimie_1013zhgt12_aac_smk, uvaimie_1013zhgt15_aac_smk,&
 !           uvaimie_800zhgt5_aac_smk, uvaimie_800zhgt6_aac_smk, uvaimie_800zhgt7_aac_smk, uvaimie_800zhgt8_aac_smk,&
 !           uvaimie_800zhgt11_aac_smk, uvaimie_800zhgt14_aac_smk, uvaimie_800zhgt17_aac_smk, STAT=status)

 DEALLOCATE(rad388_1013zhgt3_aac_dst, rad388_1013zhgt4_aac_dst, rad388_1013zhgt5_aac_dst, rad388_1013zhgt6_aac_dst,&
            rad388_1013zhgt9_aac_dst, rad388_1013zhgt12_aac_dst, rad388_1013zhgt15_aac_dst,&
            rad388_800zhgt5_aac_dst, rad388_800zhgt6_aac_dst, rad388_800zhgt7_aac_dst, rad388_800zhgt8_aac_dst,&
            rad388_800zhgt11_aac_dst, rad388_800zhgt14_aac_dst, rad388_800zhgt17_aac_dst,&
            uvaimie_1013zhgt3_aac_dst, uvaimie_1013zhgt4_aac_dst, uvaimie_1013zhgt5_aac_dst, uvaimie_1013zhgt6_aac_dst,&
            uvaimie_1013zhgt9_aac_dst, uvaimie_1013zhgt12_aac_dst, uvaimie_1013zhgt15_aac_dst,&
            uvaimie_800zhgt5_aac_dst, uvaimie_800zhgt6_aac_dst, uvaimie_800zhgt7_aac_dst, uvaimie_800zhgt8_aac_dst,&
            uvaimie_800zhgt11_aac_dst, uvaimie_800zhgt14_aac_dst, uvaimie_800zhgt17_aac_dst, STAT=status)

 DEALLOCATE(rad388_1013zhgt3_aac_dsto, rad388_1013zhgt4_aac_dsto, rad388_1013zhgt5_aac_dsto, rad388_1013zhgt6_aac_dsto,&
            rad388_1013zhgt9_aac_dsto, rad388_1013zhgt12_aac_dsto, rad388_1013zhgt15_aac_dsto,&
            rad388_800zhgt5_aac_dsto, rad388_800zhgt6_aac_dsto, rad388_800zhgt7_aac_dsto, rad388_800zhgt8_aac_dsto,&
            rad388_800zhgt11_aac_dsto, rad388_800zhgt14_aac_dsto, rad388_800zhgt17_aac_dsto,&
            uvaimie_1013zhgt3_aac_dsto, uvaimie_1013zhgt4_aac_dsto, uvaimie_1013zhgt5_aac_dsto, uvaimie_1013zhgt6_aac_dsto,&
            uvaimie_1013zhgt9_aac_dsto, uvaimie_1013zhgt12_aac_dsto, uvaimie_1013zhgt15_aac_dsto,&
            uvaimie_800zhgt5_aac_dsto, uvaimie_800zhgt6_aac_dsto, uvaimie_800zhgt7_aac_dsto, uvaimie_800zhgt8_aac_dsto,&
            uvaimie_800zhgt11_aac_dsto, uvaimie_800zhgt14_aac_dsto, uvaimie_800zhgt17_aac_dsto, STAT=status)
!

90   status = nf90_close(nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to close lut_nc4 file: ", status
      return
   end if
!
 END SUBROUTINE Read_aac_LUTparams


 SUBROUTINE Interpol_aac_LUTparams(sun_za,sat_za,phi,ocean,atype,inpterr, inzhgt, inssa, insalb354, insalb388, &
                                   fint_rad388_aac, fint_uvaimie_aac)
!
   USE LookupTableModule_nc4

   IMPLICIT NONE

  INTEGER(KIND=4) :: status, version, ocean
  INTEGER(KIND=2), INTENT(IN)  :: atype
  REAL(KIND=4),    INTENT(IN)  ::  inpterr, inzhgt, inssa, insalb354, insalb388, sun_za,sat_za,phi

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: fint_rad388_aac, fint_uvaimie_aac

  REAL(KIND=4)     ::  inpterr2
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE    :: rad_smokeout_aac, rad_dustout_aac, rad_dustoout_aac 
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE    :: uvaimie_smokeout_aac, uvaimie_dustout_aac, uvaimie_dustoout_aac

  REAL(KIND=4), DIMENSION(:,:,:,:,:), ALLOCATABLE:: intrad388_1013_aac, intrad388_800_aac,&
                                                    intuvaimie_1013_aac, intuvaimie_800_aac

  REAL(KIND=4), DIMENSION(:,:,:,:), ALLOCATABLE  :: int2_rad388_1013_aac, int2_rad388_800_aac, &
                                                    int2_uvaimie_1013_aac, int2_uvaimie_800_aac
  REAL(KIND=4), DIMENSION(:,:,:),   ALLOCATABLE  :: int3_rad388_1013_aac, int3_rad388_800_aac, &
                                                    int3_uvaimie_1013_aac, int3_uvaimie_800_aac
  REAL(KIND=4), DIMENSION(:,:),     ALLOCATABLE  :: int4_rad388_1013_aac, int4_rad388_800_aac, &
                                                    int4_uvaimie_1013_aac, int4_uvaimie_800_aac

  REAL(KIND=4), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE    :: rad_tmp

  INTEGER(KIND=4)                                     :: sza_indx1, sza_indx2
  INTEGER(KIND=4)                                     :: vza_indx1, vza_indx2
  INTEGER(KIND=4)                                     :: phi_indx1, phi_indx2


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
  REAL(KIND=4), DIMENSION(7) :: ssa388val_dst_land = (/0.779210, 0.837780, 0.876060, 0.905320, 0.948860, 0.972340, 1.000/)
  REAL(KIND=4), DIMENSION(7) :: ssa388val_dst_ocean = (/0.767990, 0.829980, 0.870520, 0.901320, 0.946900, 0.971240, 1.000/)
  REAL(KIND=4), DIMENSION(7) :: zhgt1013val = (/3.0, 4.0, 5.0, 6.0, 9.0, 12.0, 15.0/)
  REAL(KIND=4), DIMENSION(5) :: salbval = (/0.0, 0.05, 0.10, 0.15, 0.20/)

  INTEGER(KIND=4) :: zhgtlw, zhgtup, ssalw, ssaup, salblw354, salbup354, salblw388, salbup388
  INTEGER(KIND=4) ::  ialo, iz

!
   ALLOCATE(intrad388_1013_aac(nsalb_aac,7,nssa_aac,naod_aac,ncod_aac), &
             intrad388_800_aac(nsalb_aac,7,nssa_aac,naod_aac,ncod_aac) )

   ALLOCATE(intuvaimie_1013_aac(nsalb_aac,7,nssa_aac,naod_aac,ncod_aac), &
             intuvaimie_800_aac(nsalb_aac,7,nssa_aac,naod_aac,ncod_aac) )
!
   ALLOCATE(    rad_smokeout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE(     rad_dustout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE(    rad_dustoout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE(uvaimie_smokeout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE( uvaimie_dustout_aac(nssa_aac,naod_aac,ncod_aac) )
   ALLOCATE(uvaimie_dustoout_aac(nssa_aac,naod_aac,ncod_aac) )
   
   ALLOCATE(rad_tmp(nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac))


!=============================================================================
! Find the two index values for input zhgt vlaues
!=============================================================================
  status=FindTableEntry(inzhgt, zhgt1013val, 7, theta1,theta2, zhgtlw, zhgtup,&
                                                                 fractionalZhgt)
  if (zhgtlw .ge. 7) zhgtup = zhgtlw
  if (zhgtlw .lt. 1) then
      zhgtlw = 1
      zhgtup = zhgtlw
  endif

!=============================================================================
! Find the two index values for input SSA vlaues
!=============================================================================
  IF(atype .eq. 1) THEN
  status=FindTableEntry(inssa, ssa388val_smk, 7, theta1,theta2, ssalw, ssaup,&
                                                           fractionalSSA)
  ENDIF

  IF(atype .eq. 2) THEN
  status=FindTableEntry(inssa, ssa388val_dst_land, 7, theta1,theta2, ssalw, ssaup,&
                                                           fractionalSSA)
  ENDIF

  IF(atype .eq. 4) THEN
  status=FindTableEntry(inssa, ssa388val_dst_ocean, 7, theta1,theta2, ssalw, ssaup,&
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


 ! Find bounding indices of LUT for observed SZA, VZA, and PHI...
  CALL Find_Indices(nraa_aac, raatbl_aac, phi, phi_indx1, phi_indx2)
  CALL Find_Indices(nvza_aac, vzatbl_aac, sat_za, vza_indx1, vza_indx2)
  CALL Find_Indices(nsza_aac, szatbl_aac, sun_za, sza_indx1, sza_indx2)


  rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt3_aac, &
                           (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

  IF(atype .eq. 1) THEN
   !
  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))                                      
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))                                      
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))                                      
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))                                      
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))                                      
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))                                      

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_smokeout_aac)

   intrad388_1013_aac(ialo,zhgtlw,:,:,:) = rad_smokeout_aac         


   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   
   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_smokeout_aac)

   intrad388_1013_aac(ialo,zhgtup,:,:,:) = rad_smokeout_aac
   
  ENDDO   !DO ialo = salblw388, salbup388
   

  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,uvaimie_smokeout_aac)                                      

   intuvaimie_1013_aac(ialo,zhgtlw,:,:,:) = uvaimie_smokeout_aac


   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,uvaimie_smokeout_aac)
                                      
   intuvaimie_1013_aac(ialo,zhgtup,:,:,:) = uvaimie_smokeout_aac
  ENDDO   !DO ialo = salblw388, salbup388


 ! -- 800 hPa --------------------------------------------------
   !
  DO ialo = salblw388, salbup388


   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_smokeout_aac)

   intrad388_800_aac(ialo,zhgtlw,:,:,:) = rad_smokeout_aac


   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_smokelin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_smokeout_aac)
                                      
   intrad388_800_aac(ialo,zhgtup,:,:,:) = rad_smokeout_aac
  ENDDO   !DO ialo = salblw388, salbup388


  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,uvaimie_smokeout_aac)

   intuvaimie_800_aac(ialo,zhgtlw,:,:,:) = uvaimie_smokeout_aac

   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_smokelin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,uvaimie_smokeout_aac)

   intuvaimie_800_aac(ialo,zhgtup,:,:,:) = uvaimie_smokeout_aac
  ENDDO   !DO ialo = salblw388, salbup388

 ENDIF


 IF(atype .eq. 2) THEN     !over-land pixel
  !
  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)

   intrad388_1013_aac(ialo,zhgtlw,:,:,:) = rad_dustout_aac


   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)

   intrad388_1013_aac(ialo,zhgtup,:,:,:) = rad_dustout_aac
  ENDDO   !DO ialo = salblw388, salbup388

  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)

   intuvaimie_1013_aac(ialo,zhgtlw,:,:,:) = uvaimie_dustout_aac


   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)

   intuvaimie_1013_aac(ialo,zhgtup,:,:,:) = uvaimie_dustout_aac
  ENDDO   !DO ialo = salblw388, salbup388

 ! -- 800 hPa --------------------------------------------------
  !
  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)

   intrad388_800_aac(ialo,zhgtlw,:,:,:) = rad_dustout_aac


   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustlin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)
   intrad388_800_aac(ialo,zhgtup,:,:,:) = rad_dustout_aac
  ENDDO   !DO ialo = salblw388, salbup388

  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)

   intuvaimie_800_aac(ialo,zhgtlw,:,:,:) = uvaimie_dustout_aac

   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustlin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)

   intuvaimie_800_aac(ialo,zhgtup,:,:,:) = uvaimie_dustout_aac
  ENDDO   !DO ialo = salblw388, salbup388
 ENDIF



IF(atype .eq. 4) THEN     !over-ocean pixel
  !
  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustoout_aac)

   intrad388_1013_aac(ialo,zhgtlw,:,:,:) = rad_dustoout_aac

   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustoout_aac)

   intrad388_1013_aac(ialo,zhgtup,:,:,:) = rad_dustoout_aac
  ENDDO   !DO ialo = salblw388, salbup388

  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,uvaimie_dustoout_aac)

   intuvaimie_1013_aac(ialo,zhgtlw,:,:,:) = uvaimie_dustoout_aac


   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt3_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt4_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt9_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt12_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin1013zhgt15_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,uvaimie_dustoout_aac)

   intuvaimie_1013_aac(ialo,zhgtup,:,:,:) = uvaimie_dustoout_aac
  ENDDO   !DO ialo = salblw388, salbup388


! -- 800 hPa --------------------------------------------------
  !
  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)

   intrad388_800_aac(ialo,zhgtlw,:,:,:) = rad_dustoout_aac

   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(rad388_dustolin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,rad_dustout_aac)

   intrad388_800_aac(ialo,zhgtup,:,:,:) = rad_dustoout_aac
  ENDDO   !DO ialo = salblw388, salbup388

  DO ialo = salblw388, salbup388

   IF(zhgtlw .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtlw .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,uvaimie_dustoout_aac)

   intuvaimie_800_aac(ialo,zhgtlw,:,:,:) = uvaimie_dustoout_aac


   IF(zhgtup .EQ. 1) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt5_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 2) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt6_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 3) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt7_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 4) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt8_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 5) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt11_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 6) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt14_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))
   IF(zhgtup .EQ. 7) rad_tmp(:,:,:,:,:,:,:) = RESHAPE(uvaimie_dustolin800zhgt17_aac, &
                                              (/nsalb_aac,nssa_aac,naod_aac,ncod_aac,nsza_aac,nvza_aac,nraa_aac/))

   status = InterpRadiance_Linear(nssa_aac,naod_aac,ncod_aac,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,PHI_indx1, &
                 PHI_indx2,sun_za,sat_za,phi,rad_tmp,uvaimie_dustoout_aac)

   intuvaimie_800_aac(ialo,zhgtup,:,:,:) = uvaimie_dustoout_aac
  ENDDO   !DO ialo = salblw388, salbup388
 ENDIF


! --------------------------------------------------------------------------------

!!  -- Do linear interpolation on SSA --
 ALLOCATE(int2_rad388_1013_aac(nsalb_aac,7,naod_aac,ncod_aac), &
           int2_rad388_800_aac(nsalb_aac,7,naod_aac,ncod_aac), &
         int2_uvaimie_1013_aac(nsalb_aac,7,naod_aac,ncod_aac), &
	  int2_uvaimie_800_aac(nsalb_aac,7,naod_aac,ncod_aac) )

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
  wt = (LOG(1013.25) - LOG(inpterr2)) / (LOG(1013.25) - LOG(800.0))

  fint_rad388_aac(:,:) = int4_rad388_1013_aac(:,:)*(1.0 -wt) + int4_rad388_800_aac(:,:)*wt
  fint_uvaimie_aac(:,:) = int4_uvaimie_1013_aac(:,:)*(1.0 -wt) + int4_uvaimie_800_aac(:,:)*wt

!
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

  DEALLOCATE(rad_tmp, STAT=status)

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
  USE LookupTableModule_nc4

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
!               5 -> the number of surface albedo nodes

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



FUNCTION InterpRadiance_Linear(nssa_in,naod_in,ncod_in,SZA_indx1,SZA_indx2,VZA_indx1,VZA_indx2,&
          RAA_indx1, RAA_indx2,sun_za,sat_za,phi,rad_in,outRad_ai) RESULT(status)

USE LookupTableModule_nc4

IMPLICIT NONE

INTEGER(KIND=4)  :: iSSA, iAOD, iCOD

INTEGER(KIND=4)              :: iSZA, iVZA, iRAA
INTEGER(KIND=4),  INTENT(IN) :: RAA_indx1, RAA_indx2
INTEGER(KIND=4), INTENT(IN)  :: VZA_indx1, VZA_indx2
INTEGER(KIND=4), INTENT(IN)  :: SZA_indx1, SZA_indx2

REAL(KIND=4), INTENT(IN)     :: sun_za, sat_za, phi


INTEGER(KIND=4), INTENT(IN) :: nssa_in,naod_in,ncod_in
REAL(KIND=4), DIMENSION(5,7,7,10,7,14,11), INTENT(IN)    :: rad_in
REAL(KIND=4), DIMENSION(:,:,:), INTENT(OUT)              :: outRad_ai

!Local arrays
REAL(KIND=4) :: rad_PHI(nphi), val_intp
REAL(KIND=4) :: rad_dim1(nssa_in, naod_in, ncod_in, nsza, ntheta)
REAL(KIND=4) :: rad_dim2(nssa_in, naod_in, ncod_in, nsza)

  INTEGER(KIND=4)        :: status
  CHARACTER(LEN=12)     :: mytime

  status = -1

DO iSSA = 1, nssa_in
DO iAOD = 1, naod_in
DO iCOD = 1, ncod_in

    !Solar Zenith Angle Loop?~@?
    DO iSZA = SZA_indx1, SZA_indx2

    !Viewing Zenith Angle Loop?~@?
    DO iVZA = VZA_indx1, VZA_indx2

    ! Resolving RAA here...
    DO iRAA = 1, nphi
    rad_PHI(iRAA) = rad_in(1, iSSA, iAOD, iCOD, iSZA, iVZA, iRAA)
    ENDDO

    CALL Interpolation(nphi, phi_table, rad_PHI, phi, val_intp)
    rad_dim1(iSSA, iAOD, iCOD, iSZA, iVZA) = val_intp

    ENDDO  ! iVZA = VZA_indx1, VZA_indx2

    ! Resolving VZA here...
    CALL Interpolation(2, theta_table(VZA_indx1:VZA_indx2), rad_dim1(iSSA, iAOD, iCOD, iSZA, VZA_indx1:VZA_indx2), sat_za, val_intp)
    rad_dim2(iSSA, iAOD, iCOD, iSZA) = val_intp

    ENDDO  ! iSZA = SZA_indx1, SZA_indx2

! Resolving SZA here...
    CALL Interpolation(2,sza_table(SZA_indx1:SZA_indx2), &
           rad_dim2(iSSA, iAOD, iCOD, SZA_indx1:SZA_indx2), sun_za, val_intp)
    outRad_ai(iSSA, iAOD, iCOD) = val_intp

    ENDDO
    ENDDO
    ENDDO

    RETURN

END FUNCTION InterpRadiance_Linear

!**************************************************************
SUBROUTINE Find_Indices(n_nodes,val_nodes,val_true,indx1,indx2)
!**************************************************************

 IMPLICIT NONE

 INTEGER(KIND = 4)                              :: m
 INTEGER, INTENT(IN)                            :: n_nodes
 REAL, INTENT(IN), DIMENSION(n_nodes)           :: val_nodes
 REAL, INTENT(IN)                               :: val_true
 INTEGER, INTENT(OUT)                           :: indx1, indx2

 indx1 = 0
 indx2 = 0

 DO m = 1, n_nodes-1

 IF ( val_true .GE. val_nodes(m) .AND. val_true .LE. val_nodes(m+1)) THEN
      indx1 = m
      indx2 = m+1
 ENDIF

 ENDDO

END SUBROUTINE Find_Indices

SUBROUTINE INTERPOLATION(n_nodes,nodes,val_at_nodes,node_true,val_intp)

 IMPLICIT NONE

 INTEGER(KIND = 4)                              :: m
 INTEGER, INTENT(IN)                            :: n_nodes
 REAL, INTENT(IN), DIMENSION(n_nodes)           :: nodes
 REAL, INTENT(IN), DIMENSION(n_nodes)           :: val_at_nodes
 REAL(KIND=4)                                   :: node_true
 REAL, INTENT(OUT)                              :: val_intp
 INTEGER(KIND=4)                                :: indx

 val_intp = -9.99
 indx = 0

 DO m = 1, n_nodes-1

 IF ( nodes(m) .GE. 0 .AND. nodes(m+1) .GE. 0 .AND. indx .EQ. 0 ) THEN

    IF ((node_true .GE. nodes(m) .AND. node_true .LE. nodes(m+1)) .OR. &
        (node_true .LE. nodes(m) .AND. node_true .GE. nodes(m+1))) THEN
         val_intp = val_at_nodes(m) + ((val_at_nodes(m+1)-val_at_nodes(m)) &
         / (nodes(m+1)-nodes(m)) * (node_true-nodes(m)))
         indx = 1
    ENDIF

 ENDIF

 ENDDO

END SUBROUTINE Interpolation

END MODULE Get_omacaLUT7dim_H5module_nc4 
