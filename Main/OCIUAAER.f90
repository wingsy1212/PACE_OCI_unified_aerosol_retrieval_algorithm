!
! This is the main program for OCI Unified Algorithm AERosol product.
! SHORTNAME: OCIUAAER
!
PROGRAM OCIUAAER


USE OCIUAAER_Config_Module
USE OCIUAAER_L1BModule
USE Proxydata_L1BModule
USE OCI_UV1km_DataModule
USE OCI_UVAI_MieCloudModule
USE GetLUT_Miecloud_AI_H5module_nc4 
USE compute_Ocean_main
USE compute_DT_Land_main
USE compute_DB_Land_main
use calendars, only: gdatetime, gregorian_from_doy
use write_Pace_merged
use write_Pace_merged_1KM
use write_intrim_forUVtau
use read_intrim_forUVtau
USE HDF5

USE read_for_nuv_test
USE NUV_PACKAGE_MODULES 
USE  merge_land_ocean


IMPLICIT NONE

Include 'common_l1b_var.inc'
Include 'output_Variables.inc'

!==============================================================================
! Declare Global variables.
!==============================================================================
CHARACTER(LEN=256)                         :: l1b_filename, l2_file, time_str
INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE :: grn_lwmask
REAL(KIND=4),DIMENSION(:,:,:),ALLOCATABLE  :: UVtoSWIR_Reflectances
REAL(KIND=4),DIMENSION(:),ALLOCATABLE      :: UVtoSWIR_wavelengths
INTEGER(KIND=4)                            :: Day, Month, Year, ifile, doy, hr, mn
INTEGER(KIND=4)                            :: STATUS, hdferr, nXTrack, nLines
INTEGER(KIND=4), PARAMETER                 :: UVtoSWIR_nWavel=14
REAL(KIND=4),DIMENSION(UVtoSWIR_nWavel) :: waveTemp = [340.0, 354.0, 388.0, &
                         412.0, 488.0, 550.0, 670.0,680.0,688.0,860.0, &
                        1250.0,1378.0,1615.0,2260.0]

logical          ::  do_testfile
type(gdatetime)  ::  gdt1
REAL,DIMENSION(:,:,:),ALLOCATABLE  :: dtspec_test,dtaod_test
REAL,DIMENSION(:,:),ALLOCATABLE    :: dtlon,dtlat
!REAL,DIMENSION(IXRet_B,IYRet_B,5) :: uvdbdtaod
REAL,DIMENSION(:,:,:),ALLOCATABLE  :: uvdbdtaod
REAL,DIMENSION(:,:),ALLOCATABLE    :: dbdtfmf
REAL,DIMENSION(:,:,:),ALLOCATABLE  :: dbdt_refl
REAL,DIMENSION(:,:),ALLOCATABLE    :: dbdt_cld
INTEGER, DIMENSION(2)              :: cell_dims
INTEGER IX,IY,IL,IXX,IYY,jjx,jjy,clear,IJ,IK,iwav,land_pixels
real :: start, finish,cloudy
real :: lat_min,lat_max,lon_max,lon_min

                  
  INTEGER :: num_args
  CHARACTER(255) :: prog_name, par_file, l1b_file, met1_file, met2_file, out_file,&
   interm_file, out_file_1KM
  CHARACTER(255) :: metdt_file
  CHARACTER(2048) :: pycommand

  CHARACTER(len=8) :: date
  CHARACTER(len=10) :: rtime
  character(len=4) :: year_str
  character(len=2) :: month_str, day_str, hour_str, minute_str, second_str
  CHARACTER(len=5) :: zone
  integer :: values(8)

  num_args = command_argument_count()
  IF (num_args == 7) THEN
      call get_command_argument(0,prog_name)
      call get_command_argument(1,par_file)
      call get_command_argument(2,l1b_file)
      call get_command_argument(3,met1_file)
      call get_command_argument(4,met2_file)
      call get_command_argument(5,out_file)
      call get_command_argument(6,interm_file)
      call get_command_argument(7,out_file_1KM)
  ELSE
    PRINT *,'Error : Not enough arguments.'
    PRINT *, 'oci_ua_aer <Par_Filename> <L1B_Filename> <GEO_Data_start_filename> &
    <GEO_Data_End_filename> <L2_Output_Filename> Interm_L2.nc <L2_1KM_Output_Filename>'
        CALL EXIT(1)
  ENDIF

CALL h5open_f(hdferr)

!==============================================================================
! Read the Runtime paramters/LUT from the Config file.
!==============================================================================
  STATUS = read_ociuaaer_config(par_file, cfg)
  IF (STATUS < 0) THEN 
    PRINT *,'Error : Unable to Read OCIUAAER Config file..'
    CALL EXIT(1)
  ENDIF 
  doy=0
  cfg%input_l1file = l1b_file
  l1b_filename = cfg%input_l1file

  call date_and_time(date, rtime, zone, values)

  ! Convert integers to strings
  write(year_str, '(I4)') values(1)
  write(month_str, '(I2.2)') values(2)  ! Zero-padded 2-digit month
  write(day_str, '(I2.2)') values(3)    ! Zero-padded 2-digit day
  write(hour_str, '(I2.2)') values(5)   ! Zero-padded 2-digit hour
  write(minute_str, '(I2.2)') values(6) ! Zero-padded 2-digit minute
  write(second_str, '(I2.2)') values(7) ! Zero-padded 2-digit second

  cfg%history = TRIM(year_str) // '-' // TRIM(month_str) // '-' // TRIM(day_str) // 'T' // &
                TRIM(hour_str) // ':' // TRIM(minute_str) // ':' // TRIM(second_str) // ': ' // &
                TRIM(prog_name) // ' ' // TRIM(par_file) // ' ' // TRIM(l1b_file) // ' ' // &
                TRIM(met1_file) // ' ' // TRIM(met2_file) // ' ' // TRIM(out_file) // ' ' // &
                TRIM(interm_file) // ' ' // TRIM(out_file_1KM)

  ! print *, 'History = ', cfg%history

  do_testfile = .false.  !.true.  !.false.
! do_testfile = .true.
! if (do_testfile) GOTO 7007 ! For NUV-testing
 
 if (.not. do_testfile) then ! normal process
   call cpu_time(start)
   print *, 'start L1b read', start
   
   ! Allocate space for UVtoSWIR wavelengths
   ALLOCATE(UVtoSWIR_wavelengths(UVtoSWIR_nWavel), stat=STATUS) 
   UVtoSWIR_wavelengths = waveTemp(1:UVtoSWIR_nWavel)
  
   print*,'Inputfile = ',cfg%input_l1file

   IF (cfg%input_l1file /= 'NULL') THEN
    STATUS = OCI_l1b_get_meta_data()
    IF (STATUS < 0) THEN 
      PRINT *,'Error : Unable to Read L1B-Coverage Meta-Data'
      CALL EXIT(1)
    ENDIF 

    ! print *, 'Time Coverage Start = ', cfg%coverage_start
    ! print *, 'Time Coverage End = ', cfg%coverage_end

    ! Read Synthetic data file
     STATUS = OCI_l1b_get_reflect(cfg%input_l1file, l1b_nXTrack, l1b_nLines, l1b_date_time_str, &
                                Latitude, Longitude, SolarZenithAngle, ViewingZenithAngle, &
                                SolarAzimuthAngle, ViewingAzimuthAngle, TerrainHeight, &
                                UVtoSWIR_nWavel, UVtoSWIR_wavelengths, & 
                                UVtoSWIR_Reflectances)
     IF (STATUS < 0) THEN 
       PRINT *,'Error : Unable to Read L1B-Synthetic Data'
       CALL EXIT(1)
     ENDIF 

     ! Extract year, month, day, and time info. and create L2 filename.
     Read(l1b_date_time_str(1:4), '(I4)' )Year
     Read(l1b_date_time_str(5:6), '(I2)' )Month
     Read(l1b_date_time_str(7:8), '(I2)' )Day
     Read(l1b_date_time_str(10:11), '(I2)' )hr
     Read(l1b_date_time_str(12:13), '(I2)' )mn

     if (index(met1_file, 'GMAO') > 0) then
        if (mn < 30) then
            metdt_file = met1_file
        else
            metdt_file = met2_file
        endif
     else
        if (mod(hr,3) * 60 + mn < 90) then
            metdt_file = met1_file
        else
            metdt_file = met2_file
        endif
     endif

   ELSE IF (cfg%proxy_l1file /= 'NULL') THEN
     ! Read Proxydata file
     STATUS = read_l1bproxydata(cfg%proxy_l1file, l1b_nXTrack, l1b_nLines, l1b_date_time_str, &
                              Latitude, Longitude, SolarZenithAngle, ViewingZenithAngle, &
                              SolarAzimuthAngle, ViewingAzimuthAngle, TerrainHeight, &
                              UVtoSWIR_nWavel, UVtoSWIR_wavelengths, & 
                              UVtoSWIR_Reflectances)
     IF (STATUS < 0) THEN 
       PRINT *,'Error : Unable to Read L1B-Proxy Data'
       CALL EXIT(1)
     ENDIF 

     !Extract year, month, day, and time info. and create L2 filename.
     Read(l1b_date_time_str(1:4), '(I4)' )Year
     Read(l1b_date_time_str(5:7), '(I3)' )doy
     gdt1  = gregorian_from_doy(year, doy)
     Month = gdt1%month
     Day = gdt1%mday

   ELSE 

     PRINT *, 'Error : NOT a valid Inputfile (Proxy or Synthetic)'
     CALL EXIT(1)

   ENDIF

   ! Read ancillary global land-water mask data and get mask for the input granule.
   STATUS = get_gmt15arc_lwmask(cfg%landwater_mask, l1b_nXTrack, l1b_nLines, &
                                  Latitude, Longitude, grn_lwmask) 
   IF (STATUS < 0) THEN 
     PRINT *,'Error : Unable to Read landwater_mask Data'
     CALL EXIT(1)
   ENDIF 


   ! Allocate space for UVAI
   ALLOCATE(               UVAI(l1b_nXTrack, l1b_nLines), &
                    Residue_1km(l1b_nXTrack, l1b_nLines), & 
	    Reflectivity_1km(2, l1b_nXTrack, l1b_nLines), stat=STATUS) 
   IF (STATUS < 0) THEN 
     PRINT *,'Error : Unable allocate space for UVAI and diagnostics.'
     CALL EXIT(1)
   ENDIF 

! Computes Aerosol Index for 1-km pixels  Right now disabling it pending..............
   CALL OCI_UVAI_Miecloud(cfg, l1b_nXTrack, l1b_nLines, Month, &
                         Latitude, Longitude, SolarZenithAngle, ViewingZenithAngle, &
                         SolarAzimuthAngle, ViewingAzimuthAngle, TerrainHeight, &
                         UVtoSWIR_Reflectances, UVAI, Residue_1km, Reflectivity_1km)
   call cpu_time(finish)                      
   print *, 'end (L1b read + 1km UVAI)', finish
!   print '("L1b read + 1km UVAI, Time = ",f10.3," seconds.")',finish-start     
   !==============================================================================
   !==============================================================================
   !
   ! DT and DB modules begins here 
   ! Variables to work with (at 1 km resolution): 
   ! 
   ! UVtoSWIR_Reflectances** (Please note) **
   ! 340 ,354,388,412,488,550,670,680,688,860,1250.15,1378.04,1615.98,2260.64 nm
   
   ! Latitude, Longitude
   ! SolarZenithAngle, ViewingZenithAngle, SolarAzimuthAngle, ViewingAzimuthAngle
   ! TerrainHeight
   ! grn_lwmask = 0-land, 1-water

   print *, 'start DB, DT retrieval'
   call cpu_time(start)
   print *, 'start DT', start 
   call DT_package_Land_main(Latitude,Longitude,SolarZenithAngle,&
                         SolarAzimuthAngle,&
                         ViewingZenithAngle,ViewingAzimuthAngle,TerrainHeight,&
                         UVtoSWIR_Reflectances,UVAI,Year,Month,Day,&
                         grn_lwmask,l1b_nXTrack,l1b_nLines,&
                         ! Output variables     
   !   Ret_tau_land  at       0.48, 0.55, 0.67, 2.25 microns   
   !   Ret_ref_allwav_land   at 0.354,0.388,0.48, 0.55, 0.67,0.87,1.24,1.64 and 2.25 
   !   CldMsk_Native_Land ( 0=cloud 1= clear)   Ret_land_Quality_Flag( 0-3)                
          Ret_ref_allwav_land,Ret_tau_land,Ret_Lat,Ret_Lon,CldMsk_Native_Land,&
          Ret_land_Quality_Flag,Ret_Small_weighting_land,Ret_CLDFRC_land_DT,&
          Ret_Xtrack,Ret_Lines,metdt_file)
          
   call cpu_time(finish)       
   print *, 'end DT', finish
   print '("DT Time = ",f10.3," seconds.")',finish-start             
   
   print *, ' DT_Land retrieval done'   
   ! Write test output file for fast process
!     call write_pace_test_data('../L2_data/PACE.test_file.hdf',&
!           Latitude,Longitude,SolarZenithAngle,&
!           SolarAzimuthAngle,&
!           ViewingZenithAngle,ViewingAzimuthAngle,TerrainHeight,&
!           UVtoSWIR_Reflectances,Year,Month,Day,&
!           grn_lwmask,UVAI,Ret_ref_allwav_land,Ret_tau_land,Ret_Lat,Ret_Lon,&
!           CldMsk_Native_Land,Ret_land_Quality_Flag,Ret_Small_weighting_land,Ret_Xtrack,Ret_Lines)  
             
 endif !if (.not. do_testfile) then      

 nXTrack=size(UVtoSWIR_Reflectances,2)
 nLines=size(UVtoSWIR_Reflectances,3)

! cell_dims =(/size(Ret_tau_land, 1), size(Ret_tau_land, 2)/)

 call cpu_time(start)
 print *, 'start DB', start

 allocate(uvdbdtaod(Ret_Xtrack,Ret_Lines,5),dbdtfmf(Ret_Xtrack,Ret_Lines), stat=status)

 call DB_package_Land_main(Latitude,Longitude,SolarZenithAngle,SolarAzimuthAngle,&
    ViewingZenithAngle,ViewingAzimuthAngle,TerrainHeight,UVtoSWIR_Reflectances,Year,Month,Day,&
    doy,cfg%proxy_l1file,grn_lwmask,nXTrack, nLines,UVtoSWIR_nWavel,cfg%db_config,Ret_ref_allwav_land,&
    Ret_tau_land,Ret_Lat,Ret_Lon,CldMsk_Native_land,Ret_land_Quality_Flag,uvdbdtaod,dbdt_refl,&
    cfg%input_l1file,Ret_Small_weighting_land,Ret_Xtrack,Ret_Lines,dbdtfmf,dbdt_cld,met1_file,met2_file, UVAI)
      
 call cpu_time(finish)
 print *, 'end DB', finish
 print '("DB Time = ",f10.3," seconds.")',finish-start      
  ! Output variables     
  ! uvdbdtaod  at 0.354,0.388 0.48, 0.55, 0.67 microns (combined AOD, array size:[400,404,5])
  ! dbdt_refl  at 0.354,0.388,0.48, 0.55, 0.67,0.87,1.24,1.64 and 2.25  microns (combined reflectance, array size:[400,404,9])
  ! dbdtfmf  (DT FMF & DB type combined, array size:[400,404])  
  ! Ret_CldMsk_500_Land (DBDT cloud mask, overwrite on DT cloud mask array. Native resolution)
      call cpu_time(start)
     print *, 'start DT_0cean', start
       
      call DT_package_Ocean_main(Latitude,Longitude,SolarZenithAngle,&
                   SolarAzimuthAngle,&
                   ViewingZenithAngle,ViewingAzimuthAngle,TerrainHeight,&
                   UVtoSWIR_Reflectances,UVAI,Year,Month,Day,&
                   grn_lwmask,l1b_nXTrack,l1b_nLines,&
! ! Output variables                        
                   Ret_Xtrack,Ret_Lines,&
                   Ret_Lat,Ret_Lon,Ret_SolZen,Ret_View_angle,&
                   Ret_View_phi,Ret_solar_phi,Ret_tau_ocean,& 
                   Ret_Small_weighting,Ret_ref_ocean,Land_sea_flag,&
                   Ret_average_Omega_Ocean_UV,CldMsk_Native_Ocean,&
                   Ret_Index_Height,Ret_ocean_Quality,Ret_CLDFRC_Ocean,metdt_file)
                      
  
!            print *,'"Ret_Xtrack,Ret_Lines:" i4,i4 ',Ret_Xtrack,Ret_Lines           
!!!!!                        
            call cpu_time(finish)
            print *, 'end DT-Ocean', finish
            print '("DT Ocean Time = ",f10.3," seconds.")',finish-start   
             
         land_pixels = 0
          Do  IY = 1,Ret_Lines  
              Do  IX =  1,Ret_Xtrack  
                IF(Land_sea_flag(IX,IY) .ge. 1 .and. uvdbdtaod(IX,IY,3 ) .ge. -0.05)  &     
                 land_pixels = land_pixels +1 
               Enddo
          Enddo
                                             
            print *, ' land pixels ', land_pixels    
           
             if(land_pixels > 1 )  then                    
        
         Call  write_Output_forUVtau(Ret_Xtrack,Ret_Lines,Ret_Lat,Ret_Lon,&
                   Ret_SolZen,Ret_View_angle,Ret_View_phi,Ret_solar_phi,uvdbdtaod,&
                   Ret_tau_ocean,Month,interm_file)
           
                        
!          print *, ' Finished writing  Interim file'
!
!******  *********  ML  Land( 0.354.0.388 0.55)
         block
         pycommand = trim(cfg%read_nc_landonly) // " " // trim(cfg%ml_340) // " " // trim(cfg%ml_380)
         pycommand = trim(pycommand) // " " // trim(interm_file)
         call execute_command_line( pycommand ) 
         endblock 
	 
         call read_Output_forUVtau(Ret_Xtrack,Ret_Lines,uvdbdtaod,interm_file)
!          print *,'Finished reading ML output'
         
         Endif
 
!  Merge Land and Ocean Variables together for processing NEARUV variables
     
      CALL   merge_land_ocean_VAR(grn_lwmask,CldMsk_Native_Ocean,dbdt_cld,&
                   Land_sea_flag,Ret_ref_ocean,dbdt_refl,Ret_tau_ocean,uvdbdtaod,& 
                   Ret_Tau_LandOcean,ret_ref_LandOcean,Cldmask_Native_LandOcean,&
                   Cloud_Frac_LandOcean,Ret_Xtrack,Ret_Lines,&
                   Ret_ocean_Quality,Ret_Quality_LandOcean,Ret_Quality_LandOcean_W0,&
                   Ret_CLDFRC_land_DT,Ret_CLDFRC_ocean,Ret_ref_LandOceanwOutUV,&
                   Ret_ref_LandOcean_UV)
                  

   CALL NUV_package_main(cfg,Year, Month, Day, UVtoSWIR_wavelengths,UVtoSWIR_Reflectances, &
          Ret_Lat,Ret_Lon,Ret_SolZen,Ret_View_angle,&
          Ret_View_phi,Ret_solar_phi,Ret_Xtrack,Ret_Lines,& 
          Ret_Small_weighting,Land_sea_flag,Ret_ref_LandOcean,Ret_Tau_LandOcean,&
	  Ret_average_Omega_Ocean_UV,Ret_Index_Height,Cloud_Frac_LandOcean,&
	  NUV_AI, NUV_COD, NUV_CldFrac, UVReflectivity, UVResidue, NUV_SSA, NUV_ALH, &
	  NUV_ACAOD, NUV_AerCorrCOD, NUV_ACAODVsHeight, NUV_AerCorrCODVsHeight, &
	  NUV_FinalAlgorithmFlagsACA, NUV_UncertaintyACAODToSSA,NUV_UncertaintyCODToSSA)
	    

   Call write_Output_merged(Ret_Lat,Ret_Lon,Ret_SolZen,Ret_View_angle,&
     Ret_View_phi,Ret_solar_phi,Ret_Xtrack,Ret_Lines,Ret_Small_weighting,&
	  Land_sea_flag,Ret_ref_LandOceanwOutUV,&
	  Ret_ref_LandOcean_UV,Ret_Tau_LandOcean,&
      Ret_average_Omega_Ocean_UV,Ret_Index_Height,Cloud_Frac_LandOcean,&
	  Ret_Quality_LandOcean,Ret_Quality_LandOcean_W0,& 
          NUV_AI, NUV_COD, NUV_CldFrac, UVReflectivity, UVResidue, &
	  NUV_SSA, NUV_ALH, NUV_ACAOD, NUV_AerCorrCOD,& 
	  NUV_ACAODVsHeight, NUV_AerCorrCODVsHeight, &
	  NUV_FinalAlgorithmFlagsACA, NUV_UncertaintyACAODToSSA,NUV_UncertaintyCODToSSA, out_file)


  Call write_Output_merged_1KM(l1b_nXTrack, l1b_nLines, Latitude, Longitude, &
	       UVAI, Residue_1km, Reflectivity_1km, out_file_1KM)	  
            

! Close the HDF-5 interface
CALL h5close_f(hdferr)


END PROGRAM OCIUAAER
!==============================================================================
!
