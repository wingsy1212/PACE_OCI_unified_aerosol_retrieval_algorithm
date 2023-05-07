Module read_for_nuv_test

USE HDF5
USE H5Util_class
USE OCIUAAER_Config_Module
USE OCIUAAER_L1BModule
USE Proxydata_L1BModule
use calendars, only: gdatetime, gregorian_from_doy

IMPLICIT NONE

contains


SUBROUTINE load_pace_test_data_forNUV(l2_testfile,&
          testl2_Lat,testl2_Lon,testl2_SolZen,testl2_View_angle,&
          testl2_View_phi,testl2_solar_phi,IXDIM,IYDIM,& 
          testl2_Small_weighting,testl2_Land_sea_flag,testl2_Ref_LandOcean,testl2_Tau_LandOcean,&
	  testl2_Ret_average_Omega_Ocean_UV,testl2_Ret_Index_Height,testl2_Cloud_Frac_LandOcean)
  
USE H5read_module

IMPLICIT NONE

Include 'output_Variables.inc'

CHARACTER(LEN=*), INTENT(IN)  :: l2_testfile
INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id
CHARACTER(LEN=128) :: str1

INTEGER(KIND=4) :: errorStatus, STATUS, tmp_rank, index1
INTEGER(KIND=8), dimension(7) :: tmp_dims
INTEGER(KIND=4), ALLOCATABLE :: arr_2di4(:,:), arr_3di4(:,:,:)
INTEGER(KIND=8), ALLOCATABLE :: arr_2di8(:,:), arr_3di8(:,:,:)

INTEGER  IXDIM, IYDIM
REAL(KIND=4), ALLOCATABLE :: tmparr(:,:), tmpref(:,:,:), tmptau(:,:,:)
REAL(KIND=4), ALLOCATABLE :: tmpdtssa(:,:,:)

       	Real(KIND=4)  testl2_Lat(IXRet_B,IYRet_B)
       	Real(KIND=4)  testl2_Lon(IXRet_B,IYRet_B)
       	Real(KIND=4)  testl2_SolZen(IXRet_B,IYRet_B)
       	Real(KIND=4)  testl2_View_angle(IXRet_B,IYRet_B)
       	Real(KIND=4)  testl2_View_phi(IXRet_B,IYRet_B)
       	Real(KIND=4)  testl2_solar_phi(IXRet_B,IYRet_B) 
       	Real(KIND=4)  testl2_Small_weighting(IXRet_B,IYRet_B) 
       	Real(KIND=4)  testl2_Ref_LandOcean(IXRet_B,IYRet_B,Num_Waves)
       	Real(KIND=4)  testl2_Tau_LandOcean(IXRet_B,IYRet_B,Num_Waves)
       	Real(KIND=4)  testl2_Land_sea_Flag(IXRet_B,IYRet_B)
        Real(KIND=4)  testl2_Ret_average_Omega_Ocean_UV(IXRet_B,IYRet_B,(Num_Wave_Land-1)) 
        Real(KIND=4)  testl2_Ret_Index_Height(IXRet_B,IYRet_B)  
        Real(KIND=4)  testl2_Cloud_Frac_LandOcean(IXRet_B,IYRet_B)  
	
status =1

! Open HDF-5 files
 CALL h5fopen_f(l2_testfile, H5F_ACC_RDONLY_F, file_id, errorStatus)

 PRINT *, 'Loading pre-developed file : ',l2_testfile
 CALL get_rank_dims(file_id, '/geophysical_data/Mean_Reflectance', tmp_rank, tmp_dims)
 !PRINT*,'L2 TEST file = ',l2_testfile
 !PRINT*,'L2 TEST file REF rank:',tmp_rank
 !PRINT*,'L2 TEST file REF shape:',tmp_dims(1:tmp_rank)

!
!Read input variables now.
IXDIM=400
IYDIM=404
ALLOCATE( tmparr(IXDIM,IYDIM), &
          tmpref(Num_Waves,IXDIM,IYDIM), &
	  tmptau(Num_Waves,IXDIM,IYDIM), &
	  tmpdtssa(Num_Wave_Land-1,IXDIM,IYDIM) )


errorStatus = read_h5_2dr4(l2_testfile, '/geolocation_data/latitude', (/IXDIM,IYDIM/), tmparr)
testl2_lat(:, :) = -9999.
testl2_lat(1:IXDIM, 1:IYDIM) = tmparr

errorStatus = read_h5_2dr4(l2_testfile, '/geolocation_data/longitude',(/IXDIM,IYDIM/), tmparr)
testl2_lon(:, :) = -9999.
testl2_lon(1:IXDIM, 1:IYDIM) = tmparr

errorStatus = read_h5_2dr4(l2_testfile, '/geolocation_data/solar_zenith_angle', (/IXDIM,IYDIM/), tmparr)
testl2_SolZen(:, :) = -9999.
testl2_SolZen(1:IXDIM, 1:IYDIM) = tmparr 

errorStatus = read_h5_2dr4(l2_testfile, '/geolocation_data/solar_azimuth_angle', (/IXDIM,IYDIM/), tmparr)
testl2_solar_phi(:, :) = -9999.
testl2_solar_phi(1:IXDIM, 1:IYDIM) = tmparr

errorStatus = read_h5_2dr4(l2_testfile, '/geolocation_data/sensor_zenith_angle', (/IXDIM,IYDIM/), tmparr)
testl2_View_angle(:, :) = -9999.
testl2_View_angle(1:IXDIM, 1:IYDIM) = tmparr 

errorStatus = read_h5_2dr4(l2_testfile, '/geolocation_data/sensor_azimuth_angle', (/IXDIM,IYDIM/), tmparr)
testl2_View_phi(:, :) = -9999.
testl2_View_phi(1:IXDIM, 1:IYDIM) = tmparr


!geophysical data
errorStatus = read_h5_3dr4(l2_testfile, '/geophysical_data/Aerosol_Optical_Depth', &
                           (/Num_Waves,IXDIM,IYDIM/), tmptau)
testl2_Tau_LandOcean( :, :, :) = -9999.
testl2_Tau_LandOcean(1:IXDIM, 1:IYDIM, 1) = tmptau(1,:,:) 
testl2_Tau_LandOcean(1:IXDIM, 1:IYDIM, 2) = tmptau(2,:,:) 
testl2_Tau_LandOcean(1:IXDIM, 1:IYDIM, 3) = tmptau(3,:,:) 
testl2_Tau_LandOcean(1:IXDIM, 1:IYDIM, 4) = tmptau(4,:,:)
testl2_Tau_LandOcean(1:IXDIM, 1:IYDIM, 5) = tmptau(5,:,:)
testl2_Tau_LandOcean(1:IXDIM, 1:IYDIM, 6) = tmptau(6,:,:) 
testl2_Tau_LandOcean(1:IXDIM, 1:IYDIM, 7) = tmptau(7,:,:) 
testl2_Tau_LandOcean(1:IXDIM, 1:IYDIM, 8) = tmptau(8,:,:)
testl2_Tau_LandOcean(1:IXDIM, 1:IYDIM, 9) = tmptau(9,:,:)

errorStatus = read_h5_3dr4(l2_testfile, '/geophysical_data/Mean_Reflectance', &
                           (/Num_Waves,IXDIM,IYDIM/), tmpref)
testl2_Ref_LandOcean( :, :, :) = -9999.
testl2_Ref_LandOcean(1:IXDIM, 1:IYDIM, 1) = tmpref(1,:,:) 
testl2_Ref_LandOcean(1:IXDIM, 1:IYDIM, 2) = tmpref(2,:,:) 
testl2_Ref_LandOcean(1:IXDIM, 1:IYDIM, 3) = tmpref(3,:,:)
testl2_Ref_LandOcean(1:IXDIM, 1:IYDIM, 4) = tmpref(4,:,:) 
testl2_Ref_LandOcean(1:IXDIM, 1:IYDIM, 5) = tmpref(5,:,:)
testl2_Ref_LandOcean(1:IXDIM, 1:IYDIM, 6) = tmpref(6,:,:) 
testl2_Ref_LandOcean(1:IXDIM, 1:IYDIM, 7) = tmpref(7,:,:) 
testl2_Ref_LandOcean(1:IXDIM, 1:IYDIM, 8) = tmpref(8,:,:) 
testl2_Ref_LandOcean(1:IXDIM, 1:IYDIM, 9) = tmpref(9,:,:) 

errorStatus = read_h5_3dr4(l2_testfile, '/geophysical_data/DT_AerosolSingleScattAlbedo', &
                           (/Num_Wave_Land-1,IXDIM,IYDIM/), tmpdtssa)
testl2_Ret_average_Omega_Ocean_UV( :, :, :) = -9999.
testl2_Ret_average_Omega_Ocean_UV(1:IXDIM, 1:IYDIM, 1) = tmpdtssa(1,:,:) 
testl2_Ret_average_Omega_Ocean_UV(1:IXDIM, 1:IYDIM, 2) = tmpdtssa(2,:,:) 
testl2_Ret_average_Omega_Ocean_UV(1:IXDIM, 1:IYDIM, 3) = tmpdtssa(3,:,:)

errorStatus = read_h5_2dr4(l2_testfile, '/geophysical_data/DT_AerosolLayerHeight', (/IXDIM,IYDIM/), tmparr)
testl2_Ret_Index_Height(:, :) = -9999.
testl2_Ret_Index_Height(1:IXDIM, 1:IYDIM) = tmparr

errorStatus = read_h5_2dr4(l2_testfile, '/geophysical_data/Aerosol_Cld_Fraction_Land_Ocean', (/IXDIM,IYDIM/), tmparr)
testl2_Cloud_Frac_LandOcean(:, :) = -9999.
testl2_Cloud_Frac_LandOcean(1:IXDIM, 1:IYDIM) = tmparr

errorStatus = read_h5_2dr4(l2_testfile, '/geophysical_data/Land_Sea_Flag', (/IXDIM,IYDIM/), tmparr)
testl2_Land_sea_flag(:, :) = -9999.
testl2_Land_sea_flag(1:IXDIM, 1:IYDIM) = tmparr

errorStatus = read_h5_2dr4(l2_testfile, '/geophysical_data/Optical_Depth_Ratio_Small_Ocean_used', &
                        (/IXDIM,IYDIM/), tmparr)
testl2_Small_weighting(:, :) = -9999.
testl2_Small_weighting(1:IXDIM, 1:IYDIM) = tmparr


CALL h5fclose_f(file_id, errorStatus)

!PRINT *,'DiRECT READ L2_Merged = ',testl2_Lat(225,210),testl2_Lon(225,210),testl2_SolZen(225,210), &
!         testl2_View_angle(225,210),testl2_View_phi(225,210),testl2_solar_phi(225,210),&
!	 testl2_Ref_LandOcean(225,210,1:5), testl2_Tau_LandOcean(225,210,1:5)
!CALL EXIT(1)

IF (ALLOCATED(arr_2di4)) DEALLOCATE(arr_2di4)
IF (ALLOCATED(arr_3di4)) DEALLOCATE(arr_3di4)

RETURN

END SUBROUTINE load_pace_test_data_forNUV
!
!
!
SUBROUTINE read_nativeL1b_for_NUVtest(cfg, Year, Month, Day, &
				      UVtoSWIR_wavelengths, & 
				      UVtoSWIR_Reflectances)

IMPLICIT NONE

Include '../../Main/common_l1b_var.inc'
 
TYPE(ociuaaer_config_type),    INTENT(IN)  :: cfg
INTEGER(KIND=4),               INTENT(OUT) :: Day, Month, Year
REAL(KIND=4),DIMENSION(:,:,:),ALLOCATABLE, INTENT(OUT)  :: UVtoSWIR_Reflectances
REAL(KIND=4),DIMENSION(:),ALLOCATABLE,     INTENT(OUT)  :: UVtoSWIR_wavelengths

TYPE(gdatetime)  	:: gdt1
INTEGER(KIND=4)	        :: dayofyear, STATUS
INTEGER(KIND=4), PARAMETER              :: UVtoSWIR_nWavel=14
REAL(KIND=4),DIMENSION(UVtoSWIR_nWavel) :: waveTemp = [340.0, 354.0, 388.0, &
                         412.0, 488.0, 550.0, 670.0,680.0,688.0,860.0, &
                        1250.0,1378.0,1615.0,2260.0]

   ALLOCATE(UVtoSWIR_wavelengths(UVtoSWIR_nWavel), stat=STATUS) 
   UVtoSWIR_wavelengths = waveTemp(1:UVtoSWIR_nWavel)

   IF (cfg%input_l1file /= 'NULL') THEN
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

     ! Extract year, month, day.
     Read(l1b_date_time_str(1:4), '(I4)') Year
     Read(l1b_date_time_str(5:6), '(I2)') Month
     Read(l1b_date_time_str(7:8), '(I2)') Day

   ELSE 
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

     !Extract year, month, day.
     Read(l1b_date_time_str(1:4), '(I4)') Year
     Read(l1b_date_time_str(5:7), '(I3)') dayofyear
     gdt1  = gregorian_from_doy(year, dayofyear)
     Month = gdt1%month
     Day = gdt1%mday
  
   ENDIF

RETURN

END SUBROUTINE read_nativeL1b_for_NUVtest

End Module read_for_nuv_test
