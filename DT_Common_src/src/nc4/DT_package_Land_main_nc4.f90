 module compute_DT_land_main
      implicit none
contains  
!----------------------------------------------------------------------------------
! CREATION HISTORY:
!       Written  by:    Shana Mattoo
!                       shana.mattoo@nasa.gov
!-----------------------------------------------------------------------------------------------
     
     
     
      Subroutine DT_package_Land_main(Latitude_in,Longitude_in,SolarZenithAngle,&
                         SolarAzimuthAngle,&
                          ViewingZenithAngle,ViewingAzimuthAngle,TerrainHeight,&
                          UVtoSWIR_Reflectances,UVAI,Year,Month,Day,&
                          Land_Sea_Flag,l1b_nXTrack,l1b_nLines,save_ref_Land_PACE,& 
                          save_corrected,save_Lat,Save_lon,CldMsk_500_Land,&
                          Save_land_ocean_Quality,Save_dust_weighting,Save_CLDFRC_land,&
                          Ret_Xtrack,Ret_Lines,anc_file)
                         
                       
!      USE write_pace_dt
      USE gather_l1b_OCI_data
      USE read_GMAO_anc 
      USE convert_rad 
      USE compute_Temperature_from_rad
      USE compute_Gascorrection_twoways
      USE Apply_Trans_Two_way
      USE compute_water_land
      USE snow_flag
      USE SET_INDEX_L1b
      USE compute_center_GEOLOC_ANGLE
      USE compute_center_Glint_Angle  
      USE compute_Boxes_and_edges
      USE compute_theshold_water_land_coastal  
      USE cld_flags_forLand
      USE Fill_Output_SDS_land 
      
      IMPLICIT  NONE
      include 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'  
      INCLUDE 'Set_Array_dimension_Land.inc'
      INCLUDE 'Save_data_declare.inc'
      INCLUDE 'Save_data_UV_declare.inc'
      SAVE 

      CHARACTER(255) :: anc_file
         
      
                      
                 
                  
 !-------------------------------------------------------------------------------------------------------------------               
!    Reading  VIIRS    MOD35 like cloud mask and L1b data and Iband for cloud masking
!-------------------------------------------------------------------------------------------------------------------                
         
   
            CALL  open_l1b_OCI_data(Latitude_in,Longitude_in,SolarZenithAngle,SolarAzimuthAngle,&
                           ViewingZenithAngle,ViewingAzimuthAngle,TerrainHeight,&
                           UVtoSWIR_Reflectances,&
                            l1b_nLines,l1b_nXTrack,OLAT,&
                            OLON,OSOLA,OVIEW,OAZIM_View,OAZIM_Solar,&
                            OLand_ocean,OHeight,OW470,OW550,OW659,OW860,&
                            OW124,OW138,OW164,OW213,OW1100,OW659_Hres,OW354,OW388,OW412)
                              
                              
!-------------------------------------------------------------------------------------------------------------------                                
!   convert Radiance into Reflectance units. Taking this off o keep it compatable with Synthetic data reader
!-------------------------------------------------------------------------------------------------------------------                
!             Call convert_rad_ref(l1b_nLines,l1b_nXTrack,OSOLA,OW470,OW550,OW659,OW860,OW124,OW138,&
!                    OW164,OW213,OW354,OW388,OW412) 
!-------------------------------------------------------------------------------------------------------------------                             
!   Set Box for retrievals 
!-------------------------------------------------------------------------------------------------------------------                     
                  call set_Boxes_and_edges(l1b_nLines,l1b_nXTrack,Ret_Xtrack,Ret_Lines)
                                          
                         
!-------------------------------------------------------------------------------------------------------------------
! Initialize array for OUtput
!-------------------------------------------------------------------------------------------------------------------                                  
          Call init_save_array(Save_Lat,Save_Lon,Save_SolZen,Save_View_angle,Save_View_phi,Save_solar_phi,& 
       Save_Scatt_angle,Save_Glint_angle,Save_Sea_Land_Flag,Save_Cld_mask_500,Save_cloud_distance_500,Save_land_ocean_Quality, &
       Save_Image_land_ocean_tau,Save_Aerosol_Type,Save_Fitting_Error_Land,Save_Surface_Reflectance_Land,Save_corrected,&
       Save_land_ocean_tau,Save_dust_weighting,Save_NUMPIXELS_land,Save_ref_Land,Save_SDS_Land,&
       Save_SDS_mass_conc_land,Save_CLDFRC_land,Save_average_Tau_Ocean,Save_SDSTAUS_average,Save_SDSTAUL_average,&
       Save_SDS_mass_conc_ocean,Save_CLDFRC_Ocean,Save_effrad,Save_PSML003_Ocean,Save_ASSY,Save_BACK,&
       Save_SDS_angs_coeff1,Save_SDS_angs_coeff2,Save_Least_error,Save_small_weighting,Save_SDS_AOT_model,Save_NUMPIXELS_Ocean,&
       Save_ref_Ocean,Save_STD_Ocean,Save_Wspeed_Ocean,Save_Topo_Height,Save_land_ocean_Q_Ret,&
       l1b_nXTrack  ,Ret_Lines,set_Counter_Ocean, Set_Counter_Land,num_joint,num_joint2,&
        num_Sea_land,set_counter_for_Gread,set_counter_for_anc,Num_water,&
        num_land,Data_Size,START_line,END_line,l1b_nXTrack,l1b_nLines,Ret_Lines,Save_Dust_Flag_ocean,&
        Save_Lat_UV,Save_Lon_UV,&
        save_SolZen_UV,save_View_angle_UV,save_View_phi_UV,save_solar_phi_UV,Save_average_Tau_Ocean_UV,&
        Save_average_Omega_Ocean_UV,Save_mode_F_FUV,Save_mode_C_FUV,Save_Small_weighting_FUV,&
        Save_Index_Albedo,Save_Index_Height,Save_ref_allwav_uv,Save_ref_allwav_PACE,&
         Save_average_Tau_PACE)      
           set_counter_for_anc =0
 
 !-------------------------------------------------------------------------------------------------------------------
 ! Compute cloud mask for Land
 !-------------------------------------------------------------------------------------------------------------------
       CALL CldMsk_Land_updated(Data_Size,OW470,OW124, &
            OW138,CldMsk_250,CldMsk_500_Land,CldMsk_cirrus,& 
            quality_cirrus,Aerosol_Cldmask_land,START_line,END_line,Ret_Lines)
            
              
!               print*,'cloud masking done for Land '  
!  Start processing data ..............              
                 
          DO   Iscan = 1, Ret_Lines   
             Set_Counter_Ocean_cloud=0 
             IF(Iscan  .eq.1) then
                  START_line=0
                  END_line=ILINE 
             ELSE
                  START_line=START_line+Iline
                  END_line=END_line+ILINE 
            ENDIF  
            
             
 !-------------------------------------------------------------------------------------------------------------------            
 !   Moving data into proper retrieval boxes for processing
 !-------------------------------------------------------------------------------------------------------------------      
       CALL Move_Data_Reso(OW659,OW860,OW470,OW550,OW124,OW138,OW164,&
          OW213,OW1100,W470_SYN,W550_SYN,W659_SYN,W865_SYN,W124_SYN,&
          W164_SYN,W213_SYN,W138_SYN,W1100_SYN,Land_Sea_Flag, &
          OLAT,OLON,OSOLA,OVIEW,OAZIM_View,OAZIM_Solar,OHeight,&
          LandSea_Flag_new,Lat,Lon,SatZen,SatAz,SolZen,SolAz,Height,&
          START_line,END_line,CldMsk_500_Ocean,CldMsk_500_Land,&
          savecldmask,New_CldMsk_500_Ocean,New_CldMsk_500_Land,New_savecldmask,&
          l1b_nXTrack,quality_cirrus,New_Q_cirrus,OW412,W412_SYN,OW8p5,W8p5_SYN,SD_3by3,&
          SD_for_Dust,W354_SYN,W388_SYN,OW354,OW388)
          
            
                     
 !------------------------------------------------------------------------------------------------------------------            
 !   Moving data into proper retrieval boxes for processing
 !-------------------------------------------------------------------------------------------------------------------           
      
        DO   IDATA = 1,Ret_Xtrack 
!------------------------------------------------------------------------------------------------------------------       
!   this subroutine sets index to read Iline * Iline box( e.g 10* 10 for MODIS or 8*8 Viirs etc)  
!------------------------------------------------------------------------------------------------------------------              
        CALL SET_INDEX(START_500,END_500,START_250,END_250,&
              START_1KM,END_1KM,IDATA)  
!------------------------------------------------------------------------------------------------------------------                    
!this subroutine computes center og Geolocation and angles for the box used for retrieval. 
!------------------------------------------------------------------------------------------------------------------                 
       CALL GEOLOC_ANGLE(LAT,LON,SatZen,SatAz,SolZen,&
           SolAz,Height,MTHET0,MTHET,MPHI0,MPHI,MDPHI,MSCATT,MHGHT,&
           START_1KM,Lat_center,Lon_center,iscan,idata) 
           
!------------------------------------------------------------------------------------------------------------------               
!   compute glint Angle for center Angles for the  for the box used for retrieval.  
!------------------------------------------------------------------------------------------------------------------                   
          CALL COMPUTE_GLINTANGLE(MTHET0,MTHET,MDPHI,&
                  GLINT_ANGLE,QA_Flag_Ocean)  
                  
!------------------------------------------------------------------------------------------------------------------                       
 ! computing water and land and in between 
 !------------------------------------------------------------------------------------------------------------------         
         
            Call  Get_water_land(START_1KM,END_1KM,&
                 LandSea_Flag_new,land,water) 
                
                
                
!------------------------------------------------------------------------------------------------------------------                  
!  set threshold for Coastal pixels ( 20% )   
!------------------------------------------------------------------------------------------------------------------        
            call threshold_water_land_coastal(water,Water_Pixels,Pure_Land,Land_Pixels,&
            Sea_Land_Flag,save_Sea_Land_Flag,iscan,idata) 
            
 !------------------------------------------------------------------------------------------------------------------ 
 !  computing Quality cirrus and cloud fraction from land cloud mask
 !------------------------------------------------------------------------------------------------------------------                         
             Call cld_flags_land(CldMsk_cirrus,New_CldMsk_500_Land,&
             CldMsk_250,cloud_num,cloud_num_land,START_1KM,END_1KM,&
             New_Q_cirrus,Ret_Quality_cirrus,checking_total_cloud)                  
                          
 !------------------------------------------------------------------------------------------------------------------ 
 !   Threshold for Solar Zenith angle
 !------------------------------------------------------------------------------------------------------------------                        
               if  (MTHET0  .le.   MAXMTHET0 )  then    
           
!  Compute Gas correction   factor , apply to wavelengths
        
              set_counter_for_anc=set_counter_for_anc+1
             
                pwat=0
                ozone=0
                ugrd=0
                ugrd=0
                RTN_NCEP =0   
 

!------------------------------------------------------------------------------------------------------------------ 
!   Read ugrd,vgrd, Ozone and water from GMAO ancillary data. All variables are interpolated to exact Lat and Lon center
!------------------------------------------------------------------------------------------------------------------                           
           
                CALL Get_An_GMAO(Lat_center,Lon_center,ugrd,vgrd,pwat,&
                ozone,skinTemp,set_counter_for_anc,RTN_NCEP,anc_file)
  
         
 
          set_counter_for_Gread=set_counter_for_Gread+1
!------------------------------------------------------------------------------------------------------------------         
!Compute Gas correction and water vapor
!------------------------------------------------------------------------------------------------------------------      
       CALL compute_Gascorrection_nc4(pwat,ozone,&
        MTHET0,MTHET,set_counter_for_Gread,Multi_factor,RTN_NCEP) 
         
!-----------------------------------------------------------------------------------------------------------------         
! Apply Gas correction and water vapor correction  to measured wavelengths
!------------------------------------------------------------------------------------------------------------------            
    
               
                 
       CALL  Trans_Two_way(START_1KM,END_1KM,W659_SYN,&
             W865_SYN,W470_SYN,W550_SYN,W124_SYN,W164_SYN,W213_SYN,W412_SYN,&
             Multi_factor) 
         
           
             
             
!-----------------------------------------------------------------------------------------------------------------         
!   Compute Temperature for Channels 8.5 and 11.00 microns
!------------------------------------------------------------------------------------------------------------------                
            
!        call compute_temp_from_rad(START_1KM,END_1KM,&
!             W1100_SYN,W8p5_SYN,W1100_Temp,W8p5_Temp) 
!-----------------------------------------------------------------------------------------------------------------         
!  Compute snow flag 
!------------------------------------------------------------------------------------------------------------------                   
     
        call compute_snow_flag(START_1KM,END_1KM,W865_SYN,&
                    W164_SYN,skinTemp,SnowMsk_Ratio)            
             
!------------------------------------------------------------------------------------------------------------------                
!   We use PROCESS_ocean if there are 100% pixels of water  
!------------------------------------------------------------------------------------------------------------------                
                
                
                IF(Land .GE. 1) THEN   
!-----------------------------------------------------------------------------------------------------------------  
! handler_for_data_read_land sets handlers for reading LUT and other ancillary data for Land
!------------------------------------------------------------------------------------------------------------------                          
              Call  handler_for_data_read_land(HANDLE_LUT466,HANDLE_LUT553,HANDLE_LUT644,&
              HANDLE_LUT213,HANDLE_LUTMAP,handle_Urban_Table_10km)   
!-----------------------------------------------------------------------------------------------------------------                
! IF  all pixels are not detected as water pixels, proceed to land algorithm       
!----------------------------------------------------------------------------------------------------------------- 
               quality_land=1  
             Set_Counter_Land=Set_Counter_Land+1 
!-----------------------------------------------------------------------------------------------------------------  
!  Go to PROCESS_Land subroutine to Retrieve land  only
!------------------------------------------------------------------------------------------------------------------                           
        
        CALL PROCESS_Land(HANDLE_LUT466,HANDLE_LUT553,HANDLE_LUT644,&
        HANDLE_LUT213,HANDLE_LUTMAP,MONTH,ISCAN,IDATA,MTHET0,&
        MTHET,MDPHI,MHGHT,Lat_center,Lon_center,START_500,END_500,& 
        START_250,END_250,START_1KM,END_1KM,W354_SYN,W388_SYN,W470_syn,W550_SYN,& 
        W659_syn,W865_syn,W124_SYN,W164_SYN,W213_syn,&
        CldMsk_250,Set_Counter_Land,QA_Flag_Land,Success_Ret_Land,&
        Fail_Ret_Land,SDSLAT,SDSLON,SDS_MTHET0,SDS_MTHET,SDS_MPHI,&
        SDS_Tau_Land_Ocean,New_CldMsk_500_land,SDS_Tau_Land_Ocean_img,&
        SDS_Aerosol_Type,SDS_SCAT_ANGLE_land,SDS_mass_conc_land,&
        SDS_angs_coeff_land,SDS_CLDFRC_land,SDS_dust_weighting,&
        SDS_est_uncer,SDS_RefF_land,&
        SDS_TranF_land,SDS_NUMPIXELS_land,SDSTAU_corrected,&
        SDS_ref_land,SDS_ref_STD_land,SDS_QCONTROL_land,&  
        G_factor,quality_land,Ret_Quality_cirrus,cloud_num_land,&
        SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,Qcontrol_special_land,&
        SDSTAU_corrected_213,Quality_flag_forJoint,SDSTAU_small_land,&
        Save_index_land,W412_SYN,W443_SYN,Optical_depth_land,&
        SnowMsk_Ratio,l1b_nXTrack,handle_Urban_Table_10km,ETA)  
        
        
           
!-----------------------------------------------------------------------------------------------------------------  
!  Fill land  SDS's ( loop around Iscan and Idata ) Subroutine also degrades the quality if water pixels 
!   greater than 20% or coastal pixels gt 50% 
!------------------------------------------------------------------------------------------------------------------               
        
        call Fill_Output_Arrays_land(water,coastal_Pix,Quality_flag_forJoint,&
            Save_CLDFRC_land,SDS_CLDFRC_land,SDSTAU_corrected,Save_land_ocean_tau,&
            SDS_NUMPIXELS_land,save_numpix_Land_Ocean,Save_Image_land_ocean_tau,&
            Save_Topo_Height,MHGHT, Save_land_ocean_Quality,Save_Aerosol_Type,&
            SDS_Aerosol_Type,Save_Fitting_Error_Land,SDS_Fitting_Error_Land,&
            Save_dust_weighting,SDS_dust_weighting,Save_SDS_mass_conc_land,&
            SDS_mass_conc_land,SDS_ref_land,SDS_ref_STD_land,save_ref_Land,&
            save_SDS_Land,save_Surface_Reflectance_Land,SDS_Surface_Reflectance_Land,&
            save_corrected,SDSTAU_corrected_213,save_NUMPIXELS_land,idata,iscan,&
            Save_index_Land_Ocean,Save_index_land,Save_land_ocean_Q_Ret,&
            Qcontrol_special_land,save_ref_Land_PACE)
 
             
!------------------------------------------------------------------------------------------------------------------                         
! endif for  pixel Land at ocean
!------------------------------------------------------------------------------------------------------------------     
               
                    ENDIF 
!------------------------------------------------------------------------------------------------------------------                 
! Endif for solar zenith angle of 84
!------------------------------------------------------------------------------------------------------------------     
            
             ENDIF 
 
 
           Save_Lat(Idata,iscan)=Lat_center  
           Save_Lon(Idata,iscan)=Lon_center
           save_SolZen(Idata,iscan)=MTHET0
           save_View_angle(Idata,iscan)=MTHET
           save_View_phi(Idata,iscan)=MPHI 
           save_solar_phi(Idata,iscan)= MPHI0
           Save_Scatt_angle(Idata,iscan) = MSCATT
           Save_Glint_angle(Idata,iscan)=Glint_angle 
           Save_smw_land(Idata,iscan)=ETA
             
         
!------------------------------------------------------------------------------------------------------------------              
!   this enddo is for  Idata         Idata   
!------------------------------------------------------------------------------------------------------------------  
            Enddo     
 
!------------------------------------------------------------------------------------------------------------------ 
!   this enddo is for  Iscan    
!------------------------------------------------------------------------------------------------------------------         
            Enddo  
!------------------------------------------------------------------------------------------------------------------------    
! call to write  output (full granule or Full disk of data)
!------------------------------------------------------------------------------------------------------------------   
!                
               go to 8000
                 
          call  make_nc_Land(Save_Lat,Save_Lon,save_SolZen,save_View_angle,&
                  save_View_phi,save_solar_phi,Save_Scatt_Angle,&
                  Save_Glint_angle,save_Sea_Land_Flag,Save_Cld_mask_500,& 
                  Save_cloud_distance_500,Save_land_ocean_Quality,&
                  Save_land_ocean_tau,Save_Image_land_ocean_tau,& 
                  Save_Surface_Reflectance_Land,Save_corrected,&
                  Save_dust_weighting,Save_NUMPIXELS_land,Save_ref_Land,&
                  Save_SDS_Land,Save_SDS_mass_conc_land,&
                  Save_CLDFRC_land,Save_average_Tau_Ocean,&
                  Save_SDSTAUS_average,Save_SDSTAUL_average,&
                  Save_SDS_mass_conc_ocean,Save_CLDFRC_Ocean,&
                  Save_effrad,Save_PSML003_Ocean,Save_ASSY,Save_BACK,&
                  Save_SDS_angs_coeff1,Save_SDS_angs_coeff2,&
                  Save_Least_error,Save_small_weighting,&
                  Save_SDS_AOT_model,Save_NUMPIXELS_Ocean,&
                  Save_ref_Ocean,Save_STD_Ocean,Save_Wspeed_Ocean,&
                  Save_Topo_Height,Save_land_ocean_Q_Ret,&
                  Save_Aerosol_Type,Save_Fitting_Error_Land,Save_ave_clddis_land_ocean,&
                  Ret_Xtrack,Ret_Lines,ipart,num_joint)
                   
                    
                   
8000            continue   


                      
! End for programme        
        Return 
        end  subroutine DT_package_Land_main
           end module  compute_DT_land_main
           
           
           
           
               
     
 
