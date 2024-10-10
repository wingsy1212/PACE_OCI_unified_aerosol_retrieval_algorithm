module compute_Ocean_main
      implicit none
contains  
!----------------------------------------------------------------------------------
! CREATION HISTORY:
!       Written  by:    Shana Mattoo
!                       shana.mattoo@nasa.gov
!-----------------------------------------------------------------------------------------------
     
     
     
      Subroutine DT_package_Ocean_main(Latitude_in,Longitude_in,SolarZenithAngle,&
                         SolarAzimuthAngle,&
                        ViewingZenithAngle,ViewingAzimuthAngle,TerrainHeight,&
                        UVtoSWIR_Reflectances,UVAI,Year,Month,Day,&
                        Land_Sea_Flag,l1b_nXTrack,l1b_nLines,&
 !     Output Variables
                    Ret_Xtrack,Ret_Lines,Save_Lat,Save_Lon,save_SolZen,save_View_angle,&
                    save_View_phi,save_solar_phi,Save_average_Tau_Ocean_UV,& 
                    Save_Small_weighting_FUV,Save_ref_allwav_uv,save_Sea_Land_Flag,&
                    Save_average_Omega_Ocean_UV,CldMsk_500_Ocean,Save_Index_Height,&
                    Save_land_ocean_Quality,Save_CLDFRC_Ocean,anc_file)
                     
     
      USE write_pace_dt_ocean
      USE write_uv 
      USE gather_l1b_OCI_data
      USE read_GMAO_anc 
      USE convert_rad 
      USE compute_Temperature_from_rad
      USE compute_Gascorrection_twoways
      USE Apply_Trans_Two_way
      USE compute_water_land
      USE SET_INDEX_L1b
      USE compute_center_GEOLOC_ANGLE
      USE compute_center_Glint_Angle 
      USE compute_indx_wspeed_forLUT
      USE compute_wspeed_spatial 
      USE compute_Boxes_and_edges
      USE Fill_Output_SDS_Ocean 
      USE Fill_Output_SDS_Ocean_UV
      
      IMPLICIT  NONE
      include 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'  
      INCLUDE 'Set_Array_dimension.inc'
      INCLUDE 'Save_data_declare.inc'
      INCLUDE 'Save_data_UV_declare.inc'
       
      SAVE 
        
      CHARACTER(255) :: anc_file
      
                 Anc_flag ='GMAO'
                
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
!  convert Radiance into Reflectance units. Taking this off o keep it compatable with Synthetic data reader
!-------------------------------------------------------------------------------------------------------------------                
!            Call convert_rad_ref(l1b_nLines,l1b_nXTrack,OSOLA,OW470,OW550,OW659,OW860,OW124,OW138,&
!                   OW164,OW213,OW354,OW388,OW412) 
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
        l1b_nXTrack,Ret_Lines,set_Counter_Ocean, Set_Counter_Land,num_joint,num_joint2,&
        num_Sea_land,set_counter_for_Gread,set_counter_for_anc,Num_water,&
        num_land,Data_Size,START_line,END_line,l1b_nXTrack,l1b_nLines,Ret_Lines,Save_Dust_Flag_ocean,&
        Save_Lat_UV,Save_Lon_UV,&
        save_SolZen_UV,save_View_angle_UV,save_View_phi_UV,save_solar_phi_UV,Save_average_Tau_Ocean_UV,&
        Save_average_Omega_Ocean_UV,Save_mode_F_FUV,Save_mode_C_FUV,Save_Small_weighting_FUV,&
        Save_Index_Albedo,Save_Index_Height,Save_ref_allwav_uv,Save_ref_allwav_PACE,&
        Save_average_Tau_PACE)      
!-------------------------------------------------------------------------------------------------------------------          
!  Compute cloud mask for Ocean
!-------------------------------------------------------------------------------------------------------------------        

 5000    continue 
       CALL  cldMsk_Ocean_updated(OW470,OW550,OW659_Hres,&
             OW138,OW124,CldMsk_500_Ocean, &
             data_size,savecldmask,START_line,END_line,Ret_Lines) 
            
 
 !-------------------------------------------------------------------------------------------------------------------
 ! Compute cloud mask for Land
 !-------------------------------------------------------------------------------------------------------------------
       
               
!               print*,'cloud masking done for Ocean '
               
                
      
                 
         CALL Compute_Sd_for_Dust(START_line,END_line,Ret_Lines,OW860,SD_3by3,data_size) 
        
       
          
         
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
 !  SAVE FLAG FOR WRITTING
 !------------------------------------------------------------------------------------------------------------------         
                    
               IF(WATER .GE. (iline*iline)) THEN  
                save_Sea_Land_Flag (Idata,iscan) =0
               elseif(land .GE. 1) then
               save_Sea_Land_Flag(Idata,iscan) = 1
               Else
              save_Sea_Land_Flag(Idata,iscan) = -9990
              ENDIF        
 !--------- --------------------------------------------------------------------------------------------------------- 
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
                if(Iscan .eq. 11 .and. IDATA  .eq. 8)&
                print*,'Lat_center,Lon_center,ugrd,vgrd,pwat,ozone,skinTemp',&
                Lat_center,Lon_center,ugrd,vgrd,pwat,ozone,skinTemp
  
         
               
!------------------------------------------------------------------------------------------------------------------         
! choosing index for wind speed computed from GDAS data. It will be used to read LUT for Wind Speed WSPEED
!------------------------------------------------------------------------------------------------------------------               
         Call indx_wspeed(ugrd,vgrd,Indx_wspeed1,Indx_wspeed2,&
           WSPEED,Wind,RTN_NCEP) 
            
 
        set_counter_for_Gread=set_counter_for_Gread+1
!------------------------------------------------------------------------------------------------------------------         
!Compute Gas correction and water vapor
!------------------------------------------------------------------------------------------------------------------      
       CALL compute_Gascorrection(pwat,ozone,&
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
!             
!------------------------------------------------------------------------------------------------------------------                
!   We use PROCESS_ocean if there are 100% pixels of water  
!------------------------------------------------------------------------------------------------------------------                
                
                
                IF(WATER .GE. (iline*iline)) THEN  
!-----------------------------------------------------------------------------------------------------------------  
!  handler_for_data_read_ocean sets handlers for reading LUT and other ancillary data for ocean
!------------------------------------------------------------------------------------------------------------------                  
             Call handler_for_data_read_ocean(HANDLE_S,HANDLE_L,HANDLE_Ext_554_O)
         
                
             Set_Counter_Ocean=Set_Counter_Ocean+1
             Set_Counter_Ocean_cloud=Set_Counter_Ocean_cloud+1 
!-----------------------------------------------------------------------------------------------------------------  
!  Go to PROCESS_ocean subroutine to Retrieve Ocean only
!------------------------------------------------------------------------------------------------------------------                   
            CALL PROCESS_ocean(HANDLE_S,HANDLE_L,&
              ISCAN,IDATA,MTHET0,MTHET,MDPHI,START_500,&
              END_500,START_250,END_250,START_1KM,END_1KM,&
              W659_SYN,W865_SYN,W470_SYN,W550_SYN,W124_SYN,&
              W164_SYN,W213_SYN,W412_SYN,W443_SYN,W8p5_Temp,W1100_Temp,&
              Sunglint_Flag,Set_Counter_Ocean,&
              SDSTAU_best,SDSTAUS_best,SDSTAUB_best,&
              SDSTAU_average,SDSTAUS_average,SDSTAUB_average,&
              SDS_Least_error,SDS_small_weighting,SDS_sol_INDX_small,&
              SDS_sol_INDX_large,SDS_ref,SDS_ref_STD,SDSASSY_best,&
              SDSASSY_average,SDSBACK_best,SDSBACK_average,SDS_effrad,&
              SDS_RefF_best,SDS_ccn,SDS_mass_conc,&
              SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,&
              SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_SCAT_ANGLE_OCEAN,&
              SDS_AOT_model,SDS_CLDFRC_ocean,&
              Set_Counter_Ocean_cloud,QA_Flag_ocean,GLINT_ANGLE,&
              SDS_angs_coeff1,SDS_angs_coeff2,SDS_Tau_Land_Ocean,&
              Qcontrol_special,SDS_correc_small_weighting,SDS_Tau_Land_Ocean_img,&
              High_Cloud_Flag_500,W138_SYN,data_size,New_CldMsk_500_Ocean, &   
              New_savecldmask,Save_index_ocean,Quality_to_pass,Indx_wspeed1,Indx_wspeed2,&
              WSPEED,Wind,AVE_ARRAy,New_Optical_depth_Ocean,Ret_Xtrack ,&
              HANDLE_Ext_554_O,ipart,SD_for_Dust,Dust_flag_10KM,&
              Ext_554_small,Ext_554_large,WAVE,W354_SYN,W388_SYN,REFW354,REFW388) 
!-----------------------------------------------------------------------------------------------------------------  
!  Fill Ocean SDS's ( loop around Iscan and Idata  )
!------------------------------------------------------------------------------------------------------------------                   
               call Fill_Output_Arrays_Ocean(Save_land_ocean_Q_Ret,Save_CLDFRC_Ocean,&
             Save_index_Land_Ocean,save_numpix_Land_Ocean,Save_Image_land_ocean_tau,&
             Save_land_ocean_Quality,save_Wspeed_Ocean,Save_small_weighting,Save_Least_error,&
             Save_SDS_angs_coeff1,Save_SDS_angs_coeff2,Save_effrad,Save_SDS_mass_conc_ocean,&
             Save_PSML003_Ocean,save_SDS_AOT_model,save_ref_Ocean,Save_STD_Ocean,Save_SDSTAUS_average,&
             Save_ASSY,Save_BACK,Save_average_Tau_Ocean,save_NUMPIXELS_Ocean,Save_land_ocean_tau,&
             Quality_to_pass,SDSTAU_average,Save_index_ocean,SDS_NUMPIXELS,SDS_small_weighting,&
             SDS_Least_error,SDS_angs_coeff1,SDS_angs_coeff2,SDS_effrad,SDS_mass_conc,SDS_ccn,&
             SDS_AOT_model,SDS_ref,SDS_ref_STD,SDSTAUS_average,SDSTAUB_average,Qcontrol_special,&
             SDSASSY_average,SDSBACK_average,Idata,Iscan,SDS_CLDFRC_ocean,WSPEED,Save_SDSTAUL_average,&
             Dust_flag_10KM,Save_Dust_Flag_ocean)
             
!-----------------------------------------------------------------------------------------------------------------  
!    Saving  Variables computed from visible channels. Optical depth at 0.55um, Models Fine and Coarse,Extinction 
!    coefficients for 0.55 um and  small mode weighting for Optical depth at 0.55um greater than 0.2
!------------------------------------------------------------------------------------------------------------------                   
                         Mode_F  = 0  
                         Mode_C  = 0 
                 
                 IF  ((SDSTAU_average(IDATA,2)     .ge. 0.0 .and. &
                      SDS_sol_INDX_small(IDATA,2)  .gt. 0.0 .and. &
                      SDS_sol_INDX_large(IDATA,2)  .gt. 0.0).and. & 
                     (Lat_center .NE. FV_GEO .and. Lat_center.NE.FV_GEO))then
                     Tau_550           =  SDSTAU_average(IDATA,2) 
                     Mode_F            = SDS_sol_INDX_small(IDATA,2) 
                     Mode_C            = SDS_sol_INDX_large(IDATA,2) 
                     Small_m_weighting  = SDS_small_weighting(IDATA,2)  
                     EXTSMALL_550       = Ext_554_small(1,Mode_F) 
                     EXTbig_550         = Ext_554_large(1,Mode_C-4)
                      DO ii = 1,NWAV_S
                       Tau_allwave(ii)       = SDSTAU_average(Idata,ii)
                       ref_allwave_Vis(ii)   = SDS_ref(Idata,ii)
                     ENDDO
                      
                    set_flag_read_UV = set_flag_read_UV+1  
!-----------------------------------------------------------------------------------------------------------------  
!  Go to Ret_Uv_SSa subroutine to Retrieve Ocean only Single scattering albedo and other diagonostics 
!------------------------------------------------------------------------------------------------------------------                   
                       
                          
                  call Ret_Uv_SSa(MTHET0,MTHET,MDPHI,Mode_F,Mode_C,wind, &
                    WSPEED,Small_m_weighting,Tau_550,set_flag_read_UV,&
                    REFW354,REFW388,Indx_wspeed1,Indx_wspeed2,&
                    Lat_center,Lon_center,&
                    EXTSMALL_550,EXTbig_550,Tau_allwave,WAVE,Omega_new,tau_new,&
                    Index_Omega_new,Height_indx,Fitting_error,Iscan,idata,ref_allwav_uv,&
                    Ret_Lines) 
                   
                     
      500          format(i6,20f10.4)             
                     
!!-----------------------------------------------------------------------------------------------------------------  
!    Fill Ocean SDS's for UV ( loop around Iscan and Idata  )
!------------------------------------------------------------------------------------------------------------------                   
                 
                   Call Fill_Output_Arrays_Ocean_UV(Iscan,Idata,tau_new,Omega_new,&
                   Tau_allwave,Mode_F,Mode_C,Small_m_weighting,&
                   Save_average_Tau_Ocean_UV,Save_average_Omega_Ocean_UV,Save_mode_F_FUV,&
                   Save_mode_C_FUV,Save_Small_weighting_FUV,Index_Omega_new,&
                   Save_Index_Albedo,Height_indx,Save_Index_Height,Fitting_error,&
                   Save_Fitting_Error,ref_allwav_uv,Save_ref_allwav_uv,ref_allwave_Vis,&
                    Save_ref_allwav_PACE,Save_average_Tau_PACE)
                    
                     
                   
!!-----------------------------------------------------------------------------------------------------------------  
!    Fill Ocean SDS's for UV ( loop around Iscan and Idata  )
!------------------------------------------------------------------------------------------------------------------                   
                      
             Endif  
             
!------------------------------------------------------------------------------------------------------------------                         
! endif for 100 pixel water at ocean
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
                  
         call  make_nc_Ocean(Save_Lat,Save_Lon,save_SolZen,save_View_angle,&
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
                  Ret_Xtrack,Ret_Lines,ipart,num_joint,Save_Dust_Flag_ocean) 
                  
                     
                   
!------------------------------------------------------------------------------------------------------------------------    
! call to write  output (full granule or Full disk of data) for UV Channels
!------------------------------------------------------------------------------------------------------------------   
!   
    
!           call  make_UV_nc(Save_Lat,Save_Lon,save_SolZen,save_View_angle,&
!                  save_View_phi,save_solar_phi,& 
!                  Save_average_Tau_Ocean_UV,Save_average_Omega_Ocean_UV,&
!                  Ret_Xtrack,Ret_Lines,ipart,yy1,yy2,Save_mode_F_FUV,&
!                  Save_mode_C_FUV,Save_Small_weighting_FUV,Save_Index_Albedo,&
!                  Save_Index_Height,Save_Fitting_Error,Save_ref_allwav_uv)          
!                    
                  
  8000         continue     
  
! End for programme        
        Return 
        end  subroutine DT_package_Ocean_main
           end module  compute_Ocean_main
           
           
           
           
               
     
 
