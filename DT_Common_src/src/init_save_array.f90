       SUBROUTINE init_save_array(Save_Lat,Save_Lon,Save_SolZen,Save_View_angle,Save_View_phi,Save_solar_phi,& 
       Save_Scatt_angle,Save_Glint_angle,Save_Sea_Land_Flag,Save_Cld_mask_500,Save_cloud_distance_500,Save_land_ocean_Quality, &
       Save_Image_land_ocean_tau,Save_Aerosol_Type,Save_Fitting_Error_Land,Save_Surface_Reflectance_Land,Save_corrected,&
       Save_land_ocean_tau,Save_dust_weighting,Save_NUMPIXELS_land,Save_ref_Land,Save_SDS_Land,&
       Save_SDS_mass_conc_land,Save_CLDFRC_land,Save_average_Tau_Ocean,Save_SDSTAUS_average,Save_SDSTAUL_average,&
       Save_SDS_mass_conc_ocean,Save_CLDFRC_Ocean,Save_effrad,Save_PSML003_Ocean,Save_ASSY,Save_BACK,&
       Save_SDS_angs_coeff1,Save_SDS_angs_coeff2,Save_Least_error,Save_small_weighting,Save_SDS_AOT_model,Save_NUMPIXELS_Ocean,&
       Save_ref_Ocean,Save_STD_Ocean,Save_Wspeed_Ocean,Save_Topo_Height,Save_land_ocean_Q_Ret,&
      iswath,Tot_scan,set_Counter_Ocean, Set_Counter_Land,num_joint,num_joint2,&
      num_Sea_land,set_counter_for_Gread,set_counter_for_anc,Num_water,&
      num_land,Data_Size,START_line,END_line,IX1KM,IY1KM,numscan,Save_Dust_Flag_ocean,&
      Save_Lat_UV,Save_Lon_UV,&
      save_SolZen_UV,save_View_angle_UV,save_View_phi_UV,save_solar_phi_UV,Save_average_Tau_Ocean_UV,&
      Save_average_Omega_Ocean_UV,Save_mode_F_FUV,Save_mode_C_FUV,Save_Small_weighting_FUV,&
      Save_Index_Albedo,Save_Index_Height,Save_ref_allwav_uv,Save_ref_allwav_PACE,&
      Save_average_Tau_PACE) 
          IMPLICIT  NONE
      include 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'  
      INCLUDE 'Save_data_declare.inc'
      INCLUDE 'Save_data_UV_declare.inc' 
      SAVE  
      integer Iscan,Idata,Iwave_num,Imodel
       integer Set_Counter_Ocean, Set_Counter_Land,num_joint,num_joint2,&
       num_Sea_land,set_counter_for_Gread,set_counter_for_anc,Num_water,&
       num_land,START_line,END_line,Data_size(2)
             
             do Iscan =1,Tot_scan*Iline
             do Idata =1,ISWATH
            Save_Cld_mask_500(Idata,iscan) = -9999 
            Save_cloud_distance_500(Idata,iscan) = -9999 
            Enddo
            Enddo
           do Iscan =1,Tot_scan
           do Idata =1,ISWATH/ILINE
           Save_Lat(Idata,iscan)=-9999
           Save_Lon(Idata,iscan)=-9999
           save_SolZen(Idata,iscan)=-9999
           save_View_angle(Idata,iscan)=-9999
           save_View_phi(Idata,iscan)=-9999
           save_solar_phi(Idata,iscan)=-9999
           Save_Scatt_angle(Idata,iscan)=-9999
           Save_Glint_angle(Idata,iscan)=-9999
           save_Sea_Land_Flag(Idata,iscan)=-9999
           Save_land_ocean_Quality(idata,iscan)=-9999
           Save_Image_land_ocean_tau(idata,iscan)=-9999  
           Save_Fitting_Error_Land(idata,iscan)=-9999 
           Save_land_ocean_tau(idata,iscan)=-9999  
           Save_land_ocean_Q_Ret(idata,iscan,1)=-9999
           Save_land_ocean_Q_Ret(idata,iscan,2)=-9999
           Save_dust_weighting(idata,iscan)=-9999
           Save_SDS_mass_conc_land(idata,iscan)=-9999
           Save_CLDFRC_land(idata,iscan)=-9999
           Save_CLDFRC_Ocean(Idata,iscan)=-9999
           Save_Aerosol_Type(idata,iscan) = -9999
           Save_SDS_mass_conc_ocean(idata,iscan) = -9999
           save_Wspeed_Ocean(idata,iscan)=-9999
           Save_effrad(idata,iscan)=-9999
           Save_PSML003_Ocean(idata,iscan)=-9999 
           Save_Wspeed_Ocean(idata,iscan)=-9999
           Save_Least_error(Idata,iscan)=-9999 
           Save_small_weighting(Idata,iscan)=-9999 
            Save_Topo_Height(Idata,iscan) =-9999
           Save_SDS_angs_coeff1(Idata,iscan) =-9999
           Save_SDS_angs_coeff2(Idata,iscan) =-9999 
           Save_ave_clddis_land_ocean(Idata,iscan) =-9999
           Save_Dust_Flag_ocean(Idata,iscan) =-9999
           Save_Lat_UV(Idata,iscan) =-9999
          Save_Lon_UV(Idata,iscan) =-9999
          Save_SolZen_UV(Idata,iscan) =-9999
          Save_View_angle_UV(Idata,iscan) =-9999
          Save_View_phi_UV(Idata,iscan) =-9999
          Save_solar_phi_UV(Idata,iscan) =-9999
          Save_mode_F_FUV(Idata,iscan) =-9999
          Save_mode_C_FUV(Idata,iscan) =-9999
          Save_Small_weighting_FUV(Idata,iscan) =-9999 
          Save_Index_Albedo(Idata,iscan) =-9999 
          Save_Index_Height(Idata,iscan) =-9999 
          Save_smw_land(Idata,iscan) =-9999 
           Do Iwave_num =1,NWAV_uv+NWAV 
             Save_average_Tau_Ocean_UV(idata,iscan,Iwave_num)=-9999 
             save_ref_Land_PACE(idata,iscan,Iwave_num)=-9999 
          Enddo   
          Do Iwave_num =1,NWAV_uv  
             Save_average_Omega_Ocean_UV(idata,iscan,Iwave_num)=-9999  
          Enddo 
          Do Iwave_num = 1,NWAV_uv+NWAV 
             Save_ref_allwav_uv(idata,iscan,Iwave_num)=-9999
         Enddo 
         
           Do Imodel =1,num_model_index
           Save_SDS_AOT_model(idata,iscan,Imodel)=-9999
           Enddo
            Do Iwave_num = 1,NWAV_uv+3 
             Save_ref_allwav_PACE(idata,iscan,Iwave_num)=-9999
             Save_average_Tau_PACE(idata,iscan,Iwave_num)=-9999
           enddo
           Do Iwave_num =1,NWAVN  
            save_ref_Ocean(idata,iscan,Iwave_num)=-9999  
            Save_ref_Ocean(idata,iscan,IWAVE_NUm)=-9999 
            Save_STD_Ocean(idata,iscan,IWAVE_NUm)=-9999 
            Save_average_Tau_Ocean(idata,iscan,IWAVE_NUm)=-9999 
            Save_NUMPIXELS_Ocean(idata,iscan,IWAVE_NUm)=-9999 
            Save_SDSTAUS_average(idata,iscan,IWAVE_NUm)=-9999 
            Save_SDSTAUL_average(idata,iscan,IWAVE_NUm)=-9999 
            Save_ASSY(idata,iscan,IWAVE_NUm)=-9999 
            Save_BACK(idata,iscan,IWAVE_NUm)=-9999 
           ENDDO
            Do Iwave_num =1,Land_Sol3
             save_Surface_Reflectance_Land(idata,iscan,Iwave_num)=-9999 
             save_ref_Land(idata,iscan,Iwave_num)=-9999
             save_SDS_Land(idata,iscan,Iwave_num)=-9999 
             save_NUMPIXELS_land(idata,iscan,Iwave_num)=-9999 
           ENDDO 
            Do Iwave_num =1,Land_Sol4
            save_corrected(idata,iscan,Iwave_num)=-9999 
           ENDDO 
           
           Enddo
           Enddo
! Set_Counter_Ocean is set to read the table once as first entry
!     into ocean algorithm 
       
       Set_Counter_Ocean=0
       Set_Counter_Land=0
       num_joint=0
       num_joint2=0
       num_Sea_land=0
       set_counter_for_Gread=0
       set_counter_for_anc=0
         Num_water=0
         Num_land=0 
         Data_Size(1)= IX1KM
         Data_Size(2)= IY1KM
          START_line=1
         END_line= Iline*numscan  
           return
           end
           
           