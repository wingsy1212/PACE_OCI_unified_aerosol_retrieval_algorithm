 module Fill_Output_SDS_Ocean
 
      implicit none

contains   
               
         SUBROUTINE Fill_Output_Arrays_Ocean(Save_land_ocean_Q_Ret,Save_CLDFRC_Ocean,&
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
                
             include 'mod04.inc' 
             INCLUDE 'read_Sat_MODIS.inc'  
             INCLUDE 'Set_Array_dimension.inc'
             INCLUDE 'Save_data_declare.inc'
               
                 Optical_depth_Ocean =-9999 
               Save_land_ocean_Q_Ret(Idata,iscan,1)  =Real(Abs(Quality_to_pass(2)))
               Save_CLDFRC_Ocean(Idata,iscan) = SDS_CLDFRC_ocean(IDATA) 
               
            
!All#quality data    !  reporting only for Average solution    
            if(SDSTAU_average(IDATA,2) .GT.-0.01 .and.&
               SDSTAU_average(Idata,2) .LE. 5.0 .or. Qcontrol_special .gt. 0) then  
            Optical_depth_Ocean=SDSTAU_average(IDATA,2)  
                    Wave_index = 4  
               do iyy=1,INT(SDS_NUMPIXELS(idata,Wave_index))  
                 Save_index_Land_Ocean(IDATA,Iscan,IYY)= Save_index_ocean(IDATA,IYY)   
                   Enddo 
            save_numpix_Land_Ocean(IDATA,iscan)= INT(SDS_NUMPIXELS(idata,Wave_index))  
            Save_Image_land_ocean_tau(idata,iscan)= Optical_depth_Ocean
           Save_land_ocean_Quality(Idata,iscan)  = Real(Quality_to_pass(1))  
           save_Wspeed_Ocean(idata,iscan)=WSPEED 
           Save_small_weighting(Idata,iscan) = SDS_small_weighting(IDATA,2) 
           Save_Least_error(Idata,iscan) = SDS_Least_error(Idata,2) 
           Save_SDS_angs_coeff1(Idata,iscan) = SDS_angs_coeff1(IDATA,2)
           Save_SDS_angs_coeff2(Idata,iscan) = SDS_angs_coeff2(IDATA,2) 
           Save_effrad(Idata,iscan) = SDS_effrad(Idata,2)
           Save_Dust_Flag_ocean(Idata,iscan) = Real(Dust_flag_10KM(idata))
           Save_SDS_mass_conc_ocean(Idata,iscan)= SDS_mass_conc(idata,2)
           Save_PSML003_Ocean(Idata,iscan) = SDS_ccn(idata,2)
!   For ABI   0.55 um is missing so we use number of pixels for 0.86         
               
            
           Do Imodels = 1,num_model_index
            save_SDS_AOT_model(idata,iscan,Imodels)= SDS_AOT_model(idata,Imodels)
            Enddo
          Do IWAVE_NUm= 1,NWAVN  
            save_ref_Ocean(idata,iscan,IWAVE_NUm)=SDS_ref(Idata,IWAVE_NUm)
            Save_STD_Ocean(idata,iscan,IWAVE_NUm)=SDS_ref_STD(Idata,IWAVE_NUm)
            Save_SDSTAUS_average(idata,iscan,IWAVE_NUm) = SDSTAUS_average(IDATA,IWAVE_NUm) 
            Save_SDSTAUL_average(idata,iscan,IWAVE_NUm) = SDSTAUB_average(IDATA,IWAVE_NUm) 
            Save_ASSY(Idata,iscan,IWAVE_NUm) =  SDSASSY_average(idata,IWAVE_NUm) 
            Save_BACK(Idata,iscan,IWAVE_NUm) =  SDSBACK_average(idata,IWAVE_NUm) 
            Save_average_Tau_Ocean(idata,iscan,IWAVE_NUm)=SDSTAU_average(Idata,IWAVE_NUm)
            save_NUMPIXELS_Ocean(idata,iscan,IWAVE_NUm) =SDS_NUMPIXELS(IDATA,IWAVE_NUm)
          enddo  
! Quality  
      
 !          If( Quality_to_pass(1) .gt.0)then
             Save_land_ocean_tau(Idata,iscan) = Optical_depth_Ocean 
 !            Endif   
           
          Endif
      RETURN
   end  subroutine Fill_Output_Arrays_Ocean
           end module  Fill_Output_SDS_Ocean
           
           
           
              
           
           
