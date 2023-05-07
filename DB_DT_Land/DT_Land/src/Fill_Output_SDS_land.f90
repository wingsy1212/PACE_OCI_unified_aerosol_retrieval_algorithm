 module Fill_Output_SDS_land
 
      implicit none

contains   
               
         SUBROUTINE Fill_Output_Arrays_land(water,coastal_Pix,Quality_flag_forJoint,&
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
            
                
             INCLUDE 'mod04.inc' 
             INCLUDE 'read_Sat_MODIS.inc'  
             INCLUDE 'Set_Array_dimension_Land.inc'
             INCLUDE 'Save_data_declare.inc'
               
             
!-----------------------------------------------------------------------------------------------------------------           
! reduce Quality if water pixels greater than 20% or coastal pixels gt 50%
!-----------------------------------------------------------------------------------------------------------------            
             Twenty_water = (Iline*iline)  * .2
             fifty_coastal =(Iline*iline)  * .5 
          
               
        if( water .gt. Twenty_water .or. coastal_Pix .gt. Fifty_coastal) Quality_flag_forJoint(1)=0   
             
             Save_land_ocean_Q_Ret(Idata,iscan,2)=Real(Quality_flag_forJoint(2)) 
             Save_CLDFRC_land(Idata,iscan) = SDS_CLDFRC_land(IDATA)  
             Optical_depth_land= SDSTAU_corrected(IDATA,2)  
          
             if(Optical_depth_land .Ge.-0.05 .and. Optical_depth_land .le. 5 &
                   .or.Qcontrol_special_land .gt.0)then  
!   For ABI   0.55 um is missing so we use number of pixels for .66um for Leo and Geo
                   Wave_index =3 
                do iyy=1,INT(SDS_NUMPIXELS_land(idata,Wave_index))
                Save_index_Land_Ocean(IDATA,Iscan,IYY)= Save_index_land(IDATA,IYY) 
                 Enddo  
             save_numpix_Land_Ocean(IDATA,Iscan)=INT(SDS_NUMPIXELS_land(idata,Wave_index)) 
             Save_Image_land_ocean_tau (idata,iscan) = Optical_depth_land 
             Save_Topo_Height(idata,iscan) = MHGHT  
             Save_land_ocean_Quality(Idata,iscan)=Quality_flag_forJoint(1) 
             Save_Aerosol_Type(idata,iscan) = SDS_Aerosol_Type(idata)
             Save_Fitting_Error_Land(Idata,iscan) =  SDS_Fitting_Error_Land(idata) 
             Save_dust_weighting(Idata,iscan) = SDS_dust_weighting(IDATA) 
             Save_SDS_mass_conc_land(Idata,iscan) = SDS_mass_conc_land(IDATA)  
  !   save only reflectance for 0.46.0.66 and 2.1            
            Do IWAVE_NUm= 1,Land_Sol3
             if(IWAVE_NUm .eq.1)Variable=SDS_ref_land(Idata,IWAVE_NUm+2)
             if(IWAVE_NUm .eq.2)Variable=SDS_ref_land(Idata,IWAVE_NUm+3)
             if(IWAVE_NUm .eq.3)Variable=SDS_ref_land(Idata,IWAVE_NUm+6)
             if(IWAVE_NUm .eq.1)Variable1=SDS_ref_STD_land(Idata,IWAVE_NUm+2)
             if(IWAVE_NUm .eq.2)Variable1=SDS_ref_STD_land(Idata,IWAVE_NUm+3)
             if(IWAVE_NUm .eq.3)Variable1=SDS_ref_STD_land(Idata,IWAVE_NUm+6)
            if( Variable.gt.0.or.Variable.le.1) then 
            save_ref_Land(idata,iscan,IWAVE_NUm) =  Variable 
            save_SDS_Land(idata,iscan,IWAVE_NUm) =  Variable1  
            endif
         Enddo 
!   Saving for PACE
             Do IWAVE_NUm= 1,9
             save_ref_Land_PACE(idata,iscan,IWAVE_NUm)= SDS_ref_land(Idata,IWAVE_NUm)
             ENDDO
         
         
!!   save only  surface reflectance  for 0.46.0.66 and 2.1            
        Do IWAVE_NUm= 1,Land_Sol3 
          save_Surface_Reflectance_Land(idata,iscan,IWAVE_NUm)=SDS_Surface_Reflectance_Land(idata,IWAVE_NUm)  
        Enddo
 !!   save only  corrected optical depth for  for 0.46 0.55 0.66 and 2.1             
           Do IWAVE_NUm= 1,Land_Sol4 
           if(IWAVE_NUm .le.3)then
           save_corrected(idata,iscan,IWAVE_NUm)=SDSTAU_corrected(IDATA,IWAVE_NUm) 
           else
           save_corrected(idata,iscan,IWAVE_NUm)=SDSTAU_corrected_213(IDATA)
           endif  
           Enddo 
  !!   save only SDS_NUMPIXELS_land for  for 0.46.0.66 and 2.1      
            Do IWAVE_NUm= 1,Land_Sol3 
             if(IWAVE_NUm .eq.1)Variable2=SDS_NUMPIXELS_land(Idata,IWAVE_NUm+2)
             if(IWAVE_NUm .eq.2)Variable2=SDS_NUMPIXELS_land(Idata,IWAVE_NUm+3)
             if(IWAVE_NUm .eq.3)Variable2=SDS_NUMPIXELS_land(Idata,IWAVE_NUm+6)
              save_NUMPIXELS_land(idata,iscan,IWAVE_NUm)=Variable2
           Enddo 
              if(  Quality_flag_forJoint(1) .gt. 0)&
            Save_land_ocean_tau(idata,iscan) =SDSTAU_corrected(IDATA,2) 

 
            
            
            
            
! Endif for optical depth range
           Endif
           
           
           
           
           
      RETURN
   end  subroutine Fill_Output_Arrays_land
           end module  Fill_Output_SDS_land
           
           
           
           
              
           
           
