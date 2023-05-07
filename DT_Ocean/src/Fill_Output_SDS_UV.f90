 module Fill_Output_SDS_Ocean_UV
 
      implicit none

contains   
               
         SUBROUTINE Fill_Output_Arrays_Ocean_UV(Iscan,Idata,tau_new,Omega_new,&
         Tau_allwave,Mode_F,Mode_C,Small_m_weighting,&
         Save_average_Tau_Ocean_UV,Save_average_Omega_Ocean_UV,Save_mode_F_FUV,&
         Save_mode_C_FUV,Save_Small_weighting_FUV,Index_Omega_new,Save_Index_Albedo,&
         Height_indx,Save_Index_Height,Fitting_error,Save_Fitting_Error,&
         ref_allwav_uv,Save_ref_allwav_uv,ref_allwave_Vis,&
          Save_ref_allwav_PACE,Save_average_Tau_PACE)
                  
                
             include 'mod04.inc' 
             INCLUDE 'read_Sat_MODIS.inc'  
             INCLUDE 'Set_Array_dimension.inc'
             INCLUDE 'Save_data_UV_declare.inc'  
          
            Do IWAVE_NUm = 1,NWAV_uv 
            Save_average_Tau_Ocean_UV(idata,iscan,IWAVE_NUm)  =  tau_new(IWAVE_NUm) 
            Save_ref_allwav_uv(idata,iscan,IWAVE_NUm) = ref_allwav_uv(IWAVE_NUm)
            enddo 
!    0.34,0.388,0.55             
            Do IWAVE_NUm = 1,NWAV_uv+1 
            Save_average_Omega_Ocean_UV(idata,iscan,IWAVE_NUm) = Omega_new(IWAVE_NUm) 
            Enddo
 ! filling  Visible  channel data         
              Do IWAVE_NUm =  3,(NWAVN+NWAV_uv) 
            Save_ref_allwav_uv(idata,iscan,IWAVE_NUm)= ref_allwave_Vis(IWAVE_NUm-2)
            Save_average_Tau_Ocean_UV(idata,iscan,IWAVE_NUm)= Tau_allwave(IWAVE_NUm-2)
            Enddo 
 !            
            Save_mode_F_FUV(idata,iscan)=Real(Mode_F)
            Save_mode_C_FUV(idata,iscan)=Real(Mode_C-4)
            Save_Small_weighting_FUV(idata,iscan)= Small_m_weighting 
            Save_Index_Albedo(idata,iscan)= Index_Omega_new(2)
            Save_Index_Height(idata,iscan)= Height_indx(2)
            Save_Fitting_error(idata,iscan)= Fitting_error(2) 
 ! Save Variables  for PACE fused Algorithm  UV 0.354 and 0.388         
            Do IWAVE_NUm = 1,NWAV_uv 
            Save_ref_allwav_PACE(idata,iscan,IWAVE_NUm)=&
            Save_ref_allwav_uv(idata,iscan,IWAVE_NUm)
            Save_average_Tau_PACE(idata,iscan,IWAVE_NUm)=&
            Save_average_Tau_Ocean_UV(idata,iscan,IWAVE_NUm)
            ENDDO
! Save Variables  for PACE fused Algorithm 48 0.55  67Um  
            Do IWAVE_NUm = 3,5        
            Save_ref_allwav_PACE(idata,iscan,IWAVE_NUm )=&
            Save_ref_allwav_uv(idata,iscan,IWAVE_NUm)
            Save_average_Tau_PACE(idata,iscan,IWAVE_NUm )=&
            Save_average_Tau_Ocean_UV(idata,iscan,IWAVE_NUm)
            Enddo
! Quality   
      RETURN
           end   subroutine Fill_Output_Arrays_Ocean_UV
           end module  Fill_Output_SDS_Ocean_UV
           
           
                 
           
           
