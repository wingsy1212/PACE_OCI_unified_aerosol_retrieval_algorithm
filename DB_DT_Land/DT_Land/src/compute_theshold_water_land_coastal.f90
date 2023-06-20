      module compute_theshold_water_land_coastal 
      implicit none 
contains   
! ---------------------------------------------------------------------
!     
!!SUBROUTINE theshold_water_land_coastal 
        
            subroutine threshold_water_land_coastal(water,Water_Pixels,Pure_Land,Land_Pixels,&
            Sea_Land_Flag,save_Sea_Land_Flag,iscan,idata) 
      IMPLICIT  NONE
      include 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'   
      INCLUDE 'Set_Array_dimension_Land.inc'
      INCLUDE 'Save_data_declare.inc' 
          
             if(water .gt.0 .or. land .gt.0 ) Land_Pixels = 0   
              Water_Pixels =iline * iline
! 20 % of pixels are coastal then  set flag to coastal    
              Num_coastal_pixels = (iline * iline) * 0.20    
             if(WATER .ge. Water_Pixels)then
              save_Sea_Land_Flag(Idata,iscan)=0 
              else 
             if(Pure_Land .gt.0)save_Sea_Land_Flag(Idata,iscan)=1  
             if(coastal_Pix .le.Num_coastal_pixels)save_Sea_Land_Flag(Idata,iscan)=1
             if(coastal_Pix .gt.Num_coastal_pixels)save_Sea_Land_Flag(Idata,iscan)=2 
             Endif 
             Sea_Land_Flag(idata) = save_Sea_Land_Flag(Idata,iscan)
             return
            end  subroutine threshold_water_land_coastal
           end module  compute_theshold_water_land_coastal
           
           
           