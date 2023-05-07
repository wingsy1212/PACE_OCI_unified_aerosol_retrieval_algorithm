 
      module compute_water_land
 
      implicit none

contains   
! ---------------------------------------------------------------------
!     
!!SUBROUTINE Get_water_land
! ---------------------------------------------------------------------

       Subroutine Get_water_land(START_1KM,END_1KM,&
       LandSea_Flag_new,land,water)
       IMPLICIT NONE
        save 
        INCLUDE 'mod04.inc' 
        INCLUDE 'read_Sat_MODIS.inc'
      INTEGER START_1KM,END_1KM,coastal_Pix,Pure_Land   
      INTEGER  IXX,IYY,LandSea_Flag_new(ISWATH_B,ILINE)   
      Integer water,land,num
             water=0
             land=0  
           DO  IXX=START_1KM,END_1KM  
            DO  IYY = 1,IGRIDY   
          if(LandSea_Flag_new(IXX,IYY) .eq. 1) water=water+1
          if(LandSea_Flag_new(IXX,IYY) .eq. 0) land=land+1 
          Enddo
          Enddo  
         Return
          end  subroutine Get_water_land
           end module compute_water_land
           
           
