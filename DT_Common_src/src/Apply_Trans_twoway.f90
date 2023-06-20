 module Apply_Trans_Two_way
 
      implicit none

contains   
! ---------------------------------------------------------------------
!   
!         
 
   Subroutine Trans_Two_way(START_1KM,END_1KM,W659_SYN,&
       W865_SYN,W470_SYN,W550_SYN,W124_SYN,W164_SYN,W213_SYN,&
       W412_SYN,Multi_factor)   
    
       IMPLICIT NONE
        save
        
       include 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'
      INTEGER START_1KM,END_1KM 
      REAL W659_SYN(ISWATH_B,ILINE),W865_SYN(ISWATH_B,ILINE),&
           W470_SYN(ISWATH_B,ILINE),W550_SYN(ISWATH_B,ILINE),&
           W124_SYN(ISWATH_B,ILINE),W164_SYN(ISWATH_B,ILINE),&
           W213_SYN(ISWATH_B,ILINE),W1100_SYN(ISWATH_B,ILINE),&
           W412_SYN(ISWATH_B,ILINE),&
           Multi_factor(10)
          
 
      INTEGER idata,iscan,IXX,IYY  
       
      integer ICLDBLUE,ICLDRED,NUM_indx_x,J,I,IWAV,JBLUE,IBLUE
      Integer  NUM_indx_y, NUM_indx_z, NUM_indx_p,NUM_indx_Q
      INTEGER IMASK,IRED,JMASK,JRED,NUMCFR,NOWATER,II,JJ,IJ,IK
      Integer water,land,Num_Land, Num_water  
      character(len=10) :: Sat_Flag  
               
        DO  IXX=START_1KM,END_1KM  
        DO  IYY = 1,Iline
         W470_SYN(IXX,IYY)=W470_SYN(IXX,IYY)*Multi_factor(1) 
         W550_SYN(IXX,IYY)=W550_SYN(IXX,IYY)*Multi_factor(2)
         W659_SYN(IXX,IYY)=W659_SYN(IXX,IYY)*Multi_factor(3)
         W865_SYN(IXX,IYY)=W865_SYN(IXX,IYY)*Multi_factor(4)
         W124_SYN(IXX,IYY)=W124_SYN(IXX,IYY)*Multi_factor(5)
         W164_SYN(IXX,IYY)=W164_SYN(IXX,IYY)*Multi_factor(6) 
         W213_SYN(IXX,IYY)=W213_SYN(IXX,IYY)*Multi_factor(7) 
         W412_SYN(IXX,IYY)=W412_SYN (IXX,IYY)*Multi_factor(8)
        Enddo
        Enddo  
         Return
            end  subroutine Trans_Two_way
           end module  Apply_Trans_Two_way
           