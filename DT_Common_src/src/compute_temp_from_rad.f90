         
 module compute_Temperature_from_rad
 
      implicit none

contains   
! ---------------------------------------------------------------------
!  compute  snow flag  
!         
          Subroutine compute_temp_from_rad(START_1KM,END_1KM,&
             W1100_SYN,W8p5_SYN,W1100_Temp,W8p5_Temp)
              
         include 'mod04.inc' 
         include  'read_Sat_MODIS.inc'  
         INTEGER START_1KM,END_1KM,IXX,IYY ,Iwave
         Real   W8p5_SYN(ISWATH_B,ILINE),W1100_SYN(ISWATH_B,ILINE)
         Real   W1100_Temp(ISWATH_B,ILINE),W8p5_Temp(ISWATH_B,ILINE) 
!temperature conversion variables
      REAL  Planck_constant,Speed_light,Boltz_cons,wav1,wav2
      REAL   w_meter,c1,c2
      character(len=10) :: Sat_Flag 
      PARAMETER(Planck_constant=6.6260755e-34,&
         Speed_light=2.9979246e+8, Boltz_cons=1.380658e-23,&
          wav1=8.75,wav2=11.0) 
            
            
!  OCI is missing thermal channel 
! Only MODIS and VIIRS has all channels to compute snow flag
!          IF(sat_Flag .eq. "MODIS" .or. sat_Flag .eq. "VIIRS")then
!          IF(sat_Flag .eq. "MODIS" )then
         c1=2.0*Planck_constant*(Speed_light*Speed_light)
         c2=(Planck_constant*Speed_light)/Boltz_cons
! convert wavelength to meters 
         Do Iwave = 1,2  
         if(Iwave .eq.1) w_meter=(1.0e-6*wav1)
         if(Iwave .eq.2) w_meter=(1.0e-6*wav2)  
         DO  IXX=START_1KM,END_1KM  
         DO  IYY = 1,Iline 
         if(Iwave .eq.1) W8p5_Temp(IXX,IYY)=&
         c2/(w_meter*alog(c1/(1.0e+6*W8p5_SYN(IXX,IYY)*w_meter**5)+1.0)) 
         If(Iwave .eq.2) W1100_Temp(IXX,IYY)= &
         c2/(w_meter*alog(c1/(1.0e+6*W1100_SYN(IXX,IYY)*w_meter**5)+1.0)) 
         Enddo
         Enddo 
!  Enddo for Wave lengths        
         Enddo 
         Return
          end  subroutine compute_temp_from_rad
           end module  compute_Temperature_from_rad