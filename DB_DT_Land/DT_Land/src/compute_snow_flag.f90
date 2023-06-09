         
 module snow_flag
 
      implicit none

contains   
! ---------------------------------------------------------------------
!  compute  snow flag  
!         
          Subroutine compute_snow_flag(START_1KM,END_1KM,W865_SYN,&
          W164_SYN,skinTemp,SnowMsk_Ratio) 
              
         include 'mod04.inc' 
         include  'read_Sat_MODIS.inc'  
         INTEGER START_1KM,END_1KM,IXX,IYY
         Integer SnowMsk_Ratio(ISWATH_B,ILINE),cc
         Real   W865_SYN(ISWATH_B,ILINE),Ratio
         Real   W124_SYN(ISWATH_B,ILINE),W1100_SYN(ISWATH_B,ILINE)
         Real  W164_SYN(ISWATH_B,ILINE),skinTemp
!temperature conversion variables
      REAL  Planck_constant,Speed_light,Boltz_cons,wav2
      REAL   w_meter,c1,c2,W1100_Temp(ISWATH_B,ILINE) 
      character(len=10) :: Sat_Flag 
      PARAMETER(Planck_constant=6.6260755e-34,&
         Speed_light=2.9979246e+8, Boltz_cons=1.380658e-23,&
          wav2=11.0) 

! FOR ..........Snow masking        
!compute Temprature(convert from radiance to temperature 11.um channel)
! Derive constants  


! comute Rong_rond ratio. 
! 1.24 channel is replaced by 1.64 channel to compute snow flag
! Since we use standard subroutine for LEOS and GEOs we are using common channel 1.64 since
! 1.24 is missing from GEOS.
! NOTE : channel 11.2 Radiance for GOES is converted to same units as LEOS and AHI
! intialize
           cc =0
         DO  IXX=START_1KM,END_1KM  
         DO  IYY = 1,Iline 
         SnowMsk_Ratio(IXX,IYY) = 1
         enddo
         enddo 
         
          Ratio=0  
      DO  IXX=START_1KM,END_1KM  
        DO  IYY = 1,Iline 
            SnowMsk_Ratio(IXX,IYY)= 1
        if(W865_SYN(IXX,IYY) .gt.0 .and.W164_SYN(IXX,IYY).gt.0) then
             Ratio=(W865_SYN(IXX,IYY)-W164_SYN(IXX,IYY))/&
              (W865_SYN(IXX,IYY)+W164_SYN(IXX,IYY))
         Else
         Ratio =0
         Endif 
! Set snow mask as snowy pixel based on Ratio and temp(Adapted from  NOAA GOES ATBD)

        If(ratio.gt.0.3 .and. (skinTemp.gt. 0 .and. &
                                skinTemp .lt. 280))  THEN  
                                SnowMsk_Ratio(IXX,IYY)=0 
           if(SnowMsk_Ratio(IXX,IYY) .eq.0) cc =cc+1
              
         Endif 
           Enddo
           enddo 
!          print*,'Snow flag and skin temp',cc ,skinTemp  
         Return
          end  subroutine compute_snow_flag
           end module  snow_flag