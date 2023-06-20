module cld_flags_forLand
 
      implicit none

contains

! ---------------------------------------------------------------------
! compute cirrus quality flag and total clouds for computing cloud fraction for Land
!  
! --------------------------------------------------------------------- 
!*********************************************************************
      SUBROUTINE cld_flags_land(CldMsk_cirrus,New_CldMsk_500_Land,&
               CldMsk_250,cloud_num,cloud_num_land,START_1KM,END_1KM,&
                New_Q_cirrus,Ret_Quality_cirrus,checking_total_cloud) 
!---------------------------------------------------------------------
!!F90
! 
      IMPLICIT NONE
      SAVE

      include 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'

      INTEGER Water,Land,Desert,Glint,Snowice,cloud_num,cloud_num_land,&
              IYY,IXX,I,J,START_1KM,END_1KM,Det_cldmsk,Y4_offset,X4_offset,&
              Y2_offset,X2_offset,Y4,X4,Y2,X2,l,p,l11,p11,N,coast
       INTEGER  CldMsk_250(IX1KM_B,IY1KM_B),&
                New_CldMsk_500_Land(IX1KM_B,ILINE) ,&
                CldMsk_cirrus(IX1KM_B,IY1KM_B)
      INTEGER LandSea_Flag(IX1KM_B,IY1KM_B), &
              SunGlint_Flag(IX1KM_B,ILINE),SnowIce_Flag(IX1KM_B,ILINE),&
              SnowMsk_Ratio(IX1KM_B,ILINE), &    
              DET_Flag(IX1KM_B,ILINE),Shadow_Flag(IX1KM_B,ILINE),&
             High_Cloud_Flag(IX1KM_B,ILINE),&
             High_Cloud_Flag_500(IX1KM_B,ILINE),idata,iscan 
      INTEGER QA_Flag_Land(19),QA_Flag_Ocean(12),cloud_num2
      integer   Ret_Quality_cirrus, New_Q_cirrus(IX1KM_B,ILINE)
      integer  Newwater,NewLand,Nothing,cloud_num_clear, checking_total_cloud
      BYTE QA_Temp
      Ret_Quality_cirrus=1
      QA_Temp=0
      Det_cldmsk=0
      cloud_num=0
      cloud_num_land=0
      cloud_num_clear=0
      Water=0
      Land=0
      Desert=0
      Glint=0
      Snowice=0
      cloud_num2=0
     
! DT package all are in 500 meter resolution
         DO IYY = 1,ILINE
          DO IXX=START_1KM,END_1KM 
            IF(CldMsk_cirrus(IXX,IYY).EQ.0) THEN
              cloud_num=cloud_num+1 
            ENDIF 
        ENDDO
          ENDDO
! calculate number of cloudy pixels to compute cloud fraction.          
          
          DO IYY = 1,ILINE
          DO IXX=START_1KM,END_1KM   
              if( New_CldMsk_500_Land(IXX,IYY).EQ.0) THEN 
                 cloud_num_land=cloud_num_land+1  
             checking_total_cloud = checking_total_cloud+1
              else
               if( New_CldMsk_500_Land(IXX,IYY).EQ.1) THEN 
                cloud_num_clear =cloud_num_clear+1 
                endif
              endif
                ENDDO
                ENDDO 
                 N=0
!   Setting quality Flag to retrive with lower quality for cirrus or smoke?
        DO IYY = 1,ILINE
          DO IXX=START_1KM,END_1KM  
          IF(New_Q_cirrus(IXX,IYY).EQ.1  .and.&
              New_CldMsk_500_Land (IXX,IYY).eq.1) N=N+1 
             ENDDO
             ENDDO  
! if there are greater than 20% pixels that are flaged as clear but also has a cirrus flag             
!            
               IF(N .gt. ((Iline*Iline)*.20)) Ret_Quality_cirrus=0  
             Return
           end  subroutine cld_flags_land
                        
   end module  cld_flags_forLand
   
   
   