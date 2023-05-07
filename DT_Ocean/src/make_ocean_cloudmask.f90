        SUBROUTINE cldMsk_Ocean_updated(OW470,OW550,OW659_Hres,&
           OW138, OW124,CldMsk_500_Ocean,data_size,savecldmask, &
           START_line,END_line,Numscan) 
          
     
              

!----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:   This subroutine generates newcloud mask for 250*250
!                pixel resolution.
!
!!INPUT PARAMETERS:
!          CLDMSK_500     cldmask for 500 meter resolution
!          num_resol
!
!!OUTPUT PARAMETERS:
! NCLDMSK_SYN new cloud mask at 250 pixel resolution
!
!
! !REVISION HISTORY:
! $Log: Process_ocean_V2.f,v $
! 10/18/1999 fhliang
! fixed prolog;
!
! !TEAM-UNIQUE HEADER:
!
! !REFERENCES AND CREDITS:
!
! !DESIGN NOTES:
!  Editd 05 March 2012 by Leigh Munchak - lowered the lower bound of
! 1.38 reflectance threshold to apply ratio test from .01 to .005
!  
!
! !END
!-----------------------------------------------------------------------


        IMPLICIT  NONE
        SAVE

        INCLUDE   'mod04.inc'
        include   'read_Sat_MODIS.inc'

         
        Integer High_Cloud_Flag(ISWATH_B,ILINE)
        integer savecldmask(IX1KM_B*2,IY1KM_B*2) 
        character(len=10):: Sat_flag  
        INTEGER START_line,END_line
        integer CldMsk_500_Ocean(IX1KM_B*2,IY1KM_B*2)
        INTEGER  I,J,N,NUM,ISTART,IEND
        INTEGER SQNUM,M,IX,IXX,IY,IYY
        integer CldMsk_1km(IX1KM_B,IY1KM_B),kline
        integer  IBLUE,JBLUE,numaver,data_size(2)
        real AREFW550,AREFW670
        integer X2_offset,Y2_offset,l,p,LL,PP,add_clear,JX,JY
        integer oHigh_Cloud_F(IX1KM_B,IY1KM_B),ct
        integer  jjx,jjy 
        real  THRSH_3by3_Var,color_ratio,cloud_threhold_land47  
        real  THRHLD1380_1,THRHLD1380_2,Ratio_cirrus
        DATA THRSH_3by3_Var/0.0025/,color_ratio/0.75/,cloud_threhold_land47/0.4/ 
        DATA THRHLD1380_1/0.03/,THRHLD1380_2/.005/,Ratio_cirrus/0.30/
          
!-------------------------------------------------------------------------------------------------------------------
!  INTIALIZE THE ARRAY FOR NEW CLOUD MASK
!-------------------------------------------------------------------------------------------------------------------

              DO IYY=1,data_size(2)*2
                DO  IXX=1,data_size(1)*2
                 CldMsk_500_Ocean(IXX,IYY)=-9999 
                  savecldmask(IXX,IYY)=-9999
                ENDDO
             ENDDO
              START_line= 1
              END_line = data_size(2)
           
 !-------------------------------------------------------------------------------------------------------------------
 ! using channel 0.66 for cloud masking  with  (operationalMODIS  used  0.54)
 !3 x 3 spatial variability   
 !-------------------------------------------------------------------------------------------------------------------

            CALL Ocean_CldMsk_3by3_pixel_Up(OW659_Hres,THRSH_3by3_Var,&
            CldMsk_500_Ocean,data_size,START_line,END_line,numscan)
            
!            else
!            CALL Ocean_CldMsk_3by3_pixel_Up(OW550,THRSH,NROWS,&
!            RMED,RMEDSQ,CldMsk_500_Ocean,data_size,START_line,END_line,&
!            numscan)
!            Endif
!-------------------------------------------------------------------------------------------------------------------
! recover  and overwrite cldmask with clear if ratio condition is satsified to recall thick dust
!-------------------------------------------------------------------------------------------------------------------
  
                 kline=0
               DO JY=START_line,END_line
               IY =   2*jy-1
                DO JX=1,data_size(1) 
                 IX =  2*jx-1 
         do  jjy = IY,2*JY
         do  jjx = IX,2*jx 
             IF(OW470(JX,JY) .gt.0 .and. OW659_Hres(jjx,jjy) .gt.0 &
                   .and. OW659_Hres(jjx,jjy) .le. 1) then
                 IF((OW470(JX,JY)/OW659_Hres(jjx,jjy)) .le. color_ratio ) THEN
                      CldMsk_500_Ocean(jjx,jjy)=1
                  ENDIF
              ENDIF       
         enddo
         enddo
              ENDDO
             ENDDO        
   
             
!-------------------------------------------------------------------------------------------------------------------
! overwrite cldmask  with cloudy data if reflectance at 0.47 gt than the Threshold
!-------------------------------------------------------------------------------------------------------------------

          DO JY=START_line,END_line 
             IY =   2*jy-1
          DO JX=1,data_size(1) 
             IX =  2*jx-1 
           IF(OW470(JX,JY).gt.cloud_threhold_land47) THEN 
                 do  jjy = IY,2*JY
                 do  jjx = IX,2*jx
                  CldMsk_500_Ocean(jjx,jjy)=0  
                 enddo
                 enddo
            ENDIF
         ENDDO
       ENDDO
 
!-------------------------------------------------------------------------------------------------------------------
!  1.38 and 1.24 channels are used for cirrus masking  
!-------------------------------------------------------------------------------------------------------------------
 
          
           DO JY=START_line,END_line
              IY =   2*jy-1
           DO JX=1,data_size(1) 
               IX =  2*jx-1 
             do  jjy = IY,2*JY
             do  jjx = IX,2*jx
         IF(OW138(JX,JY).gt. THRHLD1380_1)then  
           CldMsk_500_Ocean(jjx,jjy)=0  
           Else
       IF(OW138(JX,JY).gt. THRHLD1380_2 .and.OW138(JX,JY).le.THRHLD1380_1)then 
! set Ratio  
      IF(OW124(JX,JY).gt.0)then  
        IF((OW138(JX,JY)/OW124(JX,JY)).GE.Ratio_cirrus)CldMsk_500_Ocean(jjx,jjy)=0  
      Endif 
        ENDIF 
             ENDIF
         ENDDO
        ENDDO 
         ENDDO
       ENDDO 
            ct   =0
           num =0
          DO JY=1,data_size(2)*2 
          DO JX=1,data_size(1)*2 
          num=num+1
         savecldmask(jx,jy)= CldMsk_500_Ocean(jx,jy) 
         if ( savecldmask(jx,jy) .gt.0)ct = ct+1
         Enddo
         Enddo 
!         print*,'Total cloud and clear',ct 
         return
         end
!***************************************************************************
        SUBROUTINE Ocean_CldMsk_3by3_pixel_Up(OW659_Hres,THRHLD,&
             CLDMSK,data_size,START_line,END_line,Numscan) 
!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!             This subroutine derive cloud mask using 3 x 3 spatial  
!              variability based upon MODIS 500 m resolution data
!
!!INPUT PARAMETERS:
!
!             Data_Size   Number of pixels in along and against scan
!                ISWATH   Number of pixels at 1 km resolution along scan
!                ILINE   Number of pixels at 1 km resolution against scan
!                refl_4   Reflectance at 500 meter resolution
!              THRHLD1   Threshold1 of 3x3 spatial variability
!              THRHLD2   Threshold2 of pixel reflectance value
!
 
 
      IMPLICIT NONE
       SAVE
      include   'read_Sat_MODIS.inc'
      INTEGER IY,JY,JX,NY,Data_Size(2),ISTART
      integer IEND 
      INTEGER NROWS(IX1KM_B*2),CLDMSK(IX1KM_B*2,IY1KM_B*2),add_numcld
       REAL  THRHLD
      REAL RMED(IX1KM_B*2),RMEDSQ(IX1KM_B*2),STD,VAR
      integer START_line,END_line,kline
     
        
! 
! Initialize work arrays 
      DO JX = 1, data_size(1)*2
        NROWS(JX) = 0
        RMED(JX) = 0
        RMEDSQ(JX) = 0
      ENDDO

!
! Checking 500 meter resolution 

      NY=0
      kline=1
       DO IY=START_line+1,END_line*2-1
         kline=kline+1
        ISTART=IY-1
        IEND=ISTART+2

        DO JY=ISTART,IEND
          NY=NY+1 

          DO JX=2,data_size(1)*2-1

             IF(OW659_Hres(JX,JY).GT.0.0.AND.OW659_Hres(JX+1,JY).GT.0.0 &
                 .AND.OW659_Hres(JX-1,JY).GT.0.0) THEN
              NROWS(JX)=NROWS(JX) + 1
              RMED(JX)=RMED(JX)+(OW659_Hres(JX,JY)+OW659_Hres(JX+1,JY)+OW659_Hres(JX-1,JY)) 
              RMEDSQ(JX)=RMEDSQ(JX)+(OW659_Hres(JX,JY)*OW659_Hres(JX,JY)&
                                   +OW659_Hres(JX+1,JY)*OW659_Hres(JX+1,JY)&
                                   + OW659_Hres(JX-1,JY)*OW659_Hres(JX-1,JY))
             
             ENDIF
          
!..........make clear determination where possible and re-initialize work array elements
            IF(NY.EQ.3) THEN
              
!..............make clear/cloud determination only when all 9 pixels in 3x3 array are valid 
               IF(NROWS(JX) .EQ. 3) THEN
                  VAR=9.0/8.0*(RMEDSQ(JX)/9.0-RMED(JX)*RMED(JX)/81.0)

                  IF(VAR.GT.0.0) THEN
                     STD=SQRT(VAR)

                     IF(STD.LT.THRHLD) THEN
                        CLDMSK(JX,kline)=1
                        Else
                        CLDMSK(JX,kline)=0
                     ENDIF
                  ENDIF
               ENDIF 
               NROWS(JX)=0
               RMED(JX)=0.0
               RMEDSQ(JX)=0.0
             ENDIF

          ENDDO
        ENDDO

        NY=0

      ENDDO 
        
   DO JX=1,data_size(1)*2
         CLDMSK(JX,1)=CLDMSK(JX,2) 
         CLDMSK(JX,data_size(2)*2)=CLDMSK(JX,data_size(2)*2-1)
      ENDDO
      
      DO JY=1,1,data_size(2)*2
          CLDMSK(1,JY)=CLDMSK(2,JY)
          CLDMSK(data_size(1)*2,JY)= CLDMSK(data_size(1)*2-1,JY)
      ENDDO   
   
     RETURN
      END
      
