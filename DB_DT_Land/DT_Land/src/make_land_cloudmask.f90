          SUBROUTINE CldMsk_Land_updated(Data_Size,OW470,OW124,&
         OW138,CldMsk_250,CldMsk_500,CldMsk_1km, &
         quality_cirrus,Aerosol_Cldmask_land,START_line,END_line,Numscan) 
   
!---------------------------------------------------------------------
!!F77
!  Cloud mask Land
!!---------------------------------------------------------------------
 
 
      IMPLICIT  NONE
      include   'read_Sat_MODIS.inc'
      include   'mod04.inc'
      SAVE
      INTEGER ISTART,IEND,Data_Size(2)
      INTEGER IX,IY,l,p,p2,ll,pp,ll2,pp2,N,iscan,NROWS(IX1KM_B)
      INTEGER X2_offset,Y2_offset,X4_offset,Y4_offset  
      REAL THRHLD1380_1,THRHLD1380_2
      REAL cloud_threhold_land47_1,cloud_threhold_land47_2
      REAL RMED(IX1KM_B),RMEDSQ(IX1KM_B)
      REAL RMED_1km(IX1KM_B),RMEDSQ_1km(IX1KM_B)
      INTEGER CldMsk_250(IX1KM_B,IY1KM_B)
      INTEGER CldMsk_500(IX1KM_B,IY1KM_B) 
      INTEGER CldMsk_1km(IX1KM_B,IY1KM_B)
      INTEGER  quality_cirrus(IX1KM_B,IY1KM_B) 
      Integer Aerosol_Cldmask_land(IX1KM_B,IY1KM_B)
       Integer START_line,END_line,kline,num_cloud,clear,cf,jx,jy
      
! Standard Threshold 
        DATA cloud_threhold_land47_1/0.0025/,cloud_threhold_land47_2/0.4/ 
         DATA THRHLD1380_1/0.003/,THRHLD1380_2/0.01/
        
          
            
          
 
!-------------------------------------------------------------------------------------------------------------------
!  INTIALIZE THE ARRAY FOR NEW CLOUD MASK
!-------------------------------------------------------------------------------------------------------------------
   
          
      DO IY=1,Data_Size(2) 
        DO  IX=1,Data_Size(1) 
          CldMsk_250(IX,IY)=-9999
        ENDDO
      ENDDO

      DO IY=1,Data_Size(2) 
        DO  IX=1,Data_Size(1)
          CldMsk_500(IX,IY)=-9999
        ENDDO
      ENDDO

       DO IY=1,Data_Size(2)
         DO  IX=1,Data_Size(1)
           CldMsk_1km(IX,IY)=-9999
            Quality_cirrus(IX,IY)=-9999
         ENDDO
       ENDDO

!-------------------------------------------------------------------------------------------------------------------
!Cloud mask based upon spatial variability of 0.47 micron reflectance  data
!-------------------------------------------------------------------------------------------------------------------


               START_line= 1
              END_line = Data_Size(2) 
      CALL CldMsk_3by3_HKM_Up(OW470,cloud_threhold_land47_1,NROWS,&
           RMED,RMEDSQ,CldMsk_500,data_size,&
           START_line,END_line,Numscan)
      
!-------------------------------------------------------------------------------------------------------------------
! reflectance test for 0.47 um
!-------------------------------------------------------------------------------------------------------------------

          kline=0
          num_cloud  = 0 
         DO IY=START_line,END_line
          kline= kline+1
         DO IX=1,Data_Size(1) 
        if(OW470(ix,iy) .gt.cloud_threhold_land47_2)CldMsk_500(IX, kline)=0 
        num_cloud  = num_cloud+1
         Enddo
        Enddo  
!         print*,'clody 0.47',num_cloud
!-------------------------------------------------------------------------------------------------------------------
! Cloud mask based upon spatial variability of 1.38  micron reflectance  data
!-------------------------------------------------------------------------------------------------------------------
 
           
        CALL CldMsk_3by3_1KM_Up(Data_Size,OW138,THRHLD1380_1,&
             THRHLD1380_2,CldMsk_1km,RMED_1km,RMEDSQ_1km,OW124,&
             quality_cirrus,NROWS,START_line,END_line,Numscan) 
 
          num_cloud  = 0 
          clear  = 0
         DO IY=1,Data_Size(2)
         DO IX=1,Data_Size(1)   
           IF(CldMsk_1km(IX,IY).EQ.1 .and.CldMsk_500(IX,IY).EQ.1)THEN
            CldMsk_500(IX,IY)=1
            Aerosol_Cldmask_land(IX,IY)=1 
            clear= clear+1 
          ELSE
            CldMsk_500(IX,IY)=0
            Aerosol_Cldmask_land(IX,IY)=0 
            num_cloud = num_cloud +1
            ENDIF  
         ENDDO
         ENDDO     
          
 20   continue
!       print*,'clear,cloudy,total land', clear,num_cloud,(Data_Size(1)*Data_Size(2))   
      RETURN
      END 
!***************************************************************************
        SUBROUTINE CldMsk_3by3_HKM_Up(OW470,THRHLD,&
       NROWS,RMED,RMEDSQ,CLDMSK,data_size,START_line,END_line,Numscan)
!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:
!              This subroutine derive cloud mask using 3 x 3 spatial 
!              variability based upon MODIS 500 m resolution data
!
!!INPUT PARAMETERS:
!
!             Data_Size   Number of pixels in along and against scan
!                ISWATH   Number of pixels at 1 km resolution along scan
!                 ILINE   Number of pixels at 1 km resolution against scan
!                REF1KM   Reflectance at 500 meter resolution
!               THRHLD1   Threshold1 of 3x3 spatial variability
!               THRHLD2   Threshold2 of pixel reflectance value
!
!OUTPUT PARAMETERS:
!
!                CLDMSK  Cloud mask
!                  RMED  Mean
!                RMEDSQ  Standard deviation
!
!!REVISION HISTORY:
!
!!TEAM-UNIQUE HEADER:
!
! Developed by MODIS Aerosol Team at NASA GSFC, Greenbelt, MD
! Edited 21 Dec 2011 by Leigh Munchak.
!Implemented overwrite of cloud mask based on standard deviation
! threshold of .0075 - only called if cloud mask is considered
! cloudy. This keeps bright but spatially homogenous smoke plumes
! from being classified as cloudy.
!
 
        IMPLICIT NONE
         include   'read_Sat_MODIS.inc'
          include   'mod04.inc'
        INTEGER IY,JY,JX,NY,Data_Size(2),ISTART
        integer IEND
        INTEGER NROWS(IX1KM_B) 
        INTEGER CLDMSK(IX1KM_B,IY1KM_B),add_numcld
         REAL  THRHLD
         REAL RMED(IX1KM_B),RMEDSQ(IX1KM_B),STD,STD2,VAR
         Integer START_line,END_line,kline,cf
        SAVE
        
     
!     Initialize work arrays 
        DO JX = 1, data_size(1)
           NROWS(JX) = 0
           RMED(JX) = 0
           RMEDSQ(JX) = 0
        ENDDO
        
     
!     Checking 500 meter resolution 
        
        NY=0 
        kline=1
       DO IY=START_line+1,END_line-1
         kline=kline+1
           ISTART=IY-1
           IEND=ISTART+2
           
         DO JY=ISTART,IEND
            NY=NY+1 
            
            DO JX=2,data_size(1)-1
               
     IF(OW470(JX,JY).GE.0.0.AND.OW470(JX+1,JY).GE.0.0.AND.OW470(JX-1,JY).GE.0.0) THEN

               NROWS(JX)=NROWS(JX) + 1
               RMED(JX)=RMED(JX)+(OW470(JX,JY)+OW470(JX+1,JY)&
                    +OW470(JX-1,JY))
               
               RMEDSQ(JX)=RMEDSQ(JX)+(OW470(JX,JY)*OW470(JX,JY)&
                    + OW470(JX+1,JY)*OW470(JX+1,JY)&
                    + OW470(JX-1,JY)*OW470(JX-1,JY))
    ENDIF
!               
!..........make clear determination where possible and re-initialize work array elements
            IF(NY.EQ.3) THEN
              
!..............make clear/cloud determination only when all 9 pixels in 3x3 array are valid 
               IF(NROWS(JX) .EQ. 3) THEN
                  VAR=9.0/8.0*(RMEDSQ(JX)/9.0-RMED(JX)*RMED(JX)/81.0)

                  IF(VAR.GT.0.0) THEN
                     STD=SQRT(VAR)
! New definitation Collection6
                    STD2 =STD* (RMED(JX)/3.0) 
                     IF(STD2.LT.THRHLD) THEN
                        CLDMSK(JX,kline)=1
                    ELSE                         
                     IF(STD .LT. .0075) THEN
                        CLDMSK(JX,IY)=1
                        Else
                        CLDMSK(JX,kline)=0
                     ENDIF
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
       
      DO JX=1,data_size(1)
         CLDMSK(JX,1)=CLDMSK(JX,2)
         CLDMSK(JX,(Data_Size(2)))=CLDMSK(JX,(Data_Size(2)-1))
      ENDDO
      
      DO JY=1,Data_Size(2)
         CLDMSK(1,JY)=CLDMSK(2,JY)
         CLDMSK(data_size(1),JY)=CLDMSK(data_size(1)-1,JY)
      ENDDO 
      RETURN
      END
!***************************************************************************
      SUBROUTINE CldMsk_3by3_1KM_Up(Data_Size,OW138,THRHLD1,&
      THRHLD2,CLDMSK,RMED,RMEDSQ,OW124,Quality_cirrus,NROWS,&
      START_line,END_line,Numscan)
     
     
!-----------------------------------------------------------------------
!!F90
!
!!DESCRIPTION:
!              This subroutine derive cloud mask using 3 x 3 spatial 
!              variability based upon MODIS 500 m resolution data
!
! 
!!TEAM-UNIQUE HEADER:
!
! Developed by MODIS Aerosol Team at NASA GSFC, Greenbelt, MD
!
!!DESIGN NOTES: 
!!END
!-----------------------------------------------------------------------
 
      IMPLICIT NONE
       include   'read_Sat_MODIS.inc'
         include   'mod04.inc'
      INTEGER IY,JY,NY,JX,Data_Size(2),ISTART,IEND
      INTEGER CLDMSK(IX1KM_B,IY1KM_B)
      INTEGER  Quality_cirrus(IX1KM_B,IY1KM_B)
      REAL THRHLD1,THRHLD2,THRHLD1380_Cirrus_1,THRHLD1380_Cirrus_2
      REAL RMED(IX1KM_B),RMEDSQ(IX1KM_B),STD,VAR 
      INTEGER Y2_offset,X2_offset,p2,l,p,kline,NROWS(IX1KM_B)
      integer START_line,END_line,cloudy,clear 
       DATA THRHLD1380_Cirrus_1/0.01/,THRHLD1380_Cirrus_2/0.025/
      SAVE
      
!    Initialize work arrays 
        DO JX = 1, data_size(1) 
           NROWS(JX) = 0
           RMED(JX) = 0
           RMEDSQ(JX) = 0
        ENDDO
        
        
!     Checking 500 meter resolution 
        
        NY=0 
        kline=1
       DO IY=START_line+1,END_line-1
         kline=kline+1
           ISTART=IY-1
           IEND=ISTART+2 
          DO JY=ISTART,IEND
            NY=NY+1  
            DO JX=2,data_size(1)-1 
            
            NROWS(JX)=NROWS(JX) + 1
       IF(OW138(JX,JY).GE.0.0.AND.OW138(JX+1,JY).GE.0.0.AND.OW138(JX-1,JY).GE.0.0) THEN
              RMED(JX)=RMED(JX)+(OW138(JX,JY)+OW138(JX+1,JY)+OW138(JX-1,JY))
              RMEDSQ(JX)=RMEDSQ(JX)+(OW138(JX,JY)*OW138(JX,JY)&
                                   + OW138(JX+1,JY)*OW138(JX+1,JY)&
                                   + OW138(JX-1,JY)*OW138(JX-1,JY))
       ENDIF
               IF(NY.EQ.3) THEN
                IF(NROWS(JX) .EQ. 3) THEN
                VAR=9.0/8.0*(RMEDSQ(JX)/9.0-RMED(JX)*RMED(JX)/81.0)
                IF(VAR.GT.0.0) THEN
                  STD=SQRT(VAR) 
                     IF(STD.LT.THRHLD1) THEN 
                    CLDMSK(JX,kline)=1
                    Else
                    CLDMSK(JX,kline)=0 
                   ENDIF
               ENDIF
                NROWS(JX)=0
                RMED(JX)=0.0
                RMEDSQ(JX)=0.0

               ENDIF

             ENDIF

          ENDDO
        ENDDO

        NY=0

      ENDDO 

    
      DO JX=1,data_size(1)
         CLDMSK(JX,1)=CLDMSK(JX,2)
         CLDMSK(JX,(Data_Size(2)))=CLDMSK(JX,(Data_Size(2)-1-1))
      ENDDO
      
      DO JY=1,Data_Size(2)
         CLDMSK(1,JY)=CLDMSK(2,JY)
         CLDMSK(data_size(1),JY)=CLDMSK(data_size(1)-1,JY)
      ENDDO 
      
      
      clear =0
          DO JY=1,Data_Size(2)
          DO JX=1,Data_Size(1)     
         IF(CLDMSK(JX,JY).EQ.1)clear= clear+1 
         enddo
         enddo
!         print*,'clear cirrus inside before',clear 
! Checking for 1.38 micron reflectance
  
      DO JY=1,Data_Size(2)
       DO JX=1,Data_Size(1)  
      IF(OW138(JX,JY).gt.THRHLD1380_Cirrus_2) then
             CLDMSK(JX,JY)=0
! save for cloud fraction 
      Else
!         
! if Variabily cloud mask says cloud free and smoke or cirrus reduce quality
         IF(CLDMSK(JX,JY) .eq.1 .and.OW138(JX,JY).gt. THRHLD1380_Cirrus_1 &
         .and.  OW138(JX,JY) .le. THRHLD1380_Cirrus_2 ) then 
! cloud free  may be cirrus  or smoke , will be retrived but put a quality value  
              Quality_cirrus(JX,JY)=1   
        ENDIF
!ENDIF   cloud free  may be cirrus  or smoke 
       ENDIF 
        ENDDO
      ENDDO 
        clear =0
         DO JY=1,Data_Size(2)
          DO JX=1,Data_Size(1)     
         IF(CLDMSK(JX,JY).EQ.1)clear= clear+1 
         enddo
         enddo
!         print*,'clear cirrus inside after',clear 
      RETURN 
      END
      
     


    