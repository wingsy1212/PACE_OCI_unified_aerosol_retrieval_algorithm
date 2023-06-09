!***************************************************************************
        SUBROUTINE Compute_Sd_for_Dust(START_line,END_line,Numscan,OW860,&
         SD_3by3,data_Size)
       

        IMPLICIT NONE
         include   'read_Sat_MODIS.inc'
          include   'mod04.inc'
        INTEGER IY,JY,JX,NY,Data_Size(2),ISTART
        integer IEND,Num
        INTEGER NROWS(IX1KM_B) 
        INTEGER add_numcld
         REAL  THRHLD,SD_3by3(IX1KM_B,IY1KM_B)
         REAL RMED(IX1KM_B),RMEDSQ(IX1KM_B),STD,STD2,VAR
         Integer START_line,END_line,kline,cf
        SAVE
       print*, (Numscan*ILINE),data_size(1)  
      
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
               
     IF(OW860(JX,JY).GT.0.0.AND.OW860(JX+1,JY).GT.0.0.AND.OW860(JX-1,JY).GT.0.0) THEN

               NROWS(JX)=NROWS(JX) + 1
               RMED(JX)=RMED(JX)+(OW860(JX,JY)+OW860(JX+1,JY)&
                    +OW860(JX-1,JY))
               
               RMEDSQ(JX)=RMEDSQ(JX)+(OW860(JX,JY)*OW860(JX,JY)&
                    + OW860(JX+1,JY)*OW860(JX+1,JY)&
                    + OW860(JX-1,JY)*OW860(JX-1,JY))
     ENDIF
               
!..........make clear determination where possible and re-initialize work array elements
            IF(NY.EQ.3) THEN
              
!..............make clear/cloud determination only when all 9 pixels in 3x3 array are valid 
               IF(NROWS(JX) .EQ. 3) THEN
                  VAR=9.0/8.0*(RMEDSQ(JX)/9.0-RMED(JX)*RMED(JX)/81.0)
                    IF(VAR.GT.0.0) THEN
                           STD=SQRT(VAR)
                          SD_3by3(JX,kline)=STD
                          
                    Else
                       SD_3by3(JX,kline)=0
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
         SD_3by3(JX,1)= SD_3by3(JX,2)
         SD_3by3(JX,data_size(2))= SD_3by3(JX,data_size(2)-1) 
      ENDDO
      
      DO JY=1,data_size(2)
          SD_3by3(1,JY)= SD_3by3(2,JY)
          SD_3by3(data_size(1),JY)= SD_3by3(data_size(1)-1,JY)
      ENDDO 
      num = 0
      DO JX=1,data_size(1)
      DO JY=1,Numscan*ILINE
      if(SD_3by3(jx,jy).gt.0)num=num+1
      ENDDO 
      ENDDO  
!      print*,'sd for dust',num
      RETURN
      END
      