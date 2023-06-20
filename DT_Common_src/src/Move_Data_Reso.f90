 
       Subroutine Move_Data_Reso(OW659,OW860,OW470,OW550,OW124,OW138,OW164,&
      OW213,OW1100,W470_SYN,W550_SYN,W659_SYN,W865_SYN,W124_SYN,&
      W164_SYN,W213_SYN,W138_SYN,W1100_SYN,Land_Sea_Flag,&
      OLAT,OLON,OSOLA,OVIEW,OAZIM_View,OAZIM_Solar,OHeight,&
      LandSea_Flag_new,Lat,Lon,SatZen,SatAz,SolZen,SolAz,Height,&
      START_line,END_line,CldMsk_500_Ocean,CldMsk_500_Land,&
      savecldmask,New_CldMsk_500_Ocean,New_CldMsk_500_Land,New_savecldmask,&
      l1b_nLines,quality_cirrus,New_Q_cirrus,OW412,W412_SYN,OW8p5,W8p5_SYN,SD_3by3,&
      SD_for_Dust,W354_SYN,W388_SYN,OW354,OW388)
      
      
      
      
       IMPLICIT NONE
       INCLUDE 'mod04.inc'
       INCLUDE 'read_Sat_MODIS.inc'
 

!----------------------------------------------------------------------
      INTEGER IL,IM,START_line,END_line,iscan
      INTEGER   numaver1,numaver2,numaver3,numaver4
      INTEGER JMASK,IMASK,IBLUE,JBLUE,IXX,IYY,IX,IY,IIJ,IIK 
      REAL W659_SYN(ISWATH_B,ILINE),W865_SYN(ISWATH_B,ILINE),&
           W470_SYN(ISWATH_B,ILINE),W550_SYN(ISWATH_B,ILINE),&
           W124_SYN(ISWATH_B,ILINE),W164_SYN(ISWATH_B,ILINE),&
           W213_SYN(ISWATH_B,ILINE),W412_SYN(ISWATH_B,ILINE)
      REAL W8p5_SYN(ISWATH_B,ILINE),W138_SYN(ISWATH_B,ILINE),&
           W1100_SYN(ISWATH_B,ILINE),W1200_SYN(ISWATH_B,ILINE),&
           SD_3by3(ISWATH_B,Numscan_B*ILINE)
      Integer  Shadow_FLAG(ISWATH_B,ILINE),DET_FLAG(ISWATH_B,ILINE),&
              SunGlint_Flag(ISWATH_B,ILINE),SnowIce_Flag(ISWATH_B,ILINE),&
              LandSea_Flag(ISWATH_B,ILINE),UFQ_Flag(ISWATH_B,ILINE),&
           Non_CloudOb_Flag(ISWATH_B,ILINE),LandSea_Flag_new(ISWATH_B,ILINE) 
       REAL LAT(ISWATH_B,ILINE),LON(ISWATH_B,ILINE)
       Real SatZen(ISWATH_B,ILINE),SatAz(ISWATH_B,ILINE)
       Real SolZen(ISWATH_B,ILINE),Height(ISWATH_B,ILINE)
       Real SD_for_Dust(ISWATH_B,ILINE)
       Real SolAz(ISWATH_B,ILINE), Wave_659(ISWATH_B,ILINE),ave
       integer CldMsk_500_Ocean(ISWATH_B*2,(ILINE*Numscan_B)*2)
       integer CldMsk_500_Land(ISWATH_B,(ILINE*Numscan_B)+1)
       Integer New_CldMsk_500_Ocean(ISWATH_B*2,ILINE*2)
       Integer New_CldMsk_500_Land(ISWATH_B,ILINE+1)
       integer savecldmask(ISWATH_B*2,(Numscan_B*ILINE)*2)
       integer New_savecldmask(ISWATH_B*2,ILINE*2)
       integer quality_cirrus(iswath_B,(Numscan_B*ILINE)+1)
       Integer New_Q_cirrus(ISWATH_B,ILINE),k,j,jx,jy,num   
       !    uv
       REAL   W354_SYN(ISWATH_B,ILINE), W388_SYN(ISWATH_B,ILINE) 
        
           
          IL=0
         DO    IYY = START_line+1,END_line
               IL=IL+1
                IM=0
         DO    IXX= 1,l1b_nLines
               IM=IM+1   
         IF( OSOLA(IXX,IYY) .le. MAXMTHET0) then  
             W470_SYN(IM,IL)=OW470(IXX,IYY)
             W550_SYN(IM,IL)=OW550(IXX,IYY) 
             W659_SYN(IM,IL)= OW659(IXX,IYY)
             W865_SYN(IM,IL)=OW860(IXX,IYY)
             W124_SYN(IM,IL)=OW124(IXX,IYY)
             W164_SYN(IM,IL)=OW164(IXX,IYY) 
             W213_SYN(IM,IL)=OW213(IXX,IYY) 
             W138_SYN(IM,IL)=OW138(IXX,IYY) 
             W1100_SYN(IM,IL)=OW1100(IXX,IYY) 
             W412_SYN(IM,IL)=OW412(IXX,IYY) 
             W8p5_SYN(IM,IL)=OW8p5(IXX,IYY) 
         LandSea_Flag_new(IM,IL)= Land_Sea_Flag(IXX,IYY)  
         LAT(IM,IL)=  OLAT(IXX,IYY)
         LON(IM,IL)=  OLON(IXX,IYY)
         SatZen(IM,IL)= OVIEW(IXX,IYY)
         SolZen(IM,IL)= OSOLA(IXX,IYY) 
         SatAz(IM,IL)=OAZIM_View(IXX,IYY)
         SolAz(IM,IL)=OAZIM_Solar(IXX,IYY)
         Height(IM,IL)=OHeight(IXX,IYY) 
         New_CldMsk_500_Land(IM,IL)  =CldMsk_500_Land(IXX,IYY) 
         New_Q_cirrus(IM,IL)  = quality_cirrus(IXX,IYY) 
         SD_for_Dust(IM,IL)  = SD_3by3(IXX,IYY)
!    UV VAriables   
          W354_SYN(IM,IL)=OW354(IXX,IYY) 
          W388_SYN(IM,IL)=OW388(IXX,IYY) 
          else 
          W470_SYN(IM,IL)=-99999.
          W550_SYN(IM,IL)=-99999.
          W659_SYN(IM,IL)=-99999.
          W865_SYN(IM,IL)=-99999. 
          W124_SYN(IM,IL)=-99999.  
          W164_SYN(IM,IL)=-99999.  
          W213_SYN(IM,IL)=-99999.   
          W138_SYN(IM,IL)=-99999.  
          W1100_SYN(IM,IL)=-99999. 
          W412_SYN(IM,IL)=-99999. 
          W8p5_SYN(IM,IL)=-99999. 
          SD_for_Dust(IM,IL)=-99999.
          LandSea_Flag_new(IM,IL)= Land_Sea_Flag(IXX,IYY)  
          LAT(IM,IL)=  OLAT(IXX,IYY)
          LON(IM,IL)=  OLON(IXX,IYY)
          SatZen(IM,IL)= OVIEW(IXX,IYY)
          SolZen(IM,IL)= OSOLA(IXX,IYY)
          SatAz(IM,IL)=OAZIM_View(IXX,IYY)
          SolAz(IM,IL)=OAZIM_Solar(IXX,IYY)
          Height(IM,IL)=OHeight(IXX,IYY) 
          New_CldMsk_500_Land(IM,IL)  =CldMsk_500_Land(IXX,IYY) 
          New_Q_cirrus(IM,IL)  = -99999 
 !    UV VAriables  
          W354_SYN(IM,IL)=-99999
          W388_SYN(IM,IL)=-99999   
       endif    
           ENDDO   
           ENDDO 
                IL = 0
           DO      IYY = START_line*2+1,END_line*2 
                     IL=IL+1
                     IM=0
            DO      IXX= 1,l1b_nLines*2  
                        IM= IM+1
             New_CldMsk_500_Ocean(IM,IL)= CldMsk_500_Ocean(IXX,IYY)
             New_savecldmask(IM,IL)  =     savecldmask(IXX,IYY)  
             enddo
             enddo    
           
          RETURN
          END  
          
          
 