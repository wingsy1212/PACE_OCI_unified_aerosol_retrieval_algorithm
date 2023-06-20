         module compute_center_GEOLOC_ANGLE
 
      implicit none

contains   
! ---------------------------------------------------------------------
! SET_INDEX    
        
        SUBROUTINE GEOLOC_ANGLE(LAT,LON,SatZen,SatAz,SolZen,&
      SolAz,Height,MTHET0,MTHET,MPHI0,MPHI,MDPHI,MSCATT,MHGHT,&
      START_1KM,Lat_center,Lon_center,iscan,idata) 
      IMPLICIT NONE
      SAVE

      include 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'
        REAL diff, sumdif
      INTEGER START_1KM,flag
      INTEGER IX1,IX2,IY1,IY2,I,iscan,idata,n1,n2
      REAL Lat_center,Lon_center,Lon_Min,Lon_Max,Lon_4(4)
      Real SolAz_4(4),SolAz_Min,SolAz_Max
      Real SatAz_4(4),SatAz_Min,SatAz_Max,MHGHT 
      REAL MTHET0,MTHET,MPHI0,MPHI,MDPHI,MSCATT,AVE1,AVE2
      REAL Lat(ISWATH_B,ILINE),Lon(ISWATH_B,ILINE),&
           SatZen(ISWATH_B,ILINE),SolZen(ISWATH_B,ILINE),&
           SatAz(ISWATH_B,ILINE),SolAz(ISWATH_B,ILINE),&
           Height(ISWATH_B,ILINE)
! set if line is even or odd. If Even pick middle points and average     
        If ( mod(Iline,2) .gt. 0)then
              Flag =1
              n1 = (iline/2)
              n2 = (iline/2)
        Else 
              Flag = 2
              n1 = (iline/2)-1
              n2 = (iline/2)
        endif
        
           
        
       if( Iline .eq.1) then
       IX1=START_1KM  
       IY1=Iline  
         Lat_center = LAT(IX1,IY1) 
         Lon_center = LON(IX1,IY1)
         MTHET0= SolZen(IX1,IY1)
         MTHET = SatZen(IX1,IY1) 
! convert to kilometers         
         MHGHT = Height(IX1,IY1)* 0.001 
          MPHI0 = SolAz(IX1,IY1) 
      if ( MPHI0  .gt. +180. ) MPHI0 = MPHI0  - 360.
      if ( MPHI0  .lt. -180. ) MPHI0 = MPHI0 + 360. 
           MPHI = SatAz(IX1,IY1)  
      if (  MPHI  .gt. +180. )  MPHI =  MPHI  - 360.
      if (  MPHI  .lt. -180. )  MPHI =  MPHI +  360.
        MDPHI = abs(MPHI0 - MPHI -180.0)
       IF(MDPHI.GT.360.0) MDPHI=amod(MDPHI,360.0)
       IF(MDPHI.GT.180.0) MDPHI=360.0-MDPHI 
         IF(MTHET0.GT.0.0.AND.MTHET.GT.0.0.AND.MDPHI.GT.0.0) THEN
        MSCATT = -COS(MTHET0*DTR)*COS(MTHET*DTR)&
                +SIN(MTHET0*DTR)*SIN(MTHET*DTR)&
                *COS(MDPHI*DTR)
        MSCATT = ACOS(MSCATT)*RTD
          ELSE
        MSCATT=-99.99
          ENDIF 
!ENDIF for iline    
       Else
        if(Flag .EQ.2) then    
        IX1=START_1KM+n1
        IX2=START_1KM+n2
         IY1=ILINE/2
         IY2=ILINE/2+1
        endif
      if(Flag.EQ.1) then 
         IX1=START_1KM+n1
         IX2=START_1KM+n2
         IY1=ILINE/2+1
         IY2=ILINE/2+1 
        Endif
         
        
!
! Finding minimum and maximum of longitudes
!

      Lon_4(1)=LON(IX1,IY1)
      Lon_4(2)=LON(IX2,IY1)
      Lon_4(3)=LON(IX1,IY2)
      Lon_4(4)=LON(IX2,IY2)
      Lon_Min=Lon_4(1)
      Lon_Max=Lon_4(1)
      
      
      SolAz_4(1)=SolAz(IX1,IY1)
      SolAz_4(2)=SolAz(IX2,IY1) 
      SolAz_4(3)=SolAz(IX1,IY2)
      SolAz_4(4)=SolAz(IX2,IY2) 
      SolAz_Min=SolAz_4(1) 
      SolAz_Max=SolAz_4(1)
      
      SatAz_4(1)=SatAz(IX1,IY1)
      SatAz_4(2)=SatAz(IX2,IY1) 
      SatAz_4(3)=SatAz(IX1,IY2)
      SatAz_4(4)=SatAz(IX2,IY2) 
      SatAz_Min=SatAz_4(1) 
      SatAz_Max=SatAz_4(1) 
      
      DO I=1,4
        IF(Lon_4(I).LE.Lon_Min) Lon_Min=Lon_4(I)
        IF(Lon_4(I).GE.Lon_Max) Lon_Max=Lon_4(I)
        IF(SolAz_4(I).LE.SolAz_Min) SolAz_Min=SolAz_4(I)
        IF(SolAz_4(I).GE.SolAz_Max) SolAz_Max=SolAz_4(I)
        IF(SatAz_4(I).LE.SatAz_Min) SatAz_Min=SatAz_4(I)
        IF(SatAz_4(I).GE.SatAz_Max) SatAz_Max=SatAz_4(I)
      ENDDO

!
! If Lon_Max <-180 means fill values are found and if Lon_Min <-180 means at least one'
! fill value is found, under these two condition fill value is set to longitude,
! latitude, and all the angles
!
      IF(Lon_Max.LT.-180.0.OR.Lon_Min.LT.-180.0) THEN

        Lon_center= FV_GEO
        Lat_center= FV_GEO
        MTHET0=-FV_GEO
        MTHET=FV_GEO
        MPHI0=FV_GEO
        MPHI=-FV_GEO
        MSCATT=FV_GEO
      ELSE

!
! Otherwise, check for other condisitons that may contain fill value
!
      sumdif = 0.

!.....take longitude difference of pixels 2 through 4 relative to pixel 1
      do i = 2,4
          diff = Lon_4(i) - Lon_4(1)

!..........we are working on a sphere, take shortest arc between points
          if ( diff  .gt.  180. ) diff = diff   - 360.
          if ( diff  .lt.  -180. ) diff = 360. + diff

           sumdif = sumdif + diff
      enddo

!.....mathematically equivalent to {Lon_4(1)+Lon_4(2)+Lon_4(3)+Lon_4(4)} / 4
!    using shortest arcs
      lon_center = Lon_4(1) + sumdif/4.

      if ( lon_center  .gt. +180. ) lon_center = lon_center  - 360.
      if ( lon_center  .lt.   -180. ) lon_center = lon_center + 360.

      Lat_center=(LAT(IX1,IY1)+LAT(IX2,IY1)+&
                 LAT(IX1,IY2)+LAT(IX2,IY2))/4.
         MTHET0=(SolZen(IX1,IY1)+SolZen(IX2,IY1)+&
              SolZen(IX1,IY2)+SolZen(IX2,IY2))/4.  
   
          MTHET =(SatZen(IX1,IY1)+SatZen(IX2,IY1)+&
              SatZen(IX1,IY2)+SatZen(IX2,IY2))/4. 
              
              if(MTHET .EQ. 0 )MTHET = MTHET + 0.001
              
             
             MHGHT =(Height(IX1,IY1)+Height(IX2,IY1)+&
                Height(IX1,IY2)+Height(IX2,IY2))/4000.
 !       print*,'iline',iline,flag,ix1,ix2,iy1,iy2,Height(IX1,IY1),Height(IX2,IY1),Height(IX1,IY2),&
!        Height(IX2,IY2)
             
!        if(Height(IX1,IY1) .ge. 0  .and. Height(IX2,IY1).ge.0 .and. &
!          Height(IX1,IY2) .ge. 0 .and. Height(IX2,IY2) .ge.0)then 
!     mhight converted to kilometers and averaged 4 pixels     
!            MHGHT =(Height(IX1,IY1)+Height(IX2,IY1)+&
!               Height(IX1,IY2)+Height(IX2,IY2))/4000.
!            else
!            if(Height(IX1,IY1) .ge. 0)MHGHT = Height(IX1,IY1)/1000.
!            if(Height(IX2,IY1) .ge. 0)MHGHT = Height(IX2,IY1)/1000.
!            if(Height(IX1,IY2) .ge. 0)MHGHT = Height(IX1,IY2)/1000.
!            if(Height(IX2,IY2) .ge. 0)MHGHT = Height(IX2,IY2)/1000.
!            endif

          

    
              
             
             
      
! Otherwise, check for other condisitons that may contain fill value
!
      sumdif = 0.

!.....take longitude difference of pixels 2 through 4 relative to pixel 1
      do i = 2,4
          diff = SolAz_4(i) - SolAz_4(1)

!..........we are working on a sphere, take shortest arc between points
          if ( diff  .gt.  180. ) diff = diff - 360.
          if ( diff  .lt.  -180. ) diff = 360. + diff

           sumdif = sumdif + diff
      enddo

!..... 
!    using shortest arcs

        MPHI0 = SolAz_4(1) + sumdif/4.

      if ( MPHI0  .gt. +180. ) MPHI0 = MPHI0  - 360.
      if ( MPHI0  .lt. -180. ) MPHI0 = MPHI0 + 360. 
        
        
      
      
!
      sumdif = 0.

 
      do i = 2,4
          diff = SatAz_4(i) - SatAz_4(1)

!..........we are working on a sphere, take shortest arc between points
          if ( diff  .gt.  180. ) diff = diff - 360.
          if ( diff  .lt.  -180. ) diff = 360. + diff

           sumdif = sumdif + diff
      enddo

 
!   using shortest arcs

        MPHI = SatAz_4(1) + sumdif/4.

      if (  MPHI  .gt. +180. )  MPHI =  MPHI  - 360.
      if (  MPHI  .lt. -180. )  MPHI =  MPHI +  360.
      
 
            
 
       MDPHI = abs(MPHI0 - MPHI -180.0)
       IF(MDPHI.GT.360.0) MDPHI=amod(MDPHI,360.0)
       IF(MDPHI.GT.180.0) MDPHI=360.0-MDPHI

      IF(MTHET0.GT.0.0.AND.MTHET.GT.0.0.AND.MDPHI.GT.0.0) THEN
        MSCATT = -COS(MTHET0*DTR)*COS(MTHET*DTR)&
                +SIN(MTHET0*DTR)*SIN(MTHET*DTR)&
               *COS(MDPHI*DTR)
        MSCATT = ACOS(MSCATT)*RTD
      ELSE
        MSCATT=-99.99
      ENDIF 
! endif for lon_max.........     
      ENDIF  
! Endif for ILINE  
       ENDIF  
         
        RETURN
      end  subroutine GEOLOC_ANGLE
           end module  compute_center_GEOLOC_ANGLE
           
