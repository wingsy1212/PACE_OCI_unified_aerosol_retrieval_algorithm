!  
 module merge_land_ocean
 
! Declare use of Fortran 90 netCDF package
! If this part doesn't compile: module load sips/gcc

use netcdf

implicit none

contains

! ---------------------------------------------------------------------
!
! Begin subroutine make_nc(ALL VARIABLES GO HERE,XL,YL)
!
! ---------------------------------------------------------------------

      subroutine  merge_land_ocean_VAR(grn_lwmask,CldMsk_Native_Ocean,dbdt_cld,&
                   Land_sea_flag,Ret_ref_ocean,dbdt_refl,Ret_tau_ocean,uvdbdtaod,& 
                   Ret_Tau_LandOcean,ret_ref_LandOcean,Cldmask_Native_LandOcean,&
                   Cloud_Frac_LandOcean,Ret_Xtrack,Ret_Lines,&
                   Ret_ocean_Quality, Ret_Quality_LandOcean,Ret_Quality_LandOcean_W0,&
                   Ret_CLDFRC_land_DT,Ret_CLDFRC_ocean)
                              
                                      
!     Define input variable dimensions, etc.
!     Define parameters from Main_Driver.f90

      include 'output_Variables.inc' 
      Integer IX,IY,IXX,IYY,clear,cloudy,jjx,jjy,IL
      INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE :: grn_lwmask
      REAL,DIMENSION(:,:,:),ALLOCATABLE  :: uvdbdtaod 
      REAL,DIMENSION(:,:,:),ALLOCATABLE  :: dbdt_refl
      REAL,DIMENSION(:,:),ALLOCATABLE  :: dbdt_cld
       
    !                 
               Ret_Tau_LandOcean(1:Ret_Xtrack,1:Ret_Lines,1:9) = -9999  
               Ret_Quality_LandOcean(1:Ret_Xtrack,1:Ret_Lines) = -9999 
               ret_ref_LandOcean(1:Ret_Xtrack,1:Ret_Lines,1:9) = -9999  
               Cldmask_Native_LandOcean(1:IX_B,1:IY_B) = -9999 
               Cloud_Frac_LandOcean(1:Ret_Xtrack,1:Ret_Lines) = -9999 
               Ret_Quality_LandOcean_W0(1:Ret_Xtrack,1:Ret_Lines) = -9999 
               
                                       
 
! computing Cloud fraction.......
             
             Do  IY = 1,Ret_Lines  
             Do  IX =  1,Ret_Xtrack  
            IF(Land_sea_flag(IX,IY) .eq. 0) &     
             Cloud_Frac_LandOcean(IX,IY)= Ret_CLDFRC_ocean(IX,IY)
            IF(Land_sea_flag(IX,IY) .eq. 1) &     
              Cloud_Frac_LandOcean(IX,IY) = Ret_CLDFRC_land_DT(IX,IY) 
               Enddo
               Enddo
                
          ! Quality Flag     reverse the order from MODIS convention (0 = best now)
          !  
               Do  IY = 1,Ret_Lines 
                 Do  IX =  1,Ret_Xtrack  
        IF(Land_sea_flag(IX,IY) .eq. 0 .and. Ret_tau_ocean(IX,IY,4) .ge. -0.01) then 
            if( Ret_ocean_Quality(IX,IY) .eq.3) Ret_Quality_LandOcean(IX,IY)= 0 
            if( Ret_ocean_Quality(IX,IY) .eq.2) Ret_Quality_LandOcean(IX,IY)= 1 
            if( Ret_ocean_Quality(IX,IY) .eq.1) Ret_Quality_LandOcean(IX,IY)= 2
            if( Ret_ocean_Quality(IX,IY) .eq.0) Ret_Quality_LandOcean(IX,IY)= 3  
        ELSE 
           IF(Land_sea_flag(IX,IY) .eq. 1 .AND. uvdbdtaod(IX,IY,4) .GE.-0.05)&
           Ret_Quality_LandOcean(IX,IY) = 0   
       Endif
                 ENDDO
            ENDDO

! converted 0= ocean 1= land   in        
                    
                 
             Do  IL = 1,9
             Do  IX = 1,Ret_Xtrack
             Do  IY = 1,Ret_Lines
!     REport optical depth Ocean        
             IF(Land_sea_flag(IX,IY) .eq. 0 )then 
             Ret_Tau_LandOcean(IX,IY,IL)= Ret_tau_ocean(IX,IY,IL)
! ! setting quality for singlescattering albedo  based on optical depth at 0.55um   
!If QA_AOD = 0 and AOD > 0.3, then QA_SSA = 0.
!If QA_AOD = 0 and AOD > 0.2 but AOD < 0.3, then QA_SSA = 1 
!If QA_AOD = 1 and AOD > 0.2, then QA_SSA = 1 
!If QA_AOD = 2 and AOD > 0.2, then QA_SSA = 2 
!If QA_AOD = 3 and AOD > 0.2, then QA_SSA = 3
            if( Ret_Quality_LandOcean(IX,IY) .eq.0.and. Ret_tau_ocean(IX,IY,4) .gt. 0.2) &
                Ret_Quality_LandOcean_W0(IX,IY)= 0
            if( Ret_Quality_LandOcean(IX,IY) .eq.0 .and.(Ret_tau_ocean(IX,IY,4) .ge.  0.2 &
               .and. Ret_tau_ocean(IX,IY,4)  .lt.  0.3)) Ret_Quality_LandOcean_W0(IX,IY)=1 
          if( Ret_Quality_LandOcean(IX,IY) .eq. 1 .and. Ret_tau_ocean(IX,IY,4).gt. 0.2) &
                Ret_Quality_LandOcean_W0(IX,IY)= 1 
            if( Ret_Quality_LandOcean(IX,IY) .eq.2 .and. Ret_tau_ocean(IX,IY,4).gt. 0.2) &
                Ret_Quality_LandOcean_W0(IX,IY)= 2   
            if( Ret_Quality_LandOcean(IX,IY) .eq.3 .and. Ret_tau_ocean(IX,IY,4).gt. 0.2) &
                Ret_Quality_LandOcean_W0(IX,IY)= 3  
             Elseif(Land_sea_flag(IX,IY) .eq. 1)then 
             If( IL .le.5)Ret_Tau_LandOcean(IX,IY,IL) = uvdbdtaod(IX,IY,IL) 
             Else
             Ret_Tau_LandOcean(IX,IY,IL) = -9999
             Endif
             Enddo 
             Enddo  
             Enddo
             
                    
                       
             Do  IL = 1,9
             Do  IX = 1,Ret_Xtrack
             Do  IY = 1,Ret_Lines
      IF(Land_sea_flag(IX,IY) .eq. 0)THEN
         Ret_ref_LandOcean(IX,IY,IL)= Ret_ref_ocean(IX,IY,IL)  
      ElseIF(Land_sea_flag(IX,IY) .eq. 1 .AND. dbdt_refl(IX,IY,IL) .GE.0)THEN 
             Ret_ref_LandOcean(IX,IY,IL) =  dbdt_refl(IX,IY,IL) 
       Endif
          ENDDO 
         ENDDO  
         ENDDO   
         Return  
end subroutine   merge_land_ocean_VAR    
end module merge_land_ocean

 



