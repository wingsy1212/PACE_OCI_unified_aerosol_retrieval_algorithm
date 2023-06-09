      Subroutine Ret_Uv_SSa(MTHET0,MTHET,MPHI,Mode_F,Mode_C,wind,&
         WSPEED,Small_m_weighting,Tau_550,set_flag_read,&
         REFW354,REFW388,IWS1,IWS2,Lat_OMI,Lon_OMI,&
         EXTSMALL_550,EXTbig_550,Tau_allwave,WAVE_MODIS,Omega_new,tau_new,&
         Index_Omega_new,Height_indx,Fitting_error,Iscan,idata,ref_allwav_uv,numscan)
         USE Linear_interpolation 
        
         implicit none  
         SAVE 
        INCLUDE 'mod04.inc'
        INCLUDE 'Set_Array_dimension_UV.inc'
        integer ii
!         
                       X(1:100)= -9999
                       Y(1:100)=  -9999
                       Z(1:100)= -9999
                       x0(1:100)= -9999    
                       Omega_new(1:NWAV_uv+1) = -9999
                       Fitting_error(1:NWAV_uv)  = -9999
                       Index_Omega_new(1:NWAV_uv)  = -9999
                       Tau_new(1:NWAV_uv)  =  -9999
                       Height_indx(1:NWAV_uv)= -9999 
                       Tau_height(1:Num_heigt,1:NWAV_uv)    =   -9999
                       Omega_height(1:Num_heigt,1:NWAV_uv)  =   -9999 
                       Error_height(1:Num_heigt,1:NWAV_uv)  =   -9999
              
                       
            ref_allwav_uv(1) =   REFW354
            ref_allwav_uv(2) =   REFW388 
                   Iopt_read =   2  
                   Mode_Coarse =  Mode_C-4  
                   
                   
!           print*,'uv',MTHET0,MTHET,MPHI,Mode_F,Mode_C,REFW354,REFW388,&
!           Lat_OMI,Lon_OMI
! reading Look up table once only                   
           If( set_flag_read .EQ.1) then  
           call READ_LOOK_Uv(PHC,THET,THET0,& 
           H1_AINTS_uv1,H1_AINTS_uv2,H2_AINTS_uv1,H2_AINTS_uv2,&
           H3_AINTS_uv1,H3_AINTS_uv2,H4_AINTS_uv1,H4_AINTS_uv2,TAUAS,WAVE,&  
           H1_AINTB_uv1,H1_AINTB_uv2,H2_AINTB_uv1,H2_AINTB_uv2,&
           H3_AINTB_uv1,H3_AINTB_uv2,H4_AINTB_uv1,H4_AINTB_uv2,&
           TAUAB,JPHI,ref_rayall_uv,Iomega,&
           EXTSMALL_save,ALBEDOSMALL_save,EXTbig_save,ALBEDOBIG_save,Iheigt,Num_iomega,&
           Num_heigt)   
!   endif for reading LUT once.           
            Endif
              
!  Followig If statements checks if measured angles are
!   out of bounds from lookup table.
               
             
     IF((MTHET0 .GE. MINMTHET0 .AND. MTHET0 .LE. MAXMTHET0.AND.&
          MTHET   .GE. MINMTHET  .AND. MTHET  .LE. MAXMTHET .AND.&
          MPHI    .GE. MINMPHI   .AND. MPHI   .LE. MAXMPHI).and. &
            (REFW354.gt.0 .AND. REFW388 .gt. 0))  THEN
          
! Interpolating for measured angles         
                
                do Iheigt = 1,Num_heigt
                do Iomega = 1,Num_iomega 
                Array_min(Iomega)= -9999
! Finding indexes for Angles                
         CALL SET_index_inter(MTHET0,MTHET,MPHI,THET0,&
             KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1) 
! Interpolating for measured angles Rayleigh             
          Call Uv_INTANGLE_ray(PHC,THET,THET0,ref_rayall,&
        MTHET0,MTHET,MPHI,Ref_ray,KSZAM1,KSZAP1,KTHEM1,KTHEP1,&
        KPHIM1,KPHIP1,IWS1,IWS2,WSPEED,Wind,Iopt_read,ref_rayall_uv,iomega,Iheigt)
     
! moving around the dimensions for set Coarse and Fine mode   
       Do Iws = 1,Nwind  
        Do  ITAU = 2,NTAU 
         DO  ITH0=1,NTH0
           DO  ITH=1,NTHET
              DO  IPHI=1,NPHI   
          If(Iheigt .eq.1 )  then
        AINTS_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,1)=H1_AINTS_uv1(iomega,IWS,IPHI,ITH,ITH0,ITAU,MODE_F) 
        AINTS_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,2)=H1_AINTS_uv2(iomega,IWS,IPHI,ITH,ITH0,ITAU,MODE_F) 
        AINTB_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,1)=H1_AINTB_uv1(iomega,IWS,IPHI,ITH,ITH0,ITAU,Mode_Coarse)
        AINTB_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,2)=H1_AINTB_uv2(iomega,IWS,IPHI,ITH,ITH0,ITAU,Mode_Coarse)
         ELseif(Iheigt .eq.2 ) then
        AINTS_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,1)=H2_AINTS_uv1(iomega,IWS,IPHI,ITH,ITH0,ITAU,MODE_F) 
        AINTS_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,2)=H2_AINTS_uv2(iomega,IWS,IPHI,ITH,ITH0,ITAU,MODE_F) 
        AINTB_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,1)=H2_AINTB_uv1(iomega,IWS,IPHI,ITH,ITH0,ITAU,Mode_Coarse)
        AINTB_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,2)=H2_AINTB_uv2(iomega,IWS,IPHI,ITH,ITH0,ITAU,Mode_Coarse)
        ELseif(Iheigt .eq.3 ) then
        AINTS_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,1)=H3_AINTS_uv1(iomega,IWS,IPHI,ITH,ITH0,ITAU,MODE_F) 
        AINTS_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,2)=H3_AINTS_uv2(iomega,IWS,IPHI,ITH,ITH0,ITAU,MODE_F) 
        AINTB_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,1)=H3_AINTB_uv1(iomega,IWS,IPHI,ITH,ITH0,ITAU,Mode_Coarse)
        AINTB_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,2)=H3_AINTB_uv2(iomega,IWS,IPHI,ITH,ITH0,ITAU,Mode_Coarse)
        ELseif(Iheigt .eq.4 ) then
        AINTS_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,1)=H4_AINTS_uv1(iomega,IWS,IPHI,ITH,ITH0,ITAU,MODE_F) 
        AINTS_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,2)=H4_AINTS_uv2(iomega,IWS,IPHI,ITH,ITH0,ITAU,MODE_F) 
        AINTB_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,1)=H4_AINTB_uv1(iomega,IWS,IPHI,ITH,ITH0,ITAU,Mode_Coarse)
        AINTB_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,2)=H4_AINTB_uv2(iomega,IWS,IPHI,ITH,ITH0,ITAU,Mode_Coarse)
        Endif
            ENDDO
             ENDDO
              ENDDO
                ENDDO
                ENDDO
        
          CALL UV_INTANGLE(PHC,THET,THET0,Adummy,TAUAS,&
        Adummy,TAUAB,MTHET0,MTHET,MPHI,REFSMALL,REFBIG,&
        Ref_ray,KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1,&
        IWS1,IWS2,WSPEED,Wind,Iopt_read,AINTS_uv,AINTB_uv,iomega,iscan,idata,Numscan)
            
!     REFSMALL and   REFBIG have been chosen for coarse and fine mode in previous steps
!     REFSMALL(1,IWAV,ITAU) and  REFBIG(1,IWAV,ITAU)  are right  in dimension 1
            
          do Iwav = 1,Nwav_Uv
              Ref_Ray_Uv(Iwav) =Ref_ray(Iwav)
               do ITAU = 1,Ntau
            REFSMALLL_Uv(IWAV,ITAU)     =    REFSMALL(1,IWAV,ITAU)     
            REFBIGL_uv(IWAV,ITAU)       =    REFBIG(1,IWAV,ITAU)
              enddo
          enddo
! Compute optical depth for UV channels from extinction coefficients
!            
         call  COMPUTE_alltau_uv(EXTSMALL_save,EXTbig_save,Small_m_weighting,tau_550,&
         Mode_F,Mode_Coarse, EXTSMALL_550,EXTbig_550,Num_iomega,Num_heigt,&
        TAU_COMPUTED_uv,iomega,ALBEDOSMALL_save,ALBEDOBIG_save,Omega_COMPUTED_uv,Iheigt) 
         
        
           Min_Uv = -9990 
          
 ! Minimise the function between measured and lookup table         
        call Retrive_Uv_minfun(REFSMALLL_uv,REFBIGL_uv,Small_m_weighting,&
        Tau_550,tau55um_LUT,Ref_Ray_Uv,ref_allwav_uv,Min_Uv,ALXWAV_uv)
 !             
               do Iwav = 1,Nwav_Uv
               luT_uv(Iomega,Iwav)= ALXWAV_uv(iwav)
               enddo
 !           
                 
              if( Min_Uv .ge.0 .and. Min_Uv .lt.Threshold_error) then
               Array_min(Iomega)= Min_Uv  
              else
              Array_min(Iomega)= -9999
              endif 
 
!     enddo for SSA           
            enddo  

                               
 !         IF((Lat_OMI .ge.lat_min .and. Lat_OMI .le. lat_max)&
 !         .and. (Lon_OMI .ge.lon_min .and. Lon_OMI .le. lon_max))then  
             
            
110        format(I5,5f9.3,2i5,2f8.3,10f7.4,4(i4,f10.3),2i4,2f8.4,2i7,f5.3)   
            
            
          
           call INDEXX(Num_iomega,Array_min,INDX )
          
           do iwav = 1,NWAV_uv 
              Tau_height(Iheigt,Iwav)   =    tau_computed_uv(INDX(1),Iwav)
              Omega_height(Iheigt,Iwav) =    Omega_COMPUTED_uv(INDX(1),iwav)
              Error_height(Iheigt,Iwav) =    Array_min(INDX(1)) 
              Index_omega(Iheigt,Iwav)  =    INDX(1)
               Enddo 
! Enddo for iheight              
             ENDDO   
               
           
141        format(i5,5f9.3,2I4,2f5.3,8i6,16f8.3)     
                     
                     
                    do iwav= 1,NWAV_uv  
                    inumber = 0 
                    do Iheigt   = 1, Num_heigt  
                    if(Error_height(Iheigt,Iwav) .ge. 0 .and.&
                       Error_height(Iheigt,Iwav).lt.Threshold_error) then
                       inumber = inumber+1
                     X(inumber)=  Error_height(Iheigt,Iwav)
                     Y(inumber)=  Omega_height(Iheigt,Iwav)
                     Z(inumber)=  Index_omega(Iheigt ,Iwav) 
                     x0(inumber)= Tau_height(Iheigt,Iwav) 
                      endif 
                     Enddo  
                      
                If(inumber .gt.0) then 
                       call INDEXX(inumber,X,INDX )  
                       Omega_new(iwav) = Y(INDX(1))
                       Fitting_error(iwav)  = X(INDX(1))
                       Index_Omega_new(iwav)  =  Z(INDX(1))
                       Tau_new(iwav)  =  X0(INDX(1))
                       Height_indx(iwav)= Height_km(INDX(1)) 
                Else 
                      Omega_new(iwav) = -9999
                      Index_omega_new(iwav)  =  -9999
                      Tau_new(iwav)  =  -9999 
                      ref_allwav_uv(iwav)= -9999
                      Fitting_error(iwav)=-9999
                      Height_indx(iwav)= -9999
                      Mode_F =-9999
                      Mode_C=-9999 
                      Small_m_weighting =-9999
                Endif 
                
                  If(Tau_550.le. 0.2) then
                    Omega_new(iwav) =  -9999
                    Height_indx(iwav) = -9999
                  Endif  
                   Enddo 
          
          If(Tau_550.ge. 0.2 .and.Tau_allwave(1).GT. 0. .AND. Tau_allwave(4) .GT.0. .AND. & 
             Omega_new(1) .ge. 0 .and. Omega_new(2) .gt. 0) THEN
             AE_44and88 =  &
            (ALOG(Tau_allwave(1)/Tau_allwave(4)))/&
                        (ALOG(WAVE_MODIS(4)/WAVE_MODIS(1)))  
          
            
           If (AE_44and88 .le. 1)   then 
               a0 = 0.9420721
               a1 = - 0.0492087
               a2 = 0.81049567
               meanSSA388 = 0.9225
           Elseif (AE_44and88 .GT. 1) then 
              a0 = 0.91439614
              a1 = - 0.0162844
              a2 = 1.7315114 
              meanSSA388 = 0.9112  
           ENDIF        
             DEL_Omega = Omega_new(2)- Omega_new(1)
             SSA550_pred = a0 + (a1*AE_44and88) + (a2*DEL_Omega)
             Omega_new(3) = SSA550_pred +(Omega_new(2) - meanSSA388) 
          ENDIF          
       
       
       
              
                     
                 
                   go to 3000 

                   
! Tumbarumba         -35.708        147.95    
!2020004 03 54 and 04 00 
               lat_min= -40.0
               lat_max= -35.0
               lon_max =150
               lon_min= 147   
        
!Xitun                        24.162        120.617                0.25
!Lulin                        23.469        120.874                0.43
!Kaohsiung                    22.676        120.292                0.29
!  2019252 0548

                   lat_min= 22.0
                   lat_max= 25
                   lon_min= 119  
                   lon_max =121    

                 
!Red Sea and Persian Gulf:
!DEWA_ResearchCentre (            24.767, 55.369
!KAUST_Campus                     22.305, 39.103 
!  2019252 1054   
                 lat_min= 22.1
                 lat_max= 24.5 
                 lon_min=  39.0 
                 lon_max = 55.5
!  2019252 0548 
!Indonesia
!Bandung ( 6.888S,107.610E) 
                  lat_min= 6.0
                  lat_max= 10
                  lon_min= 107.0
                  lon_max =110.0
                  
                  
 

 ! 
!! 2020 174  1436
!Dakar_Belair ( 14.702N, 17.426W)  
                  lat_min= 14
                  lat_max= 15
                  lon_max = -16
                  lon_min=  -20                 
                  
 
            

                        
  
                       


                   
!Eilat                 29.503                34.917    the VIIRS retrieval are pretty far away 
!Nicosia               35.141                33.381                    0.3
!Karachi               24.946                67.136                    0.55   
! 174 0936  0930 0754 0748
                  lat_min= 24.5
                  lat_max= 35.5
                  lon_max =67.5  
                  lon_min=32.5   
                  


               
 
                  

                 

!   2018 Nov 10  2018 314 1036 
!Red Sea and Persian Gulf:
!DEWA_ResearchCentre (            24.767, 55.369
!KAUST_Campus                     22.305, 39.103    
                 lat_min= 22.1
                 lat_max= 24.5 
                 lon_min=  39.0 
                 lon_max = 55.5    
                 
!India:
!   2018 Nov 10  2018 314  0854
!
!Pune                        18.537, 73.805                Pune is inland, but it is embedded in the same aerosol. You will need to get me values up to 2 degrees longitude to the west (71.8 to 73.8)
!Karunya University          10.935, 76.744                Same here.  At least 2 degrees to the west (74.74 to 76.744)
!Karachi                     24.946, 67.136

                 lat_min= 10.0
                 lat_max= 25.0
                 lon_max =74.0
                 lon_min= 73.0    
                 
                  
! Transported African smoke  2019252 1224
!Simonstown_IMT (-24.193, 18.446)
!Also, -25, 38    No AERONET stations but nice smoke.
                   lat_min= -30.0
                   lat_max= -22.0 
                   lon_min= 3
                   lon_max =35 
                    

                 
!African dust:
!Capo_Verde           16.733,   -  22.935    
! 2018314.1400   2019252.1418                                         

 
                 lat_min= 15.0
                 lat_max= 17.0
                 lon_max = -20
                 lon_min=  -30 
! Cape_San_Juan          18.384                -65.620                    1.1 to 1.5
! La_Parguera            17.970                -67.045                    1.1 to 1.5
! NEON_GUAN              17.970               -66.869                     1.1 to 1.5  
! ragged point           13.165               -59.432 
! 2020 174  1754  1800    
 
                  lat_min= 13.165
                  lat_max= 18.5
                  lon_max = -59 
                  lon_min=  -67.0   
                
               
                   
!Red Sea and Persian Gulf:
!DEWA_ResearchCentre (            24.767, 55.369
!KAUST_Campus                     22.305, 39.103 
!  2020172  11 18 11 12
                 lat_min= 22.1
                 lat_max= 24.5 
                 lon_min=  39.0 
                 lon_max = 55.5   
                   
 !Transported South American smoke 2019252 1548
 ! -27.4   -40.1   There is no AERONET station there, but there is smoke and we can see what it looks like.
                  lat_min= -30.0
                  lat_max= -27.0    
                  lon_max = -40.0
                  lon_min=  -42.0   
                   
!India:
! 2020 174  0748
!
!Pune                        18.537, 73.805                Pune is inland, but it is embedded in the same aerosol. You will need to get me values up to 2 degrees longitude to the west (71.8 to 73.8)
!Karunya University          10.935, 76.744                Same here.  At least 2 degrees to the west (74.74 to 76.744)
!Karachi                     24.946, 67.136

                 lat_min= 10.0
                 lat_max= 25.0
                 lon_max =77.0
                 lon_min= 73.0                                                                                                                 
                

 

                 
 
  
!   2018 Nov 10  2018 314  2048  
!


!California:  
!Monterey                  36.593, -121.855
!UCSB                      34.415, -119.845
!USC_SEAPRISM              33.564, -118.118

                 lat_min= 33
                 lat_max= 36.8
                 lon_max = -118
                 lon_min=  -122  
                  
                   
 ! Cape_San_Juan          18.384                -65.620                    1.1 to 1.5
! La_Parguera            17.970                -67.045                    1.1 to 1.5
! NEON_GUAN              17.970               -66.869                     1.1 to 1.5  
! ragged point           13.165               -59.432 
! 2020 174  1754  1800    

                  lat_min= 13.165
                  lat_max= 18.5
                  lon_max = -59 
                  lon_min=  -67.0                                                                                             
                IF((Lat_OMI .ge.lat_min .and. Lat_OMI .le. lat_max)&
                 .and. (Lon_OMI .ge.lon_min .and. Lon_OMI .le. lon_max))then 
                if(MTHET0 .gt. 0 .and. Omega_new(1) .gt. 0)&
                 write( 37,144)Iscan,idata,Lat_OMI,Lon_OMI,MTHET0,MTHET,MPHI,&
                  (Omega_new(iwav),iwav =1,3),Tau_new(1),Tau_new(2),&
                   (Tau_allwave(IJ),IJ=1,7),(Fitting_error(IJ),IJ=1,2)  
                Endif
 3000        continue                  
        
!    endif for refl  and angles      
              endif 
143        format(2i6,2f9.3,3f9.2,2I6,4f8.3,3f10.3,2I5,4f9.3)  
144        format(2i6,2f9.3,3f9.2,3f8.3,9f10.3,2f10.4)                        
140        format(5f9.3,3f8.3,2f12.5)  
120        format(5f9.3,f8.3,3i6,f9.5,5f8.4) 
              RETURN
                 END
         SUBROUTINE  READ_LOOK_Uv(PHC,THET,THET0,& 
           H1_AINTS_uv1,H1_AINTS_uv2,H2_AINTS_uv1,H2_AINTS_uv2,&
           H3_AINTS_uv1,H3_AINTS_uv2,H4_AINTS_uv1,H4_AINTS_uv2,TAUAS,WAVE,&  
           H1_AINTB_uv1,H1_AINTB_uv2,H2_AINTB_uv1,H2_AINTB_uv2,&
           H3_AINTB_uv1,H3_AINTB_uv2,H4_AINTB_uv1,H4_AINTB_uv2,&
           TAUAB,JPHI,ref_rayall_uv,Iomega,&
           EXTSMALL_save,ALBEDOSMALL_save,EXTbig_save,ALBEDOBIG_save,Iheigt,Num_iomega,&
           Num_heigt) 
           USE OCIUAAER_Config_Module
           IMPLICIT NONE
      SAVE 
      INCLUDE 'mod04.inc' 
      CHARACTER*132 LINE 
     INTEGER Icase,IFILE,IPHI,ITAU,ITH,ITH0,&
              IWAV,IJ,Iopt_read,NNwave,Num_iomega,Num_heigt
      REAL  EXTSMALL_save(Num_iomega,Num_heigt,Numcases,NWAV_uv)	
      REAL  ALBEDOSMALL_save(Num_iomega,Num_heigt,Numcases,NWAV_uv)  
      REAL  EXTbig_save(Num_iomega,Num_heigt,NUMCASEB,NWAV_uv) 
      REAL  ALBEDOBIG_save(Num_iomega,Num_heigt,NUMCASEB,NWAV_uv) 
      REAL  EXTSMALL(NWIND,Numcases,NWAV),EXTbig(NWIND,NUMCASEB,NWAV)
      REAL  RGSS(NUMCASES),SIGMAS(NUMCASES)
      REAL  RGSB(NUMCASEB),SIGMAB(NUMCASEB)
      REAL  MOMENTSSMALL(NWIND,numcases,4),MOMENTSBIG(NWIND,NUMCASEB,4)
      REAL  CCNSMALL(NWIND,nUMCASES),TR
      REAL  EXTNORSMALL(NUMCASES,NWAV),EXTNORBIG(NUMCASEB,NWAV)
      REAL  BACKSCTTSMALL(NWIND,NUMCASES,NWAV),BACKSCTTBIG(NWIND,NUMCASEB,NWAV)
      REAL  ASSYMSMALL(NWIND,NUMCASES,NWAV),ASSYMBIG(NWIND,NUMCASEB,NWAV)
      REAL  ALBEDOSMALL(NWIND,NUMCASES,NWAV),ALBEDOBIG(NWIND,NUMCASEB,NWAV)
      REAL  ALBEDO_R_SMALL(NWIND,NTH0,NTAU,NWAV,NUMCASES)
      REAL  ALBEDO_R_BIG(NWIND,NTH0,NTAU,NWAV,NUMCASEB)
      REAL  ALBEDO_T_SMALL(NWIND,NTH0,NTAU,NWAV,NUMCASES)
      REAL  ALBEDO_T_BIG(NWIND,NTH0,NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_R_RAY(NWAV,NTH0)
      REAL ALBEDO_T_RAY(NWAV,NTH0)

      INTEGER JPHI(NPHI),IWS,Iomega,Mode_F,Mode_C,Iheigt
      REAL  PHC(NPHI),THET(NTHET),THET0(NTH0),WAVE(NWAV)
      REAL  AINTS(NWIND, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASES)
      REAL  TAUAS(NUMCASES,NWAV,NTAU)
      REAL  AINTB(NWIND, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASEB)
      REAL  TAUAB(NUMCASEB,NWAV,NTAU)
      REAL*8   DUMMY(NTH0)
      REAL  CCNdummy,ref_rayall(NWIND,NPHI,NTHET,NTH0,NWAV)
      REAL  EFFRADSMALL(NUMCASES),EFFRADBIG(NUMCASEB),QSCT
        real H1_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
        real H2_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
        real H3_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
        real H4_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
        real H1_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
        real H2_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
        real H3_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
        real H4_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
        real H1_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
        real H2_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
        real H3_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
        real H4_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
        real H1_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
        real H2_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
        real H3_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
        real H4_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
      Real  ref_rayall_uv(Num_iomega,NWIND,NPHI,NTHET,NTH0,NWAV,4),AA
      CHARACTER*1 cWS,kws,hhh 
      CHARACTER  (len=255) :: file_name ,Extension,file_dir
      
        file_dir =  cfg%Pace_ocean   
       DO  Iheigt = 1, Num_heigt  
          Write(hhh,'(i1)' )Iheigt 
         do   Iomega = 1,4 
            Write(kws,'(i1)' )Iomega   
       	  Do IWS = 1,NWIND
          WRITE(cWS, '(i1)' )IWS  
           IFILE = 11  
           
!       '/data4/smattoo/Pace/LookupTable/PACE_LUT/Height_'//hhh//'/Omega'//kWS//'/small_v'//cWS//'c1.dat.v0', &
!        status='old')    
            
!     if( Sat_flag .eq."MODIS") then    
!                open(IFILE, FILE = trim(pkg_root) // &
!       '/Tables/MODIS_VIIRS_GLI_LUTs/small_v'//cWS//'c1.dat.v6',& 
!        status='old')            
           
           Extension = 'Height_'//hhh//'/Omega'//kWS//'/small_v'//cWS//'c1.dat.v0'
           file_name= trim(file_dir)//trim(Extension) 
           OPEN (IFILE, FILE =trim(file_name),status='old')  
           
!             print*,'trim(file_name)',trim(file_name)

!         OPEN (IFILE, FILE = trim(file_name) //'Height_'//hhh//'/Omega'&
!         //kWS//'/small_v'//cWS//'c1.dat.v0',status='old') 
                                
! read rayleigh data  
           NNwave =  NWAV_UV
        DO 200 IWAV=1,NNwave 
!
! Read aerosol information

                  read(IFILE,500)LINE
                  read(IFILE,505)WAVE(IWAV)

! Read ocean information
         DO IJ =1,4
          READ(IFILE,500)LINE
        ENDDO
! Read Atmosphere information

                  READ(IFILE,500)LINE
                  READ(IFILE,515)TR,TAUAS(1,IWAV,1)
        DO IJ = 1,2
         READ(IFILE,500)LINE
        ENDDO
!
! Read all other information
        DO IJ = 1,2
          READ(IFILE,500)LINE
        ENDDO

           READ(IFILE,525)(THET0(IJ),IJ=1,NTH0)

!
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(ALBEDO_R_RAY(IWAV,IJ),IJ=1,NTH0)
                  READ(IFILE,530)(ALBEDO_T_RAY(IWAV,IJ),IJ=1,NTH0)
!
                  READ(IFILE,500)LINE
           

            DO 300 ITH0=1,NTH0
                     READ(IFILE,500)LINE
                     IF(ITH0.EQ.1)THEN
                       READ(IFILE,535)(THET(ITH),ITH=1,NTHET)
                     ELSE
                       DO IJ=1,3
                       READ(IFILE,500)LINE
                       ENDDO
                   ENDIF
!
             DO  400 IPHI=1,NPHI
               READ(IFILE,540)JPHI(IPHI),&
              (ref_rayall(IWS,Iphi,ITH,ITH0,IWAV),ITH=1,NTHET)
             Do ITH=1,NTHET
              ref_rayall_uv(iomega,IWS,Iphi,ITH,ITH0,IWAV,Iheigt)= ref_rayall(IWS,Iphi,ITH,ITH0,IWAV)
             enddo
  400       CONTINUE
!    enddo for ith0
 300       CONTINUE
!    enddo for wav
 200       CONTINUE

! LOOP IS AROUND THET0,NEW FILE FOR EACH SMALL CASE
         DO 10 Icase =1,NUMCASES
         DO 20 IWAV=1,NNwave
! Leave itau=1 to fill up with taua=0.0 later in subroutine INTANGLE
         DO 30 ITAU = 2,NTAU
! Read aerosol information

                  read(IFILE,500)LINE
                  read(IFILE,505)WAVE(IWAV)
                   
        DO IJ = 1,3
          READ(IFILE,500)LINE 
        ENDDO

       READ(IFILE,515)RGSS(Icase),SIGMAS(Icase)
       READ(IFILE,500)LINE
       READ(IFILE,515)EFFRADSMALL(Icase)
       READ(IFILE,520)MOMENTSSMALL(IWS,Icase,1),MOMENTSSMALL(IWS,Icase,2)
       READ(IFILE,520)MOMENTSSMALL(IWS,Icase,3),MOMENTSSMALL(IWS,Icase,4)
       READ(IFILE,515)ALBEDOSMALL(IWS,Icase,IWAV),ASSYMSMALL(IWS,Icase,IWAV)
       READ(IFILE,515)CCNSMALL(IWS,Icase),BACKSCTTSMALL(IWS,Icase,IWAV)
       READ(IFILE,520)QSCT,EXTSMALL(IWS,Icase,IWAV)
! ExT coefficient  and albedos are independent of windspeed so taking one wspeed.       
      EXTSMALL_save(Iomega,Iheigt,Icase,iwav)  = EXTSMALL(1,Icase,iwav)
      ALBEDOSMALL_save(Iomega,Iheigt,Icase,iwav)= ALBEDOSMALL(1,Icase,IWAV)
          
! Read ocean information
         DO IJ =1,4
          READ(IFILE,500)LINE
        ENDDO
! Read Atmosphere information

                  READ(IFILE,500)LINE
                  READ(IFILE,515)TR,TAUAS(Icase,IWAV,ITAU)
                   
        DO IJ = 1,2
         READ(IFILE,500)LINE
        ENDDO
!
! Read all other information
        DO IJ = 1,2
          READ(IFILE,500)LINE
        ENDDO

           READ(IFILE,525)(THET0(IJ),IJ=1,NTH0)  
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)  
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
        READ(IFILE,530)(ALBEDO_R_SMALL(IWS,IJ,ITAU,IWAV,Icase),IJ=1,NTH0)
        READ(IFILE,530)(ALBEDO_T_SMALL(IWS,IJ,ITAU,IWAV,Icase),IJ=1,NTH0)
!
                  READ(IFILE,500)LINE


            DO 40 ITH0=1,NTH0
                     READ(IFILE,500)LINE
                     IF(ITH0.EQ.1)THEN
                       READ(IFILE,535)(THET(ITH),ITH=1,NTHET)
                     ELSE
                       DO IJ=1,3
                       READ(IFILE,500)LINE
                       ENDDO
                   ENDIF
!
             DO  50 IPHI=1,NPHI
               READ(IFILE,540)JPHI(IPHI),&
               (AINTS(IWS,IPHI,ITH,ITH0,ITAU,IWAV,Icase),ITH=1,NTHET)  
   50       continue
!    enddo for ith0
   40       CONTINUE
!    enddo for tau
   30       CONTINUE
! enddo for nwav
! Fill the array for albedo and transmission for all cases tau=0
      TAUAS(Icase,IWAV,1)=TAUAS(1,IWAV,1)
      DO IJ=1,NTH0
      ALBEDO_R_SMALL(IWS,IJ,1,IWAV,Icase)=ALBEDO_R_RAY(IWAV,IJ)
      ALBEDO_T_SMALL(IWS,IJ,1,IWAV,Icase)=ALBEDO_T_RAY(IWAV,IJ)
      ENDDO

   20       CONTINUE 
               
!    enddo for num of size distribution for small
   10       CONTINUE
                  DO IPHI=1,NPHI
                     PHC(IPHI)=FLOAT(JPHI(IPHI))
                  ENDDO
! ENDDO for Windspeed   
          REWIND IFILE
        ENDDO
             Do IWS = 1,NWIND 
             Do Icase =1,NUMCASES  
             Do  ITAU = 2,NTAU 
             DO  ITH0=1,NTH0
             DO  ITH=1,NTHET
             DO  IPHI=1,NPHI 
         IF(Iheigt .eq.1) then    
         H1_AINTS_uv1(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) = &
         AINTS(IWS,IPHI,ITH,ITH0,ITAU,1,Icase)
        H1_AINTS_uv2(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) = &
         AINTS(IWS,IPHI,ITH,ITH0,ITAU,2,Icase)
         ENDIF
         IF(Iheigt .eq.2) then    
         H2_AINTS_uv1(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) = &
         AINTS(IWS,IPHI,ITH,ITH0,ITAU,1,Icase)
         H2_AINTS_uv2(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) = &
         AINTS(IWS,IPHI,ITH,ITH0,ITAU,2,Icase)
         ENDIF
         IF(Iheigt .eq.3) then    
         H3_AINTS_uv1(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) = &
         AINTS(IWS,IPHI,ITH,ITH0,ITAU,1,Icase)
         H3_AINTS_uv2(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) = &
         AINTS(IWS,IPHI,ITH,ITH0,ITAU,2,Icase)
         ENDIF
         IF(Iheigt .eq.4) then    
         H4_AINTS_uv1(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) = &
         AINTS(IWS,IPHI,ITH,ITH0,ITAU,1,Icase)
         H4_AINTS_uv2(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) = &
         AINTS(IWS,IPHI,ITH,ITH0,ITAU,2,Icase)
         ENDIF
         
          enddo
            enddo
              enddo
                enddo
                enddo
                enddo 
         
!      Enddo for Omega                
        Enddo  
!      Enddo for Height              
         Enddo    
         
          
       DO  Iheigt = 1, Num_heigt  
! READ LARGE CASES OF SIZEDISTRIBUTION
          Write(hhh,'(i1)' )Iheigt 
         Do  Iomega = 1,4 
         Write(kws,'(i1)' )Iomega 
         Do IWS = 1,NWIND
         WRITE(cWS, '(i1)' )IWS
         IFILE = 12 
          Extension = 'Height_'//hhh//'/Omega'//kWS//'/big_v'//cWS//'c1.dat.v01' 
           file_name= trim(file_dir)//trim(Extension) 
           OPEN (IFILE, FILE = trim(file_name),status='old')  
          
 
          DO 60 Icase =1,NUMCASEB
         DO 70 IWAV=1,NNwave
! Leave itau=1 to fill up with taua=0.0 later in subroutine INTANGLE
         DO 80 ITAU = 2,NTAU
! Read aerosol information

                  read(IFILE,500)LINE
                  read(IFILE,505)WAVE(IWAV)
        DO IJ = 1,3
          READ(IFILE,500)LINE
        ENDDO

       READ(IFILE,515)RGSB(Icase),SIGMAB(Icase)
       READ(IFILE,500)LINE
       READ(IFILE,515)EFFRADBIG(Icase)
       READ(IFILE,520)MOMENTSBIG(IWS,Icase,1),MOMENTSBIG(IWS,Icase,2)
       READ(IFILE,520)MOMENTSBIG(IWS,Icase,3),MOMENTSBIG(IWS,Icase,4)
       READ(IFILE,515)ALBEDOBIG(IWS,Icase,IWAV),ASSYMBIG(IWS,Icase,IWAV)
       READ(IFILE,515)CCNdummy,BACKSCTTBIG(IWS,Icase,IWAV)
       READ(IFILE,520)QSCT,EXTBIG(IWS,Icase,IWAV)
       EXTbig_save(Iomega,Iheigt,Icase,iwav)  = EXTBIG(1,Icase,IWAV)
       ALBEDOBIG_save(Iomega,Iheigt,Icase,iwav)= ALBEDOBIG(1,Icase,IWAV)
        
! Read ocean information
         DO IJ =1,4
          READ(IFILE,500)LINE
        ENDDO
! Read Atmosphere information

                  READ(IFILE,500)LINE
                  READ(IFILE,515)TR,TAUAB(Icase,IWAV,ITAU)
        DO IJ = 1,2
         READ(IFILE,500)LINE
        ENDDO
!
! Read all other information
        DO IJ = 1,2
          READ(IFILE,500)LINE
        ENDDO

          READ(IFILE,525)(THET0(IJ),IJ=1,NTH0)

!
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0) 
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
       READ(IFILE,530)(ALBEDO_R_BIG(IWS,IJ,ITAU,IWAV,Icase),IJ=1,NTH0)
       READ(IFILE,530)(ALBEDO_T_BIG(IWS,IJ,ITAU,IWAV,Icase),IJ=1,NTH0)
!
                  READ(IFILE,500)LINE


            DO 90 ITH0=1,NTH0
                     READ(IFILE,500)LINE
                     IF(ITH0.EQ.1)THEN
                       READ(IFILE,535)(THET(ITH),ITH=1,NTHET)
                     ELSE
                       DO IJ=1,3
                       READ(IFILE,500)LINE
                       ENDDO
                   ENDIF
!

             DO  100 IPHI=1,NPHI
               READ(IFILE,540)JPHI(IPHI),&
            (AINTB(IWS,IPHI,ITH,ITH0,ITAU,IWAV,Icase),ITH=1,NTHET) 
 100       continue
!    enddo for ith0
  90      CONTINUE
!    enddo for tau
  80       CONTINUE
!   enddo for wav
!Fill the array for albedo and transmission for all cases tau=0
      TAUAB(Icase,IWAV,1)=TAUAS(1,IWAV,1)
      DO IJ=1,NTH0
      ALBEDO_R_BIG(IWS,IJ,1,IWAV,Icase)=ALBEDO_R_RAY(IWAV,IJ)
      ALBEDO_T_BIG(IWS,IJ,1,IWAV,Icase)=ALBEDO_T_RAY(IWAV,IJ)
      ENDDO
  70      CONTINUE 
  
!    enddo for num of size distribution for large
  60       CONTINUE
! Enddo for windspeed  close (IFILE)
           
           ENDDO 
                
             Do IWS = 1,NWIND 
             Do Icase =1,NUMCASEB  
             Do  ITAU = 2,NTAU 
             DO  ITH0=1,NTH0
             DO  ITH=1,NTHET
             DO  IPHI=1,NPHI
        IF(Iheigt .eq.1) then 
        H1_AINTB_uv1(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) =&
           AINTB(IWS,IPHI,ITH,ITH0,ITAU,1,Icase)
        H1_AINTB_uv2(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase)=&
         AINTB(IWS,IPHI,ITH,ITH0,ITAU,2,Icase) 
       ENDIF 
       IF(Iheigt .eq.2) then 
        H2_AINTB_uv1(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) =&
           AINTB(IWS,IPHI,ITH,ITH0,ITAU,1,Icase)
        H2_AINTB_uv2(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase)=&
         AINTB(IWS,IPHI,ITH,ITH0,ITAU,2,Icase) 
       ENDIF 
       IF(Iheigt .eq.3) then 
        H3_AINTB_uv1(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) =&
           AINTB(IWS,IPHI,ITH,ITH0,ITAU,1,Icase)
        H3_AINTB_uv2(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase)=&
         AINTB(IWS,IPHI,ITH,ITH0,ITAU,2,Icase) 
       ENDIF 
       IF(Iheigt .eq.4) then 
        H4_AINTB_uv1(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase) =&
           AINTB(IWS,IPHI,ITH,ITH0,ITAU,1,Icase)
        H4_AINTB_uv2(Iomega,IWS,IPHI,ITH,ITH0,ITAU,Icase)=&
         AINTB(IWS,IPHI,ITH,ITH0,ITAU,2,Icase) 
       ENDIF 
         enddo
            enddo
              enddo
                enddo
                enddo
                enddo 
  ! enddo for Omega
         ENDDO  
  ! Enddo for Heights        
          ENDDO   
               
            
!
!     format statements
!
500               format(132a1)
505               format(t32,f6.4)
510               format(t32,f6.4,t65,e11.4)
!515               format(t32,f6.4,t70,f6.4)
515               format(t32,f6.4,t69,f6.4)
520               format(t26,e12.4,t64,e12.4)
525               format(t12,f6.1,3x,6(2x,f5.1,3x)/&
                         t12,f6.1,3x,6(2x,f5.1,3x))
530               format(t10,7e10.3/t10,7e10.3)
535               format(t10,f4.1,6(5x,f5.1)/t9,7(f5.1,5x)/&
                         t9,7(f5.1,5x))
540               format(i4,1x,7e10.3/5x,7e10.3/5x,7e11.3)

       RETURN
       END
       
       Subroutine COMPUTE_alltau_uv(EXTSMALL_save,EXTbig_save,Small_m_weighting,tau_550,&
       Finemode,Mode_Coarse,EXTSMALL_550,EXTbig_550,Num_iomega,Num_heigt,&
       TAU_COMPUTED_uv,iomega,ALBEDOSMALL_save,ALBEDOBIG_save,Omega_COMPUTED_uv,Iheigt)  
       IMPLICIT NONE
      SAVE 
      INCLUDE 'mod04.inc'
      INTEGER Finemode,Mode_Coarse,IWAV 
      Integer iJ,NUM,iomega,Iheigt,Num_iomega,Num_heigt 
      REAL  EXTSMALL_save(Num_iomega,Num_heigt,Numcases,NWAV_uv),&
           ALBEDOSMALL_save(Num_iomega,Num_heigt,Numcases,NWAV_uv),&
           EXTBIG_save(Num_iomega,Num_heigt,NUMCASEB,NWAV_uv),&
           ALBEDOBIG_save(Num_iomega,Num_heigt,NUMCASEB,NWAV_uv)
      REAL EXTNORSMALL(NWAV_uv),EXTNORBIG(NWAV_uv) 
      Real TAU_COMPUTED_uv(Num_iomega,NWAV_uv),tau_550,Small_m_weighting
      REAL  EXTSMALL_550,EXTbig_550
      Real  Omega_COMPUTED_uv(Num_iomega,NWAV_uv)
      
           
!   EXTSMALL_save is taken for one wind speed because it is independent of windspeed
               
            DO  IWAV=1,NWAV_uv
              EXTNORSMALL(IWAV) = EXTSMALL_save(iomega,Iheigt,Finemode,IWAV)/EXTSMALL_550
              EXTNORBIG(IWAV)   = EXTBIG_save(iomega,Iheigt,Mode_Coarse,IWAV)/EXTBIG_550
            Enddo
              
! Changes 2/28
         
        
        
         DO IWAV = 1,NWAV_uv
         TAUX55SMALL=tau_550*Small_m_weighting
         TAUX55BIG=tau_550-TAUX55SMALL
        TAU_COMPUTED_uv(Iomega,IWAV)=TAUX55SMALL*EXTNORSMALL(IWAV)&
                           +(TAUX55BIG*EXTNORBIG(IWAV))  
       Omega_COMPUTED_uv(Iomega,IWAV)=(Small_m_weighting* &
                           ALBEDOSMALL_save(iomega,Iheigt,Finemode,IWAV))&
       +((1.-Small_m_weighting)*ALBEDOBIG_save(iomega,Iheigt,Mode_Coarse,IWAV))
       ENDDO
       
!          print*,'inside',Iomega,(TAU_COMPUTED_uv(Iomega,IWAV),IWAV =1,2)
         
! end of  Changes 2/28
        RETURN
      END
        Subroutine  Retrive_Uv_minfun(REFSMALLL_uv,REFBIGL_uv,XMIN,&
      Tau_550,tau55um_LUT,Ref_Ray_Uv,ref_allwav_uv,Min_Uv,ALXWAV_uv)
       USE Linear_interpolation 
       IMPLICIT NONE 
       SAVE
      INCLUDE 'mod04.inc' 
      INTEGER ITAU,IWAV,LOPT,NJTAU,NJWAV,NUMWAV,INDEX_wave
      REAL  RADWAV,Y1,asum1,Ref_Ray_Uv(NWAV)
      REAL XMIN,XS(100,100),YS(100,100),X(100),Y(100),Denom
      Real  REFSMALLL_Uv(NWAV,NTAU),REFBIGL_uv(NWAV,NTAU)
      Real  ref_allwav_uv(NWAV),tau55um_LUT(NTAU),Tau_550
      Real  ALXWAV_uv(NWAV),FMIN,Min_Uv 
         ASUM=0.0
         ASUM1=0.0  
         Min_Uv =0 
 
!                   ******COMPUTE FUN TO BE MINIMIZED.
!

       DO 100 IWAV = 1,NWAV_UV
       DO 101 ITAU = 1,NTAU 
         YS(ITAU,IWAV)=(REFSMALLL_uv(IWAV,ITAU)*XMIN)&
                 +((1.-XMIN)*REFBIGL_uv(IWAV,ITAU)) 
101    CONTINUE
100    CONTINUE
 
!
!                 ******* FOR TAUX55 COMPUTE Reflectance FOR
!                    ALL WAVELENGTHS.
!
      DO 104 IWAV = 1,NWAV_uv
         DO 103 ITAU = 1,NTAU
            X(ITAU)=tau55um_LUT(ITAU)
            Y(ITAU)=YS(ITAU,IWAV)  
             
103      CONTINUE
         y1=0
          CALL INTERP(NTAU,Tau_550,X,Y,Y1)
            ALXWAV_uv(IWAV)=Y1  
104   CONTINUE 
       DO 105 IWAV = 1,Nwav_UV
               INDEX_wave=IWAV 
                Denom=(ref_allwav_uv(INDEX_wave)-Ref_Ray_Uv(INDEX_wave))+0.01 
 !              Denom= ref_allwav_uv(INDEX_wave)
          if( Denom .lt. 0.01)Denom=0.01 
          ASUM1= ASUM1+(((ref_allwav_uv(INDEX_wave)-ALXWAV_uv(INDEX_wave))&
        /(Denom))**2) 
!          print*,'ref',ref_allwav_uv(INDEX_wave),ALXWAV_uv(INDEX_wave)
105   CONTINUE
              Min_Uv=SQRT(ASUM1/REAL(Nwav_UV)) 
               
      RETURN
      END

  
!*********************************************************************

      SUBROUTINE UV_INTANGLE_ray(PHC,THET,THET0,ref_rayall,&
        MTHET0,MTHET,MPHI,Ref_ray,KSZAM1,KSZAP1,KTHEM1,KTHEP1,&
        KPHIM1,KPHIP1,IWS1,IWS2,WSPEED,Wind,Iopt_read,ref_rayall_uv,iomega,Iheigt)


      USE Linear_interpolation 
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      INTEGER IJ,IPHI,ITH,ITH0,IWAV,LOPT,Iheigt
      REAL  PHC(NPHI),THET(NTHET),THET0(NTH0)
      REAL  Ref_rayall( NWIND,NPHI, NTHET, NTH0,NWAV)
      REAL  Ref_rayall_uv(4,NWIND,NPHI, NTHET, NTH0,NWAV,4)
      REAL  Ref_ray(NWAV)
      REAL  MTHET0,MTHET,MPHI
      REAL  X(100),Y(100),XX1(100),YY1(100),XX2(100),YY2(100),Y1
      REAL  XX3(100),YY3(100),WSPEED,Wind(NWIND)
      INTEGER KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1
      INTEGER LL,MM,NN,iomega
      INTEGER KK,IWS1,IWS2,IWS,NNwave,Iopt_read
     
! LOOP IS AROUND THET0,NEW FILE FOR EACH THETA0
!
          if(Iopt_read .eq. 1) NNwave =  NWAV
          if(Iopt_read .eq. 2) NNwave =  NWAV_Uv
 
          DO IJ = 1,100
          X(IJ)=0.0
          Y(IJ)=0.0
          XX1(IJ)=0.0
          YY1(IJ)=0.0
          XX2(IJ)=0.0
          YY2(IJ)=0.0
          XX3(IJ)=0.0
          YY3(IJ)=0.0
          ENDDO
 
           DO 10 IWAV=1,NNwave 
          KK=0
          DO 20 IWS = IWS1,IWS2
          LL=0
          DO 30 ITH0 = KSZAM1,KSZAP1
          MM=0
          DO 40 ITH= KTHEM1,KTHEP1
          NN=0
          DO 50 IPHI =KPHIM1,KPHIP1
          NN=NN+1
          X(NN)=PHC(IPHI) 
          Y(NN)=Ref_rayall_uv(iomega,IWS,IPHI,ITH,ITH0,IWAV,Iheigt) 
 50       CONTINUE
           Y1=0.0
           CALL INTERP(NN,MPHI,X,Y,Y1)
            MM=MM+1
           XX1(MM)=THET(ITH)
           YY1(MM)=Y1 
 40         CONTINUE
            y1=0.0
           CALL INTERP(MM,MTHET,XX1,YY1,Y1)
             LL=LL+1
            XX2(LL)=THET0(ITH0)
            YY2(LL)=Y1 
 30        CONTINUE
            y1=0.0
           CALL INTERP(LL,MTHET0,XX2,YY2,Y1)  
             KK=KK+1
             XX3(KK)=WIND(IWS)
             YY3(KK)=Y1
 20        CONTINUE
           CALL INTERP(KK,WSPEED,XX3,YY3,Y1) 
           Ref_ray(IWAV)=(Y1*PI)/(COS(DTR*MTHET0))  
 10       CONTINUE
            
            RETURN
            END



!*********************************************************************

      SUBROUTINE UV_INTANGLE(PHC,THET,THET0,AINTS,TAUAS,&
                 AINTB,TAUAB,MTHET0,MTHET,MPHI,REFSMALL,REFBIG,&
              Ref_ray,KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1,&
             IWS1,IWS2,WSPEED,Wind,Iopt_read,AINTS_uv,AINTB_uv,iomega,&
             iscan,idata,Numscan)

      USE Linear_interpolation 
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      CHARACTER*132 LINE
      CHARACTER*45  LINE1
      INTEGER Icase,IJ,IPHI,ISIZE,ITAU,ITH,ITH0,IWAV,LOPT,NUMCASE
      INTEGER ISCAN,IDATA,Numscan
      REAL  PHC(NPHI),THET(NTHET),THET0(NTH0)
      REAL  AINTS(NWIND, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASES)
      REAL  TAUAS(NUMCASES,NWAV,NTAU)
      REAL  AINTB(NWIND, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASEB)
      REAL  TAUAB(NUMCASEB,NWAV,NTAU)
      REAL  Ref_ray(NWAV)
      REAL  REFSMALL(NUMCASES,NWAV,NTAU),REFBIG(NUMCASEB,NWAV,NTAU)
      REAL  MTHET0,MTHET,MPHI, WSPEED,Wind(NWIND)
      REAL  X(100),Y(100),XX1(100),YY1(100),XX2(100),YY2(100),Y1
       REAL  XX3(100),YY3(100)
      INTEGER KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1
      INTEGER LL,MM,NN,NNwave,iomega
      INTEGER OO,IWS,IWS1,IWS2,Iopt_read
      REAL AINTS_uv(4,NWIND, NPHI, NTHET, NTH0, NTAU, NWAV)
      REAL AINTB_uv(4,NWIND, NPHI, NTHET, NTH0, NTAU, NWAV) 
    

! LOOP IS AROUND THET0,NEW FILE FOR EACH THETA0
!
           
             NNwave =  NWAV_Uv
          
             

          Y1=0.0
          DO IJ = 1,100
          X(IJ)=0.0
          Y(IJ)=0.0
          XX1(IJ)=0.0
          YY1(IJ)=0.0
          XX2(IJ)=0.0
          YY2(IJ)=0.0
          XX3(IJ)=0.0
          YY3(IJ)=0.0
          ENDDO
    
          DO 5 ISIZE = 1,2
          if(Iopt_read .eq. 1)then
          IF(ISIZE .EQ.1)NUMCASE=NUMCASES
          IF(ISIZE .EQ.2)NUMCASE=NUMCASEB
          else
          IF(ISIZE .EQ.1)NUMCASE=1
          IF(ISIZE .EQ.2)NUMCASE=1
          endif
          DO 10 Icase =1,NUMCASE  
          DO 20 IWAV=1,NNwave
! interpolate starting from optical thickess indexed 2.
! at the end of routine fill the array with interpoltaed ta=0.0
          DO 30  ITAU  = 2,NTAU
             OO=0
          DO 35 IWS = IWS1, IWS2 
             NN=0
          DO  40 ITH0  = KSZAM1,KSZAP1
             MM=0
          DO  50  ITH  = KTHEM1,KTHEP1
              LL=0
          DO 60  IPHI  = KPHIM1,KPHIP1
           LL=LL+1
          X(LL)=PHC(IPHI) 
         IF( ISIZE.EQ.1)Y(LL)=AINTS_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,IWAV)
         IF( ISIZE.EQ.2)Y(LL)=AINTB_uv(iomega,IWS,IPHI,ITH,ITH0,ITAU,IWAV)  
 60       CONTINUE
           CALL INTERP(LL,MPHI,X,Y,Y1)
           MM=MM+1
           XX1(MM)=THET(ITH)
           YY1(MM)=Y1
 50         CONTINUE
            y1=0.0
           CALL INTERP(MM,MTHET,XX1,YY1,Y1)
            NN=NN+1
            XX2(NN)=THET0(ITH0)
            YY2(NN)=Y1
 40        CONTINUE
             y1=0.0
            CALL INTERP(NN,MTHET0,XX2,YY2,Y1) 
            OO=OO+1
            XX3(OO)=WIND(IWS)
            YY3(OO)=Y1
 35        CONTINUE
             y1=0.0
            CALL INTERP(OO,WSPEED,XX3,YY3,Y1)

!                  ****** change to Reflectance units
       IF(ISIZE.EQ.1)REFSMALL(Icase,IWAV,ITAU)=(Y1*PI)/(COS(DTR*MTHET0))
       IF(ISIZE.EQ.2)REFBIG(Icase,IWAV,ITAU)=(Y1*PI)/(COS(DTR*MTHET0))
        
          
  30       CONTINUE
  20       CONTINUE
  10       CONTINUE
   5        CONTINUE
           
!
! fill empty array of itau=1 with rayleigh (taua=0.0)
!
        DO   ISIZE= 1,2
          if(Iopt_read .eq. 1)then
          IF(ISIZE .EQ.1)NUMCASE=NUMCASES
          IF(ISIZE .EQ.2)NUMCASE=NUMCASEB
          else
          IF(ISIZE .EQ.1)NUMCASE=1
          IF(ISIZE .EQ.2)NUMCASE=1
          endif
              DO   Icase =1,NUMCASE
              DO   IWAV=1,NNwave
            IF(ISIZE.EQ.1)REFSMALL(Icase,IWAV,1)=Ref_ray(IWAV)
            IF(ISIZE.EQ.2)REFBIG(Icase,IWAV,1) = Ref_ray(IWAV)
             Enddo
             Enddo
         Enddo   
         
            RETURN
            END


      

  