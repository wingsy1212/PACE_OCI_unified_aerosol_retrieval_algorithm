        SUBROUTINE PROCESS_Land(HANDLE_LUT466,HANDLE_LUT553,HANDLE_LUT644,&
        HANDLE_LUT213,HANDLE_LUTMAP,IMONTH,ISCAN,IDATA,MTHET0,&
        MTHET,MDPHI,MHGHT,Lat_center,Lon_center,START_500,END_500,&
        START_250,END_250,START_1KM,END_1KM,W354_SYN,W388_SYN,W470_syn,W550_SYN,&
        W659_syn,W865_syn,W124_SYN,W164_SYN,W213_syn,&
        CldMsk_250,Set_Counter_Land,QA_Flag_Land,Success_Ret_Land,&
        Fail_Ret_Land,SDSLAT,SDSLON,SDS_MTHET0,SDS_MTHET,SDS_MPHI,&
        SDS_Tau_Land_Ocean,CLDMSK_500,SDS_Tau_Land_Ocean_img,&
        SDS_Aerosol_Type,SDS_SCAT_ANGLE_land,SDS_mass_conc_land,&
        SDS_angs_coeff_land,SDS_CLDFRC_land,SDS_dust_weighting,&
        SDS_est_uncer,SDS_RefF_land,&
        SDS_TranF_land,SDS_NUMPIXELS_land,SDSTAU_corrected,&
        SDS_ref_land,SDS_ref_STD_land,SDS_QCONTROL_land,& 
        G_factor,quality_land,Ret_Quality_cirrus,cloud_num_land,&
        SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,Qcontrol_special_land,&
        SDSTAU_corrected_213,Quality_flag_forJoint,SDSTAU_small_land,&
        Save_index,W412_SYN,W443_SYN,Optical_depth_land,&
        SnowMsk_Ratio,Iswath,handle_Urban_Table_10km,ETA)

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'
      INCLUDE 'mod04_land.inc' 
      INCLUDE 'read_Sat_MODIS.inc'
 
 
! Define SDS array for hdf output
 
      Real SDS_MTHET0(NUMCELLS_B),SDS_MTHET(NUMCELLS_B),&
                 SDS_MPHI(NUMCELLS_B)
            REAL   SDSLAT(NUMCELLS_B),SDSLON(NUMCELLS_B)

 
! Combined Land/Ocean SDS Arrays
 

      Real SDS_Tau_Land_Ocean(NUMCELLS_B),&
                  SDS_Tau_Land_Ocean_img(NUMCELLS_B)

 
!  LAND SDS_ARRAYS........
 
      BYTE         SDS_QCONTROL_land(QA_LAND,NUMCELLS_B)
      BYTE         SDS_QCONTROL_CritRef_land(QA_LAND,NUMCELLS_B)
       Real   SDS_Aerosol_Type(NUMCELLS_B),&
                     SDS_SCAT_ANGLE_land(NUMCELLS_B),&
                     SDS_angs_coeff_land(NUMCELLS_B),&
                     SDS_CLDFRC_land(NUMCELLS_B),&
                     SDS_dust_weighting(NUMCELLS_B),&
                     SDS_NUMPIXELS_land(NUMCELLS_B,Band_land+2),&
                     SDSTAU_corrected(NUMCELLS_B,Land_Sol3),&
                     SDS_ref_land(NUMCELLS_B,Band_land+2),&
                     SDS_ref_STD_land(NUMCELLS_B,Band_land+2),&
                     SDSTAU_small_land(NUMCELLS_B,Land_Sol4)
     
       Real       SDS_Surface_Reflectance_Land(NUMCELLS_B,Land_Sol3),&
                  SDS_Fitting_Error_Land(NUMCELLS_B),&
                  SDSTAU_corrected_213(NUMCELLS_B)

       REAL       SDS_mass_conc_land(NUMCELLS_B)

 
! EXTRA LAND SDS_ARRAYS..........FOR LAND Statistics ONLY
 
       Real  SDS_Mean_Reflectance_Land_All(NUMCELLS_B,Land_Sol3),&
        SDS_SDev_Reflectance_Land_All(NUMCELLS_B,Land_Sol3),&
        SDS_Path_Radiance_Land(NUMCELLS_B,Land_Sol1),&
        SDS_Critical_Reflectance_Land(NUMCELLS_B,Land_Sol1),&
        SDS_Error_Crit_Reflectance_Land(NUMCELLS_B,Land_Sol1),&
        SDS_Error_Path_Radiance_Land(NUMCELLS_B,Land_Sol1),&
        SDS_QFlag_Critical_Ref_Land(NUMCELLS_B,Land_Sol1),&
        SDS_QFlag_Path_Radiance_Land(NUMCELLS_B,Land_Sol1)

 
! Obsolete (02/2006) Land SDS Arrays
!
      Real &
                  SDS_est_uncer(NUMCELLS_B,Land_Sol1),&
                  SDS_RefF_land(NUMCELLS_B,Land_Sol2),&
                  SDS_TranF_land(NUMCELLS_B,Land_Sol1)

      INTEGER QA_Flag_Land(19),Success_Ret_Land,Fail_Ret_Land
      INTEGER NUMSQDIM,NUMXBOX,NUMYBOX,Quality_flag_forJoint(2)
      PARAMETER (NUMSQDIM=400)
      INTEGER RTN,START_DATA,END_DATA
      CHARACTER*100 fname
      CHARACTER*20 att_n,dtype
      INTEGER NUMMODEL,NUMTAU,NUMSCATT,LINES
      PARAMETER (NUMMODEL=4,NUMTAU=12,NUMSCATT=12,LINES=100)
      INTEGER IDATA,ISCAN,IT,IERROR,IRD,NUMDATA,IWEIGHT,nms
      INTEGER Ret_Quality_cirrus 
      INTEGER IFINISH,IMONTH,IPR,IPROCE,IRPHSOM,ISMK,ISULP,LINE,ISMK_H,ISMK_M
           INTEGER  IAER
!      INTEGER HANDLE_HDF,IAER,MODFIL_MOD04(MODFILLEN),
!     &        MODFIL_Geo(MODFILLEN),MODFIL_L1B(MODFILLEN),
!     &        MODFIL_Cldmsk(MODFILLEN),HANDLE_TM(5),
!     &        MODFIL_MOD05(MODFILLEN),MODFIL_MOD07(MODFILLEN),
!     &        MODFIL_ANC(MODFILLEN)
      INTEGER Sea_landFSQ(NUMSQDIM,IGRIDX,IGRIDY),cloud_num_land
      INTEGER Sea_landF(IGRIDX,IGRIDY),Set_Counter_Land
      INTEGER START_500,END_500,START_250,END_250,START_1KM,END_1KM 
      INTEGER CLDMSK_250(ISWATH_B,ILINE)
      INTEGER CLDMSK_500(ISWATH_B,ILINE)
       REAL WEIGHT,SCAT_ANGLE_LAND
      INTEGER  Good_pixels_Land(Band_land+2) 
!  L1B Reflectance data
      REAL W470_SYN(ISWATH_B,ILINE)
      REAL W550_SYN(ISWATH_B,ILINE)
      REAL W213_SYN(ISWATH_B,ILINE)
      REAL W124_SYN(ISWATH_B,ILINE)
      REAL W164_SYN(ISWATH_B,ILINE)
      REAL w865_SYN(ISWATH_B,ILINE)
      REAL W659_SYN(ISWATH_B,ILINE)
      Integer SnowMsk_Ratio(ISWATH_B,ILINE)
! NPP channels      
      REAL W412_SYN(ISWATH_B,ILINE)
      REAL W443_SYN(ISWATH_B,ILINE)
      REAL W0P75_SYN(ISWATH_B,ILINE),ss
      
!UV channels
     REAL W354_SYN(ISWATH_B,ILINE)
     REAL W388_SYN(ISWATH_B,ILINE)     

! L2 Reflectance data
      REAL REFW466_L,REFW550_L 
      REAL REFW644_L,REFW124_L,REFW164_L,REFW866_L,REFW212_L
      REAL SDW466_L,SDW550_L
      REAL SDW644_L,SDW124_L,SDW164_L,SDW866_L,SDW212_L
!NPP channels      
      REAL REFW412_L,REFW443_L,REFW0P75_L,SDW412_L,SDW443_L,SDW0P75_L 
!UVchannels      
      REAL REFW354_L,REFW388_L,SDW354_L,SDW388_L      
     
      INTEGER Num_Pixels_Used, quality_land
      REAL SAVERATIO
      REAL THR213MIN, THR213MAX

      REAL MTHET0,MTHET,MDPHI,MHGHT
      REAL LATM,LONM,LAT_CENTER,LON_CENTER
      REAL WATER_VAPOR(NUMSQDIM),OZONE_COL(NUMSQDIM)
      CHARACTER * 5 IOFLAG

      integer ij,ik ,Qcontrol_special_land
      INTEGER NUMRED(ISWATH_B,ILINE),NUMBLUE(ISWATH_B,ILINE)

! MINMTHET0,MAXMTHET0,MINMTHET,MAXMTHET,MINMPHI,MAXMPHI,MAXTAU
 
       REAL G_factor, DEGRAD

       INTEGER NUMCLDFREE,NUMPERBLUE,NUMPERRED,NumHandles
       PARAMETER (NUMCLDFREE=20,NUMPERBLUE=7,NUMPERRED=7)
       REAL FLUXDNBLUE,FLUXDNRED,FLUXUPBLUE,FLUXUPRED

 
 
      INTEGER file_version,pgs_pc_getreference,i

 
! Percentiles used for the radiances averaged 20th and 50th percentile
 

      INTEGER NLOWPERC,NUPPERC
      PARAMETER (NLOWPERC=20,NUPPERC=50)

!  Parameters Direct from the new lookup table (NL0 = New Land Initialize)

      REAL SBAR_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL OPTH_NL0(NLTAU,NLWAV,NLTABLE)
      REAL INT_NL0(NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL FdT_NL0(NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL Fd_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL T_NL0(NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL MASSCOEF_NL0(NLTAU,NLWAV,NLTABLE)
      REAL EXTNORM_NL0(NLTAU,NLWAV,NLTABLE)

!  For Determining aerosol type from aerosol map
      INTEGER nlon, nlat
      PARAMETER (nlon = 360, nlat = 180)
      INTEGER AEROSOL(nlon,nlat)
      Real Optical_depth_land
! retrieved products
      REAL ERR644
      REAL TAULAND55,RHOSFC213
      REAL MASSCON,ANGSTROM
      REAL RHOSFC(NLWAV)
      REAL ETA
      REAL AOD(NLWAV),AODF(NLWAV), AODC(NLWAV)
      INTEGER AVE_COUNT
      INTEGER FTABLE_NL
      INTEGER ETA_FLAG

! dummy
      INTEGER IWAV,IXX,IYY
       
! Rayleigh reflectance interpolated for angle

      REAL REF_RAY_NL(NLWAV)
      REAL Rayleigh_look(2)
      Integer Save_index(NUMCELLS_L,MaxPixels_left_L)
      Real Average_Urban(3601,1801),Urban_per 
      integer handle_Urban_Table_10km,cloudy_pixels
 ! Execution only for the first time
 
          
       
      
       
       IF(Set_Counter_Land.EQ.1) THEN
   Call Read_urban_Table(Average_Urban,handle_Urban_Table_10km)        
              NUMDATA=NUMSQ
! Subroutine RNLOOKUP reads the new look-up tables for land

      CALL RNLOOKUP_NC4(&
          HANDLE_LUT466,HANDLE_LUT553,HANDLE_LUT644,HANDLE_LUT213,&
          INT_NL0,Fd_NL0,T_NL0,OPTH_NL0,&
          SBAR_NL0,MASSCOEF_NL0,EXTNORM_NL0)
!
!  Determine fine mode aeorsol location map
! 
       CALL AEROSOL_MAP(HANDLE_LUTMAP,IMONTH, AEROSOL)
        
      ENDIF
       
!     Initialize retrievals

      ETA = -9999
      ERR644 = -9999
      MASSCON = -9999
      ANGSTROM = -9999
      
      DO IWAV = 1, NLWAV
        RHOSFC(NLWAV) = -9999
        AOD(NLWAV) = -9999
        AODF(NLWAV) = -9999
        AODC(NLWAV) = -9999
      ENDDO
      AVE_COUNT = -9999
      FTABLE_NL = -9999       
      IAER = -9999
      Optical_depth_land = -99999

      START_DATA=1
      END_DATA=NUMDATA
      IFINISH=0
      Qcontrol_special_land=0
 
!-----------------------------------------------------------------------------------------------
! Checking the angle bounds
!-----------------------------------------------------------------------------------------------

        IF(MTHET0 .GE. MINMTHET0 .AND. MTHET0 .LE. MAXMTHET0.AND.& 
          MTHET  .GE. MINMTHET  .AND. MTHET  .LE. MAXMTHET .AND.&
          MDPHI  .GE. MINMPHI   .AND. MDPHI  .LE. MAXMPHI) THEN
 
        CALL  COMPUTE_SCATTANGLE_LAND(MTHET0,MTHET,MDPHI,IDATA,&
        SCAT_ANGLE_LAND)

!-----------------------------------------------------------------------------------------------
! Do path radiance and critical radiance statistics
!-----------------------------------------------------------------------------------------------


        CALL INTANGLE_RAY_NL(MTHET0,MTHET,MDPHI,&
      	  INT_NL0,REF_RAY_NL)

        Rayleigh_look(1) = REF_RAY_NL(iwave_466)  
        Rayleigh_look(2) = REF_RAY_NL(iwave_644)  
          


    
   
       LATM=Lat_center
       LONM=Lon_center
       IPROCE = 0
       Call compute_urban_Percent(Lat_center,Lon_center,Average_Urban,&
       Urban_per)
       
        
!-----------------------------------------------------------------------------------------------
!  Select criterion for detecting dark targets using 2.1 micron channel
!  (reflectance at 2.1 micron between 0.01 and THR213MAX) if the second
!  criterion fails (i.e., IFINISH=0) 
!-----------------------------------------------------------------------------------------------
             
             THR213MIN = 0.01
             THR213MAX = 0.25 
        Call Average_land(ISWATH,CLDMSK_500,W354_SYN,W388_SYN,W470_SYN,W550_SYN,&
       W659_SYN,W865_SYN,W124_SYN,W164_SYN,W213_SYN,START_500,END_500,&
       START_250,END_250,START_1KM,END_1KM,REFW466_L,REFW550_L,REFW644_L,&
       REFW124_L,REFW164_L,REFW866_L,REFW212_L,SDW466_L,SDW550_L,&
       SDW644_L,SDW124_L,SDW164_L,SDW866_L,SDW212_L,IFINISH,&
       THR213MIN,THR213MAX,Num_Pixels_Used,IERROR,NUMRED,NUMBLUE,Save_index,&
       Idata,W412_SYN,W443_SYN,W0P75_SYN,REFW412_L,REFW443_L,REFW0P75_L,&
       SDW412_L,SDW443_L,SDW0P75_L,iscan,Good_pixels_Land,SnowMsk_Ratio,&
       cloudy_pixels,REFW354_L,REFW388_L,SDW354_L,SDW388_L)
     
             If(IFINISH.EQ.1)  IPROCE = 1
      
         
!-----------------------------------------------------------------------------------------------
! only continental model Procedure =2
!-----------------------------------------------------------------------------------------------
            
            IF(IFINISH.EQ.0) THEN  
                G_factor=1./COS(DTR*MTHET)+1./SQRT(COS(DTR*MTHET0)) 
                G_factor=0.5*G_factor 
               THR213MIN= 0.25  
               THR213MAX= 0.25*G_factor 
               
         IF(THR213MAX.LE.40.0) THEN
         Call Average_land(ISWATH,CLDMSK_500,W354_SYN,W388_SYN,W470_SYN,W550_SYN,&
        W659_SYN,W865_SYN,W124_SYN,W164_SYN,W213_SYN,START_500,END_500,&
       START_250,END_250,START_1KM,END_1KM,REFW466_L,REFW550_L,REFW644_L,&
        REFW124_L,REFW164_L,REFW866_L,REFW212_L,SDW466_L,SDW550_L,&
        SDW644_L,SDW124_L,SDW164_L,SDW866_L,SDW212_L,IFINISH,&
        THR213MIN,THR213MAX,Num_Pixels_Used,IERROR,NUMRED,NUMBLUE,Save_index,&
        Idata,W412_SYN,W443_SYN,W0P75_SYN,REFW412_L,REFW443_L,REFW0P75_L,&
       SDW412_L,SDW443_L,SDW0P75_L,iscan,Good_pixels_Land,SnowMsk_Ratio,&
       cloudy_pixels,REFW354_L,REFW388_L,SDW354_L,SDW388_L) 
                   If(IFINISH.EQ.1)  IPROCE = 2
          
          ENDIF
!  Endif for  IPROCE=2
            ENDIF

           
         
!-----------------------------------------------------------------------------------------------
!
         If(IFINISH.EQ.1) Then
              IERROR=0
         CALL PROCESS_Land_Inver(IPROCE,&
         MTHET0,MTHET,MDPHI,SCAT_ANGLE_LAND,IMONTH,LATM,LONM,MHGHT,&
         REFW466_L,REFW550_L,REFW644_L,REFW866_L,&
         REFW124_L,REFW164_L,REFW212_L,&
         SDW466_L,SDW550_L,SDW644_L,SDW866_L,&
         SDW124_L,SDW164_L,SDW212_L,&
         INT_NL0,Fd_NL0,T_NL0,OPTH_NL0,SBAR_NL0,AEROSOL,REF_RAY_NL,&
         ETA,ETA_FLAG,AOD,ERR644,RHOSFC,AVE_COUNT,FTABLE_NL,&
         MASSCOEF_NL0,MASSCON,&
         EXTNORM_NL0,AODF,AODC,ANGSTROM,iscan,idata,Urban_per)
         
! Set IAER     
            IAER = FTABLE_NL
            
          
!     If 2.13 Um not met procedure and error is set     
         ELSE
           IERROR=4
           IPROCE=0
           IAER =-9999
             
! ENDIF for Ifinish     
           ENDIF 
          
             Optical_depth_land = AOD(2)
         
!
! Output parameters to HDF file
!

      CALL OUTPUT(LATM,LONM,IAER,IDATA,IERROR,&
        REFW466_L,REFW550_L,REFW644_L,REFW866_L,&
        REFW124_L,REFW164_L,REFW212_L,SDW466_L,SDW550_L,&
        SDW644_L,SDW124_L,SDW164_L,SDW866_L,SDW212_L,IPROCE,& 
        SCAT_ANGLE_LAND,ANGSTROM,MTHET0,MTHET,MDPHI,&
        QA_Flag_Land,Success_Ret_Land,Fail_Ret_Land,&
        FLUXDNBLUE,FLUXDNRED,FLUXUPBLUE,FLUXUPRED,MASSCON,&
        SDSLAT,SDSLON,SDS_MTHET0,SDS_MTHET,SDS_MPHI,&
        SDS_Tau_Land_Ocean_img,&
        SDS_Aerosol_Type,SDS_SCAT_ANGLE_land,SDS_mass_conc_land,&
        SDS_angs_coeff_land,SDS_CLDFRC_land,SDS_dust_weighting,&
       SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,SDS_NUMPIXELS_land,&
       SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,&
       SDS_QCONTROL_land,quality_land,Ret_Quality_cirrus,cloud_num_land,&
       SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,&
       ETA,ETA_FLAG,AOD,AODF,RHOSFC,&
        SDSTAU_corrected_213,SDSTAU_small_land,&
       ERR644,Num_Pixels_Used,Qcontrol_special_land,Quality_flag_forJoint,&
       REFW354_L,REFW388_L,REFW0P75_L,SDW354_L,SDW388_L,SDW0P75_L,&                         
       Good_pixels_Land,cloudy_pixels) 
           Optical_depth_land= SDSTAU_corrected(idata,2) 
             
          
 
 

 
!                       **Following else for 'ifang'

      ELSE
!                        **Write the output if the angles are out of
!                          bounds from the lookup table for different
!                         options with error indicators.
!
        IERROR=1
       
      
        CALL OUTPUT(LATM,LONM,IAER,IDATA,IERROR,&                         
       REFW466_L,REFW550_L,REFW644_L,REFW866_L,&                         
       REFW124_L,REFW164_L,REFW212_L,SDW466_L,SDW550_L,&                         
       SDW644_L,SDW124_L,SDW164_L,SDW866_L,SDW212_L,IPROCE,&                         
       SCAT_ANGLE_LAND,ANGSTROM,MTHET0,MTHET,MDPHI,&                         
       QA_Flag_Land,Success_Ret_Land,Fail_Ret_Land,&                         
       FLUXDNBLUE,FLUXDNRED,FLUXUPBLUE,FLUXUPRED,MASSCON,&                         
       SDSLAT,SDSLON,SDS_MTHET0,SDS_MTHET,SDS_MPHI,&                         
       SDS_Tau_Land_Ocean_img,&                         
       SDS_Aerosol_Type,SDS_SCAT_ANGLE_land,SDS_mass_conc_land,&                         
       SDS_angs_coeff_land,SDS_CLDFRC_land,SDS_dust_weighting,&                         
      SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,SDS_NUMPIXELS_land,&                         
      SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,&                         
      SDS_QCONTROL_land,quality_land,Ret_Quality_cirrus,cloud_num_land,&                         
      SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,&                         
      ETA,ETA_FLAG,AOD,AODF,RHOSFC,&                         
       SDSTAU_corrected_213,SDSTAU_small_land,&                         
      ERR644,Num_Pixels_Used,Qcontrol_special_land,Quality_flag_forJoint,&                         
      REFW354_L,REFW388_L,REFW0P75_L,SDW354_L,SDW388_L,SDW0P75_L,&                         
       Good_pixels_Land,cloudy_pixels) 
! Endif for checking the angle bounds
 
      ENDIF 
  
   1  FORMAT(132A1)
      RETURN
      END



!***********************************************************************
      SUBROUTINE OUTPUT(LATM,LONM,IAER,IDATA,IERROR,&                         
       REFW466_L,REFW550_L,REFW644_L,REFW866_L,&                         
       REFW124_L,REFW164_L,REFW212_L,SDW466_L,SDW550_L,&                         
       SDW644_L,SDW124_L,SDW164_L,SDW866_L,SDW212_L,IPROCE,&                         
       SCAT_ANGLE_LAND,ANGSTROM,MTHET0,MTHET,MDPHI,&                         
       QA_Flag_Land,Success_Ret_Land,Fail_Ret_Land,&                         
       FLUXDNBLUE,FLUXDNRED,FLUXUPBLUE,FLUXUPRED,MASSCON,&                         
      SDSLAT,SDSLON,SDS_MTHET0,SDS_MTHET,SDS_MPHI,&                         
      SDS_Tau_Land_Ocean_img,&                         
       SDS_Aerosol_Type,SDS_SCAT_ANGLE_land,SDS_mass_conc_land,&                         
       SDS_angs_coeff_land,SDS_CLDFRC_land,SDS_dust_weighting,&                         
      SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,SDS_NUMPIXELS_land,&                         
      SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,&                         
      SDS_QCONTROL_land,quality_land,Ret_Quality_cirrus,cloud_num_land,&                         
      SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,&                         
      ETA,ETA_FLAG,AOD,AODF,RHOSFC,&                         
      SDSTAU_corrected_213,SDSTAU_small_land,&                         
      ERR644,Num_Pixels_Used,Qcontrol_special_land,Quality_flag_forJoint,&                         
      REFW354_L,REFW388_L,REFW0P75_L,SDW354_L,SDW388_L,SDW0P75_L,&                         
       Good_pixels_Land,cloudy_pixels) 
!-----------------------------------------------------------------------
!F77
!
!DESCRIPTION:
!                 This subroutine stores all SDSe arrays to HDF file
! 
!INPUT PARAMETERS: all varaiables to be wrriten as output
!
!   LATM          Latitude of data cell
!   LONM          Longitude of data cell
!   IAER          Aerosol type
!   BARLBLUE      Average reflectance for blue channel
!   IDATA         Data cell number from 1 to NUMDATA
!   NUMDATA       Number of retrieval cells across orbit swath
!   TAUABLUE      Blue channel aerosol optical thickness (Continental)
!   BARLRED       Average reflectance for red channel
!   TAUARED       Red channel aerosol optical thickness (Continental)
!   SDBLUE        STD of blue channel reflectances
!   SDRED         STD of red channel reflectances
!   TAUABLUE!     Blue channel aerosol optical thickness (Corrected)
!   TAUARED!      Red channel aerosol optical thickness (Corrected)
!   IBLUE         Number of blue channel observations
!   IRED          Number of red channel observations
!   IERROR        Error flag (0-4)
!   WEIGHT        Relative contribution of smoke/sulfate particles to
!                 dust in the computation of the aerosol optical depth
!   ISULP         Sulfate (1)/no sulfate flag (0)
!   ISMK          Smoke (1)/no smoke flag (0)
!   IPROCE        Aerosol retrieval procedure ID (0-4)
!   SAVERATIO     Aerosol path radiance ration (red to blue channel
!                 for continental model)
!   PIXELS        Number of cells across swath (same as NUMDATA)
!   LINES         Number of cells along swath (same as NUMSCAN)
!   NUMXBOX       Not used
!   NUMYBOX       Not used
 
!OUTPUT PARAMETERS: 13 variables for output to HDF FILE
!
!  SDS1           HDF array of cell latitudes
!  SDS2           HDF array of cell longitudes
!  SDS3           HDF array of spectral reflectances
!  SDS4           HDF array of aerosol optical thicknesses for
!                     continental model
!  SDS5           HDF array of the STD of spectral reflectances
!  SDS6           HDF array of aerosol optical thicknesses for
!                     corrected model
!  SDS7           HDF array of STD for corrected optical thickness
!                     blue channel for continental model)
!  SDS8           HDF array of aerosol path radiance ratio (red to
!                     blue channel for continental model)
!  SDS9           HDF array of relative aerosol optical depth
!                     (smoke/sulfate particles to dust)
!  SDS10          HDF array of the number of blue channel clear pixels
!  SDS11          HDF array of the number of red channel clear pixels
!  SDS12          HDF array of retrieval procedure ID (0-4)
!  SDS13          HDF array of aerosol type (0-3)
!  SDS14          HDF array of Error flag (0-4)
 
      IMPLICIT NONE
      INCLUDE 'mod04.inc'
      INCLUDE 'mod04_land.inc' 
      INCLUDE 'read_Sat_MODIS.inc'
      BYTE QA_Temp
      INTEGER I,QA_Flag_Land(19),Success_Ret_Land,Fail_Ret_Land,cloud_num_land
      INTEGER  quality_land, Ret_Quality_cirrus,IWAV,Num_Pixels_Used
      INTEGER  Good_pixels_Land(Band_land+2),num_wav
      INTEGER IERROR,ICLDBLUE,ICLDRED,IAER,IPROCE,IPR,ISULP,ISMK,IDATA
      INTEGER Qcontrol_special_land,Quality_flag_forJoint(2)
      REAL  TAUABLUE,TAUARED,TAUAREDC,TAUABLUEC,TAUAREDC2,TAUABLUEC2,SDBLUE,SDRED,SD0P86
      REAL  LATM,LONM,SAVERATIO,BARLBLUE,BARLRED,BARL0P86,SDS_SCAT_ANGLE
      REAL  BARRBLUE,BARRRED,TAUABLUENN,TAUAREDNN,SDTAUCBLUE,SDTAUCRED
      REAL  MTHET0,MTHET,MDPHI,WEIGHT,SCAT_ANGLE_LAND,ANGSTROM
      REAL  TAUARAT,TAUAGREENC,TAUAGREENC2
      REAL  FLUXDNBLUE,FLUXDNRED,FLUXUPBLUE,FLUXUPRED,MASSCON
      Real New_cloud_num,Quality_Flag_for_retr
      REAL REFW466_L,REFW550_L 
      REAL REFW644_L,REFW124_L,REFW164_L,REFW866_L,REFW212_L
      REAL SDW466_L,SDW550_L
      REAL SDW644_L,SDW124_L,SDW164_L,SDW866_L,SDW212_L
!NPP channels      
      REAL REFW412_L,REFW443_L,REFW0P75_L,SDW412_L,SDW443_L,SDW0P75_L 
      REAL REFW354_L,REFW388_L,SDW354_L,SDW388_L 
       integer N20,N50, Number_pixels,cloudy_pixels
!
   Real   SDS_MTHET0(NUMCELLS_B),&
                  SDS_MTHET(NUMCELLS_B),&
                  SDS_MPHI(NUMCELLS_B),&
                  SDS_Tau_Land_Ocean(NUMCELLS_B),&
                  SDS_Tau_Land_Ocean_img(NUMCELLS_B),&
                  SDS_Aerosol_Type(NUMCELLS_B),&
                  SDS_SCAT_ANGLE_land(NUMCELLS_B),&
                  SDS_angs_coeff_land(NUMCELLS_B),&
                  SDS_CLDFRC_land(NUMCELLS_B),&
                  SDS_dust_weighting(NUMCELLS_B),&
                  SDS_NUMPIXELS_land(NUMCELLS_B,Band_land+2),&
                  SDSTAU_corrected(NUMCELLS_B,Land_Sol3),&
                  SDS_ref_land(NUMCELLS_B,Band_land+2),&
                  SDS_ref_STD_land(NUMCELLS_B,Band_land+2),&
                  SDSTAU_small_land(NUMCELLS_B,Land_Sol4)
     
       Real SDS_Surface_Reflectance_Land(NUMCELLS_B,Land_Sol3),&
                 SDS_Fitting_Error_Land(NUMCELLS_B),&
                 SDSTAU_corrected_213(NUMCELLS_B)
      REAL  SDSLAT(NUMCELLS_B),SDSLON(NUMCELLS_B),&
          SDS_mass_conc_land(NUMCELLS_B)
      BYTE  SDS_QCONTROL_land(QA_LAND,NUMCELLS_B)
 
! Obsolete (02/2006) Land SDS Arrays
!
      Real  &
                 SDS_est_uncer(NUMCELLS_B,Land_Sol1),&
                 SDS_RefF_land(NUMCELLS_B,Land_Sol2),&
                 SDS_TranF_land(NUMCELLS_B,Land_Sol1)

      REAL ETA
      REAL RHOSFC(NLWAV) 
      REAL AOD(NLWAV),AODF(NLWAV),AODC(NLWAV)
      REAL ERR644
      INTEGER ETA_FLAG
       integer Num_Q_pixel1,Num_Q_pixel2,Num_Q_pixel3,Num_Q_pixel4
      SAVE
      Quality_Flag_for_retr =-9999
      SDS_CLDFRC_land(IDATA)=-9999
      New_cloud_num= -9999
      Quality_flag_forJoint(1) =-9999
      Quality_flag_forJoint(2) = -9999
      DO I=1,19 
       QA_Flag_Land(I)=-9999
      ENDDO 
! compute fraction instead of % for outpt

        IF (cloud_num_land .ge.0) THEN
      New_cloud_num=real(cloud_num_land/Real(Iline*Iline)) 
!       IF (cloudy_pixels.ge.0) THEN
!      New_cloud_num=real(cloudy_pixels/Real(Iline*Iline)) 
      ENDIF
        
      
      IF (AOD(iwave_553) .lt. -0.10) IERROR=5
      IF (AOD(iwave_553) .gt. 5.00)  IERROR=6
      IF (IPROCE.eq.0) IERROR=4
      
!If Errors  
        
      IF (IERROR.GT.0) THEN
        Quality_Flag_for_retr=11
        Fail_Ret_Land=Fail_Ret_Land+1
        QA_Flag_Land(7)=0
        QA_Flag_Land(8)=0
        QA_Flag_Land(9)=0
        QA_Flag_Land(10)=0   
! set it to 12 to indicate No retrivals        
        QA_Flag_Land(11)=Quality_Flag_for_retr
        QA_Flag_Land(12)=IERROR
        QA_Flag_Land(15)=0 
        Qcontrol_special_land=0
        Quality_flag_forJoint(1)=QA_Flag_Land(8)
        SDS_CLDFRC_land(IDATA) = New_cloud_num
! adding 20 to set Non ret. flag so that can tell them apart from Ret. Flag        
        Quality_flag_forJoint(2)= IERROR+20
      
 !       IF ( New_cloud_num .ge.0) THEN
 !         SDS_CLDFRC_land(IDATA)=New_cloud_num
 !      else
 !        SDS_CLDFRC_land(IDATA) = -9999
 !       ENDIF
! Call for Fill values
     
       CALL FILLVALUE_LAND(IDATA,SDS_Tau_Land_Ocean_img,&
                  SDS_Aerosol_Type,SDSTAU_corrected_213,&
                  SDS_SCAT_ANGLE_land,SDS_mass_conc_land,SDS_angs_coeff_land,&
                  SDS_CLDFRC_land,SDS_dust_weighting,&
                  SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,SDS_NUMPIXELS_land,&
                  SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,&
                  SDS_QCONTROL_land,SDSTAU_small_land, &
                  SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,&
                  Qcontrol_special_land)
     
     
        
     
!    IF   retrivals Fill up the SDS's  
   
      ELSE IF (IERROR.EQ.0) THEN 
      
!   Intilize the Quality_Flag_for_retr =0  to indicate good quality  

      Quality_Flag_for_retr= 0
  
! Store SDS arrays for population to HDF file

     
        Success_Ret_Land = Success_Ret_Land+1
       
        QA_Flag_Land(7)=1
        QA_Flag_Land(8)=3 
  
        

! quality_land .eq.0 is if there is one single pixel of water 

        IF (quality_land .eq.0)then
         QA_Flag_Land(8)=0 
         Quality_Flag_for_retr= 2
         ENDIF
         

      
! If Fitting error is greater than =0.25 quality is bad  
    
        IF (ERR644  .gt. 0.25) then
         QA_Flag_Land(8)=0   
          Quality_Flag_for_retr=4
         ENDIF 
        
! Set Mass concentration & Fine optical depth to zero if tau is -ve.         
       IF (AOD(iwave_553) .LT. 0.0) then  
          Qcontrol_special_land=2 
          Quality_Flag_for_retr=5
            MASSCON=0.0 
              DO IWAV = 1, NLWAV 
                AODF(IWAV) = 0.00 
              ENDDO
        ENDIF                 
         
! If -0.10 <= AOD <= -0.05 then set AOD(all waves) = -0.05 and Quality = zero

        IF (AOD(iwave_553) .GE. -0.10 &
          .AND. AOD(iwave_553) .LE. -0.05) THEN 

  
            IF (IPROCE .GT. 1) THEN
              AOD(iwave_553) = -0.05
              AOD(iwave_466) = -0.05 
            ELSE
! Set AOD to -0.05, AODFine to 0.00         
              DO IWAV = 1, NLWAV
                AOD(IWAV) = -0.05   
              ENDDO
            ENDIF
         ENDIF
                                       
        
!         if( ILINE .eq. 10) then
!            Num_Q_pixel1= 2
!            Num_Q_pixel2= 5  
!            Num_Q_pixel3= 8 
!            Num_Q_pixel4= 12  
 
             Number_pixels =  (Iline *iline)
             Num_Q_pixel1= int(Number_pixels *.03) 
             Num_Q_pixel2= int(Number_pixels *.05) 
             Num_Q_pixel3= int(Number_pixels *.075) 
             Num_Q_pixel4= int(Number_pixels *.125)  
!        print*,'number pixels', Num_Q_pixel1,Num_Q_pixel2,Num_Q_pixel3,Num_Q_pixel4,&
!        Number_pixels,Num_Pixels_Used

! set quality for the number of pixels......  
! 3%
!       if( Num_Pixels_Used .ge. 12 .and. Num_Pixels_Used .le.20)then   
        if( Num_Pixels_Used .gt. Num_Q_pixel1 .and. Num_Pixels_Used .le.Num_Q_pixel2)then   
            QA_Flag_Land(8)=0 
            Quality_Flag_for_retr=6
        Endif
!       if( Num_Pixels_Used .ge. 21 .and. Num_Pixels_Used .le.30) then  
! 5%
        if( Num_Pixels_Used .gt. Num_Q_pixel2 .and. Num_Pixels_Used .le.Num_Q_pixel3) then  
             QA_Flag_Land(8)=1
               Quality_Flag_for_retr=7
       ENdif
! 8%
!        if( Num_Pixels_Used .ge.31 .and. Num_Pixels_Used .le.50)then   
         if( Num_Pixels_Used .gt.Num_Q_pixel3 .and. Num_Pixels_Used .le. Num_Q_pixel4)then
            QA_Flag_Land(8)=2
            Quality_Flag_for_retr=8
       ENDIF
! Above 12%
!        if( Num_Pixels_Used .ge.51)QA_Flag_Land(8)=3 
           if( Num_Pixels_Used .gt. Num_Q_pixel4)QA_Flag_Land(8)=3 
        
! Changed the Order of Quality Cirrus Ret_Quality_cirrus =0 if thin cirrus detection  
      
!  ENDIF for  ILINe 
!        ENDIF   
              
         IF (Ret_Quality_cirrus .eq.0)then
         Quality_Flag_for_retr=3
          QA_Flag_Land(8)=0  
          ENDIF

         
!         if( ILINE .eq. 1 .or. ILINE .eq. 3 ) then 
!            QA_Flag_Land(8)=0 
!         endif 
        
        SDS_Aerosol_Type(IDATA)= Real(IAER) 
        SDS_SCAT_ANGLE_land(IDATA)=SCAT_ANGLE_LAND 
        SDSTAU_corrected(IDATA,1)= AOD(iwave_466) 
        SDSTAU_corrected(IDATA,2)= AOD(iwave_553) 
        SDSTAU_corrected(IDATA,3)= AOD(iwave_644)  
        SDSTAU_corrected_213(IDATA)=AOD(iwave_212) 
        SDSTAU_small_land(IDATA,1)= AODF(iwave_466) 
        SDSTAU_small_land(IDATA,2)= AODF(iwave_553) 
        SDSTAU_small_land(IDATA,3)= AODF(iwave_644) 
        SDSTAU_small_land(IDATA,4)= AODF(iwave_212) 
        SDS_CLDFRC_land(IDATA) = New_cloud_num
 !       if( New_cloud_num .GE. 0) then
 !       SDS_CLDFRC_land(IDATA)= New_cloud_num 
 !       Else
 !       SDS_CLDFRC_land(IDATA) = -9999
 !       endif
        SDS_dust_weighting(IDATA)= ETA 
        Do num_wav=1,Band_land+2
          SDS_NUMPIXELS_land(IDATA,num_wav)=Good_pixels_Land(num_wav)
        Enddo
        SDS_Fitting_Error_Land(Idata)=  ERR644 
        SDS_mass_conc_land(IDATA)=MASSCON 
        
! Check  ANGSTROM   
! If optical depth negative then set angs_coeff_land to fill value 
! do not set the flag to zero    
    
       if(Qcontrol_special_land .eq.2)then
         SDS_angs_coeff_land(IDATA)=-9999
       ELSE
        IF (ANGSTROM.LE.5.0 .AND. ANGSTROM.GT.-1.00) THEN
          SDS_angs_coeff_land(IDATA)=ANGSTROM 
        ELSE  
          SDS_angs_coeff_land(IDATA)=-9999
             QA_Flag_Land(8)=0  
             Quality_Flag_for_retr=9
         ENDIF
         ENDIF
         
!If Procedure is 2 set quality to report optical depths at 0.47 & 55 only

        IF (IPROCE.GT.1)then
           Qcontrol_special_land=1
           Quality_Flag_for_retr=1
! If Procedure is 2  Quality is zero
           QA_Flag_Land(8)=0
        ENDIF       
        

         
         
! set quality flag 9 & 10 so that Level 3 uses the qulity flag to average        
          QA_Flag_Land(9)=QA_Flag_Land(7)
          QA_Flag_Land(10)=QA_Flag_Land(8)  
      
!  report eta only when optical depth is < 0.2  

        IF (AOD(iwave_553) .lt. 0.2) THEN 
         Quality_Flag_for_retr=10
          SDS_dust_weighting(IDATA)=-9999 
           
        ENDIF
        
        QA_Flag_Land(11)=Quality_Flag_for_retr
        QA_Flag_Land(12)=IERROR
        
!   Quality flag .....        
         Quality_flag_forJoint(1)=QA_Flag_Land(8)
!   for Retrieval we Quality_Flag_for_retr       
         Quality_flag_forJoint(2)=Quality_Flag_for_retr
          
!    changing order for PACE
        SDS_ref_land(IDATA,1) = REFW354_L 
        SDS_ref_land(IDATA,2) = REFW388_L   
        SDS_ref_land(IDATA,3) = REFW466_L 
        SDS_ref_land(IDATA,4) = REFW550_L 
        SDS_ref_land(IDATA,5) = REFW644_L  
        SDS_ref_land(IDATA,6) = REFW866_L 
        SDS_ref_land(IDATA,7) = REFW124_L  
        SDS_ref_land(IDATA,8) = REFW164_L 
        SDS_ref_land(IDATA,9) = REFW212_L 
        SDS_ref_STD_land(IDATA,1) =SDW354_L  
        SDS_ref_STD_land(IDATA,2) =SDW388_L  
        SDS_ref_STD_land(IDATA,3) =SDW466_L  
        SDS_ref_STD_land(IDATA,4) =SDW550_L 
        SDS_ref_STD_land(IDATA,5) =SDW644_L  
        SDS_ref_STD_land(IDATA,6) =SDW866_L 
        SDS_ref_STD_land(IDATA,7) =SDW124_L 
        SDS_ref_STD_land(IDATA,8) =SDW164_L 
        SDS_ref_STD_land(IDATA,9) =SDW212_L  
      
        SDS_Surface_Reflectance_Land(IDATA,1) =RHOSFC(iwave_466) 
        SDS_Surface_Reflectance_Land(IDATA,2) =RHOSFC(iwave_644) 
        SDS_Surface_Reflectance_Land(IDATA,3) =RHOSFC(iwave_212) 
       
     
       IF ( Qcontrol_special_land .eq.1 ) &
                  CALL FILLVALUE_LAND(IDATA,SDS_Tau_Land_Ocean_img,&
                  SDS_Aerosol_Type,SDSTAU_corrected_213,&
                  SDS_SCAT_ANGLE_land,SDS_mass_conc_land,SDS_angs_coeff_land,&
                  SDS_CLDFRC_land,SDS_dust_weighting,&
                  SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,SDS_NUMPIXELS_land,&
                  SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,&
                  SDS_QCONTROL_land,SDSTAU_small_land,&
                  SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,&
                  Qcontrol_special_land) 
           
! ENDIF FOR ERROR =0 or error > 0
      ENDIF
       
      
 
!Store QA flags into Quality_Assurance_Land array according to the order
! of bits in MODIS atmosphere QA plan
 
      QA_Temp=0
      CALL BYTE_SET(QA_Flag_Land(7),0,QA_Temp)
      CALL BYTE_SET(QA_Flag_Land(8),1,QA_Temp)
       CALL BYTE_SET(QA_Flag_Land(9),4,QA_Temp)
       CALL BYTE_SET(QA_Flag_Land(10),5,QA_Temp)

      SDS_QCONTROL_land(1,IDATA)=QA_Temp

      QA_Temp=0
      QA_Flag_Land(13)=0
      QA_Flag_Land(14)=0
      CALL BYTE_SET(QA_Flag_Land(11),0,QA_Temp)
      CALL BYTE_SET(QA_Flag_Land(12),4,QA_Temp)
!      CALL BYTE_SET(QA_Flag_Land(13),6,QA_Temp)
!      CALL BYTE_SET(QA_Flag_Land(14),7,QA_Temp)

      SDS_QCONTROL_land(2,IDATA)=QA_Temp

      QA_Temp=0
      QA_Flag_Land(17)=3
      QA_Flag_Land(18)=3

      CALL BYTE_SET(QA_Flag_Land(15),0,QA_Temp)
      CALL BYTE_SET(QA_Flag_Land(16),2,QA_Temp)
      CALL BYTE_SET(QA_Flag_Land(17),4,QA_Temp)
      CALL BYTE_SET(QA_Flag_Land(18),6,QA_Temp)

      SDS_QCONTROL_land(3,IDATA)=QA_Temp

      QA_Temp=0
      QA_Flag_Land(19)=1
      CALL BYTE_SET(QA_Flag_Land(19),0,QA_Temp)
 
      SDS_QCONTROL_land(4,IDATA)=QA_Temp

      SDS_QCONTROL_land(5,IDATA)=0
      
        
       
 
!   Re-initialized working variables
 
 

      DO I=1,19 
       QA_Flag_Land(I)=-9999
      ENDDO 
      IAER=-9999
      ICLDBLUE=-9999
      ICLDRED=-9999
      IPROCE=-9999 
      IERROR=-9999
      REFW466_L=-9999
      REFW550_L=-9999
      REFW644_L=-9999
      REFW124_L=-9999
      REFW164_L=-9999
      REFW866_L=-9999
      REFW212_L=-9999
      REFW412_L=-9999
      REFW412_L=-9999 
      SDW466_L=-9999
      SDW550_L =-9999
      SDW644_L=-9999
      SDW124_L=-9999
      SDW164_L=-9999
      SDW866_L=-9999
      SDW212_L=-9999 
      SDW354_L=-9999
      SDW388_L=-9999 
!      Re-Initialize for next run
      ETA = -9999
      ERR644 = -9999
      MASSCON = -9999
      ANGSTROM = -9999
      DO IWAV = 1, NLWAV
         AOD(IWAV) = -9999
         AODF(IWAV) = -9999
         AODC(IWAV) = -9999
         RHOSFC(IWAV) = -9999
      ENDDO
      
      
      
!
100   FORMAT(i4,10(10f10.4,/))
     
      RETURN
      END

!***********************************************************************
      SUBROUTINE COMPUTE_SCATTANGLE_LAND(MTHET0,MTHET,MDPHI,IDATA,&
      SCAT_ANGLE_LAND)
!----------------------------------------------------------------------
!!F90             
!
!!DES!RIPTION:
!              This subroutine !omputes s!attering angle from MODIS
!              geometry.
!
!!INPUT PARAMETERS:
!
!       MEHET0            Solar zenith angle
!       MTHET             Satellite Viewing angle
!       MDPHI             Relative azimuth angle
!       IDATA             Index of 10x10 km box
!
!!OUTPUT PARAMETERS:
!
!       ScAT_ANGLE_LAND   Scattering angle
!
 
      IMPLICIT NONE
      INCLUDE 'mod04.inc'
      REAL SCAT_ANGLE_LAND
      REAL MTHET0,MTHET,MDPHI
      INTEGER IDATA
      SAVE

!      DTR=ACOS(-1.)/180.
!      RTD=180./ACOS(-1.)
      SCAT_ANGLE_LAND=0.0

!
! Scattering angle of MODIS geometry
 

       SCAT_ANGLE_LAND= -COS(MTHET0*DTR)*COS(MTHET*DTR)&
                        +SIN(MTHET0*DTR)*SIN(MTHET*DTR)&
                        *COS(MDPHI*DTR)

       SCAT_ANGLE_LAND= ACOS(SCAT_ANGLE_LAND)*RTD

       RETURN
       END

 

 
!*********************************************************************
       SUBROUTINE PROCESS_Land_Inver(IPROCE,&
          MTHET0,MTHET,MDPHI,SCAT_ANGLE_LAND,IMONTH,LATM,LONM,MHGHT,&
          REFW466,REFW553,REFW644,REFW866,REFW123,REFW163,REFW212,&
          SDW466,SDW553,SDW644,SDW866,SDW123,SDW163,SDW212,&
          INT_NL0,Fd_NL0,T_NL0,OPTH_NL0,SBAR_NL0,AEROSOL,REF_RAY_NL,&
         ETA,ETA_FLAG,AOD,ERR644,RHOSFC,AVE_COUNT,FTABLE_NL,&
          MASSCOEF_NL0,MASSCON,&
          EXTNORM_NL0,AODF,AODC,ANGSTROM,iscan,idata,Urban_per)

!-----------------------------------------------------------------------
!!F90   
!!DESCRIPTION:
!The subroutine derives non-cloudy aerosol optical thickness from MODIS
!measured radiances over land, using the 0.47 (blue), 0.66 (red) 
!and 2.13 (IR) micron channels.
!Input data are cloud screened and averaged reflectance at at 10 x 10 km. 
!Inversion:
!A) Select fine mode model based on geography (dynamic aerosol models)
!B) Mix fine mode and dust mode with selected discrete mixing ratios
!   best matches the observed spectral reflectance. 
!The program assumes that the input data are in reflectance units
!of PI*L/(F0*cos(theta0)) for all channels.
!It also assumes that the surface reflectance in the two visible channels
!are functions of the surface reflectance in the IR. 
!Nominally, red is one half of IR, and blue is one-half of red, respectively. 
!!INPUT PARAMETERS:
!
!MTHET0                         solar zenith angle (degree)
!MTHET                          satellite viewangle angle (degree)
!MDPHI                          relative azimuth in angle (degree)
!MHGHT                          Topographic altitude (km)
!IMONTH                         Calendar month
!LATM                           latitude of 10x10 km box
!LONM                           longitude of 10x10 km box
!REFW466,REFW644,REFW866,REFW212   averaged cloud screened reflectance
!SDW470,SDW659,SDW866,SDW212       stddev cloud screened reflectance
!INT_NL0			      reflectance from lookup table
!Fd_NL0                           flux down from lookup table
!T_NL0                           transmission from lookup table
!OPTH_NL0                        optical thickness from lookup table
!SBAR_NL0                        sbar from lookup table
!MASSCOEF_NL0                    Mass Concentration from lookup table
!EXTNORM_NL0                     Normalized Extinction Coeficients from LUT
!
!!OUTPUT PARAMETERS:
!
!ETA                Fine mode (non dust) fraction
!AOD                AOD at 4 wavelengths 
!AODF/!          Fine and Coarse mode AOD at 4 wavelengths 
!RHOSF!          Surface reflectance at 4 wavelengths
!AVE_COUNT          Number of observations with error < 3%
!FTABLE_NL          Which Fine model
!MASSCON            Mass Concentration
!ETA_FLAG           Whether ETA between 0.0 and 1.0
!ANGSTROM           Angstrom Exponenent 0.47/0.66
!
!Aerosol types (FTABLE)
!   1              Continental
!   2              Generic : SSA ~ 0.9
!   3              Smoke   : SSA ~ 0.85
!   4              Urban   : SSA ~ 0.95
!   5              Dust
!
! INTERNALS:
!
!Subroutines:
!   RNLOOKUP       -  READS THE LOOK-UP TABLES
!   SELECT_NLOOKUP -  SELECT LOOK-UP TABLE
!   INTANGLE_NL   -  Interpolates LUT parameters to measured geometry
!   TOA_SIMULATE  -  Computes simulated TOA reflectance for each TAU index
!   RETRIEVAL1    -  Derives AOT, ETA and 2.1 Surface Reflectance for
!                    IPROCE=1 (0.01 < rho2p1 < 0.25)
!   RETRIEVAL2    -  Derives AOT for Continental model, IPROCE=2(0.25 < rho2p1)
!   SHLSRT_NL     -  Sorts an array into ascending order
!   INTERP_EXTRAP -  Linear interpolation/extrapolation
!   PERCENTILE    -  Computes average value/standard deviation of selected
!                    pixels given percentiles
!   ERROR         -  SETS THE ERROR CODES
!   OUTPUT        -  WRITES THE OUTPUT TO HDF FILE
!   AEROSOL_MAP   -  Reads map of aerosol model as function of season and place
!   INT_ELEV      -  Interpolates lookup table at different wavelengths to 
!                    simulate elevation 
 
      IMPLICIT NONE
      INCLUDE 'mod04.inc'
      INCLUDE 'mod04_land.inc'

      INTEGER NUMSQDIM,NUMXBOX,NUMYBOX
      PARAMETER (NUMSQDIM=400)
      INTEGER RTN,START_DATA,END_DATA
      INTEGER IDATA,ISCAN

      INTEGER IFINISH,IMONTH,IPR,IPROCE
      REAL SCAT_ANGLE_LAND

      REAL MTHET0,MTHET,MPHI0, MPHI, MDPHI,MHGHT
      REAL LATM,LONM,LAT_CENTER,LON_CENTER
!      REAL WATER_VAPOR(NUMSQDIM),OZONE_COL(NUMSQDIM)
      CHARACTER * 5 IOFLAG
      character ( len =10):: Sat_Flag

!********************************************************************
!   Parameters used from look-up tables
! ********************************************************************
!   wav_nl=wavelengths(4)(.466,.533,.644,2.13 micrometers)       4
!   nlthet0=number of theta0(solar zenith angle)              9
!   thet0_nl=(0,12,24,36,48,54,60,66,72   degrees)
!   nlthe =number of the(observation zenith angle)           13
!   the_nl=(0,6,.....72 degrees every 6 degrees)
!   nlphi =number of phi(observation azimuth angle)          16
!   phi_nl=(0,12,24,36,48,60,72,84,96,
!        108,120,132,144,156,168,180)
!   nltau =number of opth(optical thickness)                 7 
!   opth_nl=(0.0,0.25,0.50,0.75,1.0,2.0,3.0,5.0)
!   lhght = number of hght(observation heights)              1
!   hght=(80.0 km)
!**********************************************************************

  
!  For Determining aerosol type from aerosol map
      INTEGER nlon, nlat
      PARAMETER (nlon = 360, nlat = 180)
      INTEGER AEROSOL(nlon,nlat)
      INTEGER ilat, ilon
      INTEGER FTABLE,FTABLE_NL
      REAL LONM1,Urban_per
      REAL dlat, dlon
      PARAMETER (dlat = 0.5, dlon = 0.5)
      CHARACTER*5 AEROSOL_TYPEs(5), AEROSOL_TYPE




!     Data directly read from the entire set of LUTs

      REAL SBAR_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL OPTH_NL0(NLTAU,NLWAV,NLTABLE)
      REAL INT_NL0(NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL Fd_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL T_NL0(NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL MASSCOEF_NL0(NLTAU,NLWAV,NLTABLE)
      REAL EXTNORM_NL0(NLTAU,NLWAV,NLTABLE)

      integer ij,ik,set_counter_land,qcontrol_special 

!    Data for specific LUTs (specified by aerosol model)

      REAL INT_NL1(NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL Fd_NL1(NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL SBAR_NL1(NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL T_NL1(NLTHE,NLTHET0,NLTAU,NLWAV,NLSIZE)

!     Interpolated lookup table data

      REAL INT_NL(NLTAU,NLWAV,NLSIZE)
      REAL Fd_NL(NLTAU,NLWAV,NLSIZE)
      REAL T_NL(NLTAU,NLWAV,NLSIZE)
      REAL FdT_NL(NLTAU,NLWAV,NLSIZE)
      REAL SBAR_NL(NLTAU,NLWAV,NLSIZE)
      REAL OPTH_NL(NLTAU,NLWAV,NLSIZE)
      REAL EQWAV_NL(NLWAV)
      REAL MASSCOEF_NL(NLTAU,NLWAV,NLSIZE)
      REAL EXTNORM_NL(NLTAU,NLWAV,NLSIZE)

!     Rayleigh data
      REAL REF_RAY_NL(NLWAV)


!    Simulated Reflectancea
      REAL  RHO_S212(NLTAU,NLSIZE), RHOSTAR(NLTAU,NLWAV,NLSIZE)
      REAL yint_466,slope_466,yint_644,slope_644

!    Results

      REAL ERR644
      REAL TAULAND55,RHOSFC212
      REAL RHOSFC(NLWAV)
      REAL ETA
      INTEGER ETA_FLAG
      REAL AOD(NLWAV),AODF(NLWAV),AODC(NLWAV)
      INTEGER AVE_COUNT
      REAL MASSCON
      REAL ANGSTROM


!     Dummy Indices
     
      INTEGER NSOLUTION,ISIZE,IETA,JETA,IWAV
      INTEGER IPHI,ITHE,ITAU,ITHET0,J
      INTEGER NUMDATA

!     Geometry Indices

      INTEGER KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1
      REAL REFW466,REFW553,REFW644,REFW866,REFW123,REFW163,REFW212
      REAL SDW466,SDW553,SDW644,SDW866,SDW123,SDW163,SDW212

!     INPUT Surface reflectance and TAU
!     Match corresponding aerosol model with aerosol map index

      SAVE

      DATA AEROSOL_TYPEs/'CONTN','BACKG','SMOKE','URBAN','DUST'/
     

 

!     Dummy set some stuff

! Pre-processing

      REFW466L = REFW466
      REFW644L = REFW644
      REFW212L = REFW212
      REFW124L = REFW123
      REFW866L = REFW866
      REFW164L = REFW163
      REFW553L = REFW553

!   Following If statement checks for missing data in a box 10*10
!  in any one of wavelengths.
!  If there is missing data in any one wavelength
!  arrays are filled with _fillValues.

       IF (REFW466L.LT.99999 .AND. REFW644L .LT. 99999 .AND.&
          REFW866L.LT.99999 .AND. REFW212L .LT. 99999) THEN


! Determine expected fine model (from date and location)

      ilon = 181 + NINT(LONM + dlon)
      ilat = 91 - NINT(LATM + dlat)
      if(ilon .gt. 360)ilon = 360 
      if(ilat .gt. 180)ilat = 180
      FTABLE = AEROSOL(ilon,ilat) + 2
      AEROSOL_TYPE = AEROSOL_TYPEs(FTABLE)

!  Select LUT corresponding to season and location
!  Obtain fine mode LUT (choice of tables #2, #3 or #4)
! Coarse mode LUT = #5
!  Continental mode LUT = #1

      IF (IPROCE .EQ. 1) THEN
        MTABLE_NL(1) = FTABLE
        MTABLE_NL(2) = DTABLE_NL
      ELSE IF (IPROCE .EQ. 2) THEN
        MTABLE_NL(1) = 1
        MTABLE_NL(2) = 1
      ENDIF
      FTABLE_NL = MTABLE_NL(1)


      DO iSIZE = 1, NLSIZE

         CALL SELECT_NLOOKUP(INT_NL0,Fd_NL0,T_NL0,OPTH_NL0,&
            SBAR_NL0,MASSCOEF_NL0,EXTNORM_NL0,&
            ISIZE,INT_NL1,Fd_NL1,T_NL1,SBAR_NL1,&
            OPTH_NL,MASSCOEF_NL,EXTNORM_NL,MTABLE_NL(ISIZE))

      ENDDO
           
         
!Interpolate lookup tables to measured angles. 

	CALL  SET_index_inter_NL(MTHET0,MTHET,MDPHI,&
            KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1)
     
           
 
! Interpolate lookup table to measured geometry. 
     
        CALL INTANGLE_NL(MTHET0,MTHET,MDPHI,&
      	  INT_NL1,Fd_NL1,T_NL1,SBAR_NL1,&
          INT_NL,Fd_NL,T_NL,FdT_NL,SBAR_NL,&
          KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1)
            
! Interpolate lookup table to target elevation by interpolating wavelengths


          
        CALL INT_ELEV(EQWAV_NL,INT_NL,Fd_NL,T_NL,FdT_NL,OPTH_NL,&
              SBAR_NL,REF_RAY_NL,MHGHT) 
!  Simulate TOA reflectance
!
             

           CALL TOA_SIMULATE(&
              MTHET0,SCAT_ANGLE_LAND,&
              INT_NL,FdT_NL,OPTH_NL,SBAR_NL,RHO_S212,RHOSTAR,&
              yint_466,slope_466,yint_644,slope_644,Urban_per)

 	 
! DO RETRIEVAL!!!!!!!!!!!
!     Retrieval#1 for IPROCE=1, Retrieval#2 for IPROCE=2	 
!

           ETA_FLAG = 0
           IF (IPROCE .EQ. 1) THEN	
             CALL RETRIEVAL1 (&
             yint_466,slope_466,yint_644,slope_644,&
             OPTH_NL,MASSCOEF_NL,EXTNORM_NL,&
             RHO_S212,RHOSTAR,&
             AVE_COUNT,ETA,ETA_FLAG,AOD,AODF,AODC,&
             RHOSFC212,ERR644,RHOSFC,MASSCON)
           ELSE IF (IPROCE .EQ. 2) THEN
            CALL RETRIEVAL2 (&
             yint_466,slope_466,yint_644,slope_644,&
             OPTH_NL,MASSCOEF_NL,RHO_S212,RHOSTAR,&
             AVE_COUNT,ETA,AOD,&
            RHOSFC212,ERR644,RHOSFC,MASSCON)
           ENDIF
           
              
!        Calculate Angstrom Exponent

           IF (AOD(iwave_644) .GT. 0) THEN
             ANGSTROM = -1.0 * ALOG(AOD(iwave_466)/AOD(iwave_644))/&
                 ALOG(0.466/0.644)
           ENDIF


          ENDIF
       
      RETURN		
 
      END

!*****************************************************************
      SUBROUTINE SELECT_NLOOKUP(&
            INT_NL0,Fd_NL0,T_NL0,OPTH_NL0,SBAR_NL0,MASSCOEF_NL0,&
            EXTNORM_NL0,&
            ISIZE,INT_NL1,Fd_NL1,T_NL1,SBAR_NL1,&
           OPTH_NL,MASSCOEF_NL,EXTNORM_NL,ITAB)
!---------------------------------------------------------------------
!
                  
!
!DESCRIPTION:
!            This subroutine selects look-up table (using ITABLE)
!	       Used to select both the "fine" and the "coarse/dust" land model
!
!!INPUT PARAMETERS:
!
!        INT_NL0      Reflectance
!        Fd_NL0       Flux Down
!        T_NL0        Transmssion factor
!        OPTH_NL0     Optical depth
!        SBAR_NL0     Sbar
!        MASSCOEF_NL0  MASSCON coeficient
!        EXTNORM_NL0  Normalized extinction coeficient
!        ITABLE     Index for desired lookup table
!        ISIZE      Index indicating "fine" or "coarse"
!
!OUTPUT PARAMETERS:
!
!        INT_NL1        Reflectance
!        Fd_NL1       Flux Down
!        T_NL1        Transmssion factor
!        SBAR_NL1       Sbar
!        OPTH_NL       Optical depth
!        MASSCOEF_NL     Mass Conencentration coeficient
!        EXTNORM_NL     Normalized extinction coeficient
 
      IMPLICIT NONE
      INCLUDE 'mod04_land.inc'

      INTEGER IWAV,ITAB,ITAU,ITHET0,ITHE,IPHI
      INTEGER ISIZE

      REAL OPTH_NL0(NLTAU,NLWAV,NLTABLE)
      REAL SBAR_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL INT_NL0(NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL Fd_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL T_NL0(NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
      REAL MASSCOEF_NL0(NLTAU,NLWAV,NLTABLE)
      REAL EXTNORM_NL0(NLTAU,NLWAV,NLTABLE)

      REAL INT_NL1(NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL Fd_NL1(NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL SBAR_NL1(NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL T_NL1(NLTHE,NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL OPTH_NL(NLTAU,NLWAV,NLSIZE)
      REAL MASSCOEF_NL(NLTAU,NLWAV,NLSIZE)
      REAL EXTNORM_NL(NLTAU,NLWAV,NLSIZE)
      SAVE

 
! Select look-up table variable for given ITABLE
 

      DO 10 IWAV=1,NLWAV
        DO 20 ITAU=1,NLTAU
          OPTH_NL(ITAU,IWAV,ISIZE) = OPTH_NL0(ITAU,IWAV,ITAB)
          MASSCOEF_NL(ITAU,IWAV,ISIZE) = MASSCOEF_NL0(ITAU,IWAV,ITAB)
          EXTNORM_NL(ITAU,IWAV,ISIZE) = EXTNORM_NL0(ITAU,IWAV,ITAB)
          DO 30 ITHET0=1,NLTHET0
            Fd_NL1(ITHET0,ITAU,IWAV,ISIZE) =&
               Fd_NL0(ITHET0,ITAU,IWAV,ITAB) 
            SBAR_NL1(ITHET0,ITAU,IWAV,ISIZE) =&
               SBAR_NL0(ITHET0,ITAU,IWAV,ITAB) 
            DO 40 ITHE=1,NLTHE
               T_NL1(ITHE,ITHET0,ITAU,IWAV,ISIZE) = &
                 T_NL0(ITHE,ITHET0,ITAU,IWAV,ITAB) 
              DO 50 IPHI=1,NLPHI
                  INT_NL1(IPHI,ITHE,ITHET0,ITAU,IWAV,ISIZE) = &
                    INT_NL0(IPHI,ITHE,ITHET0,ITAU,IWAV,ITAB) 
50            CONTINUE
40          CONTINUE
30        CONTINUE
20      CONTINUE
10     CONTINUE

      RETURN
      END


!*****************************************************************
       SUBROUTINE SET_index_inter_NL(MTHET0,MTHET,MDPHI,&
       KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1)

!-----------------------------------------------------------------------
!F90
                  
!
!!DESCRIPTION: This subroutine sets the index for ending and starting
!            positions from the lookup table (based on input geometry)
!
!!INPUT PARAMETERS:
!              MTHET0      Measured solar Zenith angle
!              MTHET      Measured view  Zenith angle
!             MDPHI      Measured Azimuthal Angle
!
!             THET0      array of  solar Zenith angle in look_up table
!              THE        array of  satllite Zenith angle in look_up table
!              PHI        array of azimuth angle in look_up table
!
!!OUTPUT PARAMETERS:
!               KSZAM1     Starting Index for solar zenith angle
!               KSZAP1     Ending Index for solar zenith angle
!               KTHEM1     Starting Index for view angle
!               KTHEP1     Ending   Index for  view angle
!               KPHIM1     Starting Index for  azimuth angle
!               KPHIP1     Ending   Index for  azimuth angle
!
 

       IMPLICIT NONE
       SAVE

       INCLUDE 'mod04_land.inc'

       REAL MTHET0,MTHET,MDPHI,DEL
       INTEGER KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1
       INTEGER I

!     Initialize
 
       KSZAM1=0
       KSZAP1=0
       KTHEM1=0
       KTHEP1=0
       KPHIM1=0
       KPHIP1=0

!      Find THET0 indices
       DO I = 1,NLTHET0-1
          IF ((MTHET0 .GE. THET0_NL(I)) .AND. &
            (MTHET0 .LE. THET0_NL(I+1))) THEN
	    KSZAM1=I
	    KSZAP1=I+1
          ENDIF
       ENDDO

!      Find THE indices
       DO I = 1,NLTHE-1
          IF ((MTHET .GE. THE_NL(I)) .AND. &
             (MTHET .LE. THE_NL(I+1))) THEN
	    KTHEM1=I
	    KTHEP1=I+1
          ENDIF
       ENDDO

!      Find PHI indices
       DO I = 1,NLPHI-1
          IF ((MDPHI .GE. PHI_NL(I)) .AND. &
            (MDPHI .LE. PHI_NL(I+1))) THEN
	    KPHIM1=I
	    KPHIP1=I+1
          ENDIF
       ENDDO


       RETURN
       END


			
!*********************************************************************
      SUBROUTINE   INTANGLE_NL(MTHET0,MTHET,MDPHI,&
      	  INT_NL1,Fd_NL1,T_NL1,SBAR_NL1,&
          INT_NL,Fd_NL,T_NL,FdT_NL,SBAR_NL,&
          KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1)
!----------------------------------------------------------------------
!!F90
                  
!
!DESCRIPTION:  Subroutine INTANGLE_NL interpolates the lookup
!               reflectances to the measured geometry.
!
!INPUT PARAMETERS:
!	  INT_NL1		  radiance
!	  Fd_NL1		  flux down
!	  T_NL1		  transmission
!	  SBAR_NL1		  sbar
!        KSZAM1      Starting Index for solar zenith angle
!        KSZAP1      Ending Index for solar zenith angle
!        KTHEM1      Starting Index for view angle
!         KTHEP1      Ending   Index for  view angle
!         KPHIM1      Starting Index for  azimuth angle
!        KPHIP1      Ending   Index for  azimuth angle
!        MTHET0        Solar zenith angle.
!         MTHET        View angle.
!        MDPHI        Azimuth angle.
!
!        THET0      array of  solar Zenith angle in look_up table
!         THE        array of  satllite Zenith angle in look_up table
!         PHI        array of azimuth angle in look_up table
!
!!OUTPUT PARAMETERS 
!	  INT_NL		  interpolated radiance
!	  Fd_NL		  interpolated flux down
!	  T_NL		  interpolated transmission
!	  SBAR_NL         interpolated  sbar
!	  FdT_NL	  two way transmission transmission
!         REF_RAY        Reflectance for rayleigh only
!
 !----------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04_land.inc'

!     Inputs
      REAL  MTHET0,MTHET,MDPHI
      INTEGER KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1
      
      REAL INT_NL1(NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL Fd_NL1(NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL SBAR_NL1(NLTHET0,NLTAU,NLWAV,NLSIZE)
      REAL T_NL1(NLTHE,NLTHET0,NLTAU,NLWAV,NLSIZE)

!     Outputs
      REAL INT_NL(NLTAU,NLWAV,NLSIZE),Fd_NL(NLTAU,NLWAV,NLSIZE)
      REAL T_NL(NLTAU,NLWAV,NLSIZE),FdT_NL(NLTAU,NLWAV,NLSIZE)
      REAL SBAR_NL(NLTAU,NLWAV,NLSIZE)

!      CHARACTER*132 LINE
      CHARACTER*45  LINE1
      INTEGER ICASE,IJ,IPHI,ISIZE,ITAU,ITH,ITH0,IWAV,IETA,LOPT
      REAL  X(2),Y(2),XX1(2),YY1(2),XX2(2),YY2(2),Y1
      REAL  WW2(2),WW1(2),W1
      REAL  VV2(2),V1
      REAL  UU2(2),U1
      INTEGER LL,MM,NN

! LOOP IS AROUND THET0,NEW FILE FOR EACH THETA0
!

        Y1=0.0
        DO IJ = 1,2
           X(IJ)=0.0
           Y(IJ)=0.0
           XX1(IJ)=0.0
           XX2(IJ)=0.0
           YY1(IJ)=0.0
           YY2(IJ)=0.0
           WW1(IJ)=0.0
           WW2(IJ)=0.0
           VV2(IJ)=0.0
           UU2(IJ)=0.0
        ENDDO
        DO 5 ISIZE= 1,NLSIZE
          DO 15 IWAV=1,NLWAV
            DO 25  ITAU  = 1,NLTAU
              LL=0
              y1 = 0.0
              DO  40 ITH0  = KSZAM1,KSZAP1
                MM=0
                DO  50  ITH  = KTHEM1,KTHEP1
                  NN=0
                  DO 60  IPHI  = KPHIM1,KPHIP1
                    NN=NN+1
                    X(NN)=PHI_NL(IPHI)
                    Y(NN)=INT_NL1(IPHI,ITH,ITH0,ITAU,IWAV,ISIZE)
60                CONTINUE

                  CALL INTERP_EXTRAP(NN,MDPHI,X,Y,Y1,1)
                  MM=MM+1
                  XX1(MM)=THE_NL(ITH)
                  YY1(MM)=Y1
                  WW1(MM)=T_NL1(ITH,ITH0,ITAU,IWAV,ISIZE)

50              CONTINUE
                y1=0.0
                w1=0.0 
                CALL INTERP_EXTRAP(MM,MTHET,XX1,YY1,Y1,1)
                CALL INTERP_EXTRAP(MM,MTHET,XX1,WW1,W1,1)
                LL=LL+1
                XX2(LL)=THET0_NL(ITH0)
                YY2(LL)=Y1
                WW2(LL)=W1
                VV2(LL)=Fd_NL1(ITH0,ITAU,IWAV,ISIZE)
                UU2(LL)=SBAR_NL1(ITH0,ITAU,IWAV,ISIZE)

40            CONTINUE
              y1=0.0
              w1=0.0
              v1=0.0
              u1=0.0
              CALL INTERP_EXTRAP(LL,MTHET0,XX2,YY2,Y1,1)
              CALL INTERP_EXTRAP(LL,MTHET0,XX2,WW2,W1,1)
              CALL INTERP_EXTRAP(LL,MTHET0,XX2,VV2,V1,1)
              CALL INTERP_EXTRAP(LL,MTHET0,XX2,UU2,U1,1)

              INT_NL(ITAU,IWAV,ISIZE) = Y1
              T_NL(ITAU,IWAV,ISIZE) = W1
              Fd_NL(ITAU,IWAV,ISIZE) = V1
              SBAR_NL(ITAU,IWAV,ISIZE) = U1
              FdT_NL(ITAU,IWAV,ISIZE) = W1*V1

25          CONTINUE
15        CONTINUE
 5      CONTINUE


        RETURN
        END



!*********************************************************************
      SUBROUTINE  TOA_SIMULATE(&
          MTHET0,SCAT_ANGLE_LAND,&
          INT_NL,FdT_NL,OPTH_NL,SBAR_NL,RHO_S212,RHOSTAR,&
          yint_466,slope_466,yint_644,slope_644,Urban_per)

!-----------------------------------------------------------------------
!!F90
                  
!
!
!!DESCRIPTION:
!This subroutine simulates TOA reflectance
!
!
!! INPUTS
!!OUTPUTS (INTO INCLUDE FILE)
!	       RHOSTAR	  array of TOA reflectance
!-----------------------------------------------------------------------

        IMPLICIT NONE
        Save
        INCLUDE 'mod04_land.inc'
       
!     Inputs
      character ( len=10):: Sat_Flag
      REAL MTHET0,SCAT_ANGLE_LAND
      REAL INT_NL(NLTAU,NLWAV,NLSIZE),FdT_NL(NLTAU,NLWAV,NLSIZE),&
          SBAR_NL(NLTAU,NLWAV,NLSIZE),OPTH_NL(NLTAU,NLWAV,NLSIZE)
      REAL REF_RAY_NL(NLWAV)

!     for estimating surface reflectance
      REAL yint_466, slope_466
      REAL yint_644, slope_644
      REAL RHO_S466, RHO_S644
      REAL MVI,Urban_per

!     Outputs

      REAL  RHO_S212(NLTAU,NLSIZE), RHOSTAR(NLTAU,NLWAV,NLSIZE)

!     Dummy

        INTEGER IETA,ITAU,IWAV,ISIZE,I


! 	At each tau index, calculate surface reflectance at 2.1 and
!           apparent reflectance at all wavelengths
!      
       DO ISIZE = 1, NLSIZE
          DO ITAU = 1,NLTAU
	    RHO_S212(ITAU,ISIZE) = 0.0
	    RHO_S212(ITAU,ISIZE) = &
               (INT_NL(ITAU,iwave_212,ISIZE)-REFW212L) / &
               ((SBAR_NL(ITAU,iwave_212,ISIZE) * &
                (INT_NL(ITAU,iwave_212,ISIZE)-REFW212L)) - &
                (FdT_NL(ITAU,iwave_212,ISIZE))) 

            IF (RHO_S212(ITAU,ISIZE) .GT. 1.0) THEN
              RHO_S212(ITAU,ISIZE) = -9.999
            ENDIF 


 
 
! VIIRS surface Ratio   
!for simple rations, overwrite (or comment out all lines above)
!  "Current VIIRS " slope_466 = ratio of M3(.48) / M5(0.67)
!  slope_644 =ratio of M5(0.67) / M11 (2.25)

            slope_644 = 0.559
            yint_644 = 0.0
            slope_466 = 0.645
            yint_466 = 0.0                 
            
     
         
            RHO_S644 = slope_644*RHO_S212(ITAU,ISIZE) + yint_644
            RHO_S466 = slope_466*RHO_S644 + yint_466
 

!         Compute model differentiated apparent reflectance
           RHOSTAR(ITAU,iwave_466,ISIZE) = 0.0
           RHOSTAR(ITAU,iwave_466,ISIZE) = &
             INT_NL(ITAU,iwave_466,ISIZE) +&
             ((FdT_NL(ITAU,iwave_466,ISIZE) * RHO_S466) / &
             (1 - (SBAR_NL(ITAU,iwave_466,ISIZE) * RHO_S466)))

           RHOSTAR(ITAU,iwave_553,ISIZE) = 0.0 
           RHOSTAR(ITAU,iwave_644,ISIZE) = 0.0
           RHOSTAR(ITAU,iwave_644,ISIZE) =  &
             INT_NL(ITAU,iwave_644,ISIZE) + &
            ((FdT_NL(ITAU,iwave_644,ISIZE) * RHO_S644) / &
             (1 - (SBAR_NL(ITAU,iwave_644,ISIZE) * RHO_S644)))

           RHOSTAR(ITAU,iwave_212,ISIZE) = 0.0
           RHOSTAR(ITAU,iwave_212,ISIZE) = &
                INT_NL(ITAU,iwave_212,ISIZE) + &
                 ((FdT_NL(ITAU,iwave_212,ISIZE) * &
                 (RHO_S212(ITAU,ISIZE))) / &
                 (1 - (SBAR_NL(ITAU,iwave_212,ISIZE) * &
                 (RHO_S212(ITAU,ISIZE)))))

          ENDDO

        ENDDO
        RETURN

      END


!*********************************************************************

      SUBROUTINE RETRIEVAL1(&
        yint_466,slope_466,yint_644,slope_644,&
         OPTH_NL,MASSCOEF_NL,EXTNORM_NL,RHO_S212,RHOSTAR,&
         AVE_COUNT,ETA,ETA_FLAG,AOD,AODF,AODC,&
         RHOSFC212,ERR644,RHOSFC,MASSCON)
!-----------------------------------------------------------------------
!!F90 
!!DESCRIPTION:
!               This subroutine retrieves optical thickness, surface
!	 	reflectance and error parameters by comparing
!	        MODIS observations and LUT data
!
!
!! INPUTS
!              REFW466L,REFW644L,REFW212L   reflectance
!              OPTH_NL      LUT based optical thickness
!              MASSCOEF_NL  LUT based mass concentration coeficients
!              EXTNORM_NL   LUT based normalized extinction coeficients
!              RHO_S212   simulated surface reflectance
!              RHOSTAR    simulated TOA reflectance
!              yint644,yint466,slope466,slope644  
!                 characteristics of VIS/IR surface reflectance relationship
!! OUTPUTS
!	       ERR644	  retrieval fitting errors 
!	       AOD        retrieved Optical depths 
!	       AODF/C     Fine and coarse mode optical depth 
!	       RHOSFC     retrieved surface reflectance 
!	       ETA        retrieved fine mode weighting
!	       MASSCON    retrieved mass concentration coeficient
!              ETA_FLAG   0 if eta within 0.0 - 1.0
!
!-----------------------------------------------------------------------

        IMPLICIT NONE
        INCLUDE 'mod04_land.inc'

!     INPUTS
      REAL OPTH_NL(NLTAU,NLWAV,NLSIZE)
      REAL MASSCOEF_NL(NLTAU,NLWAV,NLSIZE)
      REAL EXTNORM_NL(NLTAU,NLWAV,NLSIZE)
      REAL RHO_S212(NLTAU,NLSIZE), RHOSTAR(NLTAU,NLWAV,NLSIZE)
      REAL yint_466,slope_466,yint_644,slope_644

!     Intermediate
      REAL ETA_TEMP
      REAL AOD553_TEMP(NLETA)
      REAL RHOSFC212_TEMP(NLETA),ERR644_TEMP(NLETA)
      REAL ERR644_SORT(NLETA)
      REAL NLETA_SORT(NLETA)
      REAL RHOSFC212_SORT(NLETA)
      REAL AOD553_SORT(NLETA)


!     OUTput TAU, ETA and SFCReflectance
      REAL SFC212OUT,ETAOUT,ERR644OUT,AOD553OUT
      REAL SFC212BEST,ETABEST,ERR644BEST
      REAL SFC212TOT,ETATOT,ERR644TOT
      REAL SFC212AVE,ETAAVE,ERR644AVE
      REAL AOD553AVE,AOD553BEST,AOD553TOT
      INTEGER AVE_COUNT


      REAL AODF553, AODC553, AOD553, EXTNORM
      REAL RHOSTAR_TOT(NLTAU,NLWAV)
      REAL ERRWAVE(NLWAV)
      REAL RHO_S212_TOT(NLTAU)
      REAL RHOSTAR_TOT1(NLWAV,NLETA)

!     OUTPUT
      REAL ETA,ERR644,RHOSFC(NLWAV),AOD(NLWAV)
      REAL RHOSFC212
      REAL MASSCON,MASSCONF,MASSCONC
      REAL AODF(NLWAV), AODC(NLWAV)
      INTEGER ETA_FLAG

!     Dummy

      REAL XMIN,XS(100,100),YS(100,100),X(100),Y(100),Denom
      REAL y1
      INTEGER IETA,JETA,ITAU,IWAV,ISIZE,I,NSOLUTION

      SAVE	

!     Initialize

      ETA = -9999
      ERR644 = -9999
      MASSCON = -9999
      DO IWAV = 1, NLWAV
        RHOSFC(NLWAV) = -9999
        AOD(NLWAV) = -9999
      ENDDO
   

! Now loop through ETA values (kind of like ocean)
      NSOLUTION = 0
      DO 199 IETA = 1, NLETA

! Average, best, QCONTROL, FILLVALUE, OUTPUT, etc
        ETA_TEMP = -0.2 + (IETA * 0.1)

! 	Compute total apparent reflectance

        DO ITAU = 1, NLTAU
          DO IWAV = 1, NLWAV
	     RHOSTAR_TOT(ITAU,IWAV) = &
                ETA_TEMP * RHOSTAR(ITAU,IWAV,1) + &
               ((1-ETA_TEMP) * RHOSTAR(ITAU,IWAV,2))
          ENDDO
!	Compute total surface reflectance
          RHO_S212_TOT(ITAU) = &
                ETA_TEMP * RHO_S212(ITAU,1) + &
                ((1-ETA_TEMP) * RHO_S212(ITAU,2))


        ENDDO


!       Interpolate everything based on measured reflectance at 466.


!                   ********COMPUTE FOR TAU FOR  USING
!                           RADAINCE VALUE OF 0.466 UM AND
!                           TAU VALUES OF WAV553
!

         DO ITAU = 1, NLTAU
             X(ITAU) = RHOSTAR_TOT(ITAU,iwave_466)
             Y(ITAU) = OPTH_NL(ITAU,iwave_553,1)
          ENDDO
          Y1=0
          CALL INTERP_EXTRAP(NLTAU,REFW466L,X,Y,Y1,1)
          AOD553_TEMP(IETA) = y1

 
!                 ******* FOR TAUX55 COMPUTE appararent Reflectance FOR
!                    all wavelengths
 
         DO IWAV = 1,NLWAV
            DO ITAU = 1,NLTAU
               X(ITAU)=OPTH_NL(ITAU,iwave_553,1)
               Y(ITAU)=RHOSTAR_TOT(ITAU,IWAV)
            ENDDO
            y1=0
            CALL INTERP_EXTRAP(NLTAU,AOD553_TEMP(IETA),X,Y,Y1,1)
            RHOSTAR_TOT1(IWAV,IETA)=Y1
         ENDDO

 
!                ******* FOR TAUX55 COMPUTE Surface Reflectance at 2.1
!
          DO ITAU = 1,NLTAU
               X(ITAU)=OPTH_NL(ITAU,iwave_553,1)
               Y(ITAU)=RHO_S212_TOT(ITAU)
          ENDDO
          y1=0
          CALL INTERP_EXTRAP(NLTAU,AOD553_TEMP(IETA),X,Y,Y1,1)
          RHOSFC212_TEMP(IETA)=Y1

!       Errors

         ERRWAVE(1) = ABS(REFW466L - RHOSTAR_TOT1(1,IETA)) / REFW466L
         ERRWAVE(3) = ABS(REFW644L - RHOSTAR_TOT1(3,IETA)) / REFW644L
         ERRWAVE(4) = ABS(REFW212L - RHOSTAR_TOT1(4,IETA)) / REFW212L

 
!                ******* Determine error at 644nm
 
 

         ERR644_TEMP(IETA) = ERRWAVE(iwave_644)


 199      ENDDO

! Sort results in ascending error order
          DO IETA = 1, NLETA
             NLETA_SORT(IETA) = IETA
             ERR644_SORT(IETA) = ERR644_TEMP(IETA)
             RHOSFC212_SORT(IETA) = RHOSFC212_TEMP(IETA)
             AOD553_SORT(IETA) = AOD553_TEMP(IETA)
          ENDDO
          CALL SHLSRT_NL(NLETA,ERR644_SORT,NLETA_SORT)
          DO IETA = 1, NLETA
             JETA = NLETA_SORT(IETA)
             RHOSFC212_SORT(IETA) = RHOSFC212_SORT(JETA)
             AOD553_SORT(IETA) = AOD553_SORT(JETA)
          ENDDO

! Best solution


          ETABEST = (-0.2 + (NLETA_SORT(1) * 0.1))
          ERR644BEST = ERR644_SORT(1)
          SFC212BEST = RHOSFC212_SORT(1)
          AOD553BEST = AOD553_SORT(1)

! Average solution

          AOD553TOT = 0.0
          AOD553AVE = -9999
          ETATOT = 0.0
          ERR644TOT = 0.0
          SFC212TOT = 0.0
          ETAAVE = -9999
          ERR644AVE = -9999
          SFC212AVE = -9999
          AVE_COUNT = 0

          DO IETA = 1, NLETA
             IF (ERR644_SORT(IETA) .LE. 0.03 &
                  .AND. AVE_COUNT .LT. 3) THEN
                ERR644OUT = ERR644_SORT(IETA)
                ETAOUT = (-0.2 + (NLETA_SORT(IETA) * 0.1))
                SFC212OUT = RHOSFC212_SORT(IETA)
                AOD553OUT = AOD553_SORT(IETA)
                ETATOT = ETATOT + ETAOUT
                ERR644TOT = ERR644TOT + ERR644OUT
                SFC212TOT = SFC212TOT + SFC212OUT
                AOD553TOT = AOD553TOT + AOD553OUT
                AVE_COUNT = AVE_COUNT + 1
             ENDIF
          ENDDO

          IF (AVE_COUNT .GT. 0) THEN
             ETAAVE = ETATOT / AVE_COUNT
             ERR644AVE = ERR644TOT / AVE_COUNT
             SFC212AVE = SFC212TOT / AVE_COUNT
             AOD553AVE = AOD553TOT / AVE_COUNT
          ENDIF

 
          AVE_COUNT = -1
          ETA = ETABEST
          ERR644 = ERR644BEST
          RHOSFC212 = SFC212BEST
          AOD553 = AOD553BEST
          RHOSFC(iwave_212) = RHOSFC212
          RHOSFC(iwave_644) = yint_644 + slope_644*RHOSFC212
          RHOSFC(iwave_466) = yint_466 + slope_466*RHOSFC(iwave_644)
 

! set  out of bounds Eta to 0 or 1

       IF( ETA .LT. 0) THEN
           ETA=0.0
           ETA_FLAG = -1
       ENDIF
       IF( ETA .GT. 1) THEN
           ETA=1.0
           ETA_FLAG = -1
       ENDIF

!   Compute fine and coarse AOD for mass concentration
           AODF553 = AOD553 * ETA
           AODC553 = AOD553 * (1.-ETA)
 
!                  ******* FOR TAUX55 COMPUTE Mass Concentration
!

!   Fine
          DO ITAU = 1,NLTAU
               X(ITAU)=OPTH_NL(ITAU,iwave_553,1)
               Y(ITAU)=MASSCOEF_NL(ITAU,iwave_553,1)
          ENDDO
          y1=0
          CALL INTERP_EXTRAP(NLTAU,AODF553,X,Y,Y1,1)
          MASSCONF=Y1*AODF553

!  Coarse
          DO ITAU = 1,NLTAU
               X(ITAU)=OPTH_NL(ITAU,iwave_553,2)
               Y(ITAU)=MASSCOEF_NL(ITAU,iwave_553,2)
          ENDDO
          y1=0.0
          CALL INTERP_EXTRAP(NLTAU,AODC553,X,Y,Y1,1)
          MASSCONC=Y1*AODC553

          MASSCON = MASSCONF + MASSCONC

! Compute AOD for fine mode at all wavelengths
         DO IWAV = 1,NLWAV
            DO ITAU = 1,NLTAU
               X(ITAU)=OPTH_NL(ITAU,iwave_553,1)
               Y(ITAU)=EXTNORM_NL(ITAU,IWAV,1)
            ENDDO
            y1=0.0
            CALL INTERP_EXTRAP(NLTAU,AODF553,X,Y,Y1,1)
            EXTNORM=Y1
            AODF(IWAV) = AODF553 * EXTNORM
         ENDDO

! Compute AOD for coarse mode at all wavelengths
! Compute total AOD at all wavelengths
         DO IWAV = 1,NLWAV
            DO ITAU = 1,NLTAU
               X(ITAU)=OPTH_NL(ITAU,iwave_553,2)
               Y(ITAU)=EXTNORM_NL(ITAU,IWAV,2)
            ENDDO
            y1=0.0
            CALL INTERP_EXTRAP(NLTAU,AODC553,X,Y,Y1,1)
            EXTNORM=Y1
            AODC(IWAV) = AODC553 * EXTNORM
            AOD(IWAV) = AODF(IWAV) + AODC(IWAV)
         ENDDO


      RETURN
      END

!*********************************************************************

      SUBROUTINE RETRIEVAL2(&
         yint_466,slope_466,yint_644,slope_644,&
         OPTH_NL,MASSCOEF_NL,RHO_S212,RHOSTAR,&
        AVE_COUNT,ETA,AOD,RHOSFC212,ERR644,RHOSFC,MASSCON)
!-----------------------------------------------------------------------
!!F90
                  
!
!
!!DESCRIPTION:
!               This subroutine retrieves optical thickness, surface
!	 	reflectance and error parameters by comparing
!	        MODIS observations and LUT data. The "continental"
!              aerosol model is assumed.
!
!
! INPUTS
!              REFW466L,REFW644L,REFW212L   reflectance
!              OPTH       array of LUT optical thickness
!              RHO_S212   simulated surface reflectance
!              RHOSTAR    simulated TOA reflectance
!              yint644,yint466,slope466,slope644  
!                 characteristics of VIS/IR surface reflectance relationship
! OUTPUTS
!	       ERR644	  retrieval fitting errors 
!	       AOD        retrieved optical depths 
!	       RHOSFC     retrieved surface reflectance 
!	       ETA        retrieved fine mode weighting
!
!-----------------------------------------------------------------------

        IMPLICIT NONE
        INCLUDE 'mod04_land.inc'

!     INPUTS
      REAL OPTH_NL(NLTAU,NLWAV,NLSIZE)
      REAL MASSCOEF_NL(NLTAU,NLWAV,NLSIZE)
      REAL RHO_S212(NLTAU,NLSIZE), RHOSTAR(NLTAU,NLWAV,NLSIZE)
      REAL yint_466,slope_466,yint_644,slope_644

!     Intermediate

      INTEGER AVE_COUNT

     
      REAL RHOSTAR_TOT(NLTAU,NLWAV)
      REAL ERRWAVE(NLWAV)
      REAL RHO_S212_TOT(NLTAU)
      REAL RHOSTAR_TOT1(NLWAV)

!     OUTPUT
      REAL ETA,ERR644,RHOSFC(NLWAV),AOD(NLWAV)
      REAL RHOSFC212
      REAL MASSCON

!     Dummy

      REAL XMIN,XS(100,100),YS(100,100),X(100),Y(100),Denom
      REAL y1
      INTEGER JETA,ITAU,IWAV,ISIZE,I,NSOLUTION

      SAVE	

!     Initialize

      ETA = -9999
      ERR644 = -9999
      MASSCON = -9999
      DO IWAV = 1, NLWAV
        RHOSFC(NLWAV) = -9999
        AOD(NLWAV) = -9999
      ENDDO

!Now loop through ETA values (kind of like ocean)
      NSOLUTION = 0
! Average, best, QCONTROL, FILLVALUE, OUTPUT, etc

! 	Compute total apparent reflectance

        DO ITAU = 1, NLTAU
          DO IWAV = 1, NLWAV
	     RHOSTAR_TOT(ITAU,IWAV) = RHOSTAR(ITAU,IWAV,1)
          ENDDO
! 	Compute total surface reflectance
          RHO_S212_TOT(ITAU) =  RHO_S212(ITAU,1)
        ENDDO


!       Interpolate everything based on measured reflectance at 466.

!                   ********COMPUTE FOR TAU FOR  USING
!                           RADAINCE VALUE OF 0.466 UM AND
!                           TAU VALUES OF WAV553
!

        DO IWAV = 1, NLWAV
	  DO ITAU = 1, NLTAU
	     X(ITAU) = RHOSTAR_TOT(ITAU,iwave_466)
	     Y(ITAU) = OPTH_NL(ITAU,IWAV,1)
	  ENDDO
          Y1=0
          CALL INTERP_EXTRAP(NLTAU,REFW466L,X,Y,Y1,1)
          AOD(IWAV)=y1
        ENDDO

 
!                 ******* FOR TAUX55 COMPUTE appararent Reflectance FOR
!                     all wavelengths
 
         DO IWAV = 1,NLWAV
            DO ITAU = 1,NLTAU
               X(ITAU)=OPTH_NL(ITAU,iwave_553,1)
               Y(ITAU)=RHOSTAR_TOT(ITAU,IWAV)
	    ENDDO
            y1=0
            CALL INTERP_EXTRAP(NLTAU,AOD(iwave_553),X,Y,Y1,1)
            RHOSTAR_TOT1(IWAV)=Y1
         ENDDO

 
!                 ******* FOR TAUX55 COMPUTE Surface Reflectance at 2.1
 
          DO ITAU = 1,NLTAU
               X(ITAU)=OPTH_NL(ITAU,iwave_553,1)
               Y(ITAU)=RHO_S212_TOT(ITAU)
	  ENDDO
          y1=0
          CALL INTERP_EXTRAP(NLTAU,AOD(iwave_553),X,Y,Y1,1)
          RHOSFC212=Y1

!   	Errors

	 ERRWAVE(1) = ABS(REFW466L - RHOSTAR_TOT1(1)) / REFW466L
	 ERRWAVE(3) = ABS(REFW644L - RHOSTAR_TOT1(3)) / REFW644L
	 ERRWAVE(4) = ABS(REFW212L - RHOSTAR_TOT1(4)) / REFW212L

 
!                 ******* Determine error at 644nm  
 
 

 	 ERR644 = ERRWAVE(iwave_644)

         AVE_COUNT = -1
  
         RHOSFC(iwave_212) = RHOSFC212
         RHOSFC(iwave_553) = -9999
         RHOSFC(iwave_644) = yint_644 + slope_644*RHOSFC212 
         RHOSFC(iwave_466) = yint_466 + slope_466*RHOSFC(iwave_644)

 
!                 ******* FOR TAUX55 COMPUTE Mass Concentration  
 

!   Continental only 
          DO ITAU = 1,NLTAU
               X(ITAU)=OPTH_NL(ITAU,iwave_553,1)
               Y(ITAU)=MASSCOEF_NL(ITAU,iwave_553,1)
          ENDDO
          y1=0
          CALL INTERP_EXTRAP(NLTAU,AOD(iwave_553),X,Y,Y1,1)
          MASSCON=Y1*AOD(iwave_553)


      RETURN
      END

!*********************************************************************

       subroutine SHLSRT_NL(n,arr1,arr2)

!-----------------------------------------------------------------------
! !F77
!
! !DESCRIPTION: From numerical recepies Sorts arr1ay 'arr1' of
!              length 'n' into ascending order ('arr1' and arr2
!              is replaced on OUTPUT by its sorted rearr1angement)
!
!!INPUT PARAMETERS:
!                 N         length of ARR1
!                ARR1       Unsorted array ARR1
!                ARR2       Unsorted array ARR2
!
! !OUTPUT PARAMETERS:
!                ARR1       Sorted array ARR1
!                ARR2       Sorted array ARR2
!
!
 

      implicit none

      integer i,j,k,l,lognb2,m,n,nn
      real aln2i,tiny
      parameter (aln2i = 1.4426950, tiny = 1.e-5)
      real  t,q,arr1(n),arr2(n)
 

      logical flag

      lognb2 = int(alog(float(n))*aln2i+tiny)
      m = n
      do 12 nn = 1,lognb2
        m = m/2
        k = n-m
        do 11 j = 1,k
          i = j

 

          flag=.true.

          do while(flag.eqv..true.)
          l = i+m
          flag=.false.
          if(arr1(l).lt.arr1(i)) then
            t = arr1(i)
            arr1(i) = arr1(l)
            arr1(l) = t
            q = arr2(i)
            arr2(i) = arr2(l)
            arr2(l) = q
 
            i = i-m
            if(i.ge.1) flag=.true.
          endif
          end do 

11      continue
12    continue

      return

      end

!*********************************************************************

      SUBROUTINE INTERP_EXTRAP(M,X1,X,Y,Y1,INDEX)

!----------------------------------------------------------------------
! !F77
!
! !DESCRIPTION: This subroutine is a general purpose routine and
!              interpolates linearly.  Value of y1 is interpolated
!              or extrapolated for x1
!
! !INPUT PARAMETERS:I,M,X1,X,Y
!
! !OUTPUT PARAMETERS:Y1,LOPT
!
! !REVISION HISTORY:
 

      IMPLICIT NONE

      INTEGER IL,LL,LOPT,M,INDEX
      REAL PINTEN,PPHI,SINTEN,SPHI
      REAL  X(M),Y(M),Y1,X1,Diff

      Y1=0.0
      LOPT=0
      LL=M-1
        DO 230 IL=1,LL
!        Extrapolation on lower bound
       IF(X1 .LE.X(1))THEN
          PPHI=X(1)
          SPHI=X(2)
          PINTEN=Y(1)
          SINTEN=Y(2)
           Diff=(SPHI-PPHI)
           if(Diff .Le. 0.0) Diff=1
          Y1=PINTEN+((SINTEN-PINTEN)*((X1-PPHI)/Diff))
               LOPT=1

 

           RETURN

!        Extrapolation on upper bound
       ELSEIF(X1 .GE.X(LL+1)) THEN
           PPHI=X(LL)
           SPHI=X(LL+1)
         PINTEN=Y(LL)
         SINTEN=Y(LL+1)
         Diff=(SPHI-PPHI)
           if(Diff .Le. 0.0) Diff=1
          Y1=PINTEN+((SINTEN-PINTEN)*((X1-PPHI)/Diff))
              LOPT=1

 

           RETURN

!       interpolation
       ELSEIF (X1 .GE.X(IL) .AND.X1.LE. X(IL+1)) THEN
         PPHI=X(IL)
         SPHI=X(IL+1)
         PINTEN=Y(IL)
         SINTEN=Y(IL+1)
         Diff=(SPHI-PPHI)
           if(Diff .Le. 0.0) Diff=1
          Y1=PINTEN+((SINTEN-PINTEN)*((X1-PPHI)/Diff))
          LOPT=1 

           RETURN 

         ENDIF
  230    CONTINUE

  290    RETURN
          END



!***********************************************************************
      SUBROUTINE RNLOOKUP(&
        HANDLE_LUT466,HANDLE_LUT553,HANDLE_LUT644,HANDLE_LUT213,&
        INT_NL0,Fd_NL0,T_NL0,OPTH_NL0,&
       SBAR_NL0,MASSCOEF_NL0,EXTNORM_NL0)

!-----------------------------------------------------------------------
!!F90   
!!DESCRIPTION:
!               This subroutine reads 4 lookup tables
!
!!INPUT PARAMETERS:
!
!!OUTPUT PARAMETERS:
!
!          INT_NL0           Reflectance
!          Fd_NL0           Flux down
!          T_NL0           Transmission factor
!          OPTH_NL0          Optical depth
!          SBAR_NL0          Sbar
!          MASSCOEF_NL0      Mass Concentration coeficient
!          EXTNORM_NL0       Normalized extinction coeficient
!
!!OUTPUT PARAMETERS (INTO INCLUDE FILE)
!          THET0           Array of LUT Solar Zenith Angles
!          THE             Array of LUT Satellite View angles
!          PHI             Array of LUT Azimuth angle differences
!          WAV             Array of LUT wavelengths
!
!!DESIGN NOTES:
!
!  Read from look-up table.
!   wav=wavelengths(2)(.466,.533,.644 and 2.13 micrometers)  4
!   lthet0=number of theta0(solar zenith angle)              9
!   thet0=(0,12,24,36,48,54,60,66,72   degrees)
!   lthe =number of the(observation zenith angle)           16
!   the=(0,6,.....72 degrees every 6 degrees)
!   lphi =number of phi(observation azimuth angle)          31
!   phi=(0,6,12,18,24,30,36,42,48,54,60,66,72,78,84,90,96,
!        102,108,114,120,126,132,138,144,150,156,162,168,
!       174,180)
!   ltau =number of opth(optical thickness)                  7
!   opth=(0.0,0.25,0.50,1.0,2.0,3.0,5.0)
!   lhght = number of hght(observation heights)              1
!   hght=(80.0 km) 
        USE OCIUAAER_Config_Module
        IMPLICIT NONE 
        INCLUDE 'mod04_land.inc'
        INTEGER HANDLE_LUT466,HANDLE_LUT553
        INTEGER HANDLE_LUT644,HANDLE_LUT213
        INTEGER I,IPHI,ITAU,ITAB,NFILE
        INTEGER ITHE,ITHET0,IWAV
        REAL OPTH_NL0(NLTAU,NLWAV,NLTABLE)
        REAL MASSCOEF_NL0(NLTAU,NLWAV,NLTABLE)
        REAL EXTNORM_NL0(NLTAU,NLWAV,NLTABLE)
        REAL SSA_NL0(NLTAU,NLWAV,NLTABLE)
        REAL QEXT_NL0(NLTAU,NLWAV,NLTABLE)
        REAL BEXT_NL0(NLTAU,NLWAV,NLTABLE)
        REAL VEXT_NL0(NLTAU,NLWAV,NLTABLE)
        REAL MEXT_NL0(NLTAU,NLWAV,NLTABLE)
        REAL SBAR_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
        REAL INT_NL0(NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
        REAL Fd_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
        REAL T_NL0(NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
        CHARACTER*1 LINE(132)
        character (len =10):: Sat_Flag
        REAL  Omega0(NLTAU,NLWAV,NLTABLE),ROD(NLWAV),GOD(NLWAV)
        CHARACTER  (len=255) ::  file_name,Extension 
        SAVE
          
 

             file_name = cfg%VIIRS_land  
             Extension = 'lookup_land_w0488.npp3' 
         OPEN (HANDLE_LUT466, FILE = trim(file_name)//trim(Extension),status='old')  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ! Keep this Table for 0.544 MODIS  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            Extension = 'lookup_land_w0554.v4'  
         OPEN (HANDLE_LUT553, FILE = trim(file_name)//trim(Extension),status='old')  
!  
           Extension = 'lookup_land_w0670.npp3' 
         OPEN (HANDLE_LUT644,FILE = trim(file_name)//trim(Extension),status='old')    
!  
          Extension = 'lookup_land_w2257.npp3' 
        OPEN (HANDLE_LUT213,FILE = trim(file_name)//trim(Extension),status='old')  
!        print*,'VIRRS table Land'
        
        
   

  
     
      DO 10 IWAV = 1,NLWAV
      DO 20 ITAB=1,NLTABLE

        IF (IWAV .EQ. 1) NFILE = HANDLE_LUT466        
        IF (IWAV .EQ. 2) NFILE = HANDLE_LUT553
        IF (IWAV .EQ. 3) NFILE = HANDLE_LUT644
        IF (IWAV .EQ. 4) NFILE = HANDLE_LUT213
        
        
        
 
! Reads observation zenith angles(the),and observation azimuth angles
! from look-up tables
 
         
        READ(NFILE,1111)(THE_NL(I),I=1,NLTHE)
        READ(NFILE,1000)(PHI_NL(I),I=1,NLPHI) 
         
        
 
! Reads wavelength(wav),optical thickness(opth),solar zenith angle(thet0),
! reflectance ofatmosphere(sbarw) from look-up tables.
 

        DO 100 ITAU=1,NLTAU
          READ(NFILE,*)
          READ(NFILE,1001) SSA_NL0(ITAU,IWAV,ITAB),&
              QEXT_NL0(ITAU,IWAV,ITAB),&
              BEXT_NL0(ITAU,IWAV,ITAB),&
              VEXT_NL0(ITAU,IWAV,ITAB),&
              MEXT_NL0(ITAU,IWAV,ITAB),&
             MASSCOEF_NL0(ITAU,IWAV,ITAB)
          READ(NFILE,1002) WAV_NL(IWAV),OPTH_NL0(ITAU,IWAV,ITAB),&
              ROD(IWAV),GOD(IWAV)
 
          DO 101 ITHET0=1,NLTHET0
             READ(NFILE,1003) &
               THET0_NL(ITHET0), MU0_NL(ITHET0),&
               SBAR_NL0(ITHET0,ITAU,IWAV,ITAB), &
               Fd_NL0(ITHET0,ITAU,IWAV,ITAB)

 
! Reads transmission as a function of observation zenith angle,
! optical thickness
 

             READ(NFILE,1004)&
                (T_NL0(ITHE,ITHET0,ITAU,IWAV,ITAB),ITHE=1,NLTHE)       
             READ(NFILE,*)
 
! Reads atmospheric radiance (int) as a function of solar zenith angle,
! optical thickness, height, observation zenith angle and azimuth angle
! from look-up table.
 
              DO 103 ITHE=1,NLTHE
                 READ(NFILE,1005)&
          (INT_NL0(IPHI,ITHE,ITHET0,ITAU,IWAV,ITAB),IPHI=1,NLPHI)
103           CONTINUE
101       CONTINUE
100     CONTINUE

!  Set extinction parameters for "AOD = 0.0" case

        QEXT_NL0(1,IWAV,ITAB) = QEXT_NL0(2,IWAV,ITAB)
        BEXT_NL0(1,IWAV,ITAB) = BEXT_NL0(2,IWAV,ITAB)
        VEXT_NL0(1,IWAV,ITAB) = VEXT_NL0(2,IWAV,ITAB)
        MEXT_NL0(1,IWAV,ITAB) = MEXT_NL0(2,IWAV,ITAB)
        MASSCOEF_NL0(1,IWAV,ITAB) = MASSCOEF_NL0(2,IWAV,ITAB)

20    CONTINUE
10    CONTINUE
        
       DO ITAB=1,NLTABLE
       DO ITAU=1,NLTAU
       DO IWAV=1,NLWAV
       EXTNORM_NL0(ITAU,IWAV,ITAB) = &
           QEXT_NL0(ITAU,IWAV,ITAB) / QEXT_NL0(ITAU,iwave_553,ITAB)
       ENDDO
       ENDDO
       ENDDO 
         close(HANDLE_LUT466)
         close(HANDLE_LUT553)
         close(HANDLE_LUT644)
         close(HANDLE_LUT213)
1000  FORMAT(5X,16F8.3)
1111  FORMAT(5X,15F8.3)
1001  FORMAT(6(5x,f9.4))
1002  FORMAT(10x,F9.4,3(4x,f9.4))
1003  FORMAT(5x,f8.3, 5x, f8.5, 2(5x,e11.4))
1004  FORMAT(3X,15e11.5)
1005  FORMAT(16e11.5)

      RETURN
      END


 


!*********************************************************************
      SUBROUTINE INT_ELEV(EQWAV_NL,INT_NL,Fd_NL,T_NL,FdT_NL,OPTH_NL,&
        SBAR_NL,REF_RAY_NL,MHGHT)
!----------------------------------------------------------------------
!!F90 
!!DESCRIPTION:  Subroutine INTELEV interpolates the lookup
!               reflectances to the target elevation.
!               Basically it is fudge:
!              Interpolate between wavelengths to simulate
!               elevation by a longer wavelength
!
!!INPUT PARAMETERS:
!         WAV_NL           wavelengths
!	  INT_NL		  radiance
!	  Fd_NL		  flux down
!         T_NL            transmission
!         FdT_NL          2 way transmission factor
!	  SBAR_NL	  sbar
!        OPTH_NL          optical depth
!        MGHTH           elevation
!
!!OUTPUT PARAMETERS 
!	  INT_NL		  interpolated radiance
!	  Fd_NL		          interpolated flux down
!         T_NL                    interpolated transmission
!	  FdT_NL		  interpolated transmission
!	  SBAR_NL	  interpolated sbar
!         OPTH_NL          interpolated optical depth
!         REF_RAY        Reflectance for rayleigh only
!        EQWAV_NL        Interpolated wavelengths
!
 
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04_land.inc'

!     Inputs
      REAL  MHGHT
      
!     Inputs and Outputs
      REAL INT_NL(NLTAU,NLWAV,NLSIZE),FdT_NL(NLTAU,NLWAV,NLSIZE),&
          Fd_NL(NLTAU,NLWAV,NLSIZE), T_NL(NLTAU,NLWAV,NLSIZE),&
          SBAR_NL(NLTAU,NLWAV,NLSIZE),OPTH_NL(NLTAU,NLWAV,NLSIZE)
      REAL REF_RAY_NL(NLWAV)
      REAL INT_NL9(NLTAU,NLWAV,NLSIZE),FdT_NL9(NLTAU,NLWAV,NLSIZE),&
          Fd_NL9(NLTAU,NLWAV,NLSIZE), T_NL9(NLTAU,NLWAV,NLSIZE),&
          SBAR_NL9(NLTAU,NLWAV,NLSIZE),OPTH_NL9(NLTAU,NLWAV,NLSIZE)

      REAL ROD_1013(NLWAV)
      REAL p, p0, expfactor
      PARAMETER (p0 = 1013.0)
      REAL ROD_PRES(NLWAV)
      REAL EQWAV_NL(NLWAV)

      REAL lambda0, lambda1, lambda2, diff0, exp0, exp1, exp2
      REAL tau0, tau1, tau2

!     Dummy
      INTEGER ICASE,IJ,ISIZE,ITAU,IWAV,IWAV1,IWAV2,JWAV
      REAL  X(2),Y(2),W(2),V(2),Z(2),T(2),U(2)
      REAL  Y1,W1,V1,Z1,T1,U1
      INTEGER LL


!     Estimate surface pressure (hypsometric EQ)
      p = p0 * exp(-(MHGHT/7.5))
      

!     Estimate ROD at nominal wavelenghts at p0 and at pres
      DO IWAV = 1, 3
         expfactor = -4.15 + (0.2 * WAV_NL(IWAV))
         ROD_1013(IWAV) = 0.0088 * WAV_NL(IWAV)**expfactor
         ROD_PRES(IWAV) = 0.0088*(p/p0) * WAV_NL(IWAV)**expfactor

!    Estimate equivalent wavelengths for ROD at pressure
         lambda1 = 0.1
         lambda2 = 2.0
         diff0 = 99.
         DO WHILE (diff0 .GT. 0.00001)
           lambda0 = (lambda1 + lambda2) / 2.
           exp0 = -4.15 + 0.2*lambda0
           exp1 = -4.15 + 0.2*lambda1
           exp2 = -4.15 + 0.2*lambda2
           tau0 = 0.0088*lambda0**exp0
           tau1 = 0.0088*lambda1**exp1
           tau2 = 0.0088*lambda2**exp2
           IF (tau1 .GT. ROD_PRES(IWAV) .AND. &
        tau2 .LT. ROD_PRES(IWAV)) THEN
              IF (tau0 .GT. ROD_PRES(IWAV)) THEN
                lambda1 = (lambda1 + lambda2)/2.
              ELSE 
                lambda2 = (lambda1 + lambda2)/2.
              ENDIF
           ENDIF
           diff0 = ABS(tau0 - ROD_PRES(IWAV))
         ENDDO
         EQWAV_NL(IWAV) = lambda0 
       ENDDO
           
!     INterpolate lookup tables to equiv Waves (let's start only with
!      1st two wavelengths until we derive 0.86 table)


        DO IJ = 1,2
           X(IJ)=0.0
           Y(IJ)=0.0
           V(IJ)=0.0
           W(IJ)=0.0
           Z(IJ)=0.0
           U(IJ)=0.0
           T(IJ)=0.0
        ENDDO
        DO 5 ISIZE= 1,NLSIZE
          DO 15 IWAV=1,3
            DO 25  ITAU  = 1,NLTAU
              LL=0
              Y1=0.0
              W1=0.0
              Z1=0.0
              V1=0.0
              IF (IWAV .EQ. 3) THEN
               IWAV1 = IWAV-1
               IWAV2 = IWAV
              ELSE
               IWAV1 = IWAV
               IWAV2 = IWAV+1
              ENDIF
              DO 60  JWAV = IWAV1, IWAV2
                LL=LL+1

!    Interpolate on log log scale
 
                X(LL)=ALOG(WAV_NL(JWAV))
                Y(LL)=ALOG(INT_NL(ITAU,JWAV,ISIZE))
                W(LL)=ALOG(FdT_NL(ITAU,JWAV,ISIZE))
                Z(LL)=ALOG(SBAR_NL(ITAU,JWAV,ISIZE))
                V(LL)=OPTH_NL(ITAU,JWAV,ISIZE)
                U(LL)=ALOG(Fd_NL(ITAU,JWAV,ISIZE))
                T(LL)=ALOG(T_NL(ITAU,JWAV,ISIZE))
                IF (OPTH_NL(ITAU,JWAV,ISIZE) .GT. 0.) THEN
                  V(LL)=ALOG(OPTH_NL(ITAU,JWAV,ISIZE))
                ENDIF
60            CONTINUE

              CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,Y,Y1,1)
              CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,W,W1,1)
              CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,Z,Z1,1)
              CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,V,V1,1)
              CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,U,U1,1)
              CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,T,T1,1)

              INT_NL9(ITAU,IWAV,ISIZE) = EXP(Y1)
              FdT_NL9(ITAU,IWAV,ISIZE) = EXP(W1)
              Fd_NL9(ITAU,IWAV,ISIZE) = EXP(U1)
              T_NL9(ITAU,IWAV,ISIZE) = EXP(T1)
              SBAR_NL9(ITAU,IWAV,ISIZE) = EXP(Z1)
              OPTH_NL9(ITAU,IWAV,ISIZE) = EXP(V1)
              IF (V1 .EQ. 0.) THEN
                OPTH_NL9(ITAU,IWAV,ISIZE) = V1
              ENDIF


25          CONTINUE
15        CONTINUE
 5      CONTINUE

        
        DO ISIZE= 1,NLSIZE
          DO IWAV=1,3
            DO ITAU  = 1,NLTAU
              INT_NL(ITAU,IWAV,ISIZE) = INT_NL9(ITAU,IWAV,ISIZE)
              FdT_NL(ITAU,IWAV,ISIZE) = FdT_NL9(ITAU,IWAV,ISIZE)
              Fd_NL(ITAU,IWAV,ISIZE) = Fd_NL9(ITAU,IWAV,ISIZE)
              T_NL(ITAU,IWAV,ISIZE) = T_NL9(ITAU,IWAV,ISIZE)
              SBAR_NL(ITAU,IWAV,ISIZE) = SBAR_NL9(ITAU,IWAV,ISIZE)
              OPTH_NL(ITAU,IWAV,ISIZE) = OPTH_NL9(ITAU,IWAV,ISIZE)
            ENDDO
            REF_RAY_NL(IWAV) = INT_NL(1,IWAV,ISIZE)
          ENDDO
        ENDDO
          
        
        RETURN
        END


!*************************************************
        SUBROUTINE AEROSOL_MAP(HANDLE_LUTMAP,IMONTH,AEROSOL)
!----------------------------------------------------
!!F90
                  
!
!!DESCRIPTION:  Subroutine AEROSOL_MAP reads
!               the aerosol map (1 degree resolution)
!               to determine expected aerosol type at
!               a given location and season
!
!!INPUT PARAMETERS:
!         HANDLE_LUTMAP   Logical Unit # for file
!	  IMONTH          Integer month 1-12
!
!!OUTPUT PARAMETERS 
!	  AEROSOL	  360x180 map of aerosol type for appropriate season
 
        USE OCIUAAER_Config_Module  
        IMPLICIT none

        INTEGER IMONTH, IMONTH1,iseason
        INTEGER HANDLE_LUTMAP
	INTEGER nlon, nlat
	PARAMETER (nlon = 360, nlat = 180)
 	INTEGER AEROSOL(nlon,nlat)
 	INTEGER AEROSOL_all(nlon,nlat,4)
        INTEGER ilat, ilon
         CHARACTER  (len=255) ::file_name ,Extension 
! Loops through each wavelength and aerosol model
      
        CHARACTER*3 MMMs(4), MMM
        DATA MMMs/'DJF','MAM','JJA','SON'/
        file_name = cfg%VIIRS_land  
         Extension = 'aerosol_land_map.v3'
         OPEN (HANDLE_LUTMAP, FILE = trim(file_name)//trim(Extension),status='old')  
        DO iseason = 1, 4
         READ(HANDLE_LUTMAP,*) 
         DO ilat = 1, nlat
            READ(HANDLE_LUTMAP,99) &
            (AEROSOL_all(ilon,ilat,iseason), ilon=1, nlon)
         ENDDO
        ENDDO

        IMONTH1 = IMONTH
        IF (IMONTH .EQ. 12) THEN 
          IMONTH1 = 0
        ENDIF
        iseason = IMONTH1/3 + 1
        MMM = MMMs(iseason)
        DO ilon = 1, nlon
          DO ilat = 1, nlat
            AEROSOL(ilon,ilat) = AEROSOL_all(ilon,ilat,iseason)
          ENDDO
        ENDDO
        close(HANDLE_LUTMAP)
 99     FORMAT (i1,359(1x,i1))
  
        RETURN 
        END
!********************************************************************
       SUBROUTINE Average_land(ISWATH,CLDMSK_500,W354_SYN,W388_SYN,W470_SYN,W550_SYN,&
       W659_SYN,W865_SYN,W124_SYN,W164_SYN,W213_SYN,START_500,END_500,&
       START_250,END_250,START_1KM,END_1KM,REFW466_L,REFW550_L,REFW644_L,&
       REFW124_L,REFW164_L,REFW866_L,REFW212_L,SDW466_L,SDW550_L,&
       SDW644_L,SDW124_L,SDW164_L,SDW866_L,SDW212_L,IFINISH,&
       THR213MIN,THR213MAX,Num_Pixels_Used,IERROR,NUMRED,NUMBLUE,Save_index,&
       Idata,W412_SYN,W443_SYN,W0P75_SYN,REFW412_L,REFW443_L,REFW0P75_L,&
       SDW412_L,SDW443_L,SDW0P75_L,iscan,Good_pixels_Land,SnowMsk_Ratio,&
       cloudy_pixels,REFW354_L,REFW388_L,SDW354_L,SDW388_L)
!----------------------------------------------------------------------
!!F90           
!
!!DESCRIPTION:
!
!   This subroutine processes 10*10 pixel box for cloud detection and
!   finds the average reflectance for red and blue channels. Surface
!   Reflectance from wavelength 2.13 . This surface Reflectance and
!   average Reflectance for red and blue channel are send to lookup
!   table and optical thickness is derived.
!
!!INPUT PARAMETERS:
!
!      IGRIDX        Number of pixels in x-direction
!      IGRIDY        Number of pixels in y-direction
!      W470_SYN      Reflectance for wav=0.47um
!      W213_SYN      Reflectance for wav=2.13um
!      W659_SYN      Reflectance for wav=0.66um
!      CLDMSK_250    Cloud mask at 250 m resolution (1=cloudy,o=clear)
!      MTHET0        Measured solar zenith angle in deg.
!      MTHET         Measured viewangle from ground in deg.
!      MDPHI         Measured relative azimuth angle deg.
!      THR213MIN     Minumum threshold for 2.13 micro wavelength.
!      THR213MAX     Maximum threshold for 2.13 micro wavelength.
!
 
       IMPLICIT NONE 
       save
       INCLUDE 'mod04_land.inc' 
       INCLUDE 'read_Sat_MODIS.inc'
      INTEGER IFINISH,IX,IY
      INTEGER  IERROR,IOPT,INUM
      INTEGER START_500,END_500,START_250,END_250,START_1KM,END_1KM
      INTEGER ICLDBLUE,ICLDRED,I213,indX(Max_Pixels_L),N20,N50,INN
!      W354_SYN,W388_SYN   UV channels
       REAL W354_SYN(ISWATH_B,ILINE)
       REAL W388_SYN(ISWATH_B,ILINE)

      REAL W470_SYN(ISWATH_B,ILINE)
      REAL W550_SYN(ISWATH_B,ILINE)
      REAL W213_SYN(ISWATH_B,ILINE)
      REAL W124_SYN(ISWATH_B,ILINE)
      REAL W164_SYN(ISWATH_B,ILINE)
      REAL w865_SYN(ISWATH_B,ILINE)
      REAL W659_SYN(ISWATH_B,ILINE)
! NPP channels      
      REAL W412_SYN(ISWATH_B,ILINE)
      REAL W443_SYN(ISWATH_B,ILINE)
      REAL W0P75_SYN(ISWATH_B,ILINE)
      Integer SnowMsk_Ratio(ISWATH_B,ILINE)
      
      INTEGER NUMRED(ISWATH_B,ILINE),NUMBLUE(ISWATH_B,ILINE)
      INTEGER  I,Idata,iscan 
      INTEGER CLDMSK_500(ISWATH_B,ILINE)
      REAL THR213MIN,THR213MAX
      INTEGER Num_Pixels_Used,MODIS_waves,UV_waves
       parameter( MODIS_waves=7, UV_waves=2)
      REAL   Ref_Inter_Land(MODIS_waves+UV_waves,Max_Pixels_L)
      REAL   Array(Max_Pixels_L),Array_new(Max_Pixels_L) 
      INTEGER  Good_pixels_Land(MODIS_waves+UV_waves) 
      REAL     array_interm(MODIS_waves+UV_waves,Max_Pixels_L)
       REAL     ref_allwav(MODIS_waves+UV_waves)
      REAL REFW466_L,REFW550_L,sd_allwav(MODIS_waves+UV_waves)
      REAL REFW644_L,REFW124_L,REFW164_L,REFW866_L,REFW212_L
      REAL SDW466_L,SDW550_L
      REAL SDW644_L,SDW124_L,SDW164_L,SDW866_L,SDW212_L
      REAL ref_val,sd_val
      REAL REFW412_L,REFW443_L,REFW0P75_L,SDW412_L,SDW443_L,SDW0P75_L
      REAL REFW354_L,REFW388_L,SDW354_L,SDW388_L
       Integer Save_index(NUMCELLS_L,MaxPixels_left_L)
       integer Save_Array_index(MaxPixels_left_L),Number_pixels 
      integer cloudy_pixels,Pixels_for_ret
! Initialization
 
           
         Num_Pixels_Used=0
         IFINISH=0
         IERROR=4
          REFW466_L=-9999
          REFW550_L=-9999
          REFW644_L=-9999
          REFW124_L=-9999
          REFW164_L=-9999
          REFW866_L=-9999
          REFW212_L=-9999 
          SDW466_L=-9999
          SDW550_L =-9999
          SDW644_L=-9999
          SDW124_L=-9999
          SDW164_L=-9999
          SDW866_L=-9999
          SDW212_L=-9999 
          
! NPP channels          
          REFW412_L=-9999
          REFW443_L=-9999
          REFW0P75_L=-9999
          SDW412_L=-9999
          SDW443_L=-9999
          SDW0P75_L=-9999
!UV   
          REFW354_L=-9999
          REFW388_L=-9999  
          SDW354_L=-9999
          SDW388_L=-9999        
         
         
      CALL Cldmask_NDVI(ISWATH,CLDMSK_500,W354_SYN,W388_SYN,W470_SYN,W550_SYN,&
        W659_SYN,W865_SYN,W124_SYN,W164_SYN,W213_SYN,&
        START_500,END_500,THR213MIN,THR213MAX,Ref_Inter_Land,&
        ICLDBLUE,ICLDRED,I213,NUMRED,NUMBLUE,W412_SYN,W443_SYN,W0P75_SYN,&
        START_1KM,END_1KM,MODIS_waves,UV_waves,idata,iscan,SnowMsk_Ratio,&
        cloudy_pixels)
 
!         print*,'I213,ICLDBLUE',I213,ICLDBLUE
 
! I213>0 and Icld>0 indicate that the threshold for 2.13 micron
! is met and the data set is cloud free as defined.
 

      IF(I213 .GT. 0 .AND.ICLDBLUE .GT.0)THEN
 
! Sort in asending order and get the index for Reflactances at 0.66 um
         do ix=1,ICLDBLUE
         array(ix)=Ref_Inter_Land(3,ix) 
         enddo 
           inum=ICLDBLUE  
             CALL INDEXX(ICLDBLUE,array,INDX)  
              N20=INT(REAL(inum)/REAL(5))
              N50=INT(REAL(inum)/REAL(2))
              if ( N20 .eq. 0) N20=1
              if ( N50 .eq. 0) N50=1 
! Setting Threshold for pixel selection based on pixels left after               
! Rejecting  rightest and darkest cloudfree pixels to eliminate
! the residual shadows and sub_cloud pixels. Adding 1 to make it whole number.
!  
!! thresholding is done based on10% of pixels left after throwing out dark and bright pixels             
 
           Pixels_for_ret =     (Int(Iline*Iline)/2) -(Int(Iline*Iline)/5)
              Number_pixels= (Pixels_for_ret * 0.10)   
            
        IF( (N50-N20) .gt. Number_pixels) then
          IFINISH =1
           
            DO IY=1,MODIS_waves+UV_waves
                inn=0
                Good_pixels_Land(IY)=0  
                If( Iline .eq.1) then 
                   array_new(1)=Ref_Inter_Land(IY,1)  
                IF ( array_new(1) .GT.0. .AND.  array_new(1) .LE.1)then
                Good_pixels_Land(IY)= Good_pixels_Land(IY)+1 
                array_interm(IY,Good_pixels_Land(IY))= array_new(1)  
                endif
                 Else 
                 DO ix=N20,N50  
                     inn=inn+1 
                     array_new(inn)=Ref_Inter_Land(IY,indx(ix))  
                if( Iy .eq.3)Save_Array_index(inn)=indx(ix) 
                ENDDO  
!Throw away bad pixles in all wavelengths
             DO IX=1,inn
              if ( array_new(IX) .GT.0. .AND.  array_new(IX) .LE.1)then
                Good_pixels_Land(IY)= Good_pixels_Land(IY)+1 
                array_interm(IY,Good_pixels_Land(IY))= array_new(ix) 
              if( Iy .eq.3)&
          save_index(Idata,Good_pixels_Land(IY))=Save_Array_index(Ix) 
               endif
          ENDDO 
! Endif  for line = 1          
           endif
!  ENDDO for wavelengths
          ENDDO
! Report  valid pixels used for Blue channel at 500 meter resolution

           Num_Pixels_Used=Good_pixels_Land(3) 
            DO iy=1,MODIS_waves+UV_waves
! Call to subroutine ave_std computes average and standard deviation for
! Reflectances

          IF(Good_pixels_Land(IY) .gt.0) then 
              DO IX=1,Good_pixels_Land(IY)
              array(IX)= array_interm(IY,ix)
             
              ENDDO
            call ave_std(array, Good_pixels_Land(IY),ref_val,sd_val)
             ref_allwav(iy)=ref_val
             sd_allwav(iy)=sd_val  
             
         ELSE
              ref_allwav(iy)=-99999.0
              sd_allwav(iy)=-9999.0
         ENDIF
!  ENDDO for number of wavelengths
        ENDDO
         
! Check if  3 wavelengths have valid range          
       if(  ref_allwav(1) .gt.0 .and. &
       ref_allwav(3) .gt. 0 .and. ref_allwav(7)  .gt. 0) then
     
          do iy = 1,MODIS_waves+UV_waves
   
 
!   WAVE 0.470 UM
 
          if(iy .eq.1)then
          REFW466_L=ref_allwav(iy)
          SDW466_L=sd_allwav(iy)
!
!   WAVE 0.550 UM
!
          elseif(iy .eq.2)then
          REFW550_L=ref_allwav(iy)
          SDW550_L=sd_allwav(iy)
!
!   WAVE 0.659 UM
!
         elseif (iy .eq.3)then
         REFW644_L=ref_allwav(iy)
         SDW644_L=sd_allwav(iy) 
!
!   WAVE 0.865 UM
!
         elseif(iy .eq.4)then
        REFW866_L=ref_allwav(iy)
        SDW866_L=sd_allwav(iy)

!
!   WAVE 1.24 UM
!
         elseif(iy .eq.5)then
         REFW124_L=ref_allwav(iy)
         SDW124_L=sd_allwav(iy)
         
        
 
!   WAVE 1.64 UM
!
      elseif(iy .eq. 6)then
        REFW164_L=ref_allwav(iy)
        SDW164_L=sd_allwav(iy)
 
!
!   WAVE 2.13  UM
!
         elseif(iy .eq. 7)then
        REFW212_L=ref_allwav(iy)
        SDW212_L=sd_allwav(iy) 
          
        
!   WAVE 0.354 um  
!
         elseif(iy .eq. 8)then
         REFW354_L=ref_allwav(iy)
         SDW354_L=sd_allwav(iy)  
                  
         
!   WAVE 0.388 um  
!
         elseif(iy .eq. 9)then
        REFW388_L=ref_allwav(iy)
        SDW388_L=sd_allwav(iy) 
        
        endif
 
           
! enddo for wavelengths
        ENDDO
        
!        print*,'wave',REFW466_L,REFW550_L,REFW644_L,REFW866_L,REFW124_L,REFW164_L,REFW212_L
        
! Else and endif for to Check if  3 wavelengths have valid range else put fill values        
      Else
          
         IERROR = 2
         REFW466_L=-9999
         REFW550_L=-9999
         REFW644_L=-9999
         REFW124_L=-9999
          REFW164_L=-9999
          REFW866_L=-9999
          REFW212_L=-9999
          SDW466_L=-9999
          SDW550_L =-9999
          SDW644_L=-9999
          SDW124_L=-9999
          SDW164_L=-9999
          SDW866_L=-9999
          SDW212_L=-9999 
          REFW412_L  =-9999
          REFW443_L=-9999
          REFW0P75_L=-9999
          SDW412_L=-9999
          SDW443_L=-9999
          SDW0P75_L=-9999
       ENDIF
        
!  ENDIF for N50-N20
         IERROR = 3 
         ENDIF
        
!ENDIF for  I213>0 and Icld>0 indicate that the threshold for 2.13 micron
! is met and the data set is cloud free as defined.
       
        ENDIF 
           RETURN
           END

!***********************************************************************
      SUBROUTINE Cldmask_NDVI(ISWATH,CLDMSK_500,W354_SYN,W388_SYN,W470_SYN,W550_SYN,&
        W659_SYN,W865_SYN,W124_SYN,W164_SYN,W213_SYN,&
        START_500,END_500,THR213MIN,THR213MAX,Ref_Inter_Land,&
       ICLDBLUE,ICLDRED,I213,NUMRED,NUMBLUE,W412_SYN,W443_SYN,W0P75_SYN,&
       START_1KM,END_1KM,MODIS_waves,UV_waves,idata,iscan,SnowMsk_Ratio,&
       cloudy_pixels)
!----------------------------------------------------------------------
!!F90
!!DESCRIPTION:
!              The following subroutine finds the dark targets using a
!              threshold in the 2.13 micron channel.it averages the
!              cloud free pixels in red channel to the 0.5 km resolution.
!              If too cloudy or too bright at 2.13 it leaves the value
!              as zero.
!
!!INPUT PARAMETERS:
!
!           CLDMSK_250      Cloud mask at 250 m resolution (0=cloudy,1=clear)
!             W213_SYN      Reflectance for at 2.13 um
!             W470_SYN      Reflectance for at 0.47 um
!             W659_SYN      Reflectance for at 0.66 um
!               THR213      Threshold to detect dark pixels for 2.13 um
!               IGRIDX      Number of pixels in 1x1 km resolution
!               IGRIDY      Number of pixels in 1x1 km resolution
!
!OUTPUT PARAMETERS:
!
!                 I213      number of pixels indicating cloudfree
!              ICLDRED      Number of pixels in red channel(cloudfree)
!             ICLDBLUE      Number of pixels in blue channel(cloudfree)
!             Ref_Inter_land Reflectance in 7 channels
!
 
        IMPLICIT NONE
        save
       INCLUDE 'mod04_land.inc'
       INCLUDE 'read_Sat_MODIS.inc'
      INTEGER IFINISH
      INTEGER  IERROR,IOPT,Iwave,I213
      INTEGER START_500,END_500,START_250,END_250,START_1KM,END_1KM
      REAL W470_SYN(ISWATH_B,ILINE)
      REAL W550_SYN(ISWATH_B,ILINE)
      REAL W213_SYN(ISWATH_B,ILINE)
      REAL W124_SYN(ISWATH_B,ILINE)
      REAL W164_SYN(ISWATH_B,ILINE)
      REAL w865_SYN(ISWATH_B,ILINE)
      REAL W659_SYN(ISWATH_B,ILINE)
! NPP channels      
      REal W412_SYN(ISWATH_B,ILINE)
      REAL W443_SYN(ISWATH_B,ILINE)
      REAL W0P75_SYN(ISWATH_B,ILINE)
      
! uv
      REAL W354_SYN(ISWATH_B,ILINE)
      REAL W388_SYN(ISWATH_B,ILINE)      
      
      Integer SnowMsk_Ratio(ISWATH_B,ILINE)
      INTEGER MODIS_waves,UV_waves,idata,iscan
      INTEGER IXX,IYY 
      INTEGER CLDMSK_500(ISWATH_B,ILINE)
      REAL  W865_SYN_ForAve 
      REAL  W659_SYN_ForAve
      REAL   Ref_Inter_Land(MODIS_waves+UV_waves,Max_Pixels_L)
      real THR213MIN,THR213MAX,aa
      INTEGER NUMRED(ISWATH_B,ILINE),NUMBLUE(ISWATH_B,ILINE)
      integer ICLDBLUE,ICLDRED,NUM_indx_x,J,I,IWAV,JBLUE,IBLUE,NN
      Integer  NUM_indx_y, NUM_indx_z, NUM_indx_p,NUM_indx_Q
      INTEGER IMASK,IRED,JMASK,JRED,NUMCFR,NOWATER,II,JJ,IJ,IK
       Integer  Y2_offset, X2_offset,l,p,X2,Y2      
      integer   cloudy_pixels  

 
!  Initialization of variables and arrays
 
      
 

       I213 =0
       ICLDBLUE =0
       ICLDRED  =0
       NUM_indx_x=0
       NUM_indx_y=0
       NUM_indx_Z=0
       NUM_indx_p=0
       NUM_indx_Q=0
       cloudy_pixels =0
      
          DO J = 1,ILINE 
          DO I = 1,ISWATH
          NUMBLUE(I,J)=0 
         ENDDO
         ENDDO
        
         do iwav=1,(MODIS_waves+UV_waves)
         do j=1,Max_Pixels_L
         Ref_Inter_Land(Iwav,j)=0
         ENDDO
         ENDDO
        
!Loop around 1km to add NPP wavelengths         
         DO  IXX=START_1KM,END_1KM  
         DO  IYY = 1,IGRIDY 
         if( CLDMSK_500(IXX,IYY) .EQ. 0)cloudy_pixels=cloudy_pixels+1 
! If snowy pixels set cloud mask to Zero  
!          if(SnowMsk_Ratio(IXX,IYY) .eq. 0)CLDMSK_500(IXX,IYY) =0 
            ICLDRED=0  
            NUMCFR=0
            NOWATER=0 
! Check if therhold for 2.13 um is met
! (missing pixels or noisy detectors are ignored)
          
              
             
            IF(W213_SYN(IXX,IYY) .LE. THR213MAX .AND.&
             W213_SYN(IXX,IYY) .GT. THR213MIN .and.&
             CLDMSK_500(IXX,IYY) .EQ. 1) THEN  
     
                  ICLDRED=ICLDRED+1  
                  I213=I213+1
         
           IF(W865_SYN(IXX,IYY) .GT. 0.0 .AND. &
           W659_SYN(IXX,IYY) .GT. 0.0) THEN  
 ! NDVI test according to Eric Vermote ( Not working on Synthetic data)
                IF((W865_SYN(IXX,IYY)-W659_SYN(IXX,IYY))/ &
                 (W865_SYN(IXX,IYY)+W659_SYN(IXX,IYY)).GT.0.1) THEN  
             IF(W470_SYN(IXX,IYY) .GT. 0.0) THEN
               NUM_indx_x=NUM_indx_x+1 
              Ref_Inter_Land(3,NUM_indx_x)= W659_SYN(IXX,IYY) 
              Ref_Inter_Land(4,NUM_indx_x)= W865_SYN(IXX,IYY)
              Ref_Inter_Land(1,NUM_indx_x)= W470_SYN(IXX,IYY)
              Ref_Inter_Land(2,NUM_indx_x)= W550_SYN(IXX,IYY) 
              Ref_Inter_Land(5,NUM_indx_x)= W124_SYN(IXX,IYY) 
              Ref_Inter_Land(6,NUM_indx_x)= W164_SYN(IXX,IYY)
              Ref_Inter_Land(7,NUM_indx_x)= W213_SYN(IXX,IYY)
              Ref_Inter_Land(8,NUM_indx_x)= W354_SYN(IXX,IYY)
              Ref_Inter_Land(9,NUM_indx_x)= W388_SYN(IXX,IYY)
              ICLDBLUE = ICLDBLUE + 1 
!              NUMBLUE(NUM_indx_x,J)=NUMBLUE(NUM_indx_x,J)+1   
! Endif for W470
              Endif  
! NDVI test endif
           ENDIF 
!   W865 and w65 Valid data
          ENDIF 
! 2.13  & cldmask
          ENDIF  
! Enddo for 1km data           
        Enddo
        Enddo 
        
       RETURN
       END
       
       
!*********************************************************************
      SUBROUTINE INTANGLE_RAY_NL(MTHET0,MTHET,MDPHI,&
      	     INT_NL0,REF_RAY_NL)		
!----------------------------------------------------------------------
!!F90            
!
!!DESCRIPTION:  Subroutine INTANGLE_NL interpolates the lookup
!               Rayleigh reflectances to the measured geometry.
!
!!INPUT PARAMETERS:
!	  INT_NL0	radiance from LUT
!         MTHET0        Solar zenith angle.
!         MTHET        View angle.
!         MDPHI        Azimuth angle.
!
!         THET0      array of  solar Zenith angle in look_up table
!        THE        array of  satllite Zenith angle in look_up table
!         PHI        array of azimuth angle in look_up table
!
!!OUTPUT PARAMETERS 
!         REF_RAY_NL0      Reflectance for rayleigh only
!
 

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04_land.inc'

!     Inputs
      REAL  MTHET0,MTHET,MDPHI
      INTEGER KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1
      
      REAL INT_NL0(NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)

!     Outputs
      REAL REF_RAY_NL(NLWAV)

!     Dummy
      CHARACTER*132 LINE
      CHARACTER*45  LINE1
      INTEGER ICASE,IJ,IPHI,ITABLE,ITAU,ITH,ITH0,IWAV,IETA,LOPT
      REAL  X(2),Y(2),XX1(2),YY1(2),XX2(2),YY2(2),Y1
      INTEGER LL,MM,NN

!      Initialize
 
       KSZAM1=0
       KSZAP1=0
       KTHEM1=0
       KTHEP1=0
       KPHIM1=0
       KPHIP1=0

!      Find THET0 indices
       DO ITH0 = 1,NLTHET0-1
          IF ((MTHET0 .GE. THET0_NL(ITH0)) .AND. &
             (MTHET0 .LE. THET0_NL(ITH0+1))) THEN
	    KSZAM1=ITH0
	    KSZAP1=ITH0+1
          ENDIF
       ENDDO

!      Find THE indices
       DO ITH = 1,NLTHE-1
          IF ((MTHET .GE. THE_NL(ITH)) .AND. &
             (MTHET .LE. THE_NL(ITH+1))) THEN
	    KTHEM1=ITH
	    KTHEP1=ITH+1
          ENDIF
       ENDDO

!      Find PHI indices
       DO IPHI = 1,NLPHI-1
          IF ((MDPHI .GE. PHI_NL(IPHI)) .AND. &
           (MDPHI .LE. PHI_NL(IPHI+1))) THEN
	    KPHIM1=IPHI
	    KPHIP1=IPHI+1
          ENDIF
       ENDDO

 
! Initialize
 
        Y1=0.0
        DO IJ = 1,2
           X(IJ)=0.0
           Y(IJ)=0.0
           XX1(IJ)=0.0
           XX2(IJ)=0.0
           YY1(IJ)=0.0
           YY2(IJ)=0.0
        ENDDO

 
! Interpolate Rayleigh
 


          ITABLE = 1
          ITAU = 1
          DO 15 IWAV=1,NLWAV
              LL=0
              y1 = 0.0
              DO  40 ITH0  = KSZAM1,KSZAP1
                MM=0
                DO  50  ITH  = KTHEM1,KTHEP1
                  NN=0
                  DO 60  IPHI  = KPHIM1,KPHIP1
                    NN=NN+1
                    X(NN)=PHI_NL(IPHI)
                    Y(NN)=INT_NL0(IPHI,ITH,ITH0,ITAU,IWAV,ITABLE)
60                CONTINUE

                  CALL INTERP_EXTRAP(NN,MDPHI,X,Y,Y1,1)
                  MM=MM+1
                  XX1(MM)=THE_NL(ITH)
                  YY1(MM)=Y1

50              CONTINUE
                y1=0.0
                CALL INTERP_EXTRAP(MM,MTHET,XX1,YY1,Y1,1)
                LL=LL+1
                XX2(LL)=THET0_NL(ITH0)
                YY2(LL)=Y1

40            CONTINUE
              y1=0.0
              CALL INTERP_EXTRAP(LL,MTHET0,XX2,YY2,Y1,1)

            REF_RAY_NL(IWAV) = Y1
15        CONTINUE


        RETURN
        END
        subroutine  Read_urban_Table(Average_Urban,handle_Urban_Table_10km)
        USE OCIUAAER_Config_Module
        IMPLICIT NONE
        INCLUDE 'mod04.inc'
        Real Average_Urban(3601,1801)
        integer nfile ,nn,ij,ik,handle_Urban_Table_10km
        CHARACTER  (len=255) :: file_name,Extension
        
         file_name = cfg%VIIRS_land  
         Extension = 'Urban_Table_10km'
         OPEN (handle_Urban_Table_10km, FILE = trim(file_name)//trim(Extension), &
               status = 'old', form = 'unformatted')  
         read(handle_Urban_Table_10km,end=80)Average_Urban 
80       continue  
          Return
          end
        SUBROUTINE compute_urban_Percent(Lat_center,Lon_center,Average_Urban,&
          Urban_per)
          IMPLICIT NONE
        Real Lat_center,Lon_center,Average_Urban(3601,1801)
        Real latmin,latmax,lonmin,lonmax,dx,up,Urban_per
        integer i,j,idxi,idxj 
         dx=0.1
        Urban_per=0.0 
        latmin=-90.00
        latmax=90.00
        lonmin=-180.00
        lonmax=180.00
        if(Lat_center.le.latmax.and.Lat_center.ge.latmin.and.&
        Lon_center.le.lonmax.and.Lon_center.ge.lonmin) then
        idxi=NINT((Lon_center-lonmin)/dx)
        idxj=NINT((Lat_center-latmin)/dx)
        Urban_per=Average_Urban(idxi+1,idxj+1)   
        
        endif 
        return
        end
        





