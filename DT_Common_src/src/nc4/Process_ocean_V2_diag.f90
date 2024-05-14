           SUBROUTINE PROCESS_ocean(HANDLE_S,HANDLE_L,&
              ISCAN,IDATA,MTHET0,MTHET,MPHI,START_500,&
               END_500,START_250,END_250,START_1KM,END_1KM,&
              W659_SYN,W865_SYN,W470_SYN,W550_SYN,W124_SYN,&
              W164_SYN,W213_SYN,W412_SYN,W443_SYN,W8p5_Temp,W1100_Temp,&
              Sunglint_Flag,Set_Counter_Ocean,&
              SDSTAU_best,SDSTAUS_best,SDSTAUB_best,&
              SDSTAU_average,SDSTAUS_average,SDSTAUB_average,&
              SDS_Least_error,SDS_small_weighting,SDS_sol_INDX_small,&
              SDS_sol_INDX_large,SDS_ref,SDS_ref_STD,SDSASSY_best,&
              SDSASSY_average,SDSBACK_best,SDSBACK_average,SDS_effrad,&
              SDS_RefF_best,SDS_ccn,SDS_mass_conc,&
              SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,&
              SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_SCAT_ANGLE_OCEAN,&
              SDS_AOT_model,SDS_CLDFRC_ocean, Set_Counter_Ocean_cloud,&
              QA_Flag_ocean,GLINT_ANGLE,SDS_angs_coeff1,SDS_angs_coeff2,&
               SDS_Tau_Land_Ocean, New_qcontrol_special,&
              SDS_correc_small_weighting,SDS_Tau_Land_Ocean_img,&
               High_Cloud_Flag_500,W138_SYN,data_size,&
               NCLDMSK_SYN,savecldmask,Save_index,Quality_to_pass,&
               Indx_wspped1,Indx_wspped2,WSPEED,wind,AVE_ARRAY,&
               Optical_depth,NUMSQ,HANDLE_Ext_554_O,ipart,SD_for_Dust,&
               Dust_flag_10KM,Ext_554_small,Ext_554_large,WAVE,W354_SYN,W388_SYN,&
               REFW354,REFW388)
                
     
!  ----------------------------------------------------------------------
!  !F90      
!  
!  !DES!  RIPTION:
!  The MODIS o!  ean aerosol produ!  t !  onsists of aerosol opti!  al thi!  kness
!  and size parameters estimates derived on a 10x10 (1-km) pixel spatial
!  array.  The measured radian!  es in a wide spe!  tral range (0.47-2.13
!  mi!  rons) are inverted assuming a bi-lognormal size distribution.
!  The volume and the mean parti!  le size for ea!  h log-normal mode are
!  determined.  When fully developed, the aerosol o!  ean algorithm with
!  use the seven MODIS bands: 1, 2, 3, 4, 5, 6, and 7.
!  
!  !INPUT PARAMETERS:   NONE
!  
!  !OUTPUT PARAMETERS:  NONE
!  SDS_ref_STD     Standard deviation of refle!  tan!  es at 7 bands
!  SDSTAU_best     Opti!  al thi!  kness for best solution
!  SDSTAUS_best    Opti!  al thi!  kness !  ontribution small parti!  les for best solution
!  SDSTAUB_best    Opti!  al thi!  kness !  ontribution large parti!  les for best solution
!  SDSTAU_average  Opti!  al thi!  kness for best solution
!  SDSTAUS_average Opti!  al thi!  kness !  ontribution small parti!  les for best solution
!  SDSTAUB_average Opti!  al thi!  kness !  ontribution large parti!  les for best solution
!  SDS_Least_error         Least square error estimated
!  SDS_small_weighting     Small mode weighting fa!  tor
!  SDS_sol_INDX_small      Solution Number index small parti!  les
!  SDS_sol_INDX_large      Solution Number index large parti!  les
!  SDSASSY_best      Asymmetry_Fa!  tor at 7 bands for best solution
!  SDSASSY_average   Asymmetry_Fa!  tor at 7 bands for average solution
!  SDSBA!  K_best      Ba!  ks!  attering ratio at 7 bands of best solution
!  SDSBA!  K_average   Ba!  ks!  attering ratio at 7 bands of average solution
!  SDS_effrad         Effe!  tive_Radius at 0.55 mi!  ron of both solutions
!  SDS_AOT_model Ratio of opti!  al depth of small mode vs effe!  tive opti!  al depth at 550
!  SDS_RefF_best     Normalized refle!  ted_flux at 7 bands of best solution
!  SDS_RefF_average  Normalized refle!  ted_flux at 7 bands of average solution
!  SDS_TranF_best    Normalized Transmitted_flux at 7 bands of best solution
!  SDS_TranF_average Normalized Transmitted_flux at 7 bands of average solution
!  SDS_S!  AT_ANGLE_O!  EAN S!  attering angle o!  ean
!  SDS_Q!  ONTROL         Quality !  ontrol SDS array
!  SDS_NUMPIXELS        Number of Pixels used for 0.55 mi!  ron
!  SDS_!  !  n              !  loud_Fra!  tion in per!  entage
!  SDS_mass_!  on!        Mass !  on!  entration
!  SDS_angs_!  oeff1      Angstrom Exponent for 0.550 and 0.865 miron
!  SDS_angs_!  oeff2      Angstrom exponent for 0.865 and 2.130 mi!  ron
!  
!  !REVISION HISTORY:
!  $Log: Pro!  ess_o!  ean_V2.f,v $
!  Revision 2.1  1996/05/28  15:34:05  vlin
!  Updated to write out one s!  an!  ube at a time
!  
!  
!  !TEAM-UNIQUE HEADER:
!  
!  !REFEREN!  ES AND !  REDITS:
!  
!  WRITTEN BY:
!  SHANA MATTOO                E-MAIL:mattoo@!  limate.gsf!  .nasa.gov
!  APPLIED RESEAR!  H !  ORPORATION          PHONE:  (301) 614-6214
!  NASA/GODDARD SPA!  E FLIGHT !  ENTER      FAX:    (301) 614-6307
!  !  ODE 913
!  GREENBELT, MD 20771
!  
!  !DESIGN NOTES:
!  
!  Internals Variables:
!  
!    Refle!  tan!  es of Wavelengths used modis,!  h1-!  h7  identified by REFW*
!  
!       REFW659
!       REFW865
!       REFW470
!       REFW550
!       REFW124
!       REFW164
!       REFW213
!      !  LDMSK(4*IGRIDX,4*IGRIDY)    !  loud mask 0=!  loudy,1=!  lear
!      MTHET0     Measured solar zenith angle in deg.
!      MTHET      Measured viewangle from ground in deg.
!      MPHI       Measured azimuth  in deg.
!      NWAV       Number of wavelengths.
!      NTAU       Number of opti!  al thi!  knesses.
!      NTH0       Number of solar zenith angles.
!      NTHET      Number of view angles.
!      NPHI       Number of azimuth angles.
!      PH!  ,JPHI   Azimuth angles.
!      THET       View angles.
!      THET0      Solar zenith angles.
!      TVALUE     Transmission fa!  tor.
!      AINT       Refle!  tan!  e from look-up table.
!      TAUA       opti!  al thi!  knesses.
!      Numdata    Number of input data sets.
!  
!  !END
!  -----------------------w----------------------
       USE Linear_interpolation 
       IMPLICIT  NONE  
       SAVE
       INCLUDE 'mod04.inc'  
       INCLUDE  'read_Sat_MODIS.inc' 
!     HDF array names..........
      BYTE         SDS_QCONTROL_ocean(QA_ocean,NUMCELLS_B)
      REAL SDS_ccn(NUMCELLS_B,NUM_solutions),&
          SDS_mass_conc(NUMCELLS_B,NUM_solutions)
      Real  SDS_ref(NUMCELLS_B,NWAV+2)
      Real  SDS_ref_STD(NUMCELLS_B,NWAV+2)
      Real        SDSTAU_best(NUMCELLS_B,NWAV),&
                  SDSTAU_average(NUMCELLS_B,NWAV),&
                  SDSTAUB_best(NUMCELLS_B,NWAV),&
                  SDSTAUB_average(NUMCELLS_B,NWAV),&
                  SDSTAUS_best(NUMCELLS_B,NWAV),&
                  SDSTAUS_average(NUMCELLS_B,NWAV),&
                  SDSASSY_best(NUMCELLS_B,NWAV),&
                  SDSASSY_average(NUMCELLS_B,NWAV),&
                  SDSBACK_best(NUMCELLS_B,NWAV),&
                  SDSBACK_average(NUMCELLS_B,NWAV),&
                  SDS_RefF_best(NUMCELLS_B,NWAV),&
                  SDS_RefF_average(NUMCELLS_B,NWAV),&
                  SDS_TranF_best(NUMCELLS_B,NWAV),&
                  SDS_TranF_average(NUMCELLS_B,NWAV)
       Real   &
                 SDS_small_weighting(NUMCELLS_B,NUM_solutions),&
                  SDS_correc_small_weighting(NUMCELLS_B),&
                  SDS_Least_error(NUMCELLS_B,NUM_solutions),&
                  SDS_effrad(NUMCELLS_B,NUM_solutions),& 
                  SDS_angs_coeff1(NUMCELLS_B,NUM_solutions),&
                  SDS_angs_coeff2(NUMCELLS_B,NUM_solutions),&
                  SDS_AOT_model(NUMCELLS_B,num_model_index),&
                  SDS_CLDFRC_ocean(NUMCELLS_B)
        Real  SDS_NUMPIXELS(NUMCELLS_B,NWAV+2),&
      SDS_SCAT_ANGLE_OCEAN(NUMCELLS_B),SDS_Tau_Land_Ocean(NUMCELLS_B),&
      SDS_Tau_Land_Ocean_img(NUMCELLS_B)
        Integer     SDS_sol_INDX_small(NUMCELLS_B,NUM_solutions),&
                    SDS_sol_INDX_large(NUMCELLS_B,NUM_solutions)
 
!    Extra computed quantities from Look-up table
!
      Real Wind(Lut_indx)
      REAL RGSS(NUMCASES),SIGMAS(NUMCASES),RGSB(NUMCASEB)
      REAL SIGMAB(NUMCASEB)
      REAL EXTSMALL(Lut_indx,Numcases,NWAV),EXTbig(Lut_indx,Numcaseb,NWAV)
      REAL MOMENTSSMALL(Lut_indx,numcases,4),MOMENTSBIG(Lut_indx,NUMCASEB,4)
      REAL CCNSMALL(Lut_indx,nUMCASES)
      REAL EXTNORSMALL(NUMCASES,NWAV),EXTNORBIG(NUMCASEB,NWAV)
      REAL BACKSCTTSMALL(Lut_indx,NUMCASES,NWAV),BACKSCTTBIG(Lut_indx,NUMCASEB,NWAV)
      REAL ASSYMSMALL(Lut_indx,NUMCASES,NWAV),ASSYMBIG(Lut_indx,NUMCASEB,NWAV)
      REAL NSMALL,NBIG,CCN,ASSYM(NWAV),BACKSCATT(NWAV)
      REAL TAUSMALL(NWAV),TAUBIG(NWAV)
      REAL ALBEDOSMALL(Lut_indx,NUMCASES,NWAV),ALBEDOBIG(Lut_indx,NUMCASEB,NWAV)
      REAL EFFRADIUS,EFFVARIANCE,REF_FLUX(NWAV),TRANS_FLUX(NWAV)
      Integer Dust_flag_10KM(ISWATH_B),Start_coarse,End_coarse
      INTEGER IDATA,ISMALL,IBIG
      INTEGER NUMDATA,NSOLUTION,ipart
       REAL FUNMIN,CLD_FRAC,MASS_CON_OCEAN
      INTEGER START_500,END_500,START_250,END_250,START_1KM,END_1KM
      INTEGER ISCAN,NUMDATA_550,Qcontrol_special
       REAL sd_W659,sd_W865, sd_W470,sd_W550,sd_W124,sd_W164,&
           sd_W213,sd_w354,sd_w388,sd_w0p75,sd_W869o 
!
      INTEGER NUM_ARRAY_ELEMENTS,data_size(2)
      PARAMETER(NUM_ARRAY_ELEMENTS=77)
      REAL MTHET0,MTHET,MPHI
      REAL ARRAY(NUM_ARRAY_ELEMENTS,NUMCASES*NUMCASEB)
      REAL NEW_ARRAY(NUM_ARRAY_ELEMENTS,NUMCASEB*NUMCASES)
      REAL AVE_ARRAY(NUM_ARRAY_ELEMENTS,NUM_solutions)
      real SCAT_ANGLE_OCEAN
      REAL  PHC(NPHI),THET(NTHET),THET0(NTH0)
      REAL  AINTS(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASES)
      REAL  AINTB(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASEB)
      REAL  TAUAS(NUMCASES,NWAV,NTAU)
      REAL  TAUAB(NUMCASEB,NWAV,NTAU),WAVE(NWAV)
      REAL  Ref_ray(NWAV)
      REAL ALBEDO_R_SMALL(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
      REAL ALBEDO_R_BIG(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
      REAL ALBEDO_T_BIG(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)

      REAL RGSS_nc4(NUMCASES),SIGMAS_nc4(NUMCASES),RGSB_nc4(NUMCASEB)
      REAL SIGMAB_nc4(NUMCASEB)
      REAL EXTNORSMALL_nc4(NUMCASES,NWAV),EXTNORBIG_nc4(NUMCASEB,NWAV)
      REAL  AINTS_nc4(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASES)
      REAL  AINTB_nc4(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASEB)
      REAL ref_rayall_nc4(Lut_indx,NPHI,NTHET,NTH0,NWAV)
      REAL  TAUAS_nc4(NUMCASES,NWAV,NTAU)
      REAL  TAUAB_nc4(NUMCASEB,NWAV,NTAU)
      REAL ALBEDO_R_SMALL_nc4(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
      REAL ALBEDO_R_BIG_nc4(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL_nc4(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
      REAL ALBEDO_T_BIG_nc4(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)
      REAL EXTSMALL_nc4(Lut_indx,Numcases,NWAV),EXTbig_nc4(Lut_indx,Numcaseb,NWAV)
      REAL MOMENTSSMALL_nc4(Lut_indx,numcases,4),MOMENTSBIG_nc4(Lut_indx,NUMCASEB,4)
      REAL CCNSMALL_nc4(Lut_indx,nUMCASES)
      REAL BACKSCTTSMALL_nc4(Lut_indx,NUMCASES,NWAV),BACKSCTTBIG_nc4(Lut_indx,NUMCASEB,NWAV)
      REAL ASSYMSMALL_nc4(Lut_indx,NUMCASES,NWAV),ASSYMBIG_nc4(Lut_indx,NUMCASEB,NWAV)
      REAL ALBEDOSMALL_nc4(Lut_indx,NUMCASES,NWAV),ALBEDOBIG_nc4(Lut_indx,NUMCASEB,NWAV)
      REAL Ext_554_small_nc4(Lut_indx,NUMCASES),Ext_554_large_nc4(Lut_indx,NUMCASEB)

      REAL ALBEDO_R_SMALL_tau(NTAU,NWAV,NUMCASES)
      REAL ALBEDO_R_BIG_tau(NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL_tau(NTAU,NWAV,NUMCASES)
      REAL ALBEDO_T_BIG_tau(NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_R_SMALL_final(NWAV,NUMCASES)
      REAL ALBEDO_R_BIG_final(NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL_final(NWAV,NUMCASES)
      REAL ALBEDO_T_BIG_final(NWAV,NUMCASEB)
      INTEGER JPHI(NPHI),ITAU,IWAV,TOT_cldfree_number,Quality_to_pass(2)
      REAL REFSMALL(NUMCASES,NWAV,NTAU),REFBIG(NUMCASEB,NWAV,NTAU)
      REAL ref_rayall(Lut_indx,NPHI,NTHET,NTH0,NWAV)
      INTEGER KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1
      REAL W659_SYN(ISWATH_B,ILINE),W865_SYN(ISWATH_B,ILINE),&
           W470_SYN(ISWATH_B,ILINE),W550_SYN(ISWATH_B,ILINE),&
           W124_SYN(ISWATH_B,ILINE),W164_SYN(ISWATH_B,ILINE),&
           W213_SYN(ISWATH_B,ILINE),&
           W412_SYN(ISWATH_B,ILINE),W443_SYN(ISWATH_B,ILINE),&
           W8p5_Temp(ISWATH_B,ILINE),W1100_Temp(ISWATH_B,ILINE),&
           SD_for_Dust(ISWATH_B,ILINE),&
           W354_SYN(ISWATH_B,ILINE),W388_SYN(ISWATH_B,ILINE)
       INTEGER  Set_Counter_Ocean
       INTEGER Set_Counter_Ocean_cloud,Total_coarse_used
       INTEGER NCLDMSK_SYN(ISWATH_B*2,ILINE*2)
       INTEGER Sunglint_Flag(ISWATH_B,ILINE),QA_Flag_Ocean(12)
       real W138_SYN(ISWATH_B,ILINE),WSPEED
       INTEGER   High_Cloud_Flag_500(ISWATH_B,ILINE),index_wspeed
        

       INTEGER maskoption,num_resol,savecldmask(ISWATH_B*2,ILINE*2)
       INTEGER Save_index(NUMCELLS_B,MaxPixels_left),iist
       Integer Indx_wspped1,Indx_wspped2, New_qcontrol_special
       Real   GLINT_ANGLE,Optical_depth
       integer IYY,IXX,Iblue,JBlue,Iy,IX,HANDLE_Ext_554_O, errid
       character(len=10):: Sat_Flag
       Real Ext_554_small(Lut_indx,NUMCASES),Ext_554_large(Lut_indx,NUMCASEB)   
       SDS_CLDFRC_ocean (IDATA)=-99
! start with bad QCONTROL
!
!   First time in scan cube reads the lookup table.
           New_qcontrol_special =0
           Optical_depth = -99999
           
            
       If (Set_Counter_Ocean .EQ.1) then
            NUMDATA=NUMSQ
!          Read look-up table
         

       CALL READ_LOOK_NC4(RGSS_nc4,SIGMAS_nc4,EXTSMALL_nc4,MOMENTSSMALL_nc4,&
                  CCNSMALL_nc4,EXTNORSMALL_nc4,BACKSCTTSMALL_nc4,ASSYMSMALL_nc4,&
                  RGSB_nc4,SIGMAB_nc4,EXTBIG_nc4,MOMENTSBIG_nc4,EXTNORBIG_nc4,BACKSCTTBIG_nc4,&
                  ASSYMBIG_nc4,ALBEDOSMALL_nc4,ALBEDOBIG_nc4,&
                  ALBEDO_R_SMALL_nc4,ALBEDO_R_BIG_nc4,ALBEDO_T_SMALL_nc4,ALBEDO_T_BIG_nc4,&
                  PHC,THET,THET0,AINTS_nc4,TAUAS_nc4,WAVE,&
                  AINTB_nc4,TAUAB_nc4,JPHI,ref_rayall_nc4,HANDLE_S,HANDLE_L,&
                  HANDLE_Ext_554_O,Ext_554_small_nc4,Ext_554_large_nc4)

       CALL READ_LOOK(RGSS,SIGMAS,EXTSMALL,MOMENTSSMALL,&
                  CCNSMALL,EXTNORSMALL,BACKSCTTSMALL,ASSYMSMALL,&
                  RGSB,SIGMAB,EXTBIG,MOMENTSBIG,EXTNORBIG,BACKSCTTBIG,&
                  ASSYMBIG,ALBEDOSMALL,ALBEDOBIG,&
                  ALBEDO_R_SMALL,ALBEDO_R_BIG,ALBEDO_T_SMALL,ALBEDO_T_BIG,&
                  PHC,THET,THET0,AINTS,TAUAS,WAVE,&
                  AINTB,TAUAB,JPHI,ref_rayall,HANDLE_S,HANDLE_L,&
                  HANDLE_Ext_554_O,Ext_554_small,Ext_554_large)

        errid = 0

        if (1 == 1) then
         errid = 0
         endif
        if (sum(abs(aints-aints_nc4)) > 0) then
         errid = 1
         endif
        if (sum(abs(ref_rayall-ref_rayall_nc4)) > 0) then
         errid = 2
         endif
        if (sum(abs(ALBEDO_R_SMALL-ALBEDO_R_SMALL_nc4)) > 0) then
         errid = 3
         endif
        if (sum(abs(ALBEDO_R_BIG-ALBEDO_R_BIG_nc4)) > 0) then
         errid = 4
         endif
        if (sum(abs(ALBEDO_T_SMALL-ALBEDO_T_SMALL_nc4)) > 0) then
         errid = 5
         endif
        if (sum(abs(ALBEDO_T_BIG-ALBEDO_T_BIG_nc4)) > 0) then
         errid = 6
         endif
        if (sum(abs(TAUAB-TAUAB_nc4)) > 0) then
         errid = 7
         endif
        if (sum(abs(TAUAS-TAUAS_nc4)) > 0) then
         errid = 8
         endif
        if (sum(abs(EXTSMALL-EXTSMALL_nc4)) > 0) then
         errid = 9
         endif
        if (sum(abs(EXTbig-EXTbig_nc4)) > 0) then
         errid = 10
         endif
        if (sum(abs(BACKSCTTSMALL-BACKSCTTSMALL_nc4)) > 0) then
         errid = 11
         endif
        if (sum(abs(BACKSCTTbig-BACKSCTTbig_nc4)) > 0) then
         errid = 12
         endif
        if (sum(abs(ASSYMSMALL-ASSYMSMALL_nc4)) > 0) then
         errid = 13
         endif
        if (sum(abs(ASSYMbig-ASSYMbig_nc4)) > 0) then
         errid = 14
         endif
        if (sum(abs(MOMENTSSMALL-MOMENTSSMALL_nc4)) > 0) then
         errid = 15
         endif
        if (sum(abs(MOMENTSbig-MOMENTSbig_nc4)) > 0) then
         errid = 16
         endif
        if (sum(abs(ALBEDOSMALL-ALBEDOSMALL_nc4)) > 0) then
         errid = 17
         endif
        if (sum(abs(ALBEDObig-ALBEDObig_nc4)) > 0) then
         errid = 18
         endif
        if (sum(abs(Ext_554_small-Ext_554_small_nc4)) > 0) then
         errid = 19
         endif
        if (sum(abs(Ext_554_large-Ext_554_large_nc4)) > 0) then
         errid = 20
         endif
        if (sum(abs(RGSS-RGSS_nc4)) > 0) then
         errid = 21
         endif
        if (sum(abs(SIGMAS-SIGMAS_nc4)) > 0) then
         errid = 22
         endif
        if (sum(abs(RGSB-RGSB_nc4)) > 0) then
         errid = 23
         endif
        if (sum(abs(SIGMAB-SIGMAB_nc4)) > 0) then
         errid = 24
         endif
        if (sum(abs(EXTNORSMALL-EXTNORSMALL_nc4)) > 0) then
         errid = 25
         endif
        if (sum(abs(EXTNORBIG-EXTNORBIG_nc4)) > 0) then
         errid = 26
         endif
        ENDIF

        print *, "error id: ", errid

!
!   Followig If statements checks if measured modis angles are
!   out of bounds from lookup table.

          
           
         
         
!
         IF(MTHET0 .GE. MINMTHET0 .AND. MTHET0 .LE. MAXMTHET0.AND.&
           MTHET   .GE. MINMTHET  .AND. MTHET  .LE. MAXMTHET .AND.&
           MPHI    .GE. MINMPHI   .AND. MPHI   .LE. MAXMPHI)  THEN
             
           
!   
! Call SET_index_inter to set the indexes for measured geometry to use for
!  interpolation
!
        CALL SET_index_inter(MTHET0,MTHET,MPHI,THET0,&
                      KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1)
!
! Call INTANGLE_ray to interpolate the rayleigh Reflectance for measured geo.
!
         
                   CALL INTANGLE_ray(PHC,THET,THET0,ref_rayall,& 
                  MTHET0,MTHET,MPHI,Ref_ray,KSZAM1,KSZAP1,KTHEM1,KTHEP1,& 
                  KPHIM1,KPHIP1,Indx_wspped1,Indx_wspped2,WSPEED,Wind)
     
 
          
         
!
!                **** subroutine AVERAGE averages Reflectances for 10*10 box
!
          CALL AVERAGE(W659_SYN,W865_SYN,W470_SYN,W550_SYN,& 
                  W124_SYN,W164_SYN,W213_SYN,W412_SYN,W443_SYN,& 
                  Sunglint_Flag,NCLDMSK_syn,& 
                  START_1KM,END_1KM,NUMDATA_550,sd_W659,sd_W865,sd_W470,sd_W550,& 
                  sd_W124,sd_W164,sd_W213,sd_w354,sd_w388,sd_w0p75,sd_W869o,& 
                  Ref_ray,MTHET0,WAVE,TOT_cldfree_number,GLINT_ANGLE,iscan,idata,& 
                  Qcontrol_special,High_Cloud_Flag_500,W138_SYN,CLD_FRAC,savecldmask,& 
                  save_index,W8p5_Temp,W1100_Temp,SD_for_Dust,Dust_flag_10KM,&
                  W354_SYN,W388_SYN,REFW354,REFW388)
        
       
        
!       if(CLD_FRAC .GE.0)SDS_CLDFRC_ocean (IDATA)=CLD_FRAC*100.
! Change from percent to fraction for output

          if(CLD_FRAC .GE.0)SDS_CLDFRC_ocean (IDATA)=CLD_FRAC 
        
!   Following If statement checks for missing data in a box 10*10
!  in any one of wavelengths.
!  If there is missing data in any one wavelength
!arrays are filled with _fillValues.
            New_qcontrol_special = Qcontrol_special 
          
         IF(REFW550.LT.99999 .AND. REFW659 .LT. 99999 .AND.&
                    REFW865.LT.99999  .AND. REFW213 .LT. 99999 &
             .and. (Qcontrol_special .eq. 0 .or. Qcontrol_special .eq. 2))then  
             
               
          
    
! Call INTsolarzenith_albedo_tran interpolates to measured solar Zenith angle
!
 
              CALL INTsolarzenith_albedo_tran(THET0,ALBEDO_R_SMALL,& 
                           ALBEDO_R_BIG,ALBEDO_T_SMALL,ALBEDO_T_BIG,MTHET0,& 
                           ALBEDO_R_SMALL_tau,ALBEDO_R_BIG_tau,ALBEDO_T_SMALL_tau,& 
                           ALBEDO_T_BIG_tau,KSZAM1,KSZAP1, Indx_wspped1,Indx_wspped2,& 
                            WSPEED,Wind) 
     
! Call to subroutine INTANGLE interpolates the lookup reflectances to the
! measured geometry.
!
         
          CALL INTANGLE(PHC,THET,THET0,AINTS,TAUAS,& 
                  AINTB,TAUAB,MTHET0,MTHET,MPHI,REFSMALL,REFBIG,& 
                  Ref_ray,KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1,& 
                  Indx_wspped1,Indx_wspped2,WSPEED,Wind) 
        
 
!
! Loop around the cases of small and large size-distribution.
!
         NSOLUTION=0
         DO 199 ISMALL = 1,NUMCASES
        DO ITAU = 1,NTAU
!          TAUAL(ITAU)=TAUAS(1,Index_wave_550,ITAU)
            TAUAL(ITAU)= Tau_at_55_O(Itau)
         DO IWAV=1,NWAV
            REFSMALLL(IWAV,ITAU)=REFSMALL(ISMALL,IWAV,ITAU)
          
         ENDDO
      ENDDO
      
          
       If(Dust_flag_10KM(Idata) .eq.1 .or. Quality_dust_flag_glint .eq.1)then
        Start_coarse  =   NUMCASEB
        End_coarse    =   NUMCASEB
        Total_coarse_used= 1
        else
        Start_coarse = 1
        End_coarse = NUMCASEB-1 
        Total_coarse_used= NUMCASEB-1
       Endif
    
         DO 201 IBIG   = Start_coarse,End_coarse
            DO ITAU = 1,NTAU
               DO IWAV=1,NWAV
                REFBIGL(IWAV,ITAU)=REFBIG(IBIG,IWAV,ITAU) 
               ENDDO
            ENDDO
            
            
!
!             **** subroutine SUBMIN computes the minimum of the
!                     function derived from the weighting of small and
!                     large modes. It returns optical thicknesses values.

             CALL SUBMIN(ISMALL,IBIG,FUNMIN,Ref_ray)
            
!  Call to subroutine COMPUTE_alltau computes optical thickness for
!  all wavelengths for small and coarse distributions.
!
             CALL  COMPUTE_alltau(EXTSMALL,EXTBIG,EXTNORSMALL,&
                    EXTNORBIG,ISMALL,IBIG,Indx_wspped1,Iscan,idata,&
                    Indx_wspped2,WSPEED,Wind,&
                    Ext_554_small,Ext_554_large)
            
                 
!
! Call INTtau_albedo_tran interpolates to computed tau
!
              CALL INTtau_albedo_tran(ALBEDO_R_SMALL_tau,&
                  ALBEDO_R_BIG_tau,ALBEDO_T_SMALL_tau,ALBEDO_T_BIG_tau,&
                  ALBEDO_R_SMALL_final,ALBEDO_R_BIG_final,&
                  ALBEDO_T_SMALL_final,ALBEDO_T_BIG_final,TAUAS,TAUAB,&
                  ISMALL,IBIG)

!
! Call to subroutine COMPUTE_INVERVAR computes the secondary  varaibles
! from Look-up table and computed variables
!
        CALL COMPUTE_INVERVAR(EXTSMALL,&
            EXTBIG,EXTNORSMALL,EXTNORBIG,CCNSMALL,MOMENTSSMALL,&
            MOMENTSBIG,NSMALL,NBIG,EFFRADIUS,EFFVARIANCE,TAUSMALL,&
            TAUBIG,ISMALL,IBIG,ALBEDOSMALL,ALBEDOBIG,&
            CCN,BACKSCTTSMALL,BACKSCTTBIG,ASSYMSMALL,ASSYMBIG,&
            BACKSCATT,ASSYM,REF_FLUX,TRANS_FLUX,ALBEDO_R_SMALL_final,&
            ALBEDO_R_BIG_final,ALBEDO_T_SMALL_final,ALBEDO_T_BIG_final,&
            MASS_CON_OCEAN,ALBEDO_R_SMALL_tau,MTHET0,ALBEDO_R_big_tau,&
            Indx_wspped1,Indx_wspped2,WSPEED,Wind,Ext_554_small,Ext_554_large,&
            ipart)


               NSOLUTION=NSOLUTION+1
!
! Call to store_other stores the output parameters for HDF output file
!   for each 10*10 boxes on a swath * 10 lines
!
         CALL store_other(REF_FLUX,TRANS_FLUX,FUNMIN,EFFRADIUS,&
                  BACKSCATT, ASSYM,NSOLUTION,CCN,SCAT_ANGLE_OCEAN,&
                  ARRAY,NUM_ARRAY_ELEMENTS,NUMDATA_550,ISMALL,IBIG,&
         sd_W470,sd_W550,sd_W659,sd_w865,sd_W124,sd_W164,Sd_W213,CLD_FRAC,&
         TOT_cldfree_number,MASS_CON_OCEAN,WAVE)
  201 CONTINUE
  199 CONTINUE
 
!    Call to AVERAGE_output finds the best solution and average solution.
!     

           
         CALL AVERAGE_output(ARRAY,NUM_ARRAY_ELEMENTS,AVE_ARRAY,&
                         NEW_ARRAY,Total_coarse_used) 
            
           
           New_qcontrol_special =qcontrol_special 
           
          
 
!  IF TAU55 IS IN GREATER THAN RANGE OF MAXTAU, FILL VALUES
!  Set 
      IF(AVE_ARRAY(4,2) .LE.MAXTAU )THEN
!
!    IF TAU55 IS LESS THAN -0.01 FILL WITH FILL VALUES
!
       IF(AVE_ARRAY(4,2) .GT. -0.01 )THEN
!
!  IF TAU55 IS WITHIN RANGE OF -0.01 AND 0 SET QCONTROL AND WRITE THE
!  VALUES
!
      IF( AVE_ARRAY(4,2) .GT. -0.01 .AND. AVE_ARRAY(4,2) .LT.0 )&
             QCONTROL=9
!     CALL to Fill_QAflag_ocean  uses quality control to fill the
!      quality qontrol array
!
       CALL Fill_QAflag_ocean( QA_Flag_Ocean,SDS_QCONTROL_Ocean,&
                       Idata,Quality_to_pass)

!     CALL to storeref_data   Stores Reflectance data
!     and Geolocation for HDF file.
!
 
        CALL storeref_data(IDATA,SDS_ref,SDS_ref_STD,NUMDATA_550,&
                    sd_W470,sd_W550,sd_W659,sd_w865,sd_W124,sd_W164,Sd_W213,&
                    SDS_CLDFRC_ocean,SDS_NUMPIXELS,iscan,&
                    SDS_SCAT_ANGLE_OCEAN,CLD_FRAC,SCAT_ANGLE_OCEAN,&
                    sd_w354,sd_w388,sd_w0p75)

 
!     Stores all varaibles into HDF files
 
       CALL  set_output(IDATA,AVE_ARRAY,NUM_ARRAY_ELEMENTS,&
                  SDSTAU_best,SDSTAU_average,SDSTAUB_best,&
                  SDSTAUB_average,SDSTAUS_best,SDSTAUS_average,&
                  SDS_small_weighting,SDS_Least_error,&
                  SDSASSY_best,SDSASSY_average,SDSBACK_best,SDSBACK_average,&
                  SDS_effrad,SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,&
                  SDS_TranF_average,SDS_sol_INDX_small,SDS_sol_INDX_large,&
                  SDS_ccn,SDS_mass_conc,SDS_Tau_Land_Ocean_img,&
                  SDS_angs_coeff1,SDS_angs_coeff2,SDS_AOT_model,&
                  SDS_Tau_Land_Ocean,SDS_correc_small_weighting)
   
            Optical_depth =AVE_ARRAY(4,2)  
       
      If( qcontrol_special .eq.2) then 
            
       Call Fill_QAflag_ocean( QA_Flag_Ocean,SDS_QCONTROL_Ocean,&
                       Idata,Quality_to_pass)
          
      CALL FILLVALUE_Ocean(IDATA,SDS_ref,SDS_ref_STD,&
                  SDSTAU_best,SDSTAUS_best,SDSTAUB_best,SDSTAU_average,&
                  SDSTAUS_average,SDSTAUB_average,SDS_Least_error,&
                  SDS_small_weighting,SDS_sol_INDX_small,SDS_sol_INDX_large,&
                  SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,&
                  SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,&
                  SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,&
                  SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,&
                  SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_ocean,&
                  SDS_Tau_Land_Ocean_img,Qcontrol_special,&
                  SDS_correc_small_weighting,Dust_flag_10KM) 
        Endif
 
!  ELSE AND ENDIF FOR IF TAU55 IS LESS THAN -0.01 FILL WITH FILL VALUES
 
          ELSE
        QCONTROL=-28
       Call Fill_QAflag_ocean( QA_Flag_Ocean,SDS_QCONTROL_Ocean,&
                       Idata,Quality_to_pass)
      CALL FILLVALUE_Ocean(IDATA,SDS_ref,SDS_ref_STD,&
                  SDSTAU_best,SDSTAUS_best,SDSTAUB_best,SDSTAU_average,&
                  SDSTAUS_average,SDSTAUB_average,SDS_Least_error,&
                  SDS_small_weighting,SDS_sol_INDX_small,SDS_sol_INDX_large,&
                  SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,&
                  SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,&
                  SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,&
                  SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,&
                  SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_ocean,&
                  SDS_Tau_Land_Ocean_img,Qcontrol_special,&
                  SDS_correc_small_weighting,Dust_flag_10KM)
       ENDIF
 
!  ELSE AND ENDIF IF TAU55 IS IN GREATER THAN RANGE OF MAXTAU, FILL VALUES
 

       ELSE
        QCONTROL=-29
        Call Fill_QAflag_ocean( QA_Flag_Ocean,SDS_QCONTROL_Ocean,&
                       Idata,Quality_to_pass)
         CALL FILLVALUE_Ocean(IDATA,SDS_ref,SDS_ref_STD,&
                  SDSTAU_best,SDSTAUS_best,SDSTAUB_best,SDSTAU_average,&
                  SDSTAUS_average,SDSTAUB_average,SDS_Least_error,&
                  SDS_small_weighting,SDS_sol_INDX_small,SDS_sol_INDX_large,&
                  SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,&
                  SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,&
                  SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,&
                  SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,&
                  SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_ocean,&
                  SDS_Tau_Land_Ocean_img,Qcontrol_special,&
                  SDS_correc_small_weighting,Dust_flag_10KM)
        ENDIF

 
!        This Else and endif statementis are for glint angle
!         arrays are filled with _fillValues nad Reflectance with Reflectance
!        sd and pixels
        Else
          If( qcontrol_special .eq.3) then 
!           write(36,*) 'storeref_data 2',iscan,idata,
!     *      REFW470,REFW550,REFW659,REFW865,REFW124,REFW164,REFW213,
!     *      REFW412,REFW443,REFW0p75 
        CALL storeref_data(IDATA,SDS_ref,SDS_ref_STD,NUMDATA_550,&
                    sd_W470,sd_W550,sd_W659,sd_w865,sd_W124,sd_W164,Sd_W213,&
                    SDS_CLDFRC_ocean,SDS_NUMPIXELS,iscan,&
                    SDS_SCAT_ANGLE_OCEAN,CLD_FRAC,SCAT_ANGLE_OCEAN,&
                    sd_w354,sd_w388,sd_w0p75)

 
         QCONTROL=-21
       Call Fill_QAflag_ocean( QA_Flag_Ocean,SDS_QCONTROL_Ocean,&
                       Idata,Quality_to_pass)
       CALL FILLVALUE_Ocean(IDATA,SDS_ref,SDS_ref_STD,&
                  SDSTAU_best,SDSTAUS_best,SDSTAUB_best,SDSTAU_average,&
                  SDSTAUS_average,SDSTAUB_average,SDS_Least_error,&
                  SDS_small_weighting,SDS_sol_INDX_small,SDS_sol_INDX_large,&
                  SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,&
                  SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,&
                  SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,&
                  SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,&
                  SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_ocean,&
                  SDS_Tau_Land_Ocean_img,Qcontrol_special,&
                  SDS_correc_small_weighting,Dust_flag_10KM)
         ELSE
        Call Fill_QAflag_ocean( QA_Flag_Ocean,SDS_QCONTROL_Ocean,&
                       Idata,Quality_to_pass)
       CALL FILLVALUE_Ocean(IDATA,SDS_ref,SDS_ref_STD,&
                  SDSTAU_best,SDSTAUS_best,SDSTAUB_best,SDSTAU_average,&
                  SDSTAUS_average,SDSTAUB_average,SDS_Least_error,&
                  SDS_small_weighting,SDS_sol_INDX_small,SDS_sol_INDX_large,&
                  SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,&
                  SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,&
                  SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,&
                  SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,&
                  SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_ocean,&
                  SDS_Tau_Land_Ocean_img,Qcontrol_special,&
                  SDS_correc_small_weighting,Dust_flag_10KM)
         Endif
      ENDIF

!                **** Else and endif If angles are out of bounds.
!                     fill the arrays with _FillValue
       ELSE
!  angles are out of bounds
           QCONTROL=-26
       Call Fill_QAflag_ocean( QA_Flag_Ocean,SDS_QCONTROL_Ocean,&
                       Idata,Quality_to_pass)
       CALL FILLVALUE_Ocean(IDATA,SDS_ref,SDS_ref_STD,&
                  SDSTAU_best,SDSTAUS_best,SDSTAUB_best,SDSTAU_average,&
                  SDSTAUS_average,SDSTAUB_average,SDS_Least_error,&
                  SDS_small_weighting,SDS_sol_INDX_small,SDS_sol_INDX_large,&
                  SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,&
                  SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,&
                  SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,&
                  SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,&
                  SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_ocean,&
                  SDS_Tau_Land_Ocean_img,Qcontrol_special,&
                  SDS_correc_small_weighting,Dust_flag_10KM)
      ENDIF 
              
         
        if(New_qcontrol_special .eq.1) then 
            Optical_depth =0.0 
             REFW470=99999
             REFW550=99999
             REFW659=99999
             REFW865=99999
             REFW124=99999
             REFW164=99999
             REFW213=99999 
             REFW0p75=99999
             WSPEED =99999 
             
        Endif 

       RETURN
      END

    

!*********************************************************************

       SUBROUTINE READ_LOOK(RGSS,SIGMAS,EXTSMALL,MOMENTSSMALL,&
          CCNSMALL,EXTNORSMALL,BACKSCTTSMALL,ASSYMSMALL,&
          RGSB,SIGMAB,EXTBIG,MOMENTSBIG,EXTNORBIG,BACKSCTTBIG,&
          ASSYMBIG,ALBEDOSMALL,ALBEDOBIG,&
          ALBEDO_R_SMALL,ALBEDO_R_BIG,ALBEDO_T_SMALL,ALBEDO_T_BIG,&
          PHC,THET,THET0,AINTS,TAUAS,WAVE,&
          AINTB,TAUAB,JPHI,ref_rayall,HANDLE_S,HANDLE_L,&
            HANDLE_Ext_554_O,Ext_554_small,Ext_554_large)  
          
!---------------------------------------------------------------------
!!F90
 
!
!!INPUT PARAMETERS:
!          HANDLE_S        Input read for small particles(lookup).
!          HANDLE_L        Input read for large particles(lookup).
!
!!OUTPUT PARAMETERS:
!      SMALL MODE.........
!
!                  RGSS        RADIUS
!            SIGMAS        SIGMA
! EXTSMALL       EXTINCTION COEFF
!       MOMENTSSMALL       MOMENTS SMALL MODE
! CCNSMALL       CLOUD CONDENSATION NUCLII
!       EXTNORSMALL       EXTINCTION COEFF NORMALIZED
!     BACKSCTTSMALL       BACKSCATTERING RATIO
!         ASSYMSMALL       ASSYMETRY FACTOR
!        ALBEDOSMALL       ALBEDO SINGLE SCATTERING
!     ALBEDO_R_SMALL       REFLECTED ALBEDO
!     ALBEDO_T_SMALL      TRANSMITTED ALBEDO
!     LARGE MODE..........
!
!                   RGSB       RADIUS
!   SIGMAB       SIGMA
!   EXTBIG       EXTINCTION COEFF
!         MOMENTSBIG       MOMENTS SMALL MODE
!          EXTNORBIG       EXTINCTION COEFF NORMALIZED
!        BACKSCTTBIG       BACKSCATTERING RATIO
! ASSYMBIG       ASSYMETRY FACTOR
!          ALBEDOBIG       ALBEDO SINGLE SCATTERING
!       ALBEDO_R_BIG       REFLECTED ALBEDO
!       ALBEDO_T_BIG       TRANSMITTED ALBEDO
!
! PHC       Azimuth angle.
!
!                   THET       view angle.
!
!                  THET0       Solar zenith angle.
!
!                  TFLUX       dowanward flux
!
!                   AINTS       radiance(l/fo),fo=1 small mode
!
!                   AINTB       radiance(l/fo),fo=1 large mode
!
!                   TAUAS       optical thickness SMALL MODE
!
!                   TAUAB       optical thickness LARGE MODE
!
!                   NWAV       number of wavelength
!
!                   NTAU       number of opticl thickness.
!
!                   NTH0       number of solar zenith angle.
!
!                  NTHET       number of view angle.
!
!                   NPHI       number of azimuth
!
!                   JPHI       azimuth
!        ref_rayall      radiance(l/fo),fo=1 FOR RAYLEIGH TAUA=0.0
 
      USE OCIUAAER_Config_Module
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      CHARACTER*132 LINE
      REAL EXTSMALL(Lut_indx,Numcases,NWAV),EXTbig(Lut_indx,Numcaseb,NWAV)
      REAL RGSS(NUMCASES),SIGMAS(NUMCASES)
      REAL RGSB(NUMCASEB),SIGMAB(NUMCASEB)
      REAL MOMENTSSMALL(Lut_indx,numcases,4),MOMENTSBIG(Lut_indx,NUMCASEB,4)
      REAL CCNSMALL(Lut_indx,nUMCASES),TR
      REAL EXTNORSMALL(NUMCASES,NWAV),EXTNORBIG(NUMCASEB,NWAV)
      REAL BACKSCTTSMALL(Lut_indx,NUMCASES,NWAV),BACKSCTTBIG(Lut_indx,NUMCASEB,NWAV)
      REAL ASSYMSMALL(Lut_indx,NUMCASES,NWAV),ASSYMBIG(Lut_indx,NUMCASEB,NWAV)
      REAL ALBEDOSMALL(Lut_indx,NUMCASES,NWAV),ALBEDOBIG(Lut_indx,NUMCASEB,NWAV)
      REAL ALBEDO_R_SMALL(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
      REAL ALBEDO_R_BIG(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
      REAL ALBEDO_T_BIG(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_R_RAY(Lut_indx,NWAV,NTH0)
      REAL ALBEDO_T_RAY(Lut_indx,NWAV,NTH0)
      INTEGER ICASE,IFILE,IPHI,ITAU,ITH,ITH0,IWAV,IJ
      INTEGER JPHI(NPHI),Num_lut,HANDLE_Ext_554_O,Nfile
      REAL  PHC(NPHI),THET(NTHET),THET0(NTH0),WAVE(NWAV)
      REAL  AINTS(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASES)
      REAL  TAUAS(NUMCASES,NWAV,NTAU)
      REAL  AINTB(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASEB)
      REAL  TAUAB(NUMCASEB,NWAV,NTAU),DUMMY(NTH0)
      REAL  CCNdummy,ref_rayall(Lut_indx,NPHI,NTHET,NTH0,NWAV)
      REAL EFFRADSMALL(NUMCASES),EFFRADBIG(NUMCASEB),QSCT
      Real Ext_554_small(Lut_indx,NUMCASES),Ext_554_large(Lut_indx,NUMCASEB)
      character * 1 cWS
      character ( len=10):: Sat_Flag 
      CHARACTER  (len=255) :: tbl,file_name,Extension,File  
           
      
         
!     we     are opening files to read extinction coeff for wavlength 0.544 for normalization.  
            file_name = cfg%dt_ext_coeff    
            OPEN (HANDLE_Ext_554_O, FILE = trim(file_name) ,status = 'old')  
             Nfile = HANDLE_Ext_554_O
              
            Read(Nfile,115) Line   
            Do Num_lut = 1,lut_indx
            Read(Nfile,*)(Ext_554_small(Num_lut,Icase),Icase = 1,NUMCASES)  
            enddo
            Do Num_lut = 1,lut_indx
            Read(Nfile,*)(Ext_554_large(Num_lut,Icase),Icase = 1,NUMCASEB) 
            ENDDO
115         format(132A1)  

  
             file_name = cfg%VIIRS_ocean 
             Do Num_lut = 1,lut_indx 
              WRITE(cWS, '(i1)' )Num_lut
              Extension = '/small_v'//cWS//'c1.dat.npp4' 
! read rayleigh data 
              IFILE = HANDLE_S(Num_lut) 
           OPEN (IFILE, FILE = trim(file_name)//trim(Extension),status='old')  
            
        
         
         DO 200 IWAV=1,NWAV
      
 
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
 
! Read all other information
        DO IJ = 1,2
          READ(IFILE,500)LINE
        ENDDO

           READ(IFILE,525)(THET0(IJ),IJ=1,NTH0)

 
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(DUMMY(IJ),IJ=1,NTH0)
                  READ(IFILE,530)(ALBEDO_R_RAY(Num_lut,IWAV,IJ),IJ=1,NTH0)
                  READ(IFILE,530)(ALBEDO_T_RAY(Num_lut,IWAV,IJ),IJ=1,NTH0)
 
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
 
             DO  400 IPHI=1,NPHI
               READ(IFILE,540)JPHI(IPHI),&
              (ref_rayall(Num_lut,Iphi,ITH,ITH0,IWAV),ITH=1,NTHET)
  400       CONTINUE
!    enddo for ith0
 300       CONTINUE
!    enddo for wav
 200       CONTINUE

! LOOP IS AROUND THET0,NEW FILE FOR EACH SMALL CASE
         DO 10 ICASE =1,NUMCASES
         DO 20 IWAV=1,NWAV
! Leave itau=1 to fill up with taua=0.0 later in subroutine INTANGLE
         DO 30 ITAU = 2,NTAU
! Read aerosol information

                  read(IFILE,500)LINE
                  read(IFILE,505)WAVE(IWAV)
        DO IJ = 1,3
          READ(IFILE,500)LINE
        ENDDO

       READ(IFILE,515)RGSS(ICASE),SIGMAS(ICASE)
       READ(IFILE,500)LINE
       READ(IFILE,515)EFFRADSMALL(ICASE)
       READ(IFILE,520)MOMENTSSMALL(Num_lut,ICASE,1),MOMENTSSMALL(Num_lut,ICASE,2)
       READ(IFILE,520)MOMENTSSMALL(Num_lut,ICASE,3),MOMENTSSMALL(Num_lut,ICASE,4)
       READ(IFILE,515)ALBEDOSMALL(Num_lut,ICASE,IWAV),ASSYMSMALL(Num_lut,ICASE,IWAV)
       READ(IFILE,515)CCNSMALL(Num_lut,ICASE),BACKSCTTSMALL(Num_lut,ICASE,IWAV)
       READ(IFILE,520)QSCT,EXTSMALL(Num_lut,ICASE,IWAV)

! Read ocean information
         DO IJ =1,4
          READ(IFILE,500)LINE
        ENDDO
! Read Atmosphere information

                  READ(IFILE,500)LINE
                  READ(IFILE,515)TR,TAUAS(ICASE,IWAV,ITAU)
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
        READ(IFILE,530)(ALBEDO_R_SMALL(Num_lut,IJ,ITAU,IWAV,ICASE),IJ=1,NTH0)
        READ(IFILE,530)(ALBEDO_T_SMALL(Num_lut,IJ,ITAU,IWAV,ICASE),IJ=1,NTH0)
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
              (AINTS(Num_lut,IPHI,ITH,ITH0,ITAU,IWAV,ICASE),ITH=1,NTHET)
   50       continue
!    enddo for ith0
   40       CONTINUE
!    enddo for tau
   30       CONTINUE
! enddo for nwav
! Fill the array for albedo and transmission for all cases tau=0
      TAUAS(ICASE,IWAV,1)=TAUAS(1,IWAV,1)
      DO IJ=1,NTH0
      ALBEDO_R_SMALL(Num_lut,IJ,1,IWAV,ICASE)=ALBEDO_R_RAY(Num_lut,IWAV,IJ)
      ALBEDO_T_SMALL(Num_lut,IJ,1,IWAV,ICASE)=ALBEDO_T_RAY(Num_lut,IWAV,IJ)
      ENDDO

   20       CONTINUE
!    enddo for num of size distribution for small
   10       CONTINUE
                  DO IPHI=1,NPHI
                     PHC(IPHI)=FLOAT(JPHI(IPHI))
                  ENDDO
! Endo for windspeed lookup table                  
            close (HANDLE_S(Num_lut))   
      Enddo            

! READ LARGE CASES OF SIZEDISTRIBUTION

        Do Num_lut = 1,lut_indx
             WRITE(cWS, '(i1)' )Num_lut
             Extension = '/big_v'//cWS//'c1.dat.npp4'   
             IFILE = HANDLE_L(Num_lut)
         OPEN (IFILE, FILE = trim(file_name)//trim(Extension),status='old')  
           
         DO 60 ICASE =1,NUMCASEB
         DO 70 IWAV=1,NWAV
! Leave itau=1 to fill up with taua=0.0 later in subroutine INTANGLE
         DO 80 ITAU = 2,NTAU
! Read aerosol information

                  read(IFILE,500)LINE
                  read(IFILE,505)WAVE(IWAV)
        DO IJ = 1,3
          READ(IFILE,500)LINE
        ENDDO

       READ(IFILE,515)RGSB(ICASE),SIGMAB(ICASE)
       READ(IFILE,500)LINE
       READ(IFILE,515)EFFRADBIG(ICASE)
       READ(IFILE,520)MOMENTSBIG(Num_lut,ICASE,1),MOMENTSBIG(Num_lut,ICASE,2)
       READ(IFILE,520)MOMENTSBIG(Num_lut,ICASE,3),MOMENTSBIG(Num_lut,ICASE,4)
       READ(IFILE,515)ALBEDOBIG(Num_lut,ICASE,IWAV),ASSYMBIG(Num_lut,ICASE,IWAV)
       READ(IFILE,516)CCNdummy,BACKSCTTBIG(num_lut,ICASE,IWAV)
       READ(IFILE,520)QSCT,EXTBIG(Num_lut,ICASE,IWAV)

! Read ocean information
         DO IJ =1,4
          READ(IFILE,500)LINE
        ENDDO
! Read Atmosphere information

                  READ(IFILE,500)LINE
                  READ(IFILE,515)TR,TAUAB(ICASE,IWAV,ITAU)
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
       READ(IFILE,530)(ALBEDO_R_BIG(Num_lut,IJ,ITAU,IWAV,ICASE),IJ=1,NTH0)
       READ(IFILE,530)(ALBEDO_T_BIG(Num_lut,IJ,ITAU,IWAV,ICASE),IJ=1,NTH0)
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
           (AINTB(Num_lut,IPHI,ITH,ITH0,ITAU,IWAV,ICASE),ITH=1,NTHET)
 100       continue
!    enddo for ith0
  90      CONTINUE
!    enddo for tau
  80       CONTINUE
!    enddo for wav
! Fill the array for albedo and transmission for all cases tau=0
      TAUAB(ICASE,IWAV,1)=TAUAS(1,IWAV,1)
      DO IJ=1,NTH0
      ALBEDO_R_BIG(Num_lut,IJ,1,IWAV,ICASE)=ALBEDO_R_RAY(Num_lut,IWAV,IJ)
      ALBEDO_T_BIG(Num_lut,IJ,1,IWAV,ICASE)=ALBEDO_T_RAY(Num_lut,IWAV,IJ)
      ENDDO
  70      CONTINUE
!    enddo for num of size distribution for large
  60       CONTINUE
! Enddo for windspeed  
           close (HANDLE_L(Num_lut))
           Enddo
           close (HANDLE_Ext_554_O)
!
!     format statements
!
500               format(132a1)
505               format(t32,f6.4)
510               format(t32,f6.4,t65,e11.4)
515               format(t32,f6.4,t70,f6.4)
516               format(t30,f6.4,t68,f8.4)
520               format(t26,e12.4,t64,e12.4)
525               format(t12,f6.1,3x,6(2x,f5.1,3x)/&
                   t12,f6.1,3x,6(2x,f5.1,3x))
530               format(t10,7e10.3/t10,7e10.3)
535               format(t10,f4.1,6(5x,f5.1)/t9,7(f5.1,5x)/&
                         t9,7(f5.1,5x))
540               format(i4,1x,7e10.3/5x,7e10.3/5x,7e11.3)

      RETURN
      END





!~*********************************************************************

        SUBROUTINE AVERAGE(W659_SYN,W865_SYN,W470_SYN,W550_SYN,&
        W124_SYN,W164_SYN,W213_SYN,W412_SYN,W443_SYN,&
        Sunglint_Flag,NCLDMSK_syn,&
        START_1KM,END_1KM,NUMDATA_550,sd_W659,sd_W865,sd_W470,sd_W550,&
       sd_W124,sd_W164,sd_W213,sd_w354,sd_w388,sd_w0p75,sd_W869o,&
       Ref_ray,MTHET0,WAVE,TOT_cldfree_number,GLINT_ANGLE,iscan,idata,&
       Qcontrol_special, High_Cloud_Flag_500,W138_SYN,CLD_FRAC,&
       savecldmask,save_index,W8p5_Temp,W1100_Temp,SD_for_Dust,&
       Dust_flag_10KM,W354_SYN,W388_SYN,REFW354,REFW388)
        

!-----------------------------------------------------------------------
!!F90
!
!DESCRIPTION: This subroutine averages the reflectances for each
!              10*10 pixel  square and finds standdard deviation for
!               averaged Reflectances
!
!!INPUT PARAMETERS:
!         IDATA        Index of bin number in a swath
!   Reflectances of wavelengths used modis,ch1-ch7  identified by W*_SYN
!   Extra wavelengths from ocean color ch9,ch12,ch13 and ch16 will be
!   used for substituting if Ch1 - ch7 do not have valid data.
!    These channels are identified by W*o_SYN.
!
!       W659_SYN
!       W865_SYN
!       W470_SYN
!       W550_SYN
!       W124_SYN
!       W164_SYN
!       W213_SYN
!       W443o_SYN
!       W551o_SYN
!       W667o_SYN
!       W869o_SYN
!       CLDMSK_syn   40*40 array reflectance for cloud mask.
!        Sunglint_Flag  Sunglint Flag from cloud mask
!        START_500     Starting Index for 500 meter resolution
!        END_500       Ending  Index for 500 meter resolution
!        START_250     Starting Index for 250 meter resolution
!        END_250       ending Index for 500 meter resolution
!        START_1KM     Starting Index for 1km resolution
!        END_1KM       ending Index for 1km resolution
!        QCONTROL      quality control index
!        Ref_ray       Rayleigh Reflectance
!        MTHET0        Measured Soalr zenith angle
!
!  Input from the subroutine AVERAGE_to_500meter
!
!       Reflectances averaged for 500 meter resolution
!       wavelengths used modis,ch1-ch7  identified by ref_interm_w*
!    Extra wavelengths from ocean color ch9,ch12,ch13 and ch16 will be
!    used for substituting if Ch1 - ch7 do not have valid data.
!    These channels are identified by ref_interm_w*o
!       ref_interm_w659
!       ref_interm_w865
!       ref_interm_w470
!       ref_interm_w550
!       ref_interm_w124
!       ref_interm_w164
!        ref_interm_w213
!        ref_interm_w443o
!        ref_interm_w551o
!        ref_interm_w667o
!        ref_interm_w869o
!
!        NCLDMSK_syn
! Number of data points in each 500 meter box for all wavelengths used
!        NUMDATA_659
!        NUMDATA_865
!        NUMDATA_470
!        NUMDATA_550
!        NUMDATA_124
!        NUMDATA_164
!        NUMDATA_213
!        NUMDATA_443o
!        NUMDATA_551o
!        NUMDATA_667o
!        NUMDATA_869o
!
!       Total number of valid data  500 meter resolution
!       wavelengths used modis,ch1-ch7  identified by VFlag_w*
!    Extra wavelengths from ocean color ch9,ch12,ch13 and ch16 will be
!   used for substituting if Ch1 - ch7 do not have valid data.
!    These channels are identified by VFlag_w*o
!       VFlag_w659
!       VFlag_w865
!       VFlag_w470
!       VFlag_w550
!       VFlag_w124
!       VFlag_w164
!       VFlag_w213
!       VFlag_w443o
!       VFlag_w551o
!       VFlag_w667o
!       VFlag_w869o
!!OUTPUT PARAMETERS:
!   Reflectances of Wavelengths used modis,ch1-ch7  identified by REFW*
!
!       REFW659
!       REFW865
!       REFW470
!       REFW550
!       REFW124
!       REFW164
!       REFW213
!   STANDARD DEVIATION of Wavelengths used modis,ch1-ch7  identified by sd_w*
!
!       SD_W659
!       SD_W865
!       SD_W470
!       SD_W550
!       SD_W124
!       SD_W164
!       SD_W213
!       NUMDATA_550  Number of pixels used for processing
 
      IMPLICIT NONE
       save
      INCLUDE 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'
  
      INTEGER START_1KM,END_1KM
      INTEGER IX,IY,WATER,CLDFREE
      INTEGER NUMDATA_659,NUMDATA_865,NUMDATA_470,&
       NUMDATA_550,NUMDATA_124,NUMDATA_164,NUMDATA_213,NUMDATA_354,&
       NUMDATA_388
      INTEGER NOGLINT,TOT_cldfree_number,cloudy
      REAL Ref_ray(NWAV),MTHET0,CLD_FRAC
      REAL W659_SYN(ISWATH_B,ILINE),W865_SYN(ISWATH_B,ILINE),&
           W470_SYN(ISWATH_B,ILINE),W550_SYN(ISWATH_B,ILINE),&
           W124_SYN(ISWATH_B,ILINE),W164_SYN(ISWATH_B,ILINE),&
           W213_SYN(ISWATH_B,ILINE),W138_SYN(ISWATH_B,ILINE),&
           W412_SYN(ISWATH_B,ILINE),W443_SYN(ISWATH_B,ILINE),&
           W8p5_Temp(ISWATH_B,ILINE),W1100_Temp(ISWATH_B,ILINE),&
           SD_for_Dust(ISWATH_B,ILINE),&
           W354_SYN(ISWATH_B,ILINE),W388_SYN(ISWATH_B,ILINE)
      real array(2*IGRIDX*2*IGRIDY),newarray(2*IGRIDX*2*IGRIDY)
       REAL  ref_interm(NWAV+2,2*IGRIDX*2*IGRIDY)
       REAL  array_interm(NWAV+2,2*IGRIDX*2*IGRIDY)
       REAL sd_W659,sd_W865, sd_W470,sd_W550,sd_W124,sd_W164,&
        sd_W213,sd_w354,sd_w388,sd_w0p75,sd_W869o
      REAL  sd_allwav(NWAV+2),var_allwav(NWAV+2)
      REAL del_var_wav659,del_var_wav164,del_var_waV213,del_var_wav124
      REAL del_var_wav412,del_var_wav443,del_var_wa0P75 
      real WAVE(NWAV),array_new(2*IGRIDX*2*IGRIDY)
      real ref_val,sd_val 
      INTEGER NCLDMSK_SYN(ISWATH_B*2,ILINE*2)
      integer indx(4*IGRIDX*4*IGRIDY)
      INTEGER  inum,index_set,inn,ij,Dust_flag(ISWATH_B,ILINE)
      INTEGER vFlag_w659,vFlag_w865,vFlag_w470,&
       vFlag_w550,vFlag_w124,vFlag_w164,vFlag_w213,vFlag_w354,&
       vFlag_w388,vFlag_w0P75, vFlag_w869o
      INTEGER Sunglint_Flag(ISWATH_B,ILINE)
      INTEGER n10,n40,iscan,idata,iwav,Qcontrol_special
      integer quality_wav(NWAV), Total_good_points,Dust_flag_10KM(ISWATH_B)
      real GLINT_ANGLE
      INTEGER  High_Cloud_Flag_500(ISWATH_B,ILINE)
      Integer savecldmask(ISWATH_B*2,ILINE*2)
      integer Save_index(NUMCELLS_B,MaxPixels_left)
      integer Save_Array_index(MaxPixels_left)
      integer Indx_wspped1,Indx_wspped2
      integer   Number_pixels,Tot_Number_pixels,Good_Number_pixels
      character (len =10) ::  Sat_flag    
        WATER=0
        cldfree=0
        Qcontrol1=0
        Qcontrol_special=0 
       Quality_dust_flag_glint  =0
       Quality_dust_flag_off_glint=0
        DO IWAV=1,NWAV+2
       Good_pixels(IWAV)=0
       ENDDO
 
!     Apply glint mask 
           
         IF(GLINT_ANGLE.GT.30 .AND.GLINT_ANGLE.LE.40 )QCONTROL=10
 
! Call to AVERAGE_to_500meter averages all resolutions to 500 meter resolution
 
       
        CALL   AVERAGE_to_500meter(W659_SYN,W865_SYN,W470_SYN,&   
                W550_SYN,W124_SYN,W164_SYN,W213_SYN,W412_SYN,W443_SYN,& 
                W8p5_Temp,W1100_Temp,ref_interm,NCLDMSK_syn,&
                NUMDATA_659,NUMDATA_865,NUMDATA_470,NUMDATA_550,&
                NUMDATA_124,NUMDATA_164,NUMDATA_213,NUMDATA_354,&
                NUMDATA_388,START_1KM,&
                END_1KM,VFlag_w659,VFlag_w865,VFlag_w470,VFlag_w550,VFlag_w124,&
                VFlag_w164,VFlag_w213,vFlag_w354,vFlag_w388,&
                VFlag_w869o,cloudy,iscan,High_Cloud_Flag_500,W138_SYN,Ref_ray,Wave,& 
                CLD_FRAC,savecldmask,idata,SD_for_Dust,Dust_flag_10KM,&
                W354_SYN,W388_SYN)
        
      
            
! Sort in asending order and get the index for wavelength at 0.865
         do ix=1,NUMDATA_865
         array(ix)=ref_interm(Index_wave_865,ix) 
         enddo

       CALL INDEXX(NUMDATA_865,array,INDX)
           
! Set the arrays to reject 1/4 of all bright and dark pixels based on
! Wavelength of 865 nm

          inum=NUMDATA_865
           
          
! Reject 1/4 of brightest and darkest cloudfree pixels to eliminate
! the residual shadows and sub_cloud pixels
!
       
!!!!!! I* 1 Pixel  Retere.....
        if( Iline .eq.1) then
          n10 =1
          n40 =2 
          Number_pixels = 0 
          Tot_Number_pixels=3
          Good_Number_pixels=0
         endif
         
!          go to 955
          
!!!!!!!!!! 10 * 10 Pixel Ret         
         
          n10=inum/4
          n40=(inum-n10)
!          Number_pixels =10,Tot_Number_pixels =30 Good_Number_pixels = 20 for operational  
!   used for10      
!              Number_pixels      =   3
!             Tot_Number_pixels  =    7
!             Good_Number_pixels =    5 
             if(n10 .eq.0) n10=1
             if(n40.eq.0)  n40 =1
             Number_pixels       = INT((Iline*Iline) * .025) 
             Tot_Number_pixels   = INT((Iline*Iline) * .075) 
             Good_Number_pixels  = INT((Iline*Iline) * .050) 
!   used for 3 and 8 boxes            
!           Number_pixels = 2
!          Tot_Number_pixels=3
!          Good_Number_pixels=3 
          
         
         
         
         
!   5% of total box  
!955     continue
!          
!           
!         if(Iline.gt. 1) then      
!         Number_pixels     = ((Iline*Iline) * .05)+ 2
!         Tot_Number_pixels = (((Iline*Iline) * 5) *.03) + 2
!         Good_Number_pixels =((Iline*Iline) * .1) + 2
!         endif
!            n10=(inum/4)
!           n40=(inum-n10) 
!           if(n10 .eq.0)  n10 =1
!           if(n40 .eq. 0) n40 =1
           
           
           
!         print*,'ocean',Number_pixels,Tot_Number_pixels,Good_Number_pixels,n10,n40,inum    
! THrow away 1/4 of brightest & darkest pixels in all wavelengths
 
!         IF( ( n40-n10) .gt.10) then 
!         IF( ( n40-n10) .gt.5) then
          
          IF( (n40-n10) .gt. Number_pixels) then 
             DO iy=1,NWAV+2
                inn=0
            
            If( Iline .eq.1) then 
                array_new(1)=ref_interm(IY,1) 
               if( Iy .eq.4)Save_Array_index(1)=indx(1)
                inn =1 
            Else
                
                
                DO ix=n10,n40
                     inn=inn+1 
                     array_new(inn)=ref_interm(IY,indx(ix)) 
                if( Iy .eq.4)Save_Array_index(inn)=indx(ix) 
                 ENDDO 
             Endif 
!Throw away bad pixles in all wavelengths
          
           DO IX=1,inn
              if ( array_new(IX) .GT.0. .AND. array_new(IX) .LE.1)then 
                Good_pixels(IY)= Good_pixels(IY)+1
                array_interm(IY,Good_pixels(IY))= array_new(ix)  
                 if( Iy .eq.4)&
             save_index(Idata,Good_pixels(IY))=Save_Array_index(Ix) 
              endif
               
          ENDDO 
!  ENDDO for wavelengths
          ENDDO

      
            DO iy=1,NWAV+2
! Call to subroutine ave_std computes average and standard deviation for
! Reflectances
            
        IF(Good_pixels(IY) .gt.0)then
              IF(iline .eq.1) then 
                  ref_allwav(iy) =array_interm(IY,Good_pixels(IY)) 
                  sd_allwav(iy) =-9999
                  var_allwav(iy)=-9999  
             Else
!  if Iline is greater than 1             
              DO IX=1,Good_pixels(IY)
                  array(IX)= array_interm(IY,ix)
                 ENDDO 
            call ave_std(array, Good_pixels(IY),ref_val,sd_val)
             ref_allwav(iy)=ref_val
             sd_allwav(iy)=sd_val
             var_allwav(iy)=sd_val/ref_val 
!  ENDIF for Iline             
          Endif 
!  Else if no good pixels                 
         ELSE
           sd_allwav(iy) =-9999   
           ref_allwav(iy)=-9999
           var_allwav(iy)=-9999  
!!  Endif for no good pixels               
          Endif 
!  ENDDO for number of wavelengths
         ENDDO
       
   
            

          do iy = 1,NWAV+2
 
!   WAVE 0.470 UM
 
          if(iy .eq.1)then
          REFW470=ref_allwav(iy)
          sd_w470=sd_allwav(iy)
 
!   WAVE 0.550 UM
 
          elseif(iy .eq.2)then
          REFW550=ref_allwav(iy)
          sd_w550=sd_allwav(iy)
 
!   WAVE 0.659 UM
 
         elseif (iy .eq.3)then
         REFW659=ref_allwav(iy)
         sd_w659=sd_allwav(iy)
         del_var_wav659=&
         -(1./(alog(WAVE(iy)/WAVE(Index_wave_865))))*&
         alog((1.+ var_allwav(iy))/&
             (1.+var_allwav(Index_wave_865)))

 
!   WAVE 0.865 UM
 
         elseif(iy .eq.4)then
        REFW865=ref_allwav(iy)
         sd_w865=sd_allwav(iy)

 
!   WAVE 1.24 UM
 
         elseif(iy .eq.5)then
         REFW124=ref_allwav(iy)
         sd_w124=sd_allwav(iy)
        del_var_wav124=&
         -(1./(alog(WAVE(iy)/WAVE(Index_wave_865))))*&
         alog((1.+ var_allwav(iy))/&
            (1.+var_allwav(Index_wave_865)))
 
!   WAVE 1.64 UM
 
      elseif(iy .eq. 6)then
        REFW164=ref_allwav(iy)
         sd_w164=sd_allwav(iy)
         del_var_wav164=&
        -(1./(alog(WAVE(iy)/WAVE(Index_wave_865))))*&
         alog((1.+ var_allwav(iy))/&
            (1.+var_allwav(Index_wave_865)))
 
!   WAVE 2.13  UM
 
         elseif(iy .eq. 7)then
        REFW213=ref_allwav(iy)
         sd_w213=sd_allwav(iy)
         del_var_wav213=&
         -(1./(alog(WAVE(iy)/WAVE(Index_wave_865))))*&
         alog((1.+ var_allwav(iy))/&
        (1.+var_allwav(Index_wave_865)))
      
     
 
!   WAVE 0.354  UM
 
         elseif(iy .eq. 8)then
         REFW354=ref_allwav(iy)
         sd_w354=sd_allwav(iy) 
    
!   WAVE 0.388  UM
 
         elseif(iy .eq. 9)then
        REFW388=ref_allwav(iy)
        sd_w388=sd_allwav(iy)  
         endif
! enddo for wavelengths
        enddo
!         if(ref_allwav(1) .gt. 0)&
!         print*,'ref here',(ref_allwav(iy),iy=1,9)
        
! If GLINT_ANGLE.LE.lt.glintthreshold store the Reflectances standard deviation and number of pixelsonly.
!  
!OFF GLINT..........
         If( GLINT_ANGLE.Gt.GLINT_THRESHOLD  .and.&
             ((REFW470 .gt.0  .and. REFW659 .gt.0) .and.&
             ( REFW470/ REFW659  .le. 0.75)))then 
                 Quality_dust_flag_off_glint=1
           endif
!    Recall of heavy dust in glint area.
          If( GLINT_ANGLE.le.GLINT_THRESHOLD  .and.&
            ((REFW470 .gt.0  .and. REFW659 .gt.0) .and.&
                ( REFW470/ REFW659  .le. 0.95)))then 
                  Quality_dust_flag_glint=1   
           endif
           

        If( GLINT_ANGLE.GT.GLINT_THRESHOLD  .or.&
          Quality_dust_flag_off_glint .eq.1 .or.&
           Quality_dust_flag_glint .eq.1) then  

!If there is valid data in 865 nm channel and at least one channel has valid data

      Total_good_points=Good_pixels(2)+Good_pixels(3)+Good_pixels(5)+&
      Good_pixels(6)+Good_pixels(7)
     
     
!         if(Good_pixels(Index_wave_865) .gt.10 .and. 
!    *     Total_good_points.gt.30) then


!           if(Good_pixels(Index_wave_865) .gt.5 .and. 
!     *     Total_good_points.gt.15) then

      
          if(Good_pixels(Index_wave_865) .gt.Number_pixels .and. &
          Total_good_points .gt. Tot_Number_pixels) then

             
!   If  Reflectance at wave 865 nm is gt than rhomin1
!   the aerosol content is  enough to process
!
             IF( REFW865 .GT. (Ref_ray(Index_wave_865)+&
             Ref_ray(Index_wave_865)/10.))THEN

!  Quality control flags are set for the number of cloud free pixels used
!  QCONTROL=1 means that box is less than 10% cloudfree.

 
          NUMDATA_550=good_pixels(index_wave_865)
          inn=good_pixels(index_wave_865)
!         if( inn .lt.20)QCONTROL=1
!         if( inn .gt.20)QCONTROL=0  
!          if( inn .lt.10)QCONTROL=1
!          if( inn .gt.10)QCONTROL=0


         IF( Iline .GT.1) then
           if( inn .lt.Good_Number_pixels)QCONTROL=1
           if( inn .gt.Good_Number_pixels)QCONTROL=0
         Else  
           if( inn .gt.Good_Number_pixels)QCONTROL=0
         Endif
 
!   If  Reflectance at wave 865 nm is less than   rhomin2
!   the aerosol content is not enough to derive size-distribution, but
!   enough to derive optical thickness...process
!
           IF( REFW865 .LT.(1.5*Ref_ray(Index_wave_865)))then
            QCONTROL=3
            Qcontrol_special=2
           endif

 
! Test if channel 1. 24  is good to process
 
!          if(Good_pixels(6) .ge. 10) then
           IF(REFW124 * COS(MTHET0*DTR) .LT. 3.600E-04) QCONTROL1=6
!            endif             
 
! Test if channel 2.13 and 1.64 um are good to process
  
!           if(Good_pixels(6) .ge. 10) then
           IF(REFW164 * COS(MTHET0*DTR) .LT. 3.600E-04) QCONTROL1=3
!            endif
         
!           if(Good_pixels(7) .ge. 10) then
           IF(REFW213 * COS(MTHET0*DTR) .LT. 3.110E-04) QCONTROL1=4
!            endif
          
!            if(Good_pixels(6) .ge. 10 .and. Good_pixels(7) .ge. 10)then
            IF(REFW213 * COS(MTHET0*DTR) .LT. 3.110E-04 .AND.&
            REFW164 * COS(MTHET0*DTR)  .LT. 3.600E-04)QCONTROL1=5
!            endif


! Process only for valid data
            QCONTROL=QCONTROL1
         IF( QCONTROL .GT.0) THEN
! The aerosol type and content are variable
        IF(var_allwav(Index_wave_865) .GT.0.05 .AND. &
           ABS(del_var_wav659) .GT. 0.15 .OR.&
           ABS(del_var_wav164) .GT. 0.15 .OR.&
           ABS(del_var_wav213) .GT. 0.15) QCONTROL=6
! The aerosolcontent are variable not the aerosol content
        IF(var_allwav(Index_wave_865) .GT.0.05 .AND.&
            ABS(del_var_wav659) .LT. 0.15 .OR.&
            ABS(del_var_wav164) .LT. 0.15 .OR.&
           ABS(del_var_wav213) .LT. 0.15) QCONTROL=7
         ENDIF

!  If  Reflectance at wave 865 nm is lt than rhomin1
!  the aerosol content is  not enough to process
         ELSE
!              QCONTROL=-3
               QCONTROL= 17
             Qcontrol_special=1
              
             
        ENDIF
! If there is no valid data in all channel s  
                 ELSE
                  QCONTROL=-20
                 Qcontrol_special=-2
                 ENDIF
! Glint and only Reflectances ,sd and number of pixels will bestored
                 ELSE
               Qcontrol_special=3
            NUMDATA_550=good_pixels(index_wave_865)
            QCONTROL=11
                 ENDIF
!  ELSE for number of cloudfree pixels inn box is too cloudy
        ELSE
! If cloud free pixels are lt 10 no processing is performed .
              QCONTROL=-22
              Qcontrol_special=-2
        ENDIF
         If( Iline .eq.1 .and.inn .lt.Good_Number_pixels)QCONTROL=-20
        IF (QCONTROL .LT.0 ) THEN
             REFW470=99999
             REFW550=99999
             REFW659=99999
             REFW865=99999
             REFW124=99999
             REFW164=99999
             REFW213=99999 
        ENDIF  
 !           write(36,*)'inside average',iscan,idata,&
 !           REFW470,REFW550,REFW659,REFW865,&
 !           REFW124,REFW164,REFW213,Good_Number_pixels,QCONTROL,refw354,refw388
           RETURN
           END 

!*********************************************************************

       SUBROUTINE storeref_data(IDATA,SDS_ref,SDS_ref_STD,NUMDATA_550,&
        sd_W470,sd_W550,sd_W659,sd_w865,sd_W124,sd_W164,Sd_W213,&
        SDS_CLDFRC_ocean,SDS_NUMPIXELS,iscan,&
        SDS_SCAT_ANGLE_OCEAN,CLD_FRAC,SCAT_ANGLE_OCEAN,&
        sd_w354,sd_w388,sd_w0p75)

!-----------------------------------------------------------------------
!!F90
!
!!DESCRIPTION: This subroutine stores the average Reflectance for all wavelengths
!             into arrays to be written to HDF file.
!
!!INPUT PARAMETERS:
!         IDATA          Data point Number
!
!   STANDARD DEVIATION of Wavelengths used modis,ch1-ch7  identified by sd_w
!
!       SD_W659
!       SD_W865
!       SD_W470
!       SD_W550
!       SD_W124
!       SD_W164
!       SD_W213
!       QCONTROL       quality control
!SCAT_ANGLE_OCEAN      Scattering angle
!        CLD_FRAC      cloud Fraction
!       NUMDATA_550  Number of pixels used for processing
!
!    Reflectances of Wavelengths used modis,ch1-ch7  identified by REFW*
!   passed through mod04.inc
!       REFW659
!       REFW865
!       REFW470
!       REFW550
!       REFW124
!       REFW164
!       REFW213
!!OUTPUT PARAMETERS:for HDF write
!        SDS_ref               Array for Reflectance data
!        SDS_CLDFRC_ocean             Array for Cloud Fraction
!       SDS_QCONTROL           Array for quality control
!       SDS_NUMPIXELS          Array for number of Pixels used
!       SDS_SCAT_ANGLE_OCEAN   Array for Scattering angle
!
!
 

          IMPLICIT  NONE
          SAVE

          INCLUDE 'mod04.inc' 
          INCLUDE 'read_Sat_MODIS.inc'

          INTEGER I,IDATA,NUMDATA_550,Num_Wav,iscan
          REAL sd_W659,sd_W865,sd_W470,sd_W550,sd_W124,sd_W164,&
           sd_W213,CLD_FRAC,sd_w354,sd_w388,sd_w0p75
           REAL SCAT_ANGLE_OCEAN
          REAL  SDS_ref(NUMCELLS_B,NWAV+2)
          REAL  SDS_ref_STD(NUMCELLS_B,NWAV+2)
          REAL  SDS_NUMPIXELS(NUMCELLS_B,NWAV+2),&
        SDS_CLDFRC_ocean (NUMCELLS_B),SDS_SCAT_ANGLE_OCEAN(NUMCELLS_B)

!                **** Store input reflectance data to be written for
!                     HDF output file
! NPP WAVELENGTHS   refw412,refw443,refw0p75
         SDS_ref(IDATA,1) =  REFW470 
         SDS_ref(IDATA,2) =  REFW550 
         SDS_ref(IDATA,3) =  REFW659 
         SDS_ref(IDATA,4) =  REFW865 
         SDS_ref(IDATA,5) =  REFW124 
         SDS_ref(IDATA,6) =  REFW164 
         SDS_ref(IDATA,7) =  REFW213 
         SDS_ref(IDATA,8) =  REFW354
         SDS_ref(IDATA,9) =  REFW388
         SDS_ref_STD(IDATA,1) =  sd_W470
         SDS_ref_STD(IDATA,2) =  sd_W550
         SDS_ref_STD(IDATA,3) =  sd_w659
         SDS_ref_STD(IDATA,4) =  sd_W865
         SDS_ref_STD(IDATA,5) =  sd_W124
         SDS_ref_STD(IDATA,6) =  sd_W164
         SDS_ref_STD(IDATA,7) =  sd_W213
         SDS_ref_STD(IDATA,8) =  sd_w354
         SDS_ref_STD(IDATA,9) =  sd_w388 
         Do Num_Wav=1,NWAV+2 
         SDS_NUMPIXELS(IDATA,Num_Wav)=Good_pixels(Num_Wav)
         Enddo
        
         RETURN 
         END



!*********************************************************************

      SUBROUTINE INTANGLE_ray(PHC,THET,THET0,ref_rayall,&
        MTHET0,MTHET,MPHI,Ref_ray,KSZAM1,KSZAP1,KTHEM1,KTHEP1,&
        KPHIM1,KPHIP1,Indx_wspped1,Indx_wspped2,WSPEED,&
        Wind)

!-----------------------------------------------------------------------
!!F90
!
!!DESCRIPTION: This subroutine interpolates rayleigh Reflectance for
!             measured solar zenith, view and azimuthal angle.
!!INPUT PARAMETERS:
!                PHC        Azimuth angle.
!               THET        View angle.
!               THET0       Solar zenith angle.
!          Ref_rayall      Rayleigh Reflectance
!             MTHET0      Measured solar Zenith angle
!               MTHET      Measured view  Zenith angle
!                MPHI      Measured Azimuthal Angle
!               KSZAM1     Starting Index for solar zenith angle
!               KSZAP1     Ending Index for solar zenith angle
!               KTHEM1     Starting Index for view angle
!               KTHEP1     Ending   Index for  view angle
!               KPHIM1     Starting Index for  azimuth angle
!              KPHIP1     Ending   Index for  azimuth angle
!OUTPUT PARAMETERS: interpolated Rayleigh
!             Ref_ray      interpolated Rayleigh Reflectance
!
       USE Linear_interpolation 
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      INTEGER IJ,IPHI,ITH,ITH0,IWAV,LOPT
      REAL  PHC(NPHI),THET(NTHET),THET0(NTH0)
      REAL  Ref_rayall( Lut_indx,NPHI, NTHET, NTH0,NWAV)
      REAL  Ref_ray(NWAV),Wspeed,Windspeed_Lut(Lut_indx)
      REAL  MTHET0,MTHET,MPHI,Ray(NWAV,Lut_indx)
      REAL  X(100),Y(100),XX1(100),YY1(100),XX2(100),YY2(100),Y1
      REAL XX3(100),YY3(100),Wind(Lut_indx)
      INTEGER KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1
      INTEGER LL,MM,NN,Indx_wspped1,Indx_wspped2,index_wspeed,Num


! LOOP IS AROUND THET0,NEW FILE FOR EACH THETA0
!
          
          Num=0 
          do index_wspeed = Indx_wspped1,Indx_wspped2
          Num=Num+1
          Windspeed_Lut(num)=Wind(index_wspeed)  
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
          DO 20 IWAV=1,NWAV
          LL=0
          DO 30 ITH0 = KSZAM1,KSZAP1
          MM=0
          DO 40 ITH= KTHEM1,KTHEP1
          NN=0
          DO 50 IPHI =KPHIM1,KPHIP1
          NN=NN+1
          X(NN)=PHC(IPHI)
          Y(NN)=Ref_rayall(index_wspeed,IPHI,ITH,ITH0,IWAV)
 50       CONTINUE
!
           CALL INTERP(NN,MPHI,X,Y,Y1)
            MM=MM+1
           XX1(MM)=THET(ITH)
           YY1(MM)=Y1
 40         CONTINUE
            y1=0.0
           CALL INTERP(MM,MTHET,XX1,YY1,Y1)
!
             LL=LL+1
            XX2(LL)=THET0(ITH0)
            YY2(LL)=Y1
 30        CONTINUE

             y1=0.0
           CALL INTERP(LL,MTHET0,XX2,YY2,Y1) 
            Ray(IWAV,Num)=Y1  
  20       CONTINUE
           ENDDO
           DO  IWAV=1,NWAV
           Do  IJ=1,NUM 
           XX3(Ij)=Windspeed_Lut(IJ)
           YY3(IJ)= Ray(IWAV,IJ) 
           enddo
            y1=0.0
           CALL INTERP(Num,WSPEED,XX3,YY3,Y1) 
           Ref_ray(IWAV)=(Y1*PI)/(COS(DTR*MTHET0))
           ENDDO 
           
            RETURN
            END



!*********************************************************************

         SUBROUTINE INTANGLE(PHC,THET,THET0,AINTS,TAUAS,&
                 AINTB,TAUAB,MTHET0,MTHET,MPHI,REFSMALL,REFBIG,&
                Ref_ray,KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1,&
                Indx_wspped1,Indx_wspped2,WSPEED,Wind)

!----------------------------------------------------------------------
!!F90
!
!!DESCRIPTION:  Subroutine INTANGLE interpolates the lookup
!              reflectances to the measured geometry.
!
!!INPUT PARAMETERS:
!          PHC        Azimuth angle.
!         THET        view angle.
!         THET0       Solar zenith angle.
!         AINTS       radiance(l/fo),fo=1 for small mode
!         AINTB       radiance(l/fo),fo=1 for large mode
!         TAUAS       optical thickness for small mode
!        TAUAb       optical thickness for large  mode
!         KSZAM1      Starting Index for solar zenith angle
!         KSZAP1      Ending Index for solar zenith angle
!         KTHEM1      Starting Index for view angle
!         KTHEP1      Ending   Index for  view angle
!         KPHIM1      Starting Index for  azimuth angle
!         KPHIP1      Ending   Index for  azimuth angle
!      MTHETA0        Solar zenith angle.
!        MTHET        View angle.
!        MPHI        Azimuth angle.
!
!!OUTPUT PARAMETERS:
!      REFSMALL       reflectance from input interpolated
!                   for MODIS geometry(small modes)
!     REFBIG        reflectance from input interpolated
!                    for MODIS geometry(large modes)
!     Ref_ray        Reflectance for rayleigh only
!
!-----------------------------------------------------------------------
       USE Linear_interpolation 
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      CHARACTER*132 LINE
      CHARACTER*45  LINE1
      INTEGER ICASE,IJ,IPHI,ISIZE,ITAU,ITH,ITH0,IWAV,LOPT,NUMCASE
      REAL  PHC(NPHI),THET(NTHET),THET0(NTH0)
      REAL  AINTS(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASES)
      REAL  TAUAS(NUMCASES,NWAV,NTAU)
      REAL  AINTB(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASEB)
      REAL  TAUAB(NUMCASEB,NWAV,NTAU)
      REAL  Ref_ray(NWAV)
      REAL  REFSMALL(NUMCASES,NWAV,NTAU),REFBIG(NUMCASEB,NWAV,NTAU)
      REAL  REFSMALL_inter(NUMCASES,NWAV,NTAU,Lut_indx)
      Real REFBIG_inter(NUMCASEB,NWAV,NTAU,Lut_indx)
      REAL  MTHET0,MTHET,MPHI, WSPEED,Wind(Lut_indx)
      REAL  X(100),Y(100),XX1(100),YY1(100),XX2(100),YY2(100),Y1
       REAL  XX3(100),YY3(100),Windspeed_Lut(Lut_indx)
      INTEGER KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1
      INTEGER LL,MM,NN,index_wspeed,Indx_wspped1,Indx_wspped2,num
       
    

! LOOP IS AROUND THET0,NEW FILE FOR EACH THETA0
 

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
          Num=0
         
          do index_wspeed = Indx_wspped1,Indx_wspped2
          Num=Num+1
          Windspeed_Lut(num)=Wind(index_wspeed )  
          DO 5 ISIZE= 1,2
          IF(ISIZE .EQ.1)NUMCASE=NUMCASES
          IF(ISIZE .EQ.2)NUMCASE=NUMCASEB
          DO 10 ICASE =1,NUMCASE
          DO 20 IWAV=1,NWAV
!  interpolate starting from optical thickess indexed 2.
! at the end of routine fill the array with interpoltaed ta=0.0
          DO 30  ITAU  = 2,NTAU
             NN=0
          DO  40 ITH0  = KSZAM1,KSZAP1
             MM=0
          DO  50  ITH  = KTHEM1,KTHEP1
              LL=0
          DO 60  IPHI  = KPHIM1,KPHIP1
           LL=LL+1
          X(LL)=PHC(IPHI)
         IF( ISIZE.EQ.1)Y(LL)=AINTS(index_wspeed,IPHI,ITH,ITH0,ITAU,IWAV,ICASE)
         IF( ISIZE.EQ.2)Y(LL)=AINTB(index_wspeed,IPHI,ITH,ITH0,ITAU,IWAV,ICASE)
          
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
!                  ****** change to Reflectance units
       IF(ISIZE.EQ.1)REFSMALL_inter(ICASE,IWAV,ITAU,num)=Y1
       IF(ISIZE.EQ.2)REFBIG_inter(ICASE,IWAV,ITAU,num)=Y1  
 
  30       CONTINUE
  20       CONTINUE
  10       CONTINUE
   5        CONTINUE
           ENDDO
           
          DO   ISIZE= 1,2
          IF(ISIZE .EQ.1)NUMCASE=NUMCASES
          IF(ISIZE .EQ.2)NUMCASE=NUMCASEB
          DO   ICASE =1,NUMCASE
          DO   IWAV=1,NWAV
          DO   ITAU  = 2,NTAU
          Do  index_wspeed = 1,Num
         XX3(index_wspeed)=Windspeed_Lut(index_wspeed)
          IF(ISIZE .EQ.1)&
          YY3(index_wspeed)=REFSMALL_inter(ICASE,IWAV,ITAU,index_wspeed)
         IF(ISIZE .EQ.2)&
          YY3(index_wspeed)=REFBIG_inter(ICASE,IWAV,ITAU,index_wspeed) 
          Enddo
          CALL INTERP(Num,WSPEED,XX3,YY3,Y1)
       IF(ISIZE.EQ.1)REFSMALL(ICASE,IWAV,ITAU)=(Y1*PI)/(COS(DTR*MTHET0))
       IF(ISIZE.EQ.2)REFBIG(ICASE,IWAV,ITAU)=(Y1*PI)/(COS(DTR*MTHET0))
          
        Enddo
        Enddo
        Enddo
        Enddo
        
 
!fill empty array of itau=1 with rayleigh (taua=0.0)
 
           DO   ISIZE= 1,2
           IF(ISIZE .EQ.1)NUMCASE=NUMCASES
           IF(ISIZE .EQ.2)NUMCASE=NUMCASEB
             DO   ICASE =1,NUMCASE
              DO   IWAV=1,NWAV 
            IF(ISIZE.EQ.1)REFSMALL(ICASE,IWAV,1)=Ref_ray(IWAV)
            IF(ISIZE.EQ.2)REFBIG(ICASE,IWAV,1)=Ref_ray(IWAV)
             Enddo
             Enddo
          ENddo  
            RETURN
            END 
!*********************************************************************

      SUBROUTINE SUBMIN(ISMALL,IBIG,FUNMIN,Ref_ray)

!----------------------------------------------------------------------
!!F90
!
!!DESCRIPTION:  THIS SUBROUTINE IS A GENERAL PURPOSE ROUTINE
!               AND INTERPOLATES EXPONANTIALLY. IGNORE SENDS
!               VALUE OF 0 IF SUBROUTINE INTEXP IS USED AND
!               1 IF IT FINDS THE FUNCTION TO BE LINEAR.
!
!!INPUT PARAMETERS:
!          REFSMALL    reflectance from input interpolated
!                      for MODIS geometry(small modes)
!            REFBIG    reflectance from input interpolated
!                      for MODIS geometry(large modes)
!              TAUA    Optical thickness for lookup table.
!             NSIZE    Size,big small
!           NUMCASE    Number of cases
!              NWAV    Number of wavelengths
!              NTAU    Number of optical thicknesses
!
!!OUTPUT PARAMETERS:   NONE
!
!!DESIGN NOTE:
! AX IS LOWER END-POINT OF THE INTERVAL IN WHICH MINIMUM OF FMIN
! IS TO BE LOCATED.
! BX IS THE UPPER END-POINT OF THE INTERVAL.
! TOL IS LENGTH OF THE FINAL SUBINTERVAL CONTAINING THE MINIMUM.
! FMIN IS OUTPUT , THE APPROXIMATE MINIMUM OF THE FUNCTION FMIN ON THE
! ORGINAL INTERVAL(A,B) 
!-----------------------------------------------------------------------

      IMPLICIT NONE
       Save
      INCLUDE 'mod04.inc'

! THIS SUBROUTINE IS A GENERAL PURPOSE ROUTINE  AND INTERPOLATES
! EXPONANTIALLY.IGNORE SENDS VALUE OF 0 IF SUBROUTINE INTEXP IS
! USED AND 1 IF IT FINDS THE FUNCTION TO BE LINEAR.
      INTEGER IBIG,ISMALL,ITAU,IWAV,NJTAU,NJWAV
      REAL   FB,FC,XMIN
      REAL   TOL,AX,BX,FMIN,funmin
      REAL  TAUAS(NUMCASES,NWAV,NTAU),Ref_ray(NWAV)
      real Ref_ray_save(NWAV)
      REAL  REFSMALL(NUMCASES,NWAV,NTAU),REFBIG(NUMCASEB,NWAV,NTAU)
      EXTERNAL  FMIN
      COMMON/TWO/NJTAU,NJWAV,Ref_ray_save
       
      NJTAU=NTAU
      NJWAV=NWAV
      do Iwav=1,NWAV
      Ref_ray_save(iwav)=Ref_ray(iwav)
      enddo 
      AN=0.0
! TAUAS OR TAUAB CAN BE USED HERE, THEY CARRY SAME VALUES......


!
!  SUBROUTINE MNBRAK AND GOLDEN ARE USED TO FIND THE MINIMAL OF
! FUNCTION FMIN.  LOOK FOR NUMERICAL RECIPE FOR THE EXPLANATION
!
! AX IS LOWER END-POINT OF THE INTERVAL IN WHICH MINIMUM OF FMIN
! IS TO BE LOCATED.
! BX IS THE UPPER END-POINT OF THE INTERVAL.
! TOL IS LENGTH OF THE FINAL SUBINTERVAL CONTAINING THE MINIMUM.
! FMIN IS OUTPUT , THE APPROXIMATE MINIMUM OF THE FUNCTION FMIN ON THE
! ORGINAL INTERVAL(A,B)

         AX = 0.00
         BX = 1.00

!                  ***** call to halving routine findes the minimum
!                         function fmin
!
        CALL HALVING(AX,BX,FMIN,XMIN,FUNMIN)
!                  ***** call tocompute tausmall and big for xminetc
!
              AN=XMIN
        RETURN
        END



!*****************************************************************

      REAL FUNCTION FMIN(XMIN)

!----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION:    FUNCTION FMIN IS FUNCTION TO BE MINIMIZED.
!
!!INPUT PARAMETERS:
!                XMIN    minimum value
!
!!OUTPUT PARAMETERS:
!                FMIN    value of function at minimum value
!
!
! 
!-----------------------------------------------------------------------
       USE Linear_interpolation 
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'
      INCLUDE 'read_Sat_MODIS.inc'

      INTEGER ITAU,IWAV,LOPT,NJTAU,NJWAV,NUMWAV,INDEX_wave
      REAL  RADWAV,Y1,asum1,sum_good_pixels,Ref_ray_save(NWAV)
      REAL XMIN,XS(100,100),YS(100,100),X(100),Y(100),Denom
      COMMON/TWO/NJTAU,NJWAV,Ref_ray_save

! If two long wavelengths are not available,force mode to small mode.
           if( Iline .gt.1) then
         IF( QCONTROL1 .EQ. 3 .or.QCONTROL1 .EQ. 4 .or.&
        QCONTROL1 .eq.5)XMIN=1.0
           endif
            
        ASUM=0.0
        ASUM1=0.0
         sum_good_pixels=0
        FMIN=0.0
        TAUX55=0.0
         
 
!                   ******COMPUTE FUN TO BE MINIMIZED.
!

       DO 100 IWAV = 1,NJWAV
       DO 101 ITAU = 1,NJTAU
         XS(ITAU,IWAV)=TAUAL(ITAU)
       
         YS(ITAU,IWAV)=(REFSMALLL(IWAV,ITAU)*XMIN)&
                 +((1.-XMIN)*REFBIGL(IWAV,ITAU))
101    CONTINUE
100    CONTINUE
 
!                  ********COMPUTE FOR TAU FOR  USING
!                           RADAINCE VALUE OF 0.860 UM AND
!                           TAU VALUES OF WAV550
!
         DO 102 ITAU = 1,NJTAU
            X(ITAU)=YS(ITAU,Index_wave_865)
            Y(ITAU)=XS(ITAU,Index_wave_550)
102      CONTINUE
           Y1=0
           CALL INTERP(NJTAU,REFW865,X,Y,Y1)
            TAUX55=y1

 
!                ******* FOR TAUX55 COMPUTE Reflectance FOR
!                     ALL WAVELENGTHS.
!
      DO 104 IWAV = 1,NJWAV
         DO 103 ITAU = 1,NJTAU
            X(ITAU)=TAUAL(ITAU)
            Y(ITAU)=YS(ITAU,IWAV)
103      CONTINUE
         y1=0
          CALL INTERP(NJTAU,TAUX55,X,Y,Y1)
            ALXWAV(IWAV)=Y1
104   CONTINUE
! For summation do not use 1st wavelength
          NUMWAV=NWAVN-1
!if wavelength 1.24 or 1.64 or 2.13 is small  do not use for inversion
      IF(QCONTROL1.EQ.3 .OR. QCONTROL1.EQ.4 .or. QCONTROL1.EQ.6) NUMWAV=NWAVN-2
           
        IF(QCONTROL1 .EQ.5)NUMWAV=NWAVN-3 
       
       DO 105 IWAV = 1,NUMWAV
         INDEX_wave=IWAV+1
! if 1.24 is small ignore 1.24
        IF( IWAV .ge.3 .AND.QCONTROL1.EQ.6)INDEX_wave=IWAV+2
! if 1.64 is small ignore 1.64
          IF( IWAV .EQ.5 .AND.QCONTROL1.EQ.3)INDEX_wave=IWAV+2
! if 2.13 issmall ignore 2.13
          IF( IWAV .EQ.5 .AND.QCONTROL1.EQ.4)INDEX_wave=IWAV+1
! if 1.64 and 2.13 is small  ignore 1.64 and 2.13
         IF( IWAV .EQ.4 .AND.QCONTROL1.EQ.5)INDEX_wave=IWAV+1
 
        
        Denom=(ref_allwav(INDEX_wave)-Ref_ray_save(INDEX_wave))+0.01 
        if( Denom .lt. 0.01)Denom=0.01
          ASUM1= ASUM1+(((ref_allwav( INDEX_wave)-ALXWAV(INDEX_wave))&
          /(Denom))**2)*Good_pixels(INDEX_wave)
        sum_good_pixels=sum_good_pixels+Good_pixels(INDEX_wave)
!
105   CONTINUE
               ASUM=SQRT(ASUM1/REAL(sum_good_pixels))
         FMIN=ASUM
      RETURN
      END



!********************************************************************

        SUBROUTINE COMPUTE_alltau(EXTSMALL,EXTBIG,EXTNORSMALL,&
            EXTNORBIG,ISMALL,IBIG,Indx_wspped1,Iscan,idata,&
            Indx_wspped2,WSPEED,Wind,Ext_554_small,Ext_554_large)
 
!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine computes CCN, assymetry factor,
!              backscattering ratio, effective radius, effective
!              variance, and optical thicknesses for large and small
!             mode.
!
!!INPUT PARAMETERS: all the variables are for lookup table.
!       TAUX55SMALL Optical thickness contribution small particles
!       TAUX55BIG   Optical thickness contribution large particles
!        NWAV       Number of wavelengths
!        EXTSMALL   extinction coeff small particles
!            IBIG   Index for big mode
!           ISMALL  Index for small mode.
!
!!OUTPUT PARAMETERS:   all variables are computed for large and
!      TAU_COMPUTED Optical depth for all wavelengths
!TAU_SMALL_COMPUTED Optical depth for all wavelengths for Small particles
!TAU_BIG_COMPUTED   Optical depth for all wavelengths for large particles
! 
!----------------------------------------------------------------------
       USE Linear_interpolation 
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      INTEGER IBIG,ISMALL,IWAV,index_wspeed,Indx_wspped1,Indx_wspped2 
      Integer iJ,NUM,Iscan,idata
      REAL EXTSMALL(Lut_indx,Numcases,NWAV),EXTbig(Lut_indx,Numcaseb,NWAV)
      REAL EXTNORSMALL(NUMCASES,NWAV),EXTNORBIG(NUMCASEB,NWAV)
      REAL WSPEED,Wind(Lut_indx),Windspeed_Lut(Lut_indx)
      Real EXTNORSMALL_inter(NUMCASES,NWAV,Lut_indx)
      REAL EXTNORBIG_inter(NUMCASEB,NWAV,Lut_indx)
      Real Ext_554_small(Lut_indx,NUMCASES),Ext_554_large(Lut_indx,NUMCASEB)
      Real xx1(100),yy1(100),yy2(100),y1
       DO IJ = 1,100 
          XX1(IJ)=0.0
          YY1(IJ)=0.0
          YY2(IJ)=0.0 
          ENDDO
           DO IWAV = 1,NWAV
           TAU_COMPUTED(IWAV)  =-99999
           Enddo
!      compute normalized extiction coeff
           Num=0 
       DO index_wspeed = Indx_wspped1,Indx_wspped2
          Num=Num+1
          Windspeed_Lut(num)=Wind(index_wspeed)  
           DO IWAV = 1,NWAV
        EXTNORSMALL_inter(ISMALL,IWAV,num)=EXTSMALL(index_wspeed,ISMALL,IWAV)/ &
                                           Ext_554_small(index_wspeed,ISMALL) 
!     *                              EXTSMALL(index_wspeed,ISMALL,Index_wave_550)
        EXTNORBIG_inter(IBIG,IWAV,num)= EXTBIG(index_wspeed,IBIG,IWAV)/ &
                                     Ext_554_large(index_wspeed,IBIG) 
!         EXTBIG(index_wspeed,IBIG,Index_wave_550)
          ENDDO 
       ENDDO
       
       
       
          DO  IWAV=1,NWAV
              Do index_wspeed = 1,Num
           XX1(index_wspeed)=Windspeed_Lut(index_wspeed)
           YY1(index_wspeed)= EXTNORSMALL_inter(ISMALL,IWAV,index_wspeed) 
           YY2(index_wspeed)= EXTNORBIG_inter(IBIG,IWAV,index_wspeed)
            
           enddo
            y1=0.0
           CALL INTERP(Num,WSPEED,XX1,YY1,Y1) 
           EXTNORSMALL(ISMALL,IWAV)=y1
           y1=0.0
           CALL INTERP(Num,WSPEED,XX1,YY2,Y1) 
           EXTNORBIG(IBIG,IWAV)=y1
           ENDDO 
          
!Changes 2/28
!     

        if(TAUX55 .le.0.7 .and. Quality_dust_flag_glint  .eq.1)TAUX55= -0.02 
!   if tau at 0.55um is greater than maxtau(ie 5.00) taux55 is set to maxtau for
!  glint and off glint heavy dust episodes

        if((TAUX55  .GT. MAXTAU .and. Quality_dust_flag_glint  .eq.1) &
         .or. (TAUX55 .GT. MAXTAU .and. Quality_dust_flag_off_glint .eq.1)) &
           TAUX55=MAXTAU 
     
         
       
        TAUX55SMALL=TAUX55*AN
        TAUX55BIG=TAUX55-TAUX55SMALL 
        
        
!     TAUX55 is computed at 0.544 because LUT is scaled to 0.544   report this as total for Himawari
!     Optical depths are computed for all wavelengths. We have to save  TAUX55        
            DO IWAV = 1,NWAV  
        TAU_COMPUTED(IWAV)=TAUX55SMALL*EXTNORSMALL(ISMALL,IWAV)&
                              +TAUX55BIG*EXTNORBIG(IBIG,IWAV)
       TAU_SMALL_COMPUTED(IWAV)=TAUX55SMALL*EXTNORSMALL(ISMALL,IWAV)
       TAU_BIG_COMPUTED(IWAV)=   TAUX55BIG *EXTNORBIG(IBIG,IWAV) 	
        ENDDO
           
     
      
!   Right now I will report TAUX55 in 2nd wavelength to replace computed  optical depth for 509  
              TAU_COMPUTED(2)=  TAUX55
             
! end of  Changes 2/28
       RETURN
       END


!********************************************************************

          SUBROUTINE COMPUTE_INVERVAR(EXTSMALL,&
               EXTBIG,EXTNORSMALL,EXTNORBIG,CCNSMALL,MOMENTSSMALL,&
               MOMENTSBIG,NSMALL,NBIG,EFFRADIUS,EFFVARIANCE,TAUSMALL,&
               TAUBIG,ISMALL,IBIG,ALBEDOSMALL,ALBEDOBIG,&
               CCN,BACKSCTTSMALL,BACKSCTTBIG,ASSYMSMALL,ASSYMBIG,&
              BACKSCATT,ASSYM,REF_FLUX,TRANS_FLUX,ALBEDO_R_SMALL_final,&
             ALBEDO_R_BIG_final,ALBEDO_T_SMALL_final,ALBEDO_T_BIG_final,&
             MASS_CON_OCEAN,ALBEDO_R_SMALL_tau,MTHET0,ALBEDO_R_big_tau,&
             Indx_wspped1,Indx_wspped2,WSPEED,Wind,Ext_554_small,Ext_554_large,&
             ipart)

!-----------------------------------------------------------------------
!!F90
!
!!DESCRIPTION: This subroutine computes CCN, assymetry factor,
!              backscattering ratio, effective radius, effective
!              variance, and optical thicknesses for large and small
!              mode.
!
!!INPUT PARAMETERS: all the variables are for lookup table.
!       TAUX55SMALL Optical thickness contribution small particles
!       TAUX55BIG   Optical thickness contribution large particles
!        NUMCASE    Number of large and small cases
!        NWAV       Number of wavelengths
!        EXTSMALL   extinction coeff small particles
!    MOMENTSSMALL   Moments for small particles.
!        CCNSMALL   ccn for small particles.
!     EXTNORSMALL   Normalized extinction for small particles.
!   BACKSCTTSMALL   Back scattering for small particles.
!      ASSYMSMALL   Assymetry factor for small particles.
!     ALBEDOSMALL   single scattering albedo small particles
!          EXTBIG   extinction coeff large particles
!      MOMENTSBIG   Moments for large particles.
!       EXTNORBIG   ccn for big particles.
!     BACKSCTTBIG   Back scattering for large particles.
!         ASSYMBIG  Assymetry factor for large particles.
!       ALBEDOBIG   single scattering albedo large particles
!ALBEDO_R_SMALL_final Albedo for monodirectionally incident radiation(SMALL)
!ALBEDO_R_BIG_final  Albedo for monodirectionally incident radiation(LARGE)
!ALBEDO_T_SMALL_final Trans for monodirectionally incident radiation(SMALL)
!ALBEDO_T_BIG_final   Trans for monodirectionally incident radiation(LARGE)
!            IBIG   Index for big mode
!           ISMALL  Index for small mode.
!
!!OUTPUT PARAMETERS:   all variables are computed for large and
!                     small modes for the different solutions.
!           NSMALL  ratio of taux55small to extnor(small mode)
!           NBIG    ratio of taux55small to extnor(large mode)
!         EFFRADIUS effective radius weighed for large and small mode
!       EFFVARIANCE standard deviation for effective radius
!        TAUSMALL   ratio of TAUX55SMALL/ EXTNORSMALL
!        TAUBIG     ratio of TAUX55BIG  /EXTNORBIG
!        BACKSCATT  backscattering facror for large and small
!        ASSYM      assymetry factor for small and large
!        CCN        ccn for the solution
!      TAU_COMPUTED Optical depth for all wavelengths
!TAU_SMALL_COMPUTED Optical depth for all wavelengths for Small particles
!TAU_BIG_COMPUTED   Optical depth for all wavelengths for large particles
!     REF_FLUX      Albedo for monodirectionally incident radiation
!     TRANS_FLUX    Trans for monodirectionally incident radiation
!
!----------------------------------------------------------------------
       USE Linear_interpolation 
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      INTEGER IBIG,ISMALL,IWAV,index_wspeed
      REAL EXTSMALL(Lut_indx,Numcases,NWAV),EXTbig(Lut_indx,Numcaseb,NWAV)
      REAL EXTNORSMALL(NUMCASES,NWAV),EXTNORBIG(NUMCASEB,NWAV)
      REAL MOMENTSSMALL(Lut_indx,numcases,4),MOMENTSBIG(Lut_indx,NUMCASEB,4)
      REAL CCNSMALL(Lut_indx,nUMCASES),NSMALL,NBIG
      REAL TAUSMALL(NWAV),TAUBIG(NWAV)
      REAL EFFRADIUS,EFFVARIANCE,CCN,BACKSCATT(NWAV),ASSYM(NWAV)
      REAL ALBEDOSMALL(Lut_indx,NUMCASES,NWAV),ALBEDOBIG(Lut_indx,NUMCASEB,NWAV)
      REAL BACKSCTTSMALL(Lut_indx,NUMCASES,NWAV),BACKSCTTBIG(Lut_indx,NUMCASEB,NWAV)
      REAL ASSYMSMALL(Lut_indx,NUMCASES,NWAV),ASSYMBIG(Lut_indx,NUMCASEB,NWAV)
      REAL REF_FLUX(NWAV),TRANS_FLUX(NWAV)
      REAL ALBEDO_R_SMALL_final(NWAV,NUMCASES)
      REAL ALBEDO_R_BIG_final(NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL_final(NWAV,NUMCASES)
      REAL ALBEDO_T_BIG_final(NWAV,NUMCASEB),MTHET0
      REAL MASS_CON_OCEAN,ALBEDO_R_SMALL_tau(NTAU,NWAV,NUMCASES)
      Real ALBEDO_R_big_tau(NTAU,NWAV,NUMCASEb)
      Real  EXTSMALL_inter(Lut_indx),CCN_inter(Lut_indx )
      Real  EXTBIG_inter(Lut_indx),MOMENTSSMALL3_inter(Lut_indx)
      Real  MOMENTSSMAL2_inter(Lut_indx ),MOMENTSBIG2_inter(Lut_indx )
      Real  MOMENTSBIG3_inter(Lut_indx ),MOMENTSSMAL4_inter(Lut_indx),MOMENTSBIG4_inter(Lut_indx) 
      Real  ALBEDOSMAL_inter(Lut_indx ), BACKSCTTSMAL_inter(Lut_indx ),ASSYMBIG_inter(Lut_indx)
      Real  ALBEDOBIG_inter(Lut_indx ),  BACKSCTTBIG_inter(Lut_indx ),ASSYMSMALL_inter(Lut_indx)  
      Real  y1,y2,y3,y4,y5,y6,Windspeed_Lut(Lut_indx),Wind(Lut_indx),wspeed
      Integer Indx_wspped1,Indx_wspped2,num,ipart
      Real Ext_554_small(Lut_indx,Numcases),Ext_554_large(Lut_indx,NumcaseB)
          
      NBIG=0.0
      NSMALL=0.0
      CCN=0.0
      EFFRADIUS=0.0
      MASS_CON_OCEAN=0.0
      
!........ Interpolation to be set here........temporary fix...      
           Num=0 
         DO index_wspeed = Indx_wspped1,Indx_wspped2
          Num=Num+1
          Windspeed_Lut(num)=Wind(index_wspeed)  
     
!                             *** compute  ratio extsmall
!
         
!          EXTSMALL_inter(num) = EXTSMALL(index_wspeed,ISMALL,Index_wave_550) 
          EXTSMALL_inter(num) = Ext_554_small(index_wspeed,ISMALL) 
          CCN_inter(num)      =  CCNSMALL(index_wspeed,ISMALL)
!          EXTBIG_inter(num)   = EXTBIG(index_wspeed,IBIG,Index_wave_550)
           EXTBIG_inter(num)   = Ext_554_large(index_wspeed,IBIG)
          MOMENTSSMALL3_inter(num) =  MOMENTSSMALL(index_wspeed,ISMALL,3)
          MOMENTSBIG3_inter(num)   =  MOMENTSBIG(index_wspeed,IBIG,3)
          MOMENTSSMAL2_inter(num)   = MOMENTSSMALL(index_wspeed,ISMALL,2)
          MOMENTSBIG2_inter(num)   = MOMENTSBIG(index_wspeed,IBIG,2)
          MOMENTSSMAL4_inter(num)       = MOMENTSSMALL(index_wspeed,ISMALL,4)
          MOMENTSBIG4_inter(num)         = MOMENTSBIG(index_wspeed,IBIG,4)
          enddo
          
             y1=0.0
           CALL INTERP(Num,WSPEED,Windspeed_Lut,EXTSMALL_inter,Y1)  
           NSMALL=TAUX55SMALL/y1
 
            y1=0.0
           CALL INTERP(Num,WSPEED,Windspeed_Lut,CCN_inter,Y1) 
            CCN=NSMALL* y1
 
             y1=0.0
           CALL INTERP(Num,WSPEED,Windspeed_Lut,EXTBIG_inter,Y1) 
           NBIG=TAUX55BIG/y1
 
 
         IF(NSMALL .GT. 0.0 .OR.NBIG .GT.0.0 )then
           y1=0.0
           y2=0.0
           y3=0.0
           y4=0.0
           y5=0.0
           y6=0.0
           CALL INTERP(Num,WSPEED,Windspeed_Lut,MOMENTSSMALL3_inter,Y1) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut,MOMENTSBIG3_inter,Y2) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut,MOMENTSSMAL2_inter,Y3) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut,MOMENTSBIG2_inter,Y4) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut, MOMENTSSMAL4_inter,Y5) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut,MOMENTSBIG4_inter,Y6)
           EFFRADIUS=(NSMALL*y1+ NBIG*y2) / (NSMALL*y3+NBIG*y4) 
          
 
 
! NEW definition 9/24/2003.........

        MASS_CON_OCEAN=(NSMALL*y1+ NBIG*y2)
 
!                             *** compute Effective variance

         EFFVARIANCE=((NSMALL*y5+ NBIG*y6)/((EFFRADIUS**2)*(NSMALL*y3+NBIG*4))-1.)
      
     
     
 
          ELSE
         EFFRADIUS=0.0
         EFFVARIANCE=0.0
         MASS_CON_OCEAN=0.0
         ENDIF
         
         
         
          
!                           **** compute spectral optical thickness
        DO IWAV = 1,NWAV

           Num=0 
          DO index_wspeed = Indx_wspped1,Indx_wspped2
          Num=Num+1
          Windspeed_Lut(num)=Wind(index_wspeed)  
         ALBEDOSMAL_inter(num)    =ALBEDOSMALL(index_wspeed,ISMALL,IWAV)
         BACKSCTTSMAL_inter(num)   =BACKSCTTSMALL(index_wspeed,ISMALL,IWAV)
         ALBEDOBIG_inter(num)      = ALBEDOBIG(index_wspeed,IBIG,IWAV)
         BACKSCTTBIG_inter(num)    = BACKSCTTBIG(index_wspeed,IBIG,IWAV)
         ASSYMSMALL_inter(num)     =ASSYMSMALL(index_wspeed,ISMALL,IWAV)
         ASSYMBIG_inter(num)       =ASSYMBIG(index_wspeed,IBIG,IWAV)
         Enddo
           y1=0.0
           y2=0.0
           y3=0.0
           y4=0.0
           y5=0.0
           y6=0.0
           CALL INTERP(Num,WSPEED,Windspeed_Lut,ALBEDOSMAL_inter,Y1) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut,BACKSCTTSMAL_inter,Y2) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut,ALBEDOBIG_inter,Y3) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut,BACKSCTTBIG_inter,Y4) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut,ASSYMSMALL_inter,Y5) 
           CALL INTERP(Num,WSPEED,Windspeed_Lut,ASSYMBIG_inter,Y6)
!                           **** compute backscattering ratio
          IF( TAU_SMALL_COMPUTED(IWAV) .GT.0.0 .OR. &
                 TAU_BIG_COMPUTED(IWAV) .GT. 0.0 .and.&
               (TAU_SMALL_COMPUTED(IWAV)+TAU_BIG_COMPUTED(IWAV)).GT. 0.0) THEN
     
          BACKSCATT(IWAV)= &
        (y1*TAU_SMALL_COMPUTED(IWAV)*y2 + y3*TAU_BIG_COMPUTED(IWAV)*Y4)/  &
        (y1*TAU_SMALL_COMPUTED(IWAV) + y3 *TAU_BIG_COMPUTED(IWAV))
 

         ASSYM(IWAV)=(y1*TAU_SMALL_COMPUTED(IWAV)*y5  &
                 + y3*TAU_BIG_COMPUTED(IWAV)*y6)/ &
        (y1*TAU_SMALL_COMPUTED(IWAV) + y3*TAU_BIG_COMPUTED(IWAV))
          
         REF_FLUX(IWAV)=( &
                        TAU_SMALL_COMPUTED(IWAV)* &
                  ALBEDO_R_SMALL_final(IWAV,ISMALL)&
              +   TAU_BIG_COMPUTED(IWAV)*&
                  ALBEDO_R_BIG_final(IWAV,IBIG))/&
         ( TAU_SMALL_COMPUTED(IWAV)  + &
          TAU_BIG_COMPUTED(IWAV)) 
! subtract rayleigh to compute  and * mu0 because it is divided in LUT
! subtract rayleigh to compute  and * mu0 because it is divided in LUT
        REF_FLUX(IWAV)=(REF_FLUX(IWAV)- &
             ALBEDO_R_SMALL_tau(1,IWAV,1))*(COS(MTHET0*DTR)) 

 
       
!                     **** compute transmittance
         TRANS_FLUX(IWAV)=(TAU_SMALL_COMPUTED(IWAV)* &
                  ALBEDO_T_SMALL_final(IWAV,ISMALL) + & 
                  TAU_BIG_COMPUTED(IWAV)* &
                  ALBEDO_T_BIG_final(IWAV,IBIG))/&
                  ( TAU_SMALL_COMPUTED(IWAV)  + &
                   TAU_BIG_COMPUTED(IWAV))
       ELSE
            BACKSCATT(IWAV)=0.0
            ASSYM(IWAV)=0.0
            REF_FLUX(IWAV)=0.0
            TRANS_FLUX(IWAV)=0.0
       ENDIF
!
           ENDDO

       RETURN
       END


!*********************************************************************

      SUBROUTINE  HALVING(AX,BX,F,XMIN,FUNMIN)

!----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION: This subroutine findes the minimum of function by
!              NEWTON halving concept.
!
! !INPUT PARAMETERS:
!         AX      Intial condition for halving
!         BX      end condition for the intial
!         F       function whose minimum has to be determined.
!
! !OUTPUT PARAMETERS:
!          XMIN    minimum of function
!        FUNMIN    value of the minum function at xmin
!---------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      REAL  RADWAV
      REAL AX,BX,F,XMIN,X,X1,X2,FA,FB,FX,RMIN,FX1,FX2,ARR1(5),ARR2(5)
      real funmin
      integer iter,ij,IFLAG
      EXTERNAL  F

       IFLAG=0
       XMIN=0.0
       FUNMIN=0.0
       X=(AX+BX)/2.
       FA=F(AX)
       FB=F(BX)
       FX=F(X)
 
       DO 1 ITER = 1,10
       X1=(X+AX)/2.
       X2=(X+BX)/2.
       FX1=F(X1)
       FX2=F(X2) 
       ARR1(1)=AX
       ARR1(2)=BX
       ARR1(3)=X
       ARR1(4)=X1
       ARR1(5)=X2
       ARR2(1)=FA
       ARR2(2)=FB
       ARR2(3)=FX
       ARR2(4)=FX1
       ARR2(5)=FX2
       CALL shlsrt(5,ARR2,ARR1)
       RMIN=arr1(1) 
       if( rmin .le. x1) then
        IFLAG=IFLAG+1
            ax=ax
            bx=x
           fb=fx
            x=x1
           fx=fx1 
       endif

 

       if(abs(rmin-x).lt.abs(x)*0.000001 .and. rmin.gt.x1) then
        IFLAG=IFLAG+1
         ax=x1
         fa=fx1
         bx=x2
         fb=fx2
         x=x
         fx=fx 
         endif 
 
       if(rmin.ge.x2 .and.&
        abs(rmin-x).gt.abs(x)*0.000001 .and. rmin.gt.x1) then
        IFLAG=IFLAG+1
         ax=x
         fa=fx
         bx=bx
         fb=fb
         x=x2
         fx=fx2 
         endif
 1       continue
       IF( IFLAG .GT.1) THEN
       arr1(1)=ax
       arr1(2)=bx
       arr1(3)=x
       arr2(1)=fa
       arr2(2)=fb
       arr2(3)=fx
       CALL shlsrt(3,ARR2,ARR1)
       xmin=arr1(1)
       funmin=arr2(1)
       ELSE
       XMIN=999.0
       FUNMIN=9.99
       ENDIF 
       return
        end



!********************************************************************

       SUBROUTINE store_other(REF_FLUX,TRANS_FLUX,FUNMIN,EFFRADIUS,&
       BACKSCATT, ASSYM,NSOLUTION,CCN,SCAT_ANGLE_OCEAN,&
       ARRAY,NUM_ARRAY_ELEMENTS,NUMDATA_550,ISMALL,IBIG,&
      sd_W470,sd_W550,sd_W659,sd_w865,sd_W124,sd_W164,Sd_W213,CLD_FRAC,&
      TOT_cldfree_number,MASS_CON_OCEAN,WAVE)
      
      

!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine stores all the output varaibles
!             to be wrriten into HDF format.
!
!!INPUT PARAMETERS:
!     IDATA        Index of box number
!     IBIG         Index for big mode
!     ISMALL       Index for small mode.
!     NUMCASE      Number of large and small cases
!     NWAV         Number of wavelengths
!     FUNMIN       Minimm Error function betwwen derived and
!                  computed radiances.
!     EFFRADIUS    Effective radius weighed for large and small mode.
!     EFFVARIANCE  standard deviation for effective radius.
!     CCNSMALL     CCN for small.
!     ASSYM        Assymetry factor for small and large.
!     NSMALL       Ratio of taux55small to extnor(small mode).
!     NBIG         Ratio of taux55small to extnor(large mode).
!     NSOLUTION    Index for the number of solution.
!    RGSS         Radius for small particle
!     CCN          ccn for the solution.
!     SIGMAS       Standard deviation for small mode.
!     RGSB         Radius for small particle
!     SIGMAB       Standard deviation for small mode.
!     REF_FLUX     Reflected Flux
!     TRANS_FLUX   Transmitted Flux
! NUM_ARRAY_ELEMENTS    Number of array elements.
!
!        Input through 'mod04.inc' from common
!       TAUX55SMALL Optical thickness contribution small particles
!       TAUX55BIG   Optical thickness contribution large particles
!       TAUX55      Optical thickness contribution large particles+small mode
!       AN           Weight factor for large and small mode
!!OUTPUT PARAMETERS:
!    ARRAY      Array containing all output variables for all small and
!               large solutions 
!----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      INTEGER IWAV,NSOLUTION,IJ,IK
      INTEGER NUM_ARRAY_ELEMENTS,NUMDATA_550,ISMALL,IBIG
      INTEGER TOT_cldfree_number
       REAL    sd_W470,sd_W550,sd_W659,sd_w865,sd_W124,sd_W164,Sd_W213
      REAL REF_FLUX(NWAV),TRANS_FLUX(NWAV)
      REAL CCN,ASSYM(NWAV),BACKSCATT(NWAV),ALM
      REAL VAR
      REAL EFFRADIUS,FUNMIN,ERROR(10),WAVE(NWAV)
      REAL ARRAY(NUM_ARRAY_ELEMENTS,NUMCASES*NUMCASEB)
      REAL SCAT_ANGLE_OCEAN,CLD_FRAC,MASS_CON_OCEAN
             VAR=-999.0
! Write ASCII file for scientific data only
            DO IWAV=1,NWAV
              IF( IWAV .EQ.1)ALM=REFW470
              IF( IWAV .EQ.2)ALM=REFW550
              IF( IWAV .EQ.3)ALM=REFW659
              IF( IWAV .EQ.4)ALM=REFW865
              IF( IWAV .EQ.5)ALM=REFW124
              IF( IWAV .EQ.6)ALM=REFW164
              IF( IWAV .EQ.7)ALM=REFW213
              ERROR(IWAV)=(ALM-ALXWAV(IWAV))
            ENDDO
            DO IJ =1 ,NUM_ARRAY_ELEMENTS
!  Store  INDEX FOR SMALL PARTICLES
               IF( IJ .EQ.1) VAR=ISMALL
!  Store  INDEX FOR BIG PARTICLES( NUMCASES is added to get the index given in Docu.
              IF (IJ .EQ.77)VAR=IBIG+NUMCASES

               IF( IJ .EQ.2) VAR=SCAT_ANGLE_OCEAN
!   Store values for optical thicknesses for 7 wavelengths
               IF (IJ .EQ.3)VAR=TAU_COMPUTED(1)
               IF (IJ .EQ.4)VAR=TAU_COMPUTED(2)
               IF( IJ .EQ.5)VAR=TAU_COMPUTED(3)
               IF( IJ .EQ.6)VAR=TAU_COMPUTED(4)
               IF( IJ .EQ.7)VAR=TAU_COMPUTED(5)
                IF( IJ .EQ.8)VAR=TAU_COMPUTED(6)
                IF( IJ .EQ.9)VAR=TAU_COMPUTED(7)
! Store values for optical thicknesses for 7 wavelengths for small mode
                IF( IJ .EQ.10) VAR=TAU_SMALL_COMPUTED(1)
                IF( IJ .EQ.11) VAR=TAU_SMALL_COMPUTED(2)
                IF( IJ .EQ.12) VAR=TAU_SMALL_COMPUTED(3)
                IF( IJ .EQ.13) VAR=TAU_SMALL_COMPUTED(4)
                IF( IJ .EQ.14) VAR=TAU_SMALL_COMPUTED(5)
                IF( IJ .EQ.15) VAR=TAU_SMALL_COMPUTED(6)
                IF( IJ .EQ.16) VAR=TAU_SMALL_COMPUTED(7)
! Store values for optical thicknesses for 7 wavelengths for large mode
                IF( IJ .EQ.17) VAR=TAU_BIG_COMPUTED(1)
                IF( IJ .EQ.18) VAR=TAU_BIG_COMPUTED(2)
                IF( IJ .EQ.19) VAR=TAU_BIG_COMPUTED(3)
                IF( IJ .EQ.20) VAR=TAU_BIG_COMPUTED(4)
                IF( IJ .EQ.21) VAR=TAU_BIG_COMPUTED(5)
                IF( IJ .EQ.22) VAR=TAU_BIG_COMPUTED(6)
                IF( IJ .EQ.23) VAR=TAU_BIG_COMPUTED(7)
!          1/2006 change in units  same as in Land
           IF( IJ .EQ.24)VAR=(MASS_CON_OCEAN /1.0E+6)*(4./3.* PI)
!       CCN at 0.55 um
              IF( IJ .EQ.25)VAR=CCN
! Assymetry factor for all 7 wavelengths
               IF( IJ .EQ.26) VAR= ASSYM(1)
               IF( IJ .EQ.27) VAR= ASSYM(2)
               IF( IJ .EQ.28) VAR= ASSYM(3)
               IF( IJ .EQ.29) VAR= ASSYM(4)
               IF( IJ .EQ.30) VAR= ASSYM(5)
               IF( IJ .EQ.31) VAR= ASSYM(6)
               IF( IJ .EQ.32) VAR= ASSYM(7)
! Backscattering ration for all 7 wavelengths

               IF( IJ .EQ.33) VAR=BACKSCATT(1)
               IF( IJ .EQ.34) VAR=BACKSCATT(2)
               IF( IJ .EQ.35) VAR=BACKSCATT(3)
               IF( IJ .EQ.36) VAR=BACKSCATT(4)
               IF( IJ .EQ.37) VAR=BACKSCATT(5)
               IF( IJ .EQ.38) VAR=BACKSCATT(6)
               IF( IJ .EQ.39) VAR=BACKSCATT(7)
!  save 39,40 for Angstorm coeff1 and coeff2
! Changes 7/2009
         IF(IJ .EQ.40 .AND. TAU_COMPUTED(2) .GT.0.& 
                     .AND. TAU_COMPUTED(4) .GT.0.)THEN 
          VAR=(ALOG(TAU_COMPUTED(2)/TAU_COMPUTED(4)))/&
                        (ALOG(WAVE(4)/WAVE(2)))  
          ENDIF
               IF(IJ .EQ.41 .AND. TAU_COMPUTED(4) .GT.0. &
                     .AND. TAU_COMPUTED(7) .GT.0.)THEN 
          VAR=(ALOG(TAU_COMPUTED(4)/TAU_COMPUTED(7)))/ &
                        (ALOG(WAVE(7)/WAVE(4)))  
          ENDIF
          
 
 
!  Reflected flux for 7 wavelengths
               IF( IJ .EQ.42) VAR=REF_FLUX(1)
               IF( IJ .EQ.43) VAR=REF_FLUX(2)
               IF( IJ .EQ.44) VAR=REF_FLUX(3)
               IF( IJ .EQ.45) VAR=REF_FLUX(4)
               IF( IJ .EQ.46) VAR=REF_FLUX(5)
               IF( IJ .EQ.47) VAR=REF_FLUX(6)
               IF( IJ .EQ.48) VAR=REF_FLUX(7)
!  Transmitted  flux for 7 wavelengths
               IF( IJ .EQ.49) VAR=TRANS_FLUX(1)
               IF( IJ .EQ.50) VAR=TRANS_FLUX(2)
               IF( IJ .EQ.51) VAR=TRANS_FLUX(3)
               IF( IJ .EQ.52) VAR=TRANS_FLUX(4)
               IF( IJ .EQ.53) VAR=TRANS_FLUX(5)
               IF( IJ .EQ.54) VAR=TRANS_FLUX(6)
               IF( IJ .EQ.55) VAR=TRANS_FLUX(7)
! Least sqare error in percentage
              IF (IJ .EQ.56) VAR=FUNMIN*100.
! Small mode weighting factor
               IF (IJ .EQ.57)VAR=AN
!   save 58 for Optical_depth ratio-small
            IF(IJ .EQ.58 .and. TAU_COMPUTED(2) .GT. 0.0)THEN 
                   VAR=TAU_SMALL_COMPUTED(2)/TAU_COMPUTED(2)
            ENDIF
!   save 59 for cloud fraction  
 
         IF (IJ .EQ.59) VAR=CLD_FRAC 

!   store number of pixels used
               IF (IJ .EQ.60)VAR=NUMDATA_550
! store all 7 wavelength reflectance
               IF (IJ .EQ.61)VAR=REFW470
               IF (IJ .EQ.62)VAR=REFW550
               IF (IJ .EQ.63)VAR=REFW659
               IF (IJ .EQ.64)VAR=REFW865
               IF (IJ .EQ.65)VAR=REFW124
               IF (IJ .EQ.66)VAR=REFW164
               IF (IJ .EQ.67)VAR=REFW213
! store all 7 wavelength standard deviation for Reflectance
               IF (IJ .EQ.68)VAR= sd_W470
               IF (IJ .EQ.69)VAR= sd_w550
               IF (IJ .EQ.70)VAR= sd_W659
               IF (IJ .EQ.71)VAR= sd_w865
               IF (IJ .EQ.72)VAR= sd_W124
               IF (IJ .EQ.73)VAR= Sd_W164
               IF (IJ .EQ.74)VAR= Sd_W213
!  Store quality conrol
              IF (IJ .EQ.75)VAR=QCONTROL
!  Store effective radius
              IF (IJ .EQ.76)VAR=EFFRADIUS
               ARRAY(IJ,NSOLUTION)=VAR
             ENDDO

          RETURN
          END



!********************************************************************

      SUBROUTINE AVERAGE_output(ARRAY,NUM_ARRAY_ELEMENTS,AVE_ARRAY,&
                    NEW_ARRAY,Total_coarse_used)

!----------------------------------------------------------------------
!!F90
!
!!DESCRIPTION: This subroutine sorts according to minimum
!             error and averages all variables of output
!
!!INPUT PARAMETERS:
!     ARRAY              All output variables to be printed
!     NUM_ARRAY_ELEMENTS number of variables in output array
!     AVE_ARRAY          Array with all output varaibles for best
!                        and average solutions
!     NEW_ARRAY          Array used to hold the data
!!OUTPUT PARAMETERS:
!    AVE_ARRAY          Array with all output varaibles for best
!
! 
!----------------------------------------------------------------------

        IMPLICIT NONE
        SAVE

        INCLUDE 'mod04.inc'

            integer NUM_ARRAY_ELEMENTS
            REAL ARRAY(NUM_ARRAY_ELEMENTS,NUMCASES*NUMCASEB)
            REAL ARR1(NUMCASES*NUMCASEB)
            REAL AVE_ARRAY(NUM_ARRAY_ELEMENTS,NUM_solutions),AVE,SDEV
            REAL NEW_ARRAY(NUM_ARRAY_ELEMENTS,NUMCASEB*NUMCASES)
            REAL DATA(NUMCASEB*NUMCASES)
           INTEGER IJ,INDX(NUMCASES*NUMCASEB),ICASE,IDATA,IK
           INTEGER ARRAY_ELEMENTS,NUMPOINTS,INUM,n_number,Total_coarse_used
           REAL  Tot_epsilon,Tot_data,epsilon(3)
!    order according to smallest error in minimu function
!
           INUM=0

                 DO IJ=1,NUMCASES*Total_coarse_used
                 ARR1(IJ)=ARRAY(56,IJ)
                 ENDDO
              
! Call to get ascending order index

            CALL INDEXX(NUMCASES*Total_coarse_used,ARR1,INDX)

!  Save the best possible solution with minimum error in array
!  called  AVE_ARRAY
                DO ARRAY_ELEMENTS= 1,NUM_ARRAY_ELEMENTS
       AVE_ARRAY(ARRAY_ELEMENTS,1)=ARRAY(ARRAY_ELEMENTS,INDX(1))
                ENDDO

! If the minimum error is less tahn 3% average all the observations.
      DO IJ = 1,Total_coarse_used*NUMCASES
         IF(ARRAY(56,INDX(IJ)) .LE. Threshold_LSQ_Error) THEN
            INUM=INUM+1
            DO ARRAY_ELEMENTS= 1,NUM_ARRAY_ELEMENTS
        NEW_ARRAY(ARRAY_ELEMENTS,INUM)=ARRAY(ARRAY_ELEMENTS,INDX(IJ))
           ENDDO
         ENDIF
        ENDDO
          DO ARRAY_ELEMENTS= 1,NUM_ARRAY_ELEMENTS
              IF( INUM .GT.1) THEN
                    QCONTROL=0
                  DO NUMPOINTS=1,INUM
                    DATA(NUMPOINTS)=NEW_ARRAY(ARRAY_ELEMENTS,NUMPOINTS)
                   ENDDO
!save average solution in array called  AVE_ARRAY
                     CALL  ave_std(DATA,INUM,AVE,SDEV)
                    AVE_ARRAY(ARRAY_ELEMENTS,2)=AVE
! Solution number is not averaged,solution number is taken for the best solution
             IF(ARRAY_ELEMENTS.EQ.1)&
        AVE_ARRAY(ARRAY_ELEMENTS,2)=(ARRAY(ARRAY_ELEMENTS,INDX(1)))
             IF(ARRAY_ELEMENTS.EQ.77)&
        AVE_ARRAY(ARRAY_ELEMENTS,2)=(ARRAY(ARRAY_ELEMENTS,INDX(1)))

              ELSE
!   It should atleast have 3 observations
          IF(Total_coarse_used*NUMCASES .ge.3 ) then
! If there is only one observation average atleast 3 observations

                 QCONTROL=8
                 Tot_epsilon=0.0
                  Tot_data=0.0
               n_number=0
              DO NUMPOINTS=1,3
               if(ARRAY(56,INDX(NUMPOINTS)) .gt.0.0) then
                   n_number=n_number+1
             epsilon(n_number)=(1./((ARRAY(56,INDX(NUMPOINTS))**2)))
             Tot_epsilon=tot_epsilon+epsilon(n_number)
             Tot_data=Tot_data+(ARRAY(ARRAY_ELEMENTS,INDX(NUMPOINTS))*&
                         epsilon(n_number))
              endif
              enddo

! save average solution in array called  AVE_ARRAY
                    if(Tot_epsilon .gt.0) then
                    AVE= Tot_data/Tot_epsilon
                    else
                    AVE=0.0
                    endif
            AVE_ARRAY(ARRAY_ELEMENTS,2)=AVE
! Solution number is not averaged,solution number is taken for the best solution
             IF(ARRAY_ELEMENTS.EQ.1)&
        AVE_ARRAY(ARRAY_ELEMENTS,2)=(ARRAY(ARRAY_ELEMENTS,INDX(1)))
             IF(ARRAY_ELEMENTS.EQ.77) &
       AVE_ARRAY(ARRAY_ELEMENTS,2)=(ARRAY(ARRAY_ELEMENTS,INDX(1)))
! else  and endif   for atleast 3 observation else put the best solution
         ELSE
          AVE_ARRAY(ARRAY_ELEMENTS,2)= AVE_ARRAY(ARRAY_ELEMENTS,1)
        ENDIF
! endif for inum
        ENDIF
! ENDDO for Array_elements
          ENDDO
        
         RETURN
         END



!********************************************************************

        SUBROUTINE set_output(IDATA,AVE_ARRAY,NUM_ARRAY_ELEMENTS,&
      SDSTAU_best,SDSTAU_average,SDSTAUB_best,&
      SDSTAUB_average,SDSTAUS_best,SDSTAUS_average,&
       SDS_small_weighting,SDS_Least_error,&
      SDSASSY_best,SDSASSY_average,SDSBACK_best,SDSBACK_average,&
      SDS_effrad,SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,&
      SDS_TranF_average,SDS_sol_INDX_small,SDS_sol_INDX_large,&
      SDS_ccn,SDS_mass_conc,SDS_Tau_Land_Ocean_img,&
      SDS_angs_coeff1,SDS_angs_coeff2,SDS_AOT_model,&
      SDS_Tau_Land_Ocean,SDS_correc_small_weighting)

!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine sorts according to minimum
!            error and writes  all the output varaibles
!             in ASCII format for science data only.
!
!!INPUT PARAMETERS:
!AVE_ARRAY           Two element array  for best and average solution
! NUM_ARRAY_ELEMENTS  Number of arrays for output
! IDATA               Index of box number
!!OUTPUT PARAMETERS:
! SDSTAU_best     Optical thickness for best solution
! SDSTAUS_best    Optical thickness contribution small particles for best solution
! SDSTAUB_best    Optical thickness contribution large particles for best solution
! SDSTAU_average  Optical thickness for best solution
! SDSTAUS_average Optical thickness contribution small particles for best solution
! SDSTAUB_average Optical thickness contribution large particles for best solution
!SDS_Least_error    Minimm Error function betwwen derived and computed radiances
! SDS_small_weighting Weight factor for large and small mode
! SDS_sol_INDX_small       Index for solution number small particles
!SDS_sol_INDX_large       Index for solution number large particles
! SDSASSY_best       Assymetry Factor for best solution
! SDSASSY_average    Assymetry Factor for Average solution
! SDSBACK_best       Backscattering Ratio for best solution
! SDSBACK_average    Backscattering Ratio for Average solution
! SDS_effrad         Effective Radiance
! SDS_AOT_model Ration of optical Thickess small
! SDS_RefF_best      Reflected Flux Best solution
! SDS_RefF_average   Reflected Flux Average solution
! SDS_TranF_best     Transmitted Flux Best solution
! SDS_TranF_average  Transmitted Flux Average solution
! SDS_ccn            Cloud condensation Neculii
! SDS_mass_conc      Mass concentration
! SDS_angs_coeff1    Angstrom coeff1
! SDS_angs_coeff2    Angstrom coeff2
 

        IMPLICIT NONE
        SAVE

        INCLUDE 'mod04.inc'
        INCLUDE 'read_Sat_MODIS.inc'

            integer NUM_ARRAY_ELEMENTS,IDATA,IK,ICASE
            REAL AVE_ARRAY(NUM_ARRAY_ELEMENTS,NUM_solutions)
            INTEGER IJ,NSOLUTION
            REAL  SDS_ccn(NUMCELLS_B,NUM_solutions),&
             SDS_mass_conc(NUMCELLS_B,NUM_solutions)
           INTEGER ARRAY_ELEMENTS,NUMPOINTS,INUM
           Real  SDSTAU_best(NUMCELLS_B,NWAV),&
                  SDSTAU_average(NUMCELLS_B,NWAV),&
                 SDSTAUB_best(NUMCELLS_B,NWAV),&
                  SDSTAUB_average(NUMCELLS_B,NWAV),&
                  SDSTAUS_best(NUMCELLS_B,NWAV),&
                  SDSTAUS_average(NUMCELLS_B,NWAV),&
                  SDSASSY_best(NUMCELLS_B,NWAV),&
                  SDSASSY_average(NUMCELLS_B,NWAV),&
                  SDSBACK_best(NUMCELLS_B,NWAV),&
                  SDSBACK_average(NUMCELLS_B,NWAV),&
                  SDS_RefF_best(NUMCELLS_B,NWAV),&
                  SDS_RefF_average(NUMCELLS_B,NWAV),&
                  SDS_TranF_best(NUMCELLS_B,NWAV),&
                  SDS_TranF_average(NUMCELLS_B,NWAV),&
                  SDS_Tau_Land_Ocean(NUMCELLS_B),&
                 SDS_Tau_Land_Ocean_img(NUMCELLS_B)
        Real         &
                  SDS_small_weighting(NUMCELLS_B,NUM_solutions),&
                  SDS_correc_small_weighting(NUMCELLS_B),&
                  SDS_Least_error(NUMCELLS_B,NUM_solutions),&
                  SDS_effrad(NUMCELLS_B,NUM_solutions),& 
                  SDS_angs_coeff1(NUMCELLS_B,NUM_solutions),&
                  SDS_angs_coeff2(NUMCELLS_B,NUM_solutions),&
                 SDS_AOT_model(NUMCELLS_B,num_model_index)
      Integer     SDS_sol_INDX_small(NUMCELLS_B,NUM_solutions),&
                  SDS_sol_INDX_large(NUMCELLS_B,NUM_solutions) 
!write ascii output for best and average solutions

!            IJ=1
!        WRITE(36,100)idata,(AVE_ARRAY(IK,IJ),IK=1,NUM_ARRAY_ELEMENTS)
!             IJ = 2
!         WRITE(37,100)idata,(AVE_ARRAY(IK,IJ),IK=1,NUM_ARRAY_ELEMENTS)
!
! Store into array's for HDF file format for best case(Taua )
!
      SDSTAU_best(IDATA,1)  = AVE_ARRAY(3,1)
      SDSTAU_best(IDATA,2)  = AVE_ARRAY(4,1)
      SDSTAU_best(IDATA,3)  = AVE_ARRAY(5,1)
      SDSTAU_best(IDATA,4)  = AVE_ARRAY(6,1)
      SDSTAU_best(IDATA,5)  = AVE_ARRAY(7,1)
      SDSTAU_best(IDATA,6)  = AVE_ARRAY(8,1)
      SDSTAU_best(IDATA,7)  = AVE_ARRAY(9,1)
 
!Store into array's for HDF file format for average case(Taua )
 
      SDSTAU_average(IDATA,1)  = AVE_ARRAY(3,2)
      SDSTAU_average(IDATA,2)  = AVE_ARRAY(4,2)
      SDSTAU_average(IDATA,3)  = AVE_ARRAY(5,2)
      SDSTAU_average(IDATA,4)  = AVE_ARRAY(6,2)
      SDSTAU_average(IDATA,5)  = AVE_ARRAY(7,2)
      SDSTAU_average(IDATA,6)  = AVE_ARRAY(8,2)
      SDSTAU_average(IDATA,7)  = AVE_ARRAY(9,2)
 
!Store into array's for HDF file format for best case(Taua_small)
 
      SDSTAUS_best(IDATA,1) = AVE_ARRAY(10,1)
      SDSTAUS_best(IDATA,2) = AVE_ARRAY(11,1)
      SDSTAUS_best(IDATA,3) = AVE_ARRAY(12,1)
      SDSTAUS_best(IDATA,4) = AVE_ARRAY(13,1)
      SDSTAUS_best(IDATA,5) = AVE_ARRAY(14,1)
      SDSTAUS_best(IDATA,6) = AVE_ARRAY(15,1)
      SDSTAUS_best(IDATA,7) = AVE_ARRAY(16,1)
 
!Store into array's for HDF file format for average case(Taua_small )
 
      SDSTAUS_average(IDATA,1) = AVE_ARRAY(10,2)
      SDSTAUS_average(IDATA,2) = AVE_ARRAY(11,2)
      SDSTAUS_average(IDATA,3) = AVE_ARRAY(12,2)
      SDSTAUS_average(IDATA,4) = AVE_ARRAY(13,2)
      SDSTAUS_average(IDATA,5) = AVE_ARRAY(14,2)
      SDSTAUS_average(IDATA,6) = AVE_ARRAY(15,2)
      SDSTAUS_average(IDATA,7) = AVE_ARRAY(16,2)
!
!Store into array's for HDF file format for best case(Taua_large )
 
      SDSTAUB_best(IDATA,1) = AVE_ARRAY(17,1)
      SDSTAUB_best(IDATA,2) = AVE_ARRAY(18,1)
      SDSTAUB_best(IDATA,3) = AVE_ARRAY(19,1)
      SDSTAUB_best(IDATA,4) = AVE_ARRAY(20,1)
      SDSTAUB_best(IDATA,5) = AVE_ARRAY(21,1)
      SDSTAUB_best(IDATA,6) = AVE_ARRAY(22,1)
      SDSTAUB_best(IDATA,7) = AVE_ARRAY(23,1)
 
!Store into array's for HDF file format for average case(Taua_large)
 
      SDSTAUB_average(IDATA,1) = AVE_ARRAY(17,2)
      SDSTAUB_average(IDATA,2) = AVE_ARRAY(18,2)
      SDSTAUB_average(IDATA,3) = AVE_ARRAY(19,2)
      SDSTAUB_average(IDATA,4) = AVE_ARRAY(20,2)
      SDSTAUB_average(IDATA,5) = AVE_ARRAY(21,2)
      SDSTAUB_average(IDATA,6) = AVE_ARRAY(22,2)
      SDSTAUB_average(IDATA,7) = AVE_ARRAY(23,2)
 
!Store into array's for HDF file format for best case(Assymetry)
 
      SDSASSY_best(IDATA,1) = AVE_ARRAY(26,1)
      SDSASSY_best(IDATA,2) = AVE_ARRAY(27,1)
      SDSASSY_best(IDATA,3) = AVE_ARRAY(28,1)
      SDSASSY_best(IDATA,4) = AVE_ARRAY(29,1)
      SDSASSY_best(IDATA,5) = AVE_ARRAY(30,1)
      SDSASSY_best(IDATA,6) = AVE_ARRAY(31,1)
      SDSASSY_best(IDATA,7) = AVE_ARRAY(32,1)
 
!Store into array's for HDF file format for average case(Assymetry)
 
      SDSASSY_average(IDATA,1) = AVE_ARRAY(26,2)
      SDSASSY_average(IDATA,2) = AVE_ARRAY(27,2)
      SDSASSY_average(IDATA,3) = AVE_ARRAY(28,2)
      SDSASSY_average(IDATA,4) = AVE_ARRAY(29,2)
      SDSASSY_average(IDATA,5) = AVE_ARRAY(30,2)
      SDSASSY_average(IDATA,6) = AVE_ARRAY(31,2)
      SDSASSY_average(IDATA,7) = AVE_ARRAY(32,2)

 
!Store into array's for HDF file format for best case(Backscattering)
 
      SDSBACK_best(IDATA,1) = AVE_ARRAY(33,1)
      SDSBACK_best(IDATA,2) = AVE_ARRAY(34,1)
      SDSBACK_best(IDATA,3) = AVE_ARRAY(35,1)
      SDSBACK_best(IDATA,4) = AVE_ARRAY(36,1)
      SDSBACK_best(IDATA,5) = AVE_ARRAY(37,1)
       SDSBACK_best(IDATA,6) = AVE_ARRAY(38,1)
       SDSBACK_best(IDATA,7) = AVE_ARRAY(39,1)
 
!Store into array's for HDF file format for average case(Backscattering)
 
      SDSBACK_average(IDATA,1) = AVE_ARRAY(33,2)
      SDSBACK_average(IDATA,2) = AVE_ARRAY(34,2)
      SDSBACK_average(IDATA,3) = AVE_ARRAY(35,2)
      SDSBACK_average(IDATA,4) = AVE_ARRAY(36,2)
      SDSBACK_average(IDATA,5) = AVE_ARRAY(37,2)
      SDSBACK_average(IDATA,6) = AVE_ARRAY(38,2)
      SDSBACK_average(IDATA,7) = AVE_ARRAY(39,2)
 
!Store into array's for HDF file format for best case(Reflected_flux)
 
      SDS_RefF_best(IDATA,1) = AVE_ARRAY(42,1)
      SDS_RefF_best(IDATA,2) = AVE_ARRAY(43,1)
      SDS_RefF_best(IDATA,3) = AVE_ARRAY(44,1)
      SDS_RefF_best(IDATA,4) = AVE_ARRAY(45,1)
      SDS_RefF_best(IDATA,5) = AVE_ARRAY(46,1)
      SDS_RefF_best(IDATA,6) = AVE_ARRAY(47,1)
      SDS_RefF_best(IDATA,7) = AVE_ARRAY(48,1)
 
!Store into array's for HDF file format for average case(Reflected_flux)
 
      SDS_RefF_average(IDATA,1) = AVE_ARRAY(42,2)
      SDS_RefF_average(IDATA,2) = AVE_ARRAY(43,2)
      SDS_RefF_average(IDATA,3) = AVE_ARRAY(44,2)
      SDS_RefF_average(IDATA,4) = AVE_ARRAY(45,2)
      SDS_RefF_average(IDATA,5) = AVE_ARRAY(46,2)
      SDS_RefF_average(IDATA,6) = AVE_ARRAY(47,2)
      SDS_RefF_average(IDATA,7) = AVE_ARRAY(48,2)
 
!Store into array's for HDF file format for best case(Transmitted_flux)
 
      SDS_TranF_best(IDATA,1) = AVE_ARRAY(49,1)
      SDS_TranF_best(IDATA,2) = AVE_ARRAY(50,1)
      SDS_TranF_best(IDATA,3) = AVE_ARRAY(51,1)
      SDS_TranF_best(IDATA,4) = AVE_ARRAY(52,1)
      SDS_TranF_best(IDATA,5) = AVE_ARRAY(53,1)
      SDS_TranF_best(IDATA,6) = AVE_ARRAY(54,1)
      SDS_TranF_best(IDATA,7) = AVE_ARRAY(55,1)
 
!Store into array's for HDF file format for average case(Transmitted_flux)
 
      SDS_TranF_average(IDATA,1) = AVE_ARRAY(49,2)
      SDS_TranF_average(IDATA,2) = AVE_ARRAY(50,2)
      SDS_TranF_average(IDATA,3) = AVE_ARRAY(51,2)
      SDS_TranF_average(IDATA,4) = AVE_ARRAY(52,2)
      SDS_TranF_average(IDATA,5) = AVE_ARRAY(53,2)
      SDS_TranF_average(IDATA,6) = AVE_ARRAY(54,2)
      SDS_TranF_average(IDATA,7) = AVE_ARRAY(55,2)
 
      DO ICASE = 1, NUM_solutions

      SDS_mass_conc(IDATA,ICASE)=AVE_ARRAY(24,ICASE) 
      SDS_effrad(IDATA,ICASE)=AVE_ARRAY(76,ICASE)
      SDS_ccn(IDATA,ICASE)=AVE_ARRAY(25,ICASE) 
      SDS_angs_coeff1(IDATA,ICASE)=AVE_ARRAY(40,ICASE) 
      SDS_angs_coeff2(IDATA,ICASE)=AVE_ARRAY(41,ICASE)
      SDS_Least_error(IDATA,ICASE) = (AVE_ARRAY(56,ICASE)/100)
      SDS_small_weighting(IDATA,ICASE) = AVE_ARRAY(57,ICASE) 
      SDS_sol_INDX_small(IDATA,ICASE)=AVE_ARRAY(1,ICASE)
      SDS_sol_INDX_large(IDATA,ICASE)=AVE_ARRAY(77,ICASE)

       ENDDO
        DO ICASE = 1, num_model_index
!  fill every index of models with zero
         SDS_AOT_model(IDATA,ICASE)=0.0
!  fill only small and large index with optical depth of average values for small and large mode at 0.55 um
        SDS_AOT_model(IDATA,SDS_sol_INDX_small(IDATA,2))=&
            AVE_ARRAY(11,2)
        SDS_AOT_model(IDATA,SDS_sol_INDX_large(IDATA,2))=&
            AVE_ARRAY(18,2)
       enddo
!        SDS_correc_small_weighting(IDATA)=SDS_small_weighting(IDATA,2)
100      FORMAT(i4,10(10f10.4,/))

         RETURN
         END



!********************************************************************

       subroutine shlsrt(n,arr1,arr2)

!-----------------------------------------------------------------------
! !F90       
! !DESCRIPTION: From numerical recepies Sorts arr1ay 'arr1' of
!              length 'n' into ascending order ('arr1' and arr2
!              is replaced on OUTPUT by its sorted rearr1angement)
!
! !INPUT PARAMETERS:
!                 N         length of ARR1
!                ARR1       Unsorted array ARR1
!               ARR2       Unsorted array ARR2
!
!!OUTPUT PARAMETERS:
!                ARR1       Sorted array ARR1
!                ARR2       Sorted array ARR2
!
! 
!-----------------------------------------------------------------------

      implicit none
       save
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



!********************************************************************

      SUBROUTINE INDEXX( N, ARR, INDX )

!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine findes the  index of array ARR
!
!!INPUT PARAMETERS:
!         N     Number of points in array ARR
!        ARR     Data Array
!!OUTPUT PARAMETERS:
!       INDX    INDX returns index of array ARR arranged in ascending order 
!----------------------------------------------------------------------

      IMPLICIT NONE

!       Numerical Recipes ver. 2.06
!     .. Parameters ..
      INTEGER   M, NSTACK
      PARAMETER ( M = 7, NSTACK = 50 )
!     ..
!     .. Scalar Arguments ..
      INTEGER   N
!     ..
!     .. Array Arguments ..
      INTEGER   INDX( N )
      REAL      ARR( N )
!     .. Local Scalars ..
      INTEGER   I, INDXT, IR, ITEMP, J, JSTACK, K, L
      REAL      A
!     ..
!     .. Local Arrays ..
      INTEGER   ISTACK( NSTACK )
      DO 10 J = 1, N
         INDX( J ) = J
   10 CONTINUE
      JSTACK = 0
      L      = 1
      IR     = N

   20 CONTINUE
      IF( IR - L .LT. M ) THEN
         DO 50 J = L + 1, IR
            INDXT  = INDX( J )
            A      = ARR( INDXT )
            DO 30  I = J - 1, L, - 1
               IF( ARR( INDX(I) ) .LE. A ) GO TO 40
               INDX( I + 1 ) = INDX( I )
   30       CONTINUE
            I = L - 1
   40       CONTINUE
            INDX( I + 1 ) = INDXT
   50    CONTINUE
         IF( JSTACK.EQ.0 ) RETURN
         IR     = ISTACK( JSTACK )
         L      = ISTACK( JSTACK - 1 )
         JSTACK = JSTACK - 2
      ELSE
         K      = ( L + IR ) / 2
         ITEMP  = INDX( K )
         INDX( K ) = INDX( L + 1 )
         INDX( L + 1 ) = ITEMP
         IF( ARR( INDX(L) ) .GT. ARR( INDX(IR) ) ) THEN
            ITEMP  = INDX( L )
            INDX( L ) = INDX( IR )
            INDX( IR ) = ITEMP
         END IF
         IF( ARR( INDX(L + 1) ) .GT. ARR( INDX(IR) ) ) THEN
            ITEMP  = INDX( L + 1 )
            INDX( L + 1 ) = INDX( IR )
            INDX( IR ) = ITEMP
         END IF
         IF( ARR( INDX(L) ) .GT. ARR( INDX(L + 1) ) ) THEN
            ITEMP  = INDX( L )
            INDX( L ) = INDX( L + 1 )
            INDX( L + 1 ) = ITEMP
         END IF
         I      = L + 1
         J      = IR
         INDXT  = INDX( L + 1 )
         A      = ARR( INDXT )
   60    CONTINUE
         I = I + 1
         IF( ARR( INDX(I) ) .LT. A ) GO TO 60
   70    CONTINUE
         J = J - 1
         IF( ARR( INDX(J) ) .GT. A ) GO TO 70
         IF( J.LT.I ) GO TO 80
         ITEMP  = INDX( I )
         INDX( I ) = INDX( J )
         INDX( J ) = ITEMP
         GO TO 60
   80    CONTINUE
         INDX( L + 1 ) = INDX( J )
         INDX( J ) = INDXT
         JSTACK = JSTACK + 2
         IF( JSTACK.GT.NSTACK )STOP 'NSTACK too small in INDEXX'
         IF( IR - I + 1 .GE. J - L ) THEN
            ISTACK( JSTACK ) = IR
            ISTACK( JSTACK - 1 ) = I
            IR  = J - 1
         ELSE
            ISTACK( JSTACK ) = J - 1
            ISTACK( JSTACK - 1 ) = L
            L = I
         END IF
      END IF

      GO TO 20
      END



!*****************************************************************

      SUBROUTINE ave_std(DATA,N,AVE,SDEV)

!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION:  This subroutine is borrowed from "numerical Recipes"
!               page 458.  This subroutine returns the values of mean
!               average and standard deviation
!
!!INPUT PARAMETERS:
!    INTEGER   N           Number of elements in input array, DATA
!     REAL      DATA(N)     Array of N elements
!
!!OUTPUT PARAMETERS:
!     REAL      AVE         Arithmetic average of the N values of the
!                           array DATA
!     REAL      SDEV        Standard Deviation of the N values of the
!                          array DATA.
!---------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER J,N
      REAL DATA(N),S,P,AVE,VAR,SDEV

! This line was commented out by SDST vlin
!      IF(N.LE.1)PAUSE 'N must be at least 2'
      S=0.
      DO 11 J=1,N
        S=S+DATA(J)
11    CONTINUE
      AVE=S/N
      VAR=0.0
      DO 12 J=1,N
        S=DATA(J)-AVE
        P=S*S
        VAR=VAR+P
12    CONTINUE
        IF(VAR.GT.0)THEN
           VAR=VAR/(N-1)
           SDEV=SQRT(VAR)
        ELSE
          SDEV=0.0
      ENDIF

       RETURN
       END



!*****************************************************************

      SUBROUTINE COMPUTE_SCATTANGLE_OCEAN(MTHET0,MTHET,MPHI,IDATA,&
                                    SCAT_ANGLE_OCEAN)

!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION:  This subroutine computes scattering angle from modis measured
!               geometry.
!
!!INPUT PARAMETERS:
!MTHET0     Measured MODIS solar Zenith angle
!MTHET      Measured MODIS view angle
!MPHI       Measured MODIS Azimuthal angle
!IDATA      Index of 120* 10 bin in a swath
!!OUTPUT PARAMETERS:
!         SCAT_ANGLE_OCEAN Scattering angle 
!----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

       INTEGER IDATA
       REAL SCAT_ANGLE_OCEAN
       REAL  MTHET0,MTHET,MPHI

      SCAT_ANGLE_OCEAN=0.0
!
!                 *** Change modis geometry to scattering angle.
!
       SCAT_ANGLE_OCEAN=-(COS(MTHET0*DTR))*(COS(MTHET*DTR))&
                        -((SIN(MTHET0*DTR))*(SIN(MTHET*DTR))&
                        *( COS(MPHI*DTR)))
       SCAT_ANGLE_OCEAN= (ACOS(SCAT_ANGLE_OCEAN))*RTD

             RETURN
             END 

!*****************************************************************

        SUBROUTINE GLINT_IKM(START_1KM,END_1KM,Sunglint_Flag,NOGLINT)

!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION:  This subroutine computes total number of glint free pixels
!              in a 10*10 box
!
!!INPUT PARAMETERS:
!     START_1KM    starting index for 10 * 10 box
!       END_1KM    ending  index for  10 * 10 box
! Sunglint_Flag   sungleint flag from cloud mask
!
!!OUTPUT PARAMETERS:
!
!        NOGLINT tatal number of Glint_free pixels
! 

        IMPLICIT NONE
        SAVE

        INCLUDE 'mod04.inc'
        INCLUDE 'read_Sat_MODIS.inc'

        INTEGER IYY,IXX,START_1KM,END_1KM,NOGLINT
        INTEGER Sunglint_Flag(ISWATH_B,ILINE)

        NOGLINT=0
        DO   IYY = 1,IGRIDY
        DO   IXX =START_1KM,END_1KM
             IF(Sunglint_Flag(IXX,IYY).EQ.1)THEN
                  NOGLINT=NOGLINT+1
              ENDIF
        ENDDO
       ENDDO

       RETURN
       END



!*****************************************************************  

       SUBROUTINE  AVERAGE_to_500meter(W659_SYN,W865_SYN,W470_SYN,&
       W550_SYN,W124_SYN,W164_SYN,W213_SYN,W412_SYN,W443_SYN,&
       W8p5_Temp,W1100_Temp,ref_interm,NCLDMSK_syn,&
        NUMDATA_659,NUMDATA_865,NUMDATA_470,NUMDATA_550,&
       NUMDATA_124,NUMDATA_164,NUMDATA_213,NUMDATA_354,&
       NUMDATA_388,START_1KM,&
       END_1KM,VFlag_w659,VFlag_w865,VFlag_w470,VFlag_w550,VFlag_w124,&
       VFlag_w164,VFlag_w213,vFlag_w354,vFlag_w388,&
       VFlag_w869o,cloudy,iscan, High_Cloud_Flag_500,W138_SYN,Ref_ray,Wave,&
       CLD_FRAC,savecldmask,idata,SD_for_Dust,Dust_flag_10KM,&
       W354_SYN,W388_SYN)
      
!     ----------------------------------------------------------------------
!     F90      
!     
!     DES!     RIPTION: This subroutine averages the refle!     tan!     es for ea!     h
!              10*10 pixel  square and finds standdard deviation for
!               averaged Refle!     tan!     es
!     
!     INPUT PARAMETERS:
!     IDATA        Index of bin number in a swath
!     
!        START_1KM     Begining index of the swath in 10*10 km BOX
!        END_1KM      Ending   Index of the swath in 10*10 km BOX
!     Refle!     tan!     es of wavelengths used modis,!     h1-!     h7  identified by W*_SYN
!     Extra wavelengths from o!     ean !     olor !     h9,!     h12,!     h13 and !     h16 will be
!     used for substituting if !     h1 - !     h7 do not have valid data.
!     These !     hannels are identified by W*o_SYN.
!       W659_SYN
!       W865_SYN
!       W470_SYN
!       W550_SYN
!       W124_SYN
!       W164_SYN
!       W213_SYN
!       W443o_SYN
!       W551o_SYN
!       W667o_SYN
!       W869o_SYN
!     OUTPUT PARAMETERS:
!     Refle!     tan!     es averaged for 500 meter resolution
!     wavelengths used modis,!     h1-!     h7  identified by ref_interm_w*
!     Extra wavelengths from o!     ean !     olor !     h9,!     h12,!     h13 and !     h16 will be
!     used for substituting if !     h1 - !     h7 do not have valid data.
!     These !     hannels are identified by ref_interm_w*o
!       ref_interm_w659
!       ref_interm_w865
!       ref_interm_w470
!       ref_interm_w550
!       ref_interm_w124
!       ref_interm_w164
!       ref_interm_w213
!       ref_interm_w443o
!       ref_interm_w551o
!       ref_interm_w667o
!       ref_interm_w869o
!      Number of data points in ea!     h 500 meter box for all wavelengths used
!      NUMDATA_659
!      NUMDATA_865
!      NUMDATA_470
!      NUMDATA_550
!      NUMDATA_124
!      NUMDATA_164
!      NUMDATA_213
!      NUMDATA_443o
!      NUMDATA_551o
!      NUMDATA_667o
!      NUMDATA_869o
!     Total number of valid data  500 meter resolution
!     wavelengths used modis,!     h1-!     h7  identified by VFlag_w*
!     Extra wavelengths from o!     ean !     olor !     h9,!     h12,!     h13 and !     h16 will be
!     used for substituting if !     h1 - !     h7 do not have valid data.
!     These !     hannels are identified by VFlag_w*o
!     VFlag_w659
!     VFlag_w865
!     VFlag_w470
!     VFlag_w550
!     VFlag_w124
!     VFlag_w164
!     VFlag_w213
!     VFlag_w443o
!     VFlag_w551o
!     VFlag_w667o
!     VFlag_w869o 
!     ---------------------------------------------------------------------
       IMPLICIT NONE
       SAVE

       INCLUDE 'mod04.inc' 
       INCLUDE 'read_Sat_MODIS.inc'
 
       REAL W659_SYN(ISWATH_B,ILINE),W865_SYN(ISWATH_B,ILINE),&
          W470_SYN(ISWATH_B,ILINE),W550_SYN(ISWATH_B,ILINE),&
           W124_SYN(ISWATH_B,ILINE),W164_SYN(ISWATH_B,ILINE),&
           W213_SYN(ISWATH_B,ILINE),W138_SYN(ISWATH_B,ILINE),&
           W412_SYN(ISWATH_B,ILINE),W443_SYN(ISWATH_B,ILINE),&
           W8p5_Temp(ISWATH_B,ILINE),W1100_Temp(ISWATH_B,ILINE),&
           SD_for_Dust(ISWATH_B,ILINE),&
            W354_SYN(ISWATH_B,ILINE),W388_SYN(ISWATH_B,ILINE)
      REAL arefw659,arefw865,Ref_ray(NWAV),wave(Nwav)
      REAL ref_interm(NWAV+2,2*IGRIDX*2*IGRIDY),CLD_FRAC
      INTEGER NCLDMSK_SYN(ISWATH_B*2,ILINE*2),Dust_flag
      INTEGER START_1KM,END_1KM,JJ,II
      INTEGER   numaver1,numaver2,numaver3,numaver4
      INTEGER vFlag_w659,vFlag_w865,vFlag_w470,NumDust,&
       vFlag_w550,vFlag_w124,vFlag_w164,vFlag_w213,vFlag_w412,&
       vFlag_w354,vFlag_w388, vFlag_w869o
       INTEGER NUMDATA_659,NUMDATA_865,NUMDATA_470,&
       NUMDATA_550,NUMDATA_124,NUMDATA_164,NUMDATA_213,NUMDATA_412,&
       NUMDATA_354,NUMDATA_388
      INTEGER JMASK,IMASK,IBLUE,JBLUE,IXX,IYY,IX,IY,IIJ,IIK,cloudy,iscan
      INTEGER  High_Cloud_Flag_500(ISWATH_B,ILINE),idata,Dust_flag_10KM(ISWATH_B)
      INTEGER  sed_mask(ISWATH_B,ILINE),savecldmask(ISWATH_B*2,ILINE*2)
      Integer num_cloudy_pixels, Total_pixels,IWAVE,jy,jx,clear_pixels
      character (len =10) :: Sat_flag
       CLD_FRAC=-99
       cloudy=0
      Qcontrol_cirrus=0
       aREFW659=0.0
       aREFW865=0.0
        NUMDATA_470=0
       NUMDATA_550=0
        NUMDATA_659=0
        NUMDATA_865=0
         NUMDATA_124=0
       NUMDATA_164=0
       NUMDATA_213=0
       NUMDATA_354=0
       NUMDATA_388=0 
         numaver1=0
         numaver2=0
        vFlag_w470=0
        vFlag_w550=0
        vFlag_w659=0
        vFlag_w865=0
        vFlag_w124=0
        vFlag_w164=0
        vFlag_w213=0  
        vFlag_w354=0 
        vFlag_w388=0   
        
         Do IWAVE = 1,NWAV+2
         Do IXX = 1,(2*IGRIDX*2*IGRIDY)   
         ref_interm(IWAVE,IXX) =-99999
         enddo
         enddo 
        
         num_cloudy_pixels=0
         Total_pixels=0
         NumDust =0
        DO   IYY = 1,IGRIDY 
            IY=2*IYY-1
        DO  IXX=START_1KM,END_1KM   
            IX=2*IXX-1
            clear_pixels    = 0
            dust_flag  = 0
            
! Since cloud mask has 0 and 1 values finding product of
! rellectance of channels and cldmask helps to reduce CPU time than
! using  IF_ELSE_THEN statement in cloud mask to determine clear pixels.
              
 
            
!  Sediment Mask  
              call make_sediment(W470_syn,W550_SYN,&
          W124_SYN,W164_SYN,W213_SYN,sed_mask,IXX,IYY,wave,Sat_flag)
            do  jy = IY,2*IYY
            do  jx = IX,2*IXX  
! lower the quality for retrivals 
        if( W124_syn(IXX,IYY) .gt.0) then
          IF (NCLDMSK_SYN(jx,jy) .eq.1 .and.&
        W659_syn(IXX,IYY) .gt. (1.50*Ref_ray(3)).and. &
       ((W138_syn(IXX,IYY)/W124_syn(IXX,IYY)) .gt.0.10 .and.&
       (W138_syn(IXX,IYY)/W124_syn(IXX,IYY)) .lT.0.30).and.&
       (W138_syn(IXX,IYY) .gt.0.01 .and.W138_syn(IXX,IYY).le.0.03))Qcontrol_cirrus=1  
        ENDIF   
      IF(savecldmask(jx,jy) .EQ.0)num_cloudy_pixels=num_cloudy_pixels+1
         Total_pixels=Total_pixels+1  
       IF(NCLDMSK_SYN(jx,jy)*sed_mask (IXX,IYY) .GT. 0)clear_pixels = clear_pixels+1 
! enddo for 500meter          
          enddo
         enddo  
!*********************************************************************************************** 
! Applying Dust Flag ( Yaping  before cloudmask and sediment mask)
!***********************************************************************************************          
!          Call Dust_Detection_Ocean(IXX,IYY,SD_for_Dust, &
!          W659_syn,W865_syn,W470_syn,W550_SYN,W213_SYN,W412_SYN,W8p5_Temp,&
!          W1100_Temp,dust_flag) 
!******************************************************"#F2F2F2"********************************* 
! Applying Sediment mask to Dust Flag ( Yaping  before cloudmask )
!***********************************************************************************************      
!          if(dust_flag * sed_mask(IXX,IYY) .GT. 0)NumDust = NumDust +1  
          
!*********************************************************************************************** 
! Applying Dust Flag ( Yaping  after cloudmask and sediment mask)
!*********************************************************************************************** 
         if(clear_pixels .ge. 4) then
           Call Dust_Detection_Ocean(IXX,IYY,SD_for_Dust, &
            W659_syn,W865_syn,W470_syn,W550_SYN,W213_SYN,W412_SYN,W8p5_Temp,&
            W1100_Temp,dust_flag) 
          if(dust_flag .eq.1)NumDust = NumDust +1  
         ENDIF
         if( clear_pixels .ge. 4)   then
!     WAVE 0. 659 Um 
          NUMDATA_659= NUMDATA_659+1
         ref_interm(3,NUMDATA_659)= W659_syn(IXX,IYY) 
 
!   channel 3-channel 7 and cloud mask 
!     WAVE 0.869 Um  
     
              NUMDATA_865 =NUMDATA_865+1
        ref_interm(4,NUMDATA_865)= W865_syn(IXX,IYY) 
         

  
!    Foe clear pixel find product all Reflectances for
!    channel 3-channel 7 and cloud mask
!
!     WAVE 0.470 UM
!
          NUMDATA_470= NUMDATA_470+1
         ref_interm(1,NUMDATA_470)= W470_syn(IXX,IYY) 
!     WAVE 0.550 UM
!
           
   
         vFlag_W550=vFlag_W550+1
          NUMDATA_550= NUMDATA_550+1
         ref_interm(2,NUMDATA_550)= W550_syn(IXX,IYY) 
           
        
!     WAVE 1.24 UM
 
 
         vFlag_W124=vFlag_W124+1
         NUMDATA_124=NUMDATA_124+1
         ref_interm(5,NUMDATA_124)= W124_syn(IXX,IYY) 

!     WAVE 1.64 UM
 
 
          vFlag_W164=vFlag_W164+1
          NUMDATA_164=NUMDATA_164+1
          ref_interm(6,NUMDATA_164)=W164_syn(IXX,IYY)
!    WAVE 2.13 UM 
 
         vFlag_W213=vFlag_W213+1
         NUMDATA_213= NUMDATA_213+1
         ref_interm(7,NUMDATA_213)=W213_syn(IXX,IYY)
  
   !  WAVE 354 UM 
 
         vFlag_W354=vFlag_W354+1
         NUMDATA_354= NUMDATA_354+1
         ref_interm(8,NUMDATA_354)=W354_syn(IXX,IYY)
         
  !  WAVE 388nUM 
 
         vFlag_W388=vFlag_W388+1
         NUMDATA_388= NUMDATA_388+1
         ref_interm(9,NUMDATA_388)=W388_syn(IXX,IYY)    
         
! Endif for  clear_pixels      
          ENDIF  
!      ENDDO's for 1km pixels for 10*10 box
         ENDDO
         ENDDO 
         
        if( NumDust .GE. 3) then
        Dust_flag_10KM(Idata) = 1
        Else
        Dust_flag_10KM(Idata) = 0
        ENDIF
          if( Total_pixels .gt.0)&  
         CLD_FRAC=real(num_cloudy_pixels)/Real(Total_pixels)   
         RETURN
         END
      

!*****************************************************************

      SUBROUTINE INTsolarzenith_albedo_tran(THET0,ALBEDO_R_SMALL,&
                 ALBEDO_R_BIG,ALBEDO_T_SMALL,ALBEDO_T_BIG,MTHET0,&
                 ALBEDO_R_SMALL_tau,ALBEDO_R_BIG_tau,ALBEDO_T_SMALL_tau,&
                 ALBEDO_T_BIG_tau,KSZAM1,KSZAP1,Indx_wspped1,Indx_wspped2,&
                 WSPEED,Wind) 

!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine interpolates albedo and transmission on
!             solar zenith anfgle
!
!!INPUT PARAMETERS:
! THET0          Solar Zenith angle
!  Albedo for monodirectionally incident radiation Small and large
! ALBEDO_R_SMALL
! ALBEDO_R_BIG
! Transmission for monodirectionally incident radiation Small and large
! ALBEDO_T_SMALL
! ALBEDO_T_BIG
! measured solar zenith angle
! MTHET0
!!OUTPUT PARAMETERS:
! Interpolated Albedo for monodirectionally incident radiation Small and large
! ALBEDO_R_SMALL_tau
!  ALBEDO_R_BIG_tau
! InterpolatedTransmission for monodirectionally incident radiation Small and
! large
!ALBEDO_T_SMALL_tau
! ALBEDO_T_BIG_tau) 
!----------------------------------------------------------------------
       USE Linear_interpolation 
      IMPLICIT NONE

      INCLUDE 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'

      REAL ALBEDO_R_SMALL(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
      REAL ALBEDO_R_BIG(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
      REAL ALBEDO_T_BIG(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_R_SMALL_tau(NTAU,NWAV,NUMCASES)
      REAL ALBEDO_R_BIG_tau(NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL_tau(NTAU,NWAV,NUMCASES)
      REAL ALBEDO_T_BIG_tau(NTAU,NWAV,NUMCASEB)
      Real ALBEDO_R_int(NTAU,NWAV,(NUMCASEB+NUMCASES),Lut_indx)
      Real ALBEDO_T_int(NTAU,NWAV,(NUMCASEB+NUMCASES),Lut_indx)
      INTEGER IJ,ITH0,IWAV,LOPT,ITAU,ICASE,index_wspeed,num
      REAL  THET0(NTH0),wind(Lut_indx),WSPEED
      REAL  MTHET0,Windspeed_Lut(Lut_indx)
      REAL  X(100),Y(100),YY1(100),Y1
      REAL  XX(100),YY(100),YYY(100)
      INTEGER KSZAM1,KSZAP1,LL, Indx_wspped1,Indx_wspped2
      SAVE
! LOOP IS AROUND THET0,NEW FILE FOR EACH THETA0
           
          DO IJ = 1,100
          X(IJ)=0.0
          Y(IJ)=0.0
          YY1(IJ)=0.0
          XX(Ij)=0
          YY(Ij)=0
          YYY(Ij)=0
          ENDDO
          Num=0
          do index_wspeed = Indx_wspped1,Indx_wspped2
          Num=Num+1
          Windspeed_Lut(num)=Wind(index_wspeed) 
          Y1=0.0
!    INTEROOLATE FOR SMALL CASES
          DO 5  ICASE=1,NUMCASES
          DO 10 IWAV=1,NWAV
          DO 20 ITAU=1,NTAU
          LL=0
          DO 30 ITH0 = KSZAM1,KSZAP1
             LL=LL+1
            X(LL)=THET0(ITH0)
            Y(LL)=ALBEDO_R_SMALL(index_wspeed,ITH0,ITAU,IWAV,ICASE)
            YY1(LL)=ALBEDO_T_SMALL(index_wspeed,ITH0,ITAU,IWAV,ICASE)
 30        CONTINUE
             y1=0.0
            CALL INTERP(LL,MTHET0,X,Y,Y1)
            ALBEDO_R_int(ITAU,IWAV,ICASE,num)=y1 
             y1=0.0
             CALL INTERP(LL,MTHET0,X,YY1,Y1)
            ALBEDO_T_int(ITAU,IWAV,ICASE,num)=y1 
  20       CONTINUE
  10       CONTINUE
   5       CONTINUE
           Enddo
            
          DO   ICASE=1,NUMCASES
          DO   IWAV=1,NWAV
          DO   ITAU=1,NTAU
          Do   index_wspeed = 1,num
           XX(index_wspeed)=Windspeed_Lut(index_wspeed)
           YY(index_wspeed)=ALBEDO_R_int(ITAU,IWAV,ICASE,index_wspeed)
           YYY(index_wspeed)=ALBEDO_T_int(ITAU,IWAV,ICASE,index_wspeed)
          Enddo
            Y1=0
           CALL INTERP(Num,WSPEED,XX,YY,Y1)
           ALBEDO_R_SMALL_tau(ITAU,IWAV,ICASE)=Y1
            Y1=0
            CALL INTERP(Num,WSPEED,XX,YYY,Y1)
            ALBEDO_T_SMALL_tau(ITAU,IWAV,ICASE)=Y1
         
          Enddo
          Enddo
          Enddo 
          
           Num=0
          do index_wspeed = Indx_wspped1,Indx_wspped2
          Num=Num+1
          Windspeed_Lut(num)=Wind(index_wspeed) 
          Y1=0.0   
!    INTEROOLATE FOR LARGE CASES
          DO 15  ICASE=1,NUMCASEB
          DO 110 IWAV=1,NWAV
          DO 120 ITAU=1,NTAU
           LL=0
          DO 130 ITH0 =KSZAM1,KSZAP1
             LL=LL+1
            X(LL)=THET0(ITH0)
            Y(LL)=ALBEDO_R_BIG(index_wspeed,ITH0,ITAU,IWAV,ICASE)
            YY1(LL)=ALBEDO_T_BIG(index_wspeed,ITH0,ITAU,IWAV,ICASE)
 130        CONTINUE
             y1=0.0
             CALL INTERP(LL,MTHET0,X,Y,Y1)
              ALBEDO_R_int(ITAU,IWAV,ICASE,num)=Y1 
            y1=0.0
             CALL INTERP(LL,MTHET0,X,YY1,Y1)
          ALBEDO_T_int(ITAU,IWAV,ICASE,num)=y1 
  120       CONTINUE
  110       CONTINUE
   15       CONTINUE
           ENDDO
           
          DO   ICASE=1,NUMCASEB
          DO   IWAV=1,NWAV
          DO   ITAU=1,NTAU
          Do   index_wspeed = 1,num
           XX(index_wspeed)=Windspeed_Lut(index_wspeed)
           YY(index_wspeed)=ALBEDO_R_int(ITAU,IWAV,ICASE,index_wspeed)
           YYY(index_wspeed)=ALBEDO_T_int(ITAU,IWAV,ICASE,index_wspeed)
          Enddo
            Y1=0
           CALL INTERP(Num,WSPEED,XX,YY,Y1)
          ALBEDO_R_BIG_tau(ITAU,IWAV,ICASE)=Y1
            Y1=0
            CALL INTERP(Num,WSPEED,XX,YYY,Y1)
            ALBEDO_T_BIG_tau(ITAU,IWAV,ICASE)=Y1 
          Enddo
          Enddo
          Enddo 
            RETURN
            END 
!*****************************************************************

       SUBROUTINE INTtau_albedo_tran(ALBEDO_R_SMALL_tau,&
             ALBEDO_R_BIG_tau,ALBEDO_T_SMALL_tau,ALBEDO_T_BIG_tau,&
             ALBEDO_R_SMALL_final,ALBEDO_R_BIG_final,&
             ALBEDO_T_SMALL_final,ALBEDO_T_BIG_final,TAUAS,TAUAB,&
             ISMALL,IBIG)

!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine interpolates albedo and transmission on
!            optical computed
!
!!INPUT PARAMETERS:
!  Albedo for monodirectionally incident radiation Small and large
!  ALBEDO_R_SMALL_tau
! ALBEDO_R_BIG_tau
! InterpolatedTransmission for monodirectionally incident radiation Small and
! large
! ALBEDO_T_SMALL_tau
! ALBEDO_T_BIG_tau)
! ISMALL            Index for small mode
! IBIG              Index for Large mode
! Optical thicknesses for small and large mode
!TAUAS
!TAUAB
!!OUTPUT PARAMETERS:
! Interpolated Albedo for monodirectionally incident radiation Small and large
!   ALBEDO_R_SMALL_final
!   ALBEDO_R_BIG_final
! InterpolatedTransmission for monodirectionally incident radiation Small and
! large
!  ALBEDO_T_SMALL_final
! ALBEDO_T_BIG_final 
!----------------------------------------------------------------------
       USE Linear_interpolation 
      IMPLICIT NONE
       Save
      INCLUDE 'mod04.inc'

      REAL ALBEDO_R_SMALL_tau(NTAU,NWAV,NUMCASES)
      REAL ALBEDO_R_BIG_tau(NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL_tau(NTAU,NWAV,NUMCASES)
      REAL ALBEDO_T_BIG_tau(NTAU,NWAV,NUMCASEB)
      REAL ALBEDO_R_SMALL_final(NWAV,NUMCASES)
      REAL ALBEDO_R_BIG_final(NWAV,NUMCASEB)
      REAL ALBEDO_T_SMALL_final(NWAV,NUMCASES)
      REAL ALBEDO_T_BIG_final(NWAV,NUMCASEB)
      REAL TAUAS(NUMCASES,NWAV,NTAU)
      REAL TAUAB(NUMCASEB,NWAV,NTAU)
      INTEGER IJ,IWAV,LOPT,ITAU,ISMALL,IBIG
      REAL  X(100),Y(100),yy1(100),Y1,XBAR
       
! LOOP IS AROUND THET0,NEW FILE FOR EACH THETA0
          Y1=0.0
          DO IJ = 1,100
          X(IJ)=0.0
          Y(IJ)=0.0
          YY1(IJ)=0.0
          ENDDO
!    INTEROOLATE FOR SMALL CASES for total optical depth at 0.55
        DO 10 IWAV=1,NWAV
         XBAR=TAU_COMPUTED(Index_wave_550)
          DO 20 ITAU=1,NTAU
            X(ITAU)=TAUAS(ISMALL,Index_wave_550,ITAU)
            Y(ITAU)=ALBEDO_R_SMALL_tau(ITAU,IWAV,ISMALL)
            YY1(ITAU)=ALBEDO_T_SMALL_tau(ITAU,IWAV,ISMALL) 
  20      CONTINUE
                y1=0.0
             CALL INTERP(NTAU,XBAR,X,Y,Y1)
!           CALL RATINT(X,Y,NTAU,XBAR,y1)
             ALBEDO_R_SMALL_final(IWAV,ISMALL)=Y1
            y1=0.0
              CALL INTERP(NTAU,XBAR,X,YY1,Y1)
!               CALL RATINT(X,YY1,NTAU,XBAR,y1)
           ALBEDO_T_SMALL_final(IWAV,ISMALL)=Y1
  10      CONTINUE
!    INTEROOLATE FOR LARGE CASES
          DO 110 IWAV=1,NWAV
        XBAR=TAU_COMPUTED(Index_wave_550)
          DO 120 ITAU=1,NTAU
            X(ITAU)=TAUAB(IBIG,Index_wave_550,ITAU)
            Y(ITAU)=ALBEDO_R_BIG_tau(ITAU,IWAV,IBIG)
            YY1(ITAU)=ALBEDO_T_BIG_tau(ITAU,IWAV,IBIG)
 120        CONTINUE
             y1=0.0
             CALL INTERP(NTAU,XBAR,X,Y,Y1)
!               CALL RATINT(X,Y,NTAU,XBAR,y1)
             ALBEDO_R_BIG_final(IWAV,IBIG)=Y1
            y1=0.0
              CALL INTERP(NTAU,XBAR,X,YY1,Y1)
!            CALL RATINT(X,YY1,NTAU,XBAR,y1)
           ALBEDO_T_BIG_final(IWAV,IBIG)=Y1
  110       CONTINUE

            RETURN
            END 
!*****************************************************************

       SUBROUTINE SET_index_inter(MTHET0,MTHET,MPHI,THET0,&
       KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1)

!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine sets the index for ending and starting
!             positions from the lookup table
!
!!INPUT PARAMETERS:
!              MTHET0      Measured solar Zenith angle
!               MTHET      Measured view  Zenith angle
!                MPHI      Measured Azimuthal Angle
!               THET0      array of  solar Zenith angle in look_up table
!
!!OUTPUT PARAMETERS:
!               KSZAM1     Starting Index for solar zenith angle
!               KSZAP1     Ending Index for solar zenith angle
!               KTHEM1     Starting Index for view angle
!               KTHEP1     Ending   Index for  view angle
!               KPHIM1     Starting Index for  azimuth angle
!               KPHIP1     Ending   Index for  azimuth angle
! 
!----------------------------------------------------------------------

       IMPLICIT NONE
       SAVE

       INCLUDE 'mod04.inc'

       REAL MTHET0,MTHET,MPHI,THET0(NTH0),DEL
       INTEGER KSZA,KSZAM1,KSZAP1,KTHE,KTHEM1,KTHEP1,KPHI,KPHIM1,KPHIP1

!  set theta0
       KSZA=0
       KSZAM1=0
       KSZAP1=0
       KTHE=0
       KTHEM1=0
       KTHEP1=0
       KPHI=0
       KPHIM1=0
       KPHIP1=0
       IF(MTHET0 .LE.THET0(5))THEN
          DEL=12.0
          KSZA=MTHET0/DEL+1
       ELSEIF(MTHET0 .GT.THET0(5))THEN
           DEL=6.0
           KSZA=MTHET0/DEL-3
       ENDIF
        KSZAM1=KSZA
        KSZAP1=KSZA+1
       IF(KSZAM1.LE.0)THEN
          KSZAM1=1
           KSZAP1=2
       ENDIF
       IF(KSZAP1.GE.NTH0)THEN
          KSZAM1=KSZA-1
          KSZAP1=KSZAM1+1
       ENDIF
!  set theta
       KTHE=MTHET/6.0+1
        KTHEM1=KTHE
       KTHEP1=KTHE+1
       IF(KTHEM1.LE.0)THEN
         KTHEM1=1
         KTHEP1=2
       ENDIF
       IF(KTHEP1.GE.NTHET)THEN
         KTHEM1=KTHE-1
         KTHEP1= KTHEM1+1
       ENDIF
!  set phi
        KPHI=MPHI/12.0+1
        KPHIM1=KPHI
        KPHIP1=KPHI+1
       IF(KPHIM1.LE.0)THEN
         KPHIM1=1
         KPHIP1=2
       ENDIF
       IF(KPHIP1.GE.NPHI)THEN
         KPHIM1=KPHI-1
         KPHIP1=KPHIM1+1
       ENDIF
100    format(f5.1,3i6)

       RETURN
       END



!*****************************************************************

       Subroutine Fill_QAflag_ocean( QA_Flag_Ocean,SDS_QCONTROL_Ocean,&
          Idata,Quality_to_pass)
!-----------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine sets the array for quality control flags
!
!!INPUT PARAMETERS:
!               QA_Flag_Ocean           Quality flag
!               Idata                   Index of 10*10 box
!
!!OUTPUT PARAMETERS:
!               SDS_QCONTROL_Ocean      HDF array for quality control
! ---------------------------------------------------------------------

       IMPLICIT NONE
       SAVE

       INCLUDE 'mod04.inc'
       INCLUDE 'read_Sat_MODIS.inc'

       Integer QA_Flag_Ocean(12),idata,ii,Quality_to_pass(2)
       BYTE    SDS_QCONTROL_Ocean(QA_Ocean,NUMCELLS_B),QA_Temp

       Quality_to_pass(1) = -99999
       Quality_to_pass(2) = -99999
! Product Quality & Retrieval Processing QA flags over ocean
! Quality of 3 and 4 ( Average solution)are repeated tion have
! as best and average solution have same quality.
!

!        QA_Flag_Ocean(1)=0     NO retrieval ( Not useful)
!        QA_Flag_Ocean(1)=1     retrieval    (  useful)
!
!    For all non retrieval boxes QA_Flag_Ocean is set to zero to indicate
!    not useful and QA_Flag_Ocean(5) is set to values of QCONTROL to
!    indicate the cause of Retrieval. QA_Flag_Ocean(6) is set to 11 for 
!    no retrieval
        
         If ( QCONTROL .lt. 0) then
           QA_Flag_Ocean(1)=0
           QA_Flag_Ocean(2)=0
           QA_Flag_Ocean(3)=0
           QA_Flag_Ocean(4)=0
         QA_Flag_Ocean(5)=abs(QCONTROL)
         QA_Flag_Ocean(6)=15
         Else
            QA_Flag_Ocean(1)=1
            QA_Flag_Ocean(3)=1
!  All retrieval boxes
! Estimated Quality of aerosol parameters
!        QA_Flag_Ocean(2)=3     very Good
!        QA_Flag_Ocean(2)=2     Good
!        QA_Flag_Ocean(2)=1     Marginal
!        QA_Flag_Ocean(2)=0     Bad
!        QA_Flag_Ocean(4)=3     very Good
!        QA_Flag_Ocean(4)=2     Good
!        QA_Flag_Ocean(4)=1     Marginal
!        QA_Flag_Ocean(4)=0     Bad


          if( Quality_dust_flag_glint .eq.1)QCONTROL=12
           if(Qcontrol_cirrus.eq.1)QCONTROL=13
          if( Quality_dust_flag_off_glint .eq.1)QCONTROL=14

 
!       IF QCONTROL is 0( see doc.)quality of retrieval  is very good
 
       If( QCONTROL .eq. 0 )QA_Flag_Ocean(2)=3
       If( QCONTROL .eq. 0 )QA_Flag_Ocean(4)=3
 
!       IF QCONTROL is 7 or 14( see doc.)quality of retrieval  is  good
 
       If( QCONTROL .eq. 7 .or. QCONTROL .eq. 14)QA_Flag_Ocean(2)=2
       If( QCONTROL .eq. 7 .or. QCONTROL .eq. 14)QA_Flag_Ocean(4)=2
 
!       IF QCONTROL is 1,3,4,6,8 or 10  quality of retrieval  is Average
 
       If( QCONTROL .eq. 1  .or. QCONTROL .eq. 3 .or. QCONTROL .eq. 4 &
        .or. QCONTROL .eq. 6 .or. QCONTROL .eq. 8 &
        .or. QCONTROL .eq. 10)QA_Flag_Ocean(2)=1
       If( QCONTROL .eq. 1  .or. QCONTROL .eq. 3 .or. QCONTROL .eq. 4 &
         .or. QCONTROL .eq. 6 .or. QCONTROL .eq. 8 &
         .or. QCONTROL .eq. 10)QA_Flag_Ocean(4)=1
 
!       IF QCONTROL is 2,5,9,12,13 quality of retrieval  is poor
        
 
       If( QCONTROL .eq. 2  .or. QCONTROL .eq. 5 .or. QCONTROL .eq. 9 &
                  .or. QCONTROL .eq. 12 .or. QCONTROL.eq.13) &
          QA_Flag_Ocean(2)=0
     
       If( QCONTROL .eq. 2  .or. QCONTROL .eq. 5 .or. QCONTROL .eq. 9 &
               .or. QCONTROL .eq. 12 .or. QCONTROL.eq.13) &
       QA_Flag_Ocean(4)=0

        
          
         QA_Flag_Ocean(5)=0
         QA_Flag_Ocean(6)=QCONTROL
                IF(QCONTROL .eq. 17)then 
                     QA_Flag_Ocean(2)=3
                     QA_Flag_Ocean(4)=3 
                 Endif
         endif
         
          Quality_to_pass(1)= QA_Flag_Ocean(2)
          Quality_to_pass(2)= QCONTROL
 
! Store QA flags into Quality_Assurance_Ocean array according to the order
! of bits in MODIS atmosphere QA plan
 
      QA_Temp=0
      CALL BYTE_SET(QA_Flag_Ocean(1),0,QA_Temp)
      CALL BYTE_SET(QA_Flag_Ocean(2),1,QA_Temp)
      CALL BYTE_SET(QA_Flag_Ocean(3),4,QA_Temp)
      CALL BYTE_SET(QA_Flag_Ocean(4),5,QA_Temp)
!  End of 1 byte

      SDS_QCONTROL_Ocean(1,IDATA)=QA_Temp
 
 
! Retrieval processing QA flags_Processing path flags
! For all nonRetrieval boxes change -ve numbers to +ve to repreasent
! the bit array.
! For all Retrieval boxes fill another array
      QA_Temp=0
      CALL BYTE_SET(QA_Flag_Ocean(5),0,QA_Temp)
      CALL BYTE_SET(QA_Flag_Ocean(6),4,QA_Temp)
!  End of 2 byte

      SDS_QCONTROL_Ocean(2,IDATA)=QA_Temp
!       write(33,*)'sec',idata, (QA_Flag_Ocean(ii),ii=5,6),
!     *  QA_Temp,SDS_QCONTROL_Ocean(2,IDATA)

          
!10     format(12I7)

       Return
       end



!*****************************************************************

      Subroutine Total_retrieval_ocean(Success_Ret_Ocean, Fail_Ret_Ocean) 

!----------------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine counts the total number of boxes based on the
! retreival or non-retreival.
!
!!INPUT PARAMETERS:
!
!     QCONTROL           Through the mod04.inc
!
!!OUTPUT PARAMETERS:
!
!      Success_Ret_Ocean  Total number of boxes where retreival was made
!      Fail_Ret_Ocean     Total number of boxes where no retreival was made
! 
!----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      integer Fail_Ret_Ocean,Success_Ret_Ocean

      if( QCONTROL .lt.0.0)then
         Fail_Ret_Ocean=Fail_Ret_Ocean+1
      else
        Success_Ret_Ocean=Success_Ret_Ocean+1
      endif
      return
      end
         Subroutine make_sediment(W470_syn,W550_SYN,&
      W124_SYN,W164_SYN,W213_SYN,sed_mask,IX,IY,wave,Sat_flag)
!----------------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION: This subroutine identifies the sediments 
! and sets a mask.
!
!!INPUT PARAMETERS:
!      IX        index of 10*10 box
!      IY        index of 10*10 box
!    Reflectances of wavelengths used modis,ch1-ch7  identified by W*_SYN
!       W470_SYN
!       W550_SYN
!       W124_SYN
!       W164_SYN
!       W213_SYN
!        WAVE    wavelength index
!!OUTPUT PARAMETERS:
!
!      sed_mask   Sediment mask 0 = no sediment 1= sediment
! 
!----------------------------------------------------------------------
           IMPLICIT NONE
            SAVE
          INCLUDE 'mod04.inc'
          INCLUDE 'read_Sat_MODIS.inc'
         REAL W865_SYN(ISWATH_B,ILINE),&
           W470_SYN(ISWATH_B,ILINE),W550_SYN(ISWATH_B,ILINE),&
           W124_SYN(ISWATH_B,ILINE),W164_SYN(ISWATH_B,ILINE),&
           W213_SYN(ISWATH_B,ILINE),Sed_W659(ISWATH_B,ILINE),&
           Sed_W865(ISWATH_B,ILINE),a,b,wave(nwav),xa
           Real ref_55_inter_extra,a1,a2,a3,b1,b2,b3,ave_a,ave_b
           Integer sed_mask(ISWATH_B,ILINE),ij,ix,iy,iscan,idata 
           integer waveused,opt
           real x(nwav),y(nwav),sig(nwav),ref_55_inter,Del_wav55
           character (len =10):: Sat_flag
             sed_mask(ix,iy)=1
             ref_55_inter =0
              Del_wav55 = 0 
!            print*,IX,IY, W470_syn(IX,IY),W213_SYN(IX,IY),&
!                    W164_SYN(IX,IY),W124_SYN(IX,IY) 
                    
          IF( W470_syn(IX,IY) .gt. 0 .and. W213_SYN(IX,IY) .gt.0  &
             .and. W164_SYN(IX,IY) .gt.0 .and. W124_SYN(IX,IY) .gt.0) THEN   
                 opt=1
                waveused=4
         elseif( W470_syn(IX,IY) .gt. 0 .and. W213_SYN(IX,IY) .gt.0  &
               .and. W124_SYN(IX,IY) .gt.0) THEN  
                  opt=0
                waveused=3  
         elseif( W470_syn(IX,IY) .gt. 0 .and. W213_SYN(IX,IY) .gt.0  &
                 .and. W164_SYN(IX,IY) .gt.0) THEN   
                 opt=2
                waveused=3 
          else 
!    Atleast 1.24 or 1.64  have to have valid data ( 0.47 and 2.13 have to have valid data)
!    If all wavelengths data is missing we call it bad data by declaring it sediment mask.
!    this data point will be discarded in further processing          
           sed_mask(ix,iy)=0
            waveused = 0
          endif
             
             
           
 ! compute sediment mask if no fill values          
           If ( waveused .gt.0) then
           IF( opt .eq.0 .and. waveused .eq.3)then 
!  Numerical version of  fit
              do ij=1,waveused 
             if(ij.eq.1)x(ij)=alog(wave(1))
             if(ij.eq.2)x(ij)=alog(wave(5)) 
             if(ij.eq.3)x(ij)=alog(wave(7))
             if( ij .eq.1)y(ij)=alog(W470_syn(IX,IY))
             if( ij .eq.2)y(ij)=alog(W124_SYN(IX,IY)) 
             if( ij .eq.3)y(ij)=alog(W213_SYN(IX,IY))
               enddo  
          Elseif( opt .eq.1 .and. waveused .eq.4)then
                do ij=1,waveused 
             if(ij.eq.1)x(ij)=alog(wave(1))
             if(ij.eq.2)x(ij)=alog(wave(5))
             if(ij.eq.3)x(ij)=alog(wave(6))
             if(ij.eq.4)x(ij)=alog(wave(7))
             if( ij .eq.1)y(ij)=alog(W470_syn(IX,IY))
             if( ij .eq.2)y(ij)=alog(W124_SYN(IX,IY))
             if( ij .eq.3)y(ij)=alog(W164_SYN(IX,IY))
             if( ij .eq.4)y(ij)=alog(W213_SYN(IX,IY)) 
             enddo
         Elseif(opt .eq.2 .and. waveused .eq.3)then 
               do ij=1,waveused 
             if(ij.eq.1)x(ij)=alog(wave(1))
             if(ij.eq.2)x(ij)=alog(wave(6)) 
             if(ij.eq.3)x(ij)=alog(wave(7))
             if( ij .eq.1)y(ij)=alog(W470_syn(IX,IY))
             if( ij .eq.2)y(ij)=alog(W164_SYN(IX,IY))
             if( ij .eq.3)y(ij)=alog(W213_SYN(IX,IY))
               enddo  
          ENDIF   
              call  FIT_line(X,Y,SIG,waveused,A,B)
              XA=ALOG(wave(2))
              ref_55_inter=EXP(A+(B*XA)) 
             Del_wav55=W550_SYN(IX,IY)-ref_55_inter
     
! mask1 true it is sedement+smoke+dust set mask to sedements
             if(W213_SYN(IX,IY) .lt.0.10 .and. Del_wav55 .ge. 0.015)then  
               sed_mask(ix,iy)=0
! mask2 true it is  smoke+dust and no sedements set previuos sedement mask
! to no sedement because it smoke or dust.
               if(W470_syn(IX,IY) .ge. 0.25 .and. Del_wav55 .ge. 0.015)then  
               sed_mask(ix,iy)=1
                endif
              endif 
! endif for all valid data              
         Endif  
         return 
         end  
 
      SUBROUTINE FIT_line(X,Y,SIG,NDATA,A,B)
!----------------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION:
! THIS SUBROUTINE IS TAKEN FROM NUMERICAL RECIPES PAGE 508
!GIVEN A SET OF NDATA POINTS X(I),Y(I) WITH STANDARD DEVIATIONS
! SIG(I),FIT THEM TO A STRIGHT LINE Y=A+BX BY MINIMIZING X**2
!!RETURNED ARE A,B AND THEIR RESPECTIVE PROBABLE UNCERTAINTIES.
! SIGA AND SIGB ,THE CHI-SQUARE AND THE GOODNESS OF FIT PROBABILITY Q
! (THAT THE FIT WOULD HAVE X**2 THIS LARGE OR LARGER). IF MWT=O ON INPUT
! THEN STANDARD DEVIATIONS ARE ASSUMED TO BE UNAVAILABLE:Q IS RETURNED
! AS 1.0 AND THE NORMALIZATION OF CHI-SQARE IS TO UNIT STANDARD DEVATION
! ON ALL POINTS
!
!!INPUT PARAMETERS:
!     X        data set x
!     Y        data set y
!    NDATA  number of data points
!
!!OUTPUT PARAMETERS:
!
!     A      intercept
!     B       slope
!----------------------------------------------------------------------

 
      Integer ndata
      REAL      X(NDATA),Y(NDATA),SIG(NDATA)
      SX=0.
      SY=0.
      ST2=0.
      B=0.
       MWT=0
      IF(MWT.NE.0) THEN
        SS=0.
        DO 11 I=1,NDATA
          WT=1./(SIG(I)**2)
          SS=SS+WT
          SX=SX+X(I)*WT
          SY=SY+Y(I)*WT
11      CONTINUE
      ELSE
        DO 12 I=1,NDATA
          SX=SX+X(I)
          SY=SY+Y(I)
12      CONTINUE
        SS=FLOAT(NDATA)
      ENDIF
      SXOSS=SX/SS
      IF(MWT.NE.0) THEN
        DO 13 I=1,NDATA
          T=(X(I)-SXOSS)/SIG(I)
          ST2=ST2+T*T
          B=B+T*Y(I)/SIG(I)
13      CONTINUE
      ELSE
        DO 14 I=1,NDATA
          T=X(I)-SXOSS
          ST2=ST2+T*T
          B=B+T*Y(I)
14      CONTINUE
      ENDIF
      B=B/ST2
      A=(SY-SX*B)/SS
       SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)
       SIGB=SQRT(1./ST2)
      RETURN
      END
       FUNCTION GAMMQ(A,X)
       real GAMMQ
      IF(X.LT.0..OR.A.LE.0.)STOP
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      RETURN
      END
       
      SUBROUTINE GSER(GAMSER,A,X,GLN)
!----------------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION:
! Call routine from Fit_line
!THIS SUBROUTINE IS TAKEN FROM NUMERICAL RECIPES PAGE 508
! 
!!INPUT PARAMETERS:
!     
!    
!
!!OUTPUT PARAMETERS: 
!--------------------------------------------------------------------

      integer  ITMAX
      real      EPS,GAMMLN
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)STOP
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      STOP 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END

      SUBROUTINE GCF(GAMMCF,A,X,GLN)
!----------------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION:
! Call routine from Fit_line
!THIS SUBROUTINE IS TAKEN FROM NUMERICAL RECIPES PAGE 508
!
!!INPUT PARAMETERS:
!     
!    
!
!!OUTPUT PARAMETERS:
!
!     
!
!!REVISION HISTORY:
!  
!
!!TEAM-UNIQUE HEADER:
!
!!REFERENCES AND CREDITS  NUMERICAL RECEPIES 
!----------------------------------------------------------------------
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      STOP 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
      FUNCTION GAMMLN(XX)
!----------------------------------------------------------------------------
!!F90      
!
!!DESCRIPTION:
! Call routine from Fit_line
!THIS SUBROUTINE IS TAKEN FROM NUMERICAL RECIPES PAGE 508
! 
!!INPUT PARAMETERS:
!   
!    
!
!!OUTPUT PARAMETERS:
!
!      
! 
!----------------------------------------------------------------------
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
        -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
      
!*********************************************************************


        SUBROUTINE Dust_Detection_Ocean(IXX,IYY,SD_for_Dust, &
        W659_syn,W865_syn,W470_syn,W550_SYN,W213_SYN,W412_SYN,W8p5_Temp,&
        W1100_Temp,dust_flag)
!----------------------------------------------------------------------
!!F90
!
!!DESCRIPTION:  subroutine Dust_Detection_Ocean detect dusty pixels in the
!                10x10 aerosol retrieval grid and reports a dust flag if
!		more than 3 dusty pixels are detected. 
!
!!INPUT PARAMETERS:
!	MOD04_cloudmask(10x10) in 1km resolution
!       landwater flag(10x10)
!       stdev86(10x10)   standard deviation of 0.86um reflectance in 9x9 1km pixels
!       ref65(10x10)     band1 reflectance 
!	ref86(10x10)     band2
!	ref47(10x10)     band3
!       ref55(10x10)     band4
!       ref21(10x10)     band7
!       ref41(10x10)     band8
!       BT87(10x10)      band29 brightness temperature
!       BT11(10x10)      band31 brightness temperature
!!OUTPUT PARAMETERS:
!	dust_flag   1 for dusty, 0 for non-dusty
!	dust_counts  number of dusty pixels detected
      IMPLICIT NONE
       SAVE 
       INCLUDE 'mod04.inc' 
       INCLUDE 'read_Sat_MODIS.inc'
     Real 	tdiff3,rat31,rat14,rat43,NDAI,bad_data
	 integer dust_flag, IXX,IYY
	 REAL W659_SYN(ISWATH_B,ILINE),W865_SYN(ISWATH_B,ILINE),&
           W470_SYN(ISWATH_B,ILINE),W550_SYN(ISWATH_B,ILINE),&
           W124_SYN(ISWATH_B,ILINE),W164_SYN(ISWATH_B,ILINE),&
           W213_SYN(ISWATH_B,ILINE),&
           W412_SYN(ISWATH_B,ILINE),W443_SYN(ISWATH_B,ILINE),&
           W8p5_Temp(ISWATH_B,ILINE),W1100_Temp(ISWATH_B,ILINE),&
           SD_for_Dust(ISWATH_B,ILINE)
      

	bad_data = 0

!  No Temp files are available for PACE
           
	
	  IF (W550_SYN(IXX,IYY)     .gt. bad_data .and. W865_syn(IXX,IYY) .gt. bad_data &
	    .and. W470_syn(IXX,IYY) .gt. bad_data .and. W550_SYN(IXX,IYY) .gt. bad_data &
		.and. W213_SYN(IXX,IYY) .gt. bad_data .and. W412_SYN(IXX,IYY) .gt. bad_data)THEN
!		 tdiff3 = W8p5_Temp(IXX,IYY) - W1100_Temp(IXX,IYY)
		 rat31 =  W470_syn(IXX,IYY)/W659_syn(IXX,IYY)
		 rat14 =  W659_syn(IXX,IYY)/W550_SYN(IXX,IYY)
		 rat43 =  W550_SYN(IXX,IYY)/W470_syn(IXX,IYY)
	 	  NDAI = -10.*log10(W412_SYN(IXX,IYY)/W213_SYN(IXX,IYY)) 
         IF (SD_for_Dust(IXX,IYY) .le. 0.01 .and. W412_SYN(IXX,IYY) .LE. 0.35 &
 			       .and. W412_SYN(IXX,IYY) .ge. 0.1) THEN
	              IF (rat31 .lt. 0.9) dust_flag = 1
	         IF (rat31 .ge. 0.9 .and.rat31 .lt. 1.8 ) THEN 
       			 IF (rat14 .ge. 0.70 .and. rat43 .ge. 0.68) THEN  
      		        IF (NDAI .GT. -7.0 )dust_flag = 1  
		       ENDIF
	       ENDIF	 
        ENDIF 
      ENDIF    
	   RETURN
	   
	   END 
