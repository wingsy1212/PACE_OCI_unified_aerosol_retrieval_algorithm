         
        !
! LAND SDS_ARRAYS..........FOR LAND ONLY
!
      BYTE         SDS_QCONTROL_land(QA_LAND,NUMCELLS_B)
      BYTE         SDS_QCONTROL_CritRef_land(QA_LAND,NUMCELLS_B)
      BYTE         SDS_QCONTROL_ocean(QA_Ocean,NUMCELLS_B)
      Real   SDS_Aerosol_Type(NUMCELLS_B),&
                   SDS_SCAT_ANGLE_land(NUMCELLS_B),&
                  SDS_angs_coeff_land(NUMCELLS_B),&
                   SDS_CLDFRC_land(NUMCELLS_B),&
                   SDS_dust_weighting(NUMCELLS_B),&
                   SDS_NUMPIXELS_land(NUMCELLS_B,Band_land+2),&
                   SDSTAU_corrected(NUMCELLS_B,Land_Sol3),&
                   SDS_ref_land(NUMCELLS_B,(Band_land+2)),&
                   SDS_ref_STD_land(NUMCELLS_B,(Band_land+2)) 
                   
! 9/2005 ( two new SDS's)......     
        Real SDS_Surface_Reflectance_Land(NUMCELLS_B,Land_Sol3),&
                  SDS_Fitting_Error_Land(NUMCELLS_B),&
                  SDSTAU_corrected_213(NUMCELLS_B),&
                  SDSTAU_small_land(NUMCELLS_B,Land_Sol4)

 
       Real  &
                  SDS_est_uncer(NUMCELLS_B,Land_Sol1),&
                  SDS_RefF_land(NUMCELLS_B,Land_Sol2),&
                  SDS_TranF_land(NUMCELLS_B,Land_Sol1)
! 
! OCEAN SDS_ARRAYS.........FOR OCEAN ONLY
!
      REAL         SDS_ccn(NUMCELLS_B,NUM_solutions),&
                   SDS_mass_conc(NUMCELLS_B,NUM_solutions)
      BYTE         SDS_QCONTROL(QA_OCEAN,NUMCELLS_B)
      Real   SDS_ref(NUMCELLS_B,(NWAV_S+3)),&
                  SDS_ref_STD(NUMCELLS_B,(NWAV_S+3)),&
                  SDSTAU_best(NUMCELLS_B,NWAV_S),&
                   SDSTAU_average(NUMCELLS_B,NWAV_S),&
                   SDSTAUB_best(NUMCELLS_B,NWAV_S),&
                   SDSTAUB_average(NUMCELLS_B,NWAV_S),&
                   SDSTAUS_best(NUMCELLS_B,NWAV_S),&
                   SDSTAUS_average(NUMCELLS_B,NWAV_S),&
                   SDSASSY_best(NUMCELLS_B,NWAV_S),&
                   SDSASSY_average(NUMCELLS_B,NWAV_S),&
                   SDSBACK_best(NUMCELLS_B,NWAV_S),&
                  SDSBACK_average(NUMCELLS_B,NWAV_S),&
                  SDS_RefF_best(NUMCELLS_B,NWAV_S) 
        Real       SDS_RefF_average(NUMCELLS_B,NWAV_S),&
                   SDS_TranF_best(NUMCELLS_B,NWAV_S),&
                   SDS_TranF_average(NUMCELLS_B,NWAV_S),&
                   SDS_small_weighting(NUMCELLS_B,NUM_solutions),&
                   SDS_correc_small_weighting(NUMCELLS_B),&
                   SDS_Least_error(NUMCELLS_B,NUM_solutions),&
                   SDS_effrad(NUMCELLS_B,NUM_solutions),& 
                   SDS_angs_coeff1(NUMCELLS_B,NUM_solutions),&
                   SDS_angs_coeff2(NUMCELLS_B,NUM_solutions),&
                   SDS_AOT_model(NUMCELLS_B,num_model_index),&
                   SDS_NUMPIXELS(NUMCELLS_B,NWAV_S+3),SDS_CLDFRC(NUMCELLS_B),&
                   SDS_SCAT_ANGLE_OCEAN(NUMCELLS_B),&
                   SDS_CLDFRC_ocean(NUMCELLS_B)
                   
    integer        SDS_sol_INDX_small(NUMCELLS_B,NUM_solutions),&
                   SDS_sol_INDX_large(NUMCELLS_B,NUM_solutions)
!
! Input and Output file names for PGSTK etc.......
!
      INTEGER flag_tarfox
!
! Position of image in the processing of 10x10 km box
!
      INTEGER START_500,END_500,START_250,END_250,START_1KM,END_1KM
      INTEGER Iscan,Buf_Size1,Buf_Size2,Data_Size(2),i,j
      INTEGER IMONTH(maxnum_scans),IDATA,Water,Land,cloud_num,&
             cloud_num_land,Glint,Snowice,ij,ik,Qcontrol_special,&
              Aerosol_cloud_land(iswath_B*2,ILINE*2),&
               Qcontrol_special_land,Quality_to_pass(2)
 
! Define arrays of Geoloation, cloud mask and solar and satellite angles
 
      REAL MTHET0,MTHET,MPHI0,MPHI,MDPHI,MSCATT,Lat_center,Lon_center
      REAL sfctmp,ugrd,vgrd,pwat,ozone,skinTemp
      REAL Lat(iswath_B,ILINE),Lon(iswath_B,ILINE),&
           SatZen(iswath_B,ILINE),SolZen(iswath_B,ILINE),&
           SatAz(iswath_B,ILINE),SolAz(iswath_B,ILINE),&
           RelAz(iswath_B,ILINE),Height(iswath_B,ILINE),G_factor
      Real  ref_allwave_Vis(NWAV+NWAV_uv)
!!!!! UV Variables  
     REAL Lat_UV(iswath_B,ILINE),Lon_UV(iswath_B,ILINE),&
           SatZen_UV(iswath_B,ILINE),SolZen_UV(iswath_B,ILINE),&
           SatAz_UV(iswath_B,ILINE),SolAz_UV(iswath_B,ILINE),& 
            EXTbig_550,EXTSMALL_550 
       REAL MTHET0_UV,MTHET_UV,MPHI0_UV,MPHI_UV,MDPHI_UV,MSCATT_UV ,Lat_center_UV,&
             Lon_center_UV,Tau_550,ref_allwav_uv(NWAV_uv)
       INTEGER set_flag_read_UV
       REAL tau_new(NWAV_uv),Omega_new(NWAV_uv+1),Small_m_weighting
          integer Mode_F,Mode_C  
       Real Tau_allwave(NWAV_S),WAVE(NWAV_S) 
       Real Ext_554_small(Lut_indx,NUMCASES),Ext_554_large(Lut_indx,NUMCASEB) 
       INTEGER CldMsk_500_Land(IX1KM_B,IY1KM_B)
! New
       
       INTEGER CldMsk_250(IX1KM_B,IY1KM_B),&
               CldMsk_500(IX1KM_B,IY1KM_B) ,&
               CldMsk_Cirrus(IX1KM_B,IY1KM_B),&
               DET_Flag(iswath_B,ILINE),&
                UFQ_Flag(iswath_B,ILINE),&
               DayNight_Flag(iswath_B,ILINE),&
               SunGlint_Flag(iswath_B,ILINE),&
               SnowIce_Flag(iswath_B,ILINE),&
               SnowMsk_Ratio(iswath_B,ILINE),&
               SnowMsk_500m(2*iswath_B,2*ILINE) 
     
         

           INTEGER      LandSea_Flag(iswath_B,ILINE),&
              Non_CloudOb_Flag(iswath_B,ILINE),&
              Thin_CirNIR_Flag(iswath_B,ILINE),&
              Shadow_Flag(iswath_B,ILINE),&
              Thin_CirIR_Flag(iswath_B,ILINE),&
               Cloud_SimpIR_Flag(iswath_B,ILINE),&
               High_Cloud_Flag(iswath_B,ILINE),&
               Cloud_IRTemp_Flag(iswath_B,ILINE),&
               Cloud_3p75_11_Flag(iswath_B,ILINE),&
               Cloud_VisRat_Flag(iswath_B,ILINE),&
              Cloud_SpatVar_Flag(iswath_B,ILINE),&
               High_Cloud_Flag_500(2*iswath_B,2*ILINE),&
               NROWS(iswath_B*2),MONTH
                

      

        REAL W659_SYN(iswath_B,ILINE),W865_SYN(iswath_B,ILINE),&
          W470_SYN(iswath_B,ILINE),W550_SYN(iswath_B,ILINE),&
          W124_SYN(iswath_B,ILINE),W164_SYN(iswath_B,ILINE),&
          W213_SYN(iswath_B,ILINE),&
          W412_SYN(iswath_B,ILINE),W443_SYN(iswath_B,ILINE),&
          W8P5_SYN(iswath_B,ILINE),W869o_SYN(iswath_B,ILINE),&
          W395_SYN(iswath_B,ILINE),W138_SYN(iswath_B,ILINE),&
          W1100_SYN(iswath_B,ILINE),& 
          W1200_SYN(iswath_B,ILINE),&
          W1100_Temp(iswath_B,ILINE),W8p5_Temp(iswath_B,ILINE),&
          W354_SYN(iswath_B,ILINE),W388_SYN(iswath_B,ILINE)
 
! Ancillary Data of Total precipitable water and ozone 
 
      REAL Total_H2O,Total_O3
!
! QA arrays (see MODIS Atmosphere QA Plan)
!
       INTEGER QA_Flag_Land(19),QA_Flag_Ocean(12),Success_Ret_Land,&
             Fail_Ret_Land,NO_Ret_Land,Total_Grids,&
             Success_Ret_Ocean,Fail_Ret_Ocean,NO_Ret_Ocean,&
             quality_cirrus(iswath_B,(Numscan_B*ILINE)+1),Ret_Quality_cirrus,&
             Quality_flag_forJoint(2)
        Integer Jday,Hour,Minutes
        INTEGER No_PSA,JX,JY,start_line,END_line
        PARAMETER(No_PSA=15)
        REAL QA_Metadata_MOD04(No_PSA)
         Real  Uwind(360*2,181*2),vwind(360*2,181*2),water_Vapor(360*2,181*2),Tozone(360*2,181*2) 
 
! Arrays for miscellaneous use
!
      INTEGER Set_Counter_Land,Set_Counter_Ocean,GetModisDat_MOD04
      INTEGER maskoption,num_resol,rtn,Set_Counter_Ocean_cloud
      INTEGER create_mod04,idebug,Indx_wspeed1,Indx_wspeed2
      CHARACTER*4 choice_flag
      LOGICAL error_flag,modis_flag
      Real    Glint_angle
      integer index_wave,quality_land,ix,iy,IL,IS
!    New cldmask
      REAL RMED(iswath_B ),RMEDSQ(iswath_B )
      REAL RMED_1km(iswath_B),RMEDSQ_1km(iswath_B)
      integer CldMsk_500_Ocean(IX1KM_B*2,IY1KM_B*2)
      integer Joint_Cld_mask_500(iswath_B,iline),ii,num,jj
      integer  Num_coastal_pixels 
      integer savecldmask(IX1KM_B*2,IY1KM_B*2) 
      integer num_joint 
      integer num_Sea_land,Full_Sea_Land_mask(iswath_B,ILINE*Tot_scan_B)
      integer xsize,ysize,Aerosol_Cldmask_land(iswath_B,(Numscan_B*ILINE)+1)
      integer xsize1,ysize1  
      integer Save_index_ocean(NUMCELLS_B,Max_Pixels_O),index_wspeed
      integer Save_index_land(NUMCELLS_B,Max_Pixels_O)
      integer Save_index_Land_Ocean(NUMCELLS_B,Tot_scan_B,Max_Pixels_O)
      integer save_numpix_Land_Ocean(NUMCELLS_B,Tot_scan_B),Iyy
       Real  new_Average_dis(NUMCELLS_B,Tot_scan_B),Wind(Lut_indx)  
       Integer    set_counter_for_Gread,set_counter_for_anc  
      integer  IXX,Y2_offset,X2_offset,l,p,y2,x2
      integer cloud_num_clear,IWAVE_NUm,rtn_ncep
      integer oHigh_Cloud_F(IX1KM_B,IY1KM_B),num_joint2
       Real SD_3by3(IX1KM_B,IY1KM_B)
       REal center_height,Multi_factor(10), Variable,Variable1,variable2
      integer Year,day,ifile,Water_Pixels,Land_Pixels,Kmonth
      Integer  Num_water,Num_land,coastal_Pix,Pure_Land
      integer LandSea_Flag_new(iswath_B,ILINE) 
       Real Sea_Land_Flag(iswath_B)
       Real SD_for_Dust(ISWATH_B,ILINE),ETA
       REAL AVE_ARRAY(77,2),Optical_depth_Ocean,Optical_depth_land
        Real   MHGHT,New_Optical_depth_Ocean, ss
        Integer New_savecldmask(iswath_B*2,ILINE*2)
        Integer New_CldMsk_500_Ocean(iswath_B*2,ILINE*2)
        Integer New_CldMsk_500_Land(iswath_B,ILINE),Imodels 
        integer New_Q_cirrus(iswath_B,ILINE),Return_flag 
        Real ERR644,AOD(10),cossza 
        Integer  Twenty_water, Fifty_coastal,k,HANDLE_Ext_554_O
        Integer Dust_flag_10KM(iswath_B/ILINE)
        Real*8   A_Height(IX1KM_B,IY1KM_B),G_Time(IX1KM_B,IY1KM_B)
        integer  read_once_option,yy1,yy2,xx1,Ipart,XX_ipart,YY_ipart
       integer dim_sizes(32), start(32), edges(32), stride(32),LinesPerRead 
       integer  checking_total_cloud,Status,Wave_index,ct
       character(len=10) :: SSat_Flag,Sat_Flag
       character(len=4)  ::  Anc_flag 
       
         integer Index_Omega_new(2)
         Real Fitting_error(2),Height_indx(2) 
         integer handle_Urban_Table_10km
            
          
         
       BYTE       Cloud(Buf_cldmsk,iswath_B,ILINE),&
                 QA_Cloud(Buf_cldmsk_QA,iswath_B,ILINE)
        Real  SDS_MTHET0(NUMCELLS_B),SDS_MTHET(NUMCELLS_B),&
                 SDS_MPHI0(NUMCELLS_B),SDS_MPHI(NUMCELLS_B),&
                 SDS_Scattering_Angle(NUMCELLS_B),&
                 SDS_Tau_Land_Ocean(NUMCELLS_B),&
                 SDS_Tau_Land_Ocean_img(NUMCELLS_B),&
                 SDS_ratio_small_Land_Ocean(NUMCELLS_B) ,&
                 SDS_Reflected_flux_Land_Ocean(NUMCELLS_B),&
                 SDS_Sea_Land_Flag(iswath_B,Tot_scan_B),&
                 SDS_Sea_Sunglint_Flag(iswath_B,Tot_scan_B),&
                 SDS_NCEP_Wspeed(iswath_B,Tot_scan_B), & 
                  SDS_Topo_Altitude(iswath_B,Tot_scan_B)
      REAL       SDSLAT(NUMCELLS_B),SDSLON(NUMCELLS_B),&
                 SDS_mass_conc_land(NUMCELLS_B),&
                 SLOPE_MEAN_LAND(3),SLOPE_MEAN_OCEAN(3)
      REAL        WSPEED 
      REAL*8     SDS_SST(NUMCELLS_B),SST
      BYTE       SDS_CldMskQA(iswath_B*2),QA_Temp
