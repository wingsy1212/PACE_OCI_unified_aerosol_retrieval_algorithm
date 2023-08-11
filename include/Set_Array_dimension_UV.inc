   CHARACTER*132 LINE
      Integer   Num_heigt,Num_iomega
      real      Threshold_error
      Parameter  (Num_heigt = 4 ,Num_iomega =4,Threshold_error = 5.0)
      Real Wind(NWIND)
      Real tau55um_LUT(NTAU),Height_km(Num_heigt) 
      DATA tau55um_LUT(1),tau55um_LUT(2),tau55um_LUT(3),&
           tau55um_LUT(4),tau55um_LUT(5),tau55um_LUT(6)&
           /0.0,0.20,0.5,1.0,2.0,3.0/
      DATA  Height_km(1),Height_km(2),Height_km(3),Height_km(4)&
          /1.5,3.0,6.0,10.0/
        
      Real EXTSMALL(NWIND,Numcases,NWAV),EXTbig(NWIND,NUMCASEB,NWAV)
      Real RGSS(NUMCASES),SIGMAS(NUMCASES)
      Real RGSB(NUMCASEB),SIGMAB(NUMCASEB)
      Real MOMENTSSMALL(NWIND,numcases,4),MOMENTSBIG(NWIND,NUMCASEB,4)
      Real CCNSMALL(NWIND,nUMCASES),TR
      Real EXTNORSMALL(NUMCASES,NWAV),EXTNORBIG(NUMCASEB,NWAV)
      Real BACKSCTTSMALL(NWIND,NUMCASES,NWAV),BACKSCTTBIG(NWIND,NUMCASEB,NWAV)
      Real ASSYMSMALL(NWIND,NUMCASES,NWAV),ASSYMBIG(NWIND,NUMCASEB,NWAV)
      Real ALBEDOSMALL(NWIND,NUMCASES,NWAV),ALBEDOBIG(NWIND,NUMCASEB,NWAV)
      Real ALBEDO_R_SMALL(NWIND,NTH0,NTAU,NWAV,NUMCASES)
      Real ALBEDO_R_BIG(NWIND,NTH0,NTAU,NWAV,NUMCASEB)
      Real ALBEDO_T_SMALL(NWIND,NTH0,NTAU,NWAV,NUMCASES)
      Real ALBEDO_T_BIG(NWIND,NTH0,NTAU,NWAV,NUMCASEB)
      Real ALBEDO_R_RAY(NWAV,NTH0),Tau_550_Mode_F,Tau_550_Mode_C
      Real ALBEDO_T_RAY(NWAV,NTH0),Icoeff(7),WAVE_MODIS(NWAV)
      Real  REFW354_input, REFW388_input,Ref_ray(NWAV)
      Real  PHC(NPHI),THET(NTHET),THET0(NTH0),WAVE(NWAV)
      Real Adummy(NWIND, NPHI, NTHET, NTH0, NTAU, NWAV,NUMCASEB) 
      Real  TAUAS(NUMCASES,NWAV,NTAU) 
      Real  TAUAB(NUMCASEB,NWAV,NTAU),DUMMY(NTH0),WSPEED
      Real  CCNdummy,ref_rayall_uv(4,NWIND,NPHI,NTHET,NTH0,NWAV,Num_heigt)
      Real  ref_rayall(NWIND,NPHI,NTHET,NTH0,NWAV) 
      Real  EFFRADSMALL(NUMCASES),EFFRADBIG(NUMCASEB),QSCT 
      Real MTHET0,MTHET,MPHI,MSCATT,Small_m_weighting,Tau_550,fmin 
      Real REFSMALL(NUMCASES,NWAV,NTAU),REFBIG(NUMCASEB,NWAV,NTAU)  
      Real  REFSMALLL_Uv(NWAV,NTAU),REFBIGL_uv(NWAV,NTAU),frac_fine_mode
      Real Ref_combined(NWIND,NWAV,NTAU),Ref_Ray_Uv(NWAV)
      Real  ref_allwav_uv(NWAV),Min_Uv,Array_min(4),Lat_OMI,Lon_OMI
      Real Taua_omi_354,Taua_omi_388
      Real AINTS_uv(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU, NWAV)
      Real AINTB_uv(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU, NWAV) 
      Real ALXWAV_uv(NWAV),luT_uv(4,NWAV),Tau_allwave(NWAV)
      Real  EXTSMALL_550,EXTbig_550
      Real  EXTSMALL_save(Num_iomega,Num_heigt,Numcases,NWAV_uv)	
      Real  ALBEDOSMALL_save(Num_iomega,Num_heigt,Numcases,NWAV_uv)  
      Real  EXTbig_save(Num_iomega,Num_heigt,NUMCASEB,NWAV_uv) 
      Real  ALBEDOBIG_save(Num_iomega,Num_heigt,NUMCASEB,NWAV_uv) 
      Real Omega_COMPUTED_uv(Num_iomega,NWAV_uv),tau_COMPUTED_uv(Num_iomega,NWAV_uv)
      Real  AAOD(4,NWAV_uv)
     Real slope,intercept,XA(10),YA(10),Uv2,A,B,Xi
     Real  Del_AOT388,AOT388,AOT344,Tau388,Uv1  
     Real H1_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
     Real H2_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
     Real H3_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
     Real H4_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
     Real H1_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
     Real H2_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
     Real H3_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
     Real H4_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES) 
     Real H1_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
     Real H2_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
     Real H3_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
     Real H4_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
     Real H1_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
     Real H2_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
     Real H3_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
     Real H4_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
     Real  Tau_height(Num_heigt,NWAV_uv),Omega_height(Num_heigt,NWAV_uv)
     Real Error_height(Num_heigt ,NWAV_uv)
      Real TAU_UV,Omega_new(NWAV_uv+1),X(100),Y(100),y1,Z(100),X0(100) 
      Real Extra_error_all(4*Num_heigt),&
           Extra_omega_all(4*Num_heigt,NWAV_uv),error_error,AVE,SDEV
       Real Extra_tau_all(4*Num_heigt,NWAV_uv),X1(100),Tau_new(NWAV_uv) 
       Real Fitting_error(NWAV_uv),Height_indx(NWAV_uv)
       Real SIG(10),maxsolz
       Real tau_H(Num_heigt,NWAV_uv),min_error(Num_heigt),Omega_H(Num_heigt,NWAV_uv)  
       Real Average_height,Height(Num_heigt) 
       Real lat_min,lat_max,lon_max,lon_min 
          Integer Icase,IFILE,IPHI,ITAU,ITH,ITH0,NN,jfile,Iscan,Idata,&
                  IWAV,IJ,Iopt_read,NNwave,ktau,Iomega, L_outFile,Inum_files
          Integer  KSZAM1,KSZAP1,KTHEM1,KTHEP1,KPHIM1,KPHIP1,Isens_test,Ret_option
          Integer JPHI(NPHI),IWS,iopt,set_flag_read,M_outFile
          Integer Mode_F,Mode_C,IWS1,IWS2,INDX(16) 
          Integer ifine,icoarse,itest,Iwind,Iheigt,LOPT,numscan
          Integer index_mode,Num_num,minsolz,Iheight_saved,Mode_Coarse
          Integer Index_omega(Num_heigt ,NWAV_uv),inumber,Index_omega_new(NWAV_uv)
          Real AE_44and88,a0,a1,a2,meanSSA388,DEL_Omega,SSA550_pred
          