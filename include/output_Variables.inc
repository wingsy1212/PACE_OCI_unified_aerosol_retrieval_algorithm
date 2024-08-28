     Integer       Ret_Xtrack,Ret_Lines
     Integer       IXRet_B,IYRet_B,Ibox_Ret,Num_Waves,Num_Wave_Land
     Integer       IX_B,IY_B              
     PARAMETER(IX_B=1272,IY_B=1722,Ibox_Ret =7,Num_Waves=9,Num_Wave_Land=4) 
     PARAMETER(IXRet_B=IX_B/Ibox_Ret,IYRet_B=IY_B/Ibox_Ret) 

! Define SDS data scales & offsets 
      DOUBLE PRECISION SCALE1,SCALE2,SCALE3,SCALE4,SCALE_10,&
                       OFFSET1,OFFSET2,OFFSET3,OFFSET4,OFFSET_10
      PARAMETER (SCALE1=1.0,SCALE2=100.0,SCALE3=1000.0,SCALE4=10000.0, SCALE_10=10.0,&
                 OFFSET1=0.0,OFFSET2=0.0,OFFSET3=0.0,OFFSET4=0.0, OFFSET_10=0.0) 
		 
!******************************************************************
        Real  Ret_Lat(IXRet_B,IYRet_B)
        Real  Ret_Lon(IXRet_B,IYRet_B)
        Real  Ret_SolZen(IXRet_B,IYRet_B)
        Real  Ret_View_angle(IXRet_B,IYRet_B)
        Real  Ret_View_phi(IXRet_B,IYRet_B)
        Real  Ret_solar_phi(IXRet_B,IYRet_B)  
        Real  Ret_Small_weighting(IXRet_B,IYRet_B)         
        Real  Ret_ref_LandOcean(IXRet_B,IYRet_B,Num_Waves)	
 !  not corrected  mean reflactance for 2 UV channels       
        Real  Ret_ref_LandOcean_UV(IXRet_B,IYRet_B,Num_Waves-7)
        Real  Ret_ref_LandOceanwOutUV(IXRet_B,IYRet_B,Num_Waves-2)
        Real  Ret_Tau_LandOcean(IXRet_B,IYRet_B,Num_Waves) 
        Real  Ret_tau_ocean(IXRet_B,IYRet_B,Num_Waves)
        Real  Ret_ref_ocean(IXRet_B,IYRet_B,Num_Waves)
        Real  Ret_ref_allwav_land(IXRet_B,IYRet_B,Num_Waves) 
        Real  Land_sea_Flag(IXRet_B,IYRet_B)
        Real     Ret_tau_land(IXRet_B,IYRet_B,Num_Wave_Land)
        Integer  Ret_CldMsk_500_Land(IX_B,IY_B)
        Integer   Cldmask_Native_LandOcean(IX_B,IY_B)
        Integer   CldMsk_Native_Land(IX_B,IY_B)
        Integer  CldMsk_Native_Ocean(IX_B*2,IY_B*2) 
        Real     Ret_land_Quality_Flag(IXRet_B,IYRet_B) 
        REAL     Ret_CLDFRC_land_DT(IXRet_B,IYRet_B)
        Real     Ret_ocean_Quality(IXRet_B,IYRet_B)
        Real     Ret_CLDFRC_Ocean(IXRet_B,IYRet_B)
        Real     Ret_Small_weighting_land(IXRet_B,IYRet_B) 
        Real     Ret_average_Omega_Ocean_UV(IXRet_B,IYRet_B,(Num_Wave_Land-1)) 
        Real     Ret_Index_Height(IXRet_B,IYRet_B)  
        Real     Cloud_Frac_LandOcean(IXRet_B,IYRet_B) 
        Real     Ret_Quality_LandOcean(IXRet_B,IYRet_B)
        Real     Ret_Quality_LandOcean_W0(IXRet_B,IYRet_B)
       
!******************************************************************
! Added NUV Output variables here - Vinay Kayetha       
        Real     NUV_AI(IXRet_B,IYRet_B)  
        Real     NUV_COD(IXRet_B,IYRet_B) 
        Real     NUV_CldFrac(IXRet_B,IYRet_B)     
	Real     UVResidue(IXRet_B,IYRet_B)
	Real     UVReflectivity(IXRet_B,IYRet_B,2)
	! Aerosol retreival output
        Real     NUV_ALH(IXRet_B,IYRet_B)  
        Real     NUV_SSA(IXRet_B,IYRet_B,5)   ! 354,388,480,550,670   
	! Above-cloud output
	Real     NUV_ACAOD(IXRet_B,IYRet_B,3) ! 354,388,550
	Real     NUV_ACAODVsHeight(IXRet_B,IYRet_B,3,5) ! 354,388,550
	Real     NUV_AerCorrCOD(IXRet_B,IYRet_B) ! 388
	Real     NUV_AerCorrCODVsHeight(IXRet_B,IYRet_B,5) ! 388
	Real 	 NUV_UncertaintyACAODToSSA(IXRet_B,IYRet_B,2)
	Real	 NUV_UncertaintyCODToSSA(IXRet_B,IYRet_B,2)
	Real     NUV_FinalAlgorithmFlagsACA(IXRet_B,IYRet_B)
