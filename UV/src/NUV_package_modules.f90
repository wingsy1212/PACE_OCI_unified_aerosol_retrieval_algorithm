Module NUV_PACKAGE_MODULES

USE HDF5
USE MyConstants
USE H5Util_class
USE H5write_module

USE LookupTableModule                  
USE InterpolationModule
USE OCIUAAER_Config_Module
USE OCIUAAER_L1BModule
USE GetSurfAlbLUT_H5module
USE GetLUT_LER_AI_H5module
USE GetLUT_Miecloud_AI_H5module 
USE Get_SnowIce_module
USE Get_OceanLUT_H5module 
USE Get_TerrainPressure_H5module
USE GetLUT_H5module
USE OCI_UVAI_MieCloudModule 
USE OCI_UVOutput_DataModule
USE Nearuv_alhssa_Module
USE NUV_AerosolModule

USE Get_ssaclimLUT2_H5module
USE Get_omacaLUT7dim_H5module 
USE NUV_ACAerosolModule


IMPLICIT NONE

  CHARACTER(LEN=255) :: NUV_outfile
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: SurfaceType
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: AerosolType
  REAL(KIND=4), DIMENSION(:,:)  ,  ALLOCATABLE :: RelAzAng, TerPress, Dstsmk_Index
  REAL(KIND=4), DIMENSION(:,:,:),  ALLOCATABLE :: UVSurfaceAlbedo, UVIR_AllRef, NUV_Reflectivity
  REAL(KIND=4), DIMENSION(:,:)  ,  ALLOCATABLE :: NUV_Residue, SnowIce_frac

! -- ACA data fields --------------------
  REAL (KIND=4),   DIMENSION(:,:,:),   ALLOCATABLE :: InputAerosolSingleScattAlbACA
  REAL (KIND=4),   DIMENSION(:,:),     ALLOCATABLE :: ApparentCloudOpticalDepth

CONTAINS

!!
!==============================================================================
!==============================================================================
!!
SUBROUTINE NUV_package_main(cfg, Year, Month, Day, UVtoSWIR_wavelengths, UVtoSWIR_Reflectances, &
	  Ret_Lat,Ret_Lon,Ret_SolZen,Ret_View_angle,&
          Ret_View_phi,Ret_solar_phi,Ret_Xtrack,Ret_Lines,& 
          Ret_Small_weighting,Land_sea_flag,Ret_ref_LandOcean,Ret_Tau_LandOcean,&
	  Ret_average_Omega_Ocean_UV,Ret_Index_Height,Cloud_Frac_LandOcean,&
	  NUV_AI, NUV_COD, NUV_CldFrac, NUV_SSA, NUV_ALH, &
	  NUV_ACAOD, NUV_AerCorrCOD, &
	  NUV_FinalAlgorithmFlagsACA, NUV_UncertaintyACAODToSSA, NUV_UncertaintyCODToSSA)

 
IMPLICIT NONE


Include 'output_Variables.inc'

TYPE(ociuaaer_config_type),    INTENT(IN)  :: cfg
INTEGER(KIND=4),               INTENT(IN)  :: Day, Month, Year
REAL(KIND=4),DIMENSION(:),     INTENT(IN)  :: UVtoSWIR_wavelengths
REAL(KIND=4),DIMENSION(:,:,:), INTENT(IN)  :: UVtoSWIR_Reflectances

  
 INTEGER(KIND=4) :: iPix, jPix, STATUS, SkipProcess, l1b_XDim, l1b_YDim
 REAL(KIND=4), PARAMETER :: DTOR = (PI/180.)
 REAL(KIND=4) :: cossza, airscocm
 INTEGER(KIND=4), PARAMETER :: uvWavel=2
 REAL(KIND=4),DIMENSION(1:uvWavel),PARAMETER :: uvWaves=(/354.0,388.0/)
 REAL(KIND=4),DIMENSION(1:uvWavel) :: NormRad, AllNormRad

 INTEGER(KIND=4)   :: index1
 CHARACTER(LEN=22) :: l1b_date_time_str
 CHARACTER(LEN=256)  :: NUV_outfile
  
!------------------------------------------------------------------------------
!  Variables used to reshape the radiance and transmittance matrices
!------------------------------------------------------------------------------
  REAL (KIND=4), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: aux
  REAL (KIND=4), DIMENSION(:,:,:,:,:,:),   ALLOCATABLE :: aux2
  INTEGER (KIND=4) :: linIndex(3), ialo

!------------------------------------------------------------------------------
! Read HE4 AI Mie-Atmosphere look-up tables
! GetLUT_Miecloud_AI_H5module 
! nplev_mie, nsalb_mie, ncod_mie, nsza_mie, nvza_mie, nraa_mie
! rad354, rad388, rad_lin354_ai, rad_lin388_ai
  CALL ReadLUTAIparams(cfg%uv_ai_mielut)


!------------------------------------------------------------------------------
! Read HE4 AI LER-Atmosphere look-up tables
! GetLUT_LER_AI_H5module 
! nsalb_ler, nplev_ler, nsza_ler, SURFALBSET_ler, PRESSURESET_ler
! rad354_ler, rad388_ler, rad_lin354_ai_ler, rad_lin388_ai_ler
  CALL ReadLUT_LER_AIparams(cfg%uv_ai_lerlut)

!------------------------------------------------------------------------------
! Read OCI Surface AlbeDO HE4 look-up tables (Version 5 : ocean-corrected LER)
  CALL ReadSurfAlbLUTparams(cfg%uv_surfalbfile,SRFLER354, SRFLER388, GRDLAT, GRDLON)


!------------------------------------------------------------------------------
! Read Terrain pressureH5 database [3600 longitude x 1800 latitude] bins
  CALL ReadTerrainPressparams(cfg%uv_surfprsfile,tpres, terrain_Longitude, terrain_Latitude)


!------------------------------------------------------------------------------
! Read Snow_Ice database [360 longitude x 180 latitude x 12 month] bins
  CALL snowice_Reader(cfg%uv_snowicefile,snowice, swice_Longitude, swice_Latitude)


!------------------------------------------------------------------------------
  CALL Read_OceanLUTparams(cfg%uv_ocncorr_lut, oceanler, nwave_ocean, &
                           nsza_ocean, nvza_ocean, nraa_ocean,&
                       sza_oceanset, vza_oceanset, raa_oceanset) !Fresnel_OceanLUT_v1.he4

!==============================================================================
! Read OMI Single Scattering Albedo Climatology HE4 look-up tables
!==============================================================================
  CALL Read_ssaclimLUTparams(cfg%uv_ssa388_clim_file, ssa388_smk,ssa388_dst, &
                   ssa388_clim_smk,ssa388_clim_dst,ssalon_table,ssalat_table)


!==============================================================================
! Read Aerosol Above Clouds HE4 look-up tables
!==============================================================================
  CALL Read_aac_LUTparams(cfg%uv_omacalut_file, &
  			  rad388_smokelin1013zhgt3_aac,rad388_smokelin1013zhgt4_aac,&
                          rad388_smokelin1013zhgt5_aac,rad388_smokelin1013zhgt6_aac,&
                          rad388_smokelin800zhgt5_aac,rad388_smokelin800zhgt6_aac,&
                          rad388_smokelin800zhgt7_aac,rad388_smokelin800zhgt8_aac,&
                          uvaimie_smokelin1013zhgt3_aac,uvaimie_smokelin1013zhgt4_aac,&
                          uvaimie_smokelin1013zhgt5_aac,uvaimie_smokelin1013zhgt6_aac,&
                          uvaimie_smokelin800zhgt5_aac,uvaimie_smokelin800zhgt6_aac,&
                          uvaimie_smokelin800zhgt7_aac,uvaimie_smokelin800zhgt8_aac,&

!                        ;DUST LUT Parameters...
                          rad388_dustlin1013zhgt3_aac,rad388_dustlin1013zhgt4_aac,&
                          rad388_dustlin1013zhgt5_aac,rad388_dustlin1013zhgt6_aac,&
                          rad388_dustlin800zhgt5_aac,rad388_dustlin800zhgt6_aac,&
                          rad388_dustlin800zhgt7_aac,rad388_dustlin800zhgt8_aac,&
                          uvaimie_dustlin1013zhgt3_aac,uvaimie_dustlin1013zhgt4_aac,&
                          uvaimie_dustlin1013zhgt5_aac,uvaimie_dustlin1013zhgt6_aac,&
                          uvaimie_dustlin800zhgt5_aac,uvaimie_dustlin800zhgt6_aac,&
                          uvaimie_dustlin800zhgt7_aac,uvaimie_dustlin800zhgt8_aac,&

!                         DUST OVER-OCEAN LUT Parameters...
                          rad388_dustolin1013zhgt3_aac,rad388_dustolin1013zhgt4_aac,&
                          rad388_dustolin1013zhgt5_aac,rad388_dustolin1013zhgt6_aac,&
                          rad388_dustolin800zhgt5_aac,rad388_dustolin800zhgt6_aac,&
                          rad388_dustolin800zhgt7_aac,rad388_dustolin800zhgt8_aac,&
                          uvaimie_dustolin1013zhgt3_aac,uvaimie_dustolin1013zhgt4_aac,&
                          uvaimie_dustolin1013zhgt5_aac,uvaimie_dustolin1013zhgt6_aac,&
                          uvaimie_dustolin800zhgt5_aac,uvaimie_dustolin800zhgt6_aac,&
                          uvaimie_dustolin800zhgt7_aac,uvaimie_dustolin800zhgt8_aac,&
                          wavetbl_aac, aodtbl_aac, codtbl_aac,szatbl_aac,vzatbl_aac,&
                          raatbl_aac,ssatbl_aac, salbtbl_aac)

!------------------------------------------------------------------------------
!Read Aerosol and ancillary LUTs
  status = allocateLUT()
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error allocating memory for Aerosol LUT params"
     	  CALL EXIT(1)
  	ENDIF 

  CALL ReadLUTparams(cfg%uv_aer_lut) !LUT_v5_2wave.he4
  CALL Read_ancillaryLUTs(cfg%uv_sfc_file, cfg%uv_AIRSCO_clm_file)

!------------------------------------------------------------------------------
! Redefines radiance and transmittance lookup tables matrices as vectors 

   ALLOCATE(radlin_p10(nwave*nzae*nw0model*ntau*nsza*nphi*ntheta), &
            radlin_p06(nwave*nzae*nw0model*ntau*nsza*nphi*ntheta))
   ALLOCATE(tralin_p10(nwave*nzae*nw0model*ntau*nsza*ntheta), &
            tralin_p06(nwave*nzae*nw0model*ntau*nsza*ntheta))
   ALLOCATE(       aux(nwave,nzae,nw0model,ntau,nsza,nphi,ntheta), &
                  aux2(nwave,nzae,nw0model,ntau,nsza,ntheta) )

  ! ------- Radiance -------------
  DO ialo = 1,2
     linIndex(ialo) = ialo
     aux(ialo,:,:,:,:,:,:) = radp10(linIndex(ialo),:,:,:,:,:,:)
  ENDDO
  radlin_p10(:)  = RESHAPE(aux, (/nwave*nzae*nw0model*ntau*nsza*nphi*ntheta/) )

  DO ialo = 1,2
     linIndex(ialo) = ialo
     aux(ialo,:,:,:,:,:,:) = radp6(linIndex(ialo),:,:,:,:,:,:)
  ENDDO
  radlin_p06(:)  = RESHAPE(aux, (/nwave*nzae*nw0model*ntau*nsza*nphi*ntheta/) )

  ! ------- Transmittance -------------    
  DO ialo = 1,2
     linIndex(ialo) = ialo
     aux2(ialo,:,:,:,:,:) = trp10(linIndex(ialo),:,:,:,:,:)
  ENDDO
  tralin_p10(:)  = RESHAPE(aux2, (/nwave*nzae*nw0model*ntau*nsza*ntheta/) )

  DO ialo = 1,2
     linIndex(ialo) = ialo
     aux2(ialo,:,:,:,:,:) = trp6(linIndex(ialo),:,:,:,:,:)
  ENDDO
  tralin_p06(:) = RESHAPE(aux2,(/nwave*nzae*nw0model*ntau*nsza*ntheta/) )


!------------------------------------------------------------------------------
! Allocate for variables to work with.
	ALLOCATE( RelAzAng(Ret_Xtrack, Ret_Lines), &
	          TerPress(Ret_Xtrack, Ret_Lines), &
	      SnowIce_frac(Ret_Xtrack, Ret_Lines), &
	       SurfaceType(Ret_Xtrack, Ret_Lines), &
	       AerosolType(Ret_Xtrack, Ret_Lines), &
	       NUV_Residue(Ret_Xtrack, Ret_Lines), &
          NUV_Reflectivity(Ret_Xtrack, Ret_Lines, 2), &	       
           UVSurfaceAlbedo(Ret_Xtrack, Ret_Lines, 2), &
	       UVIR_AllRef(Ret_Xtrack, Ret_Lines, 3), &
	      Dstsmk_Index(Ret_Xtrack, Ret_Lines) )


!------------------------------------------------------------------------------
! Allocate memory to ACA Fields....
    ALLOCATE(InputAerosolSingleScattAlbACA(Ret_Xtrack, Ret_Lines, 3), &
    	         ApparentCloudOpticalDepth(Ret_Xtrack, Ret_Lines) )


!------------------------------------------------------------------------------
! Intialize all variables
NUV_AI(:,:) = -999.
NUV_COD(:,:) = -999.
NUV_CldFrac(:,:) = -999.
NUV_ALH(:,:) = -999.
NUV_SSA(:,:,:) = -999.
!
RelAzAng(:,:) = -999.
TerPress(:,:) = -999.
SnowIce_frac(:,:) = -999.
SurfaceType(:,:) = 255
AerosolType(:,:) = 255
NUV_Residue(:,:) = -999.
NUV_Reflectivity(:,:,:) = -999.
UVSurfaceAlbedo(:,:,:) = -999.
UVIR_AllRef(:,:,:) = -999.
Dstsmk_Index(:,:) = -999.

! ACA-fields
NUV_ACAOD(:,:,:) = -999.
NUV_AerCorrCOD(:,:) = -999.
NUV_FinalAlgorithmFlagsACA(:,:) = -999.
ApparentCloudOpticalDepth(:,:) = -999.
InputAerosolSingleScattAlbACA(:,:,:) = -999.
NUV_UncertaintyACAODToSSA(:,:,:) = -999.
NUV_UncertaintyCODToSSA(:,:,:) = -999.


! Makes average all reflectances at aggregate resolution for 354, 388, 2230 nm.
 l1b_XDim = Size(UVtoSWIR_Reflectances, 2)
 l1b_YDim = Size(UVtoSWIR_Reflectances, 3)
 CALL make_aggres_allref(Ret_Xtrack, Ret_Lines, Ibox_Ret, &
                         UVtoSWIR_wavelengths, UVtoSWIR_Reflectances, &
                         UVIR_AllRef, Dstsmk_Index)


!------------------------------------------------------------------------------
! Start Pixels and Scan line loop
! Loop over Scan lines
DO jPix = 1, Ret_Lines 

! Loop over Pixels
    DO iPix = 1, Ret_Xtrack 

      IF (Ret_Lon(iPix, jPix) .GE. -180.0 .AND. Ret_Lon(iPix, jPix) .LE. 180.0 .AND. &
          Ret_Lat(iPix, jPix) .GE. -90.0  .AND. Ret_Lat(iPix, jPix) .LE. 90.0 .AND. &
          Ret_SolZen(iPix, jPix) .GE. 0.001 .AND. Ret_SolZen(iPix, jPix) .LE. 90.0 .AND. &
	  Ret_View_angle(iPix, jPix) .GE. 0.001 .AND. Ret_View_angle(iPix, jPix) .LE. 90.0 .AND. &
	  Ret_View_phi(iPix, jPix) .GE. -180.001 .AND. Ret_View_phi(iPix, jPix) .LE. 180.0 .AND. &
	  Ret_solar_phi(iPix, jPix) .GE. -180.001 .AND. Ret_solar_phi(iPix, jPix) .LE. 180.0) THEN
         SkipProcess = 0
      ELSE
         SkipProcess = 1
	 CYCLE
      ENDIF


      CALL compute_relazimuth_ang(Ret_solar_phi(iPix,jPix), &
                                   Ret_View_phi(iPix,jPix), &
                                       RelAzAng(iPix,jPix))


      STATUS = Get_TerrainPressure(Ret_Lat(iPix,jPix), &
                                   Ret_Lon(iPix,jPix), &
                                  TerPress(iPix,jPix))
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting Terrain Pressure"
     	  CALL EXIT(1)
  	ENDIF 


      STATUS = GetSurfaceType(Ret_Lat(iPix,jPix), Ret_Lon(iPix,jPix), &
                                              SurfaceType(iPix,jPix))
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting SurfaceType"
     	  CALL EXIT(1)
  	ENDIF 


      STATUS = Get_snowice_fraction(Month, Ret_Lat(iPix,jPix), &
                                           Ret_Lon(iPix,jPix), &
                                      SnowIce_frac(iPix,jPix))
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting Snow/Ice Fraction"
     	  CALL EXIT(1)
  	ENDIF 
                                   

      STATUS = GetAIRSCO_Clm(Month, Ret_Lat(iPix,jPix), &
                                    Ret_Lon(iPix,jPix), &
				    airscocm)
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting AIRS Clm CO"
     	  CALL EXIT(1)
  	ENDIF 

                                   
      STATUS = Get_UVSurfaceAlbedo(Month, Ret_Lat(iPix,jPix), &
                                          Ret_Lon(iPix,jPix), &
                                  UVSurfaceAlbedo(iPix,jPix,:))
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting UV Surface albedo"
     	  CALL EXIT(1)
  	ENDIF 


      STATUS = OCEAN_FRESNEL_CORRECTION(SurfaceType(iPix,jPix), Ret_SolZen(iPix,jPix), &
                                     Ret_View_angle(iPix,JPix),   RelAzAng(iPix,jPix), &
				  			   UVSurfaceAlbedo(iPix,jPix,:)) 
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting Fresnel correction for UV Surface albedo"
     	  CALL EXIT(1)
  	ENDIF 
      
      
      cossza = COS(DTOR * Ret_SolZen(iPix,jPix))
      NormRad(1) = (Ret_ref_LandOcean(iPix,jPix,1)*cossza)/PI
      NormRad(2) = (Ret_ref_LandOcean(iPix,jPix,2)*cossza)/PI
      AllNormRad(1) = (UVIR_AllRef(iPix,jPix,1)*cossza)/PI
      AllNormRad(2) = (UVIR_AllRef(iPix,jPix,2)*cossza)/PI


      ! Compute UVAerosol Index if ref > 0 (All reflectance, cloud+aerosol).
      IF( AllNormRad(1) .GT. 0 .AND. AllNormRad(2) .GT. 0 ) THEN 
          CALL calc_coeff(theta_table,sza_table,phi_table,ntheta,nsza,nphi, &
                          			 Ret_View_angle(iPix,jPix), &
                            			     Ret_SolZen(iPix,jPix), &
                        		               RelAzAng(iPix,jPix) )

          CALL calc_coeff_2dim(theta_table,sza_table,ntheta,nsza, &
                          	       Ret_View_angle(iPix,jPix), &
                            	           Ret_SolZen(iPix,jPix) )

          ! --- Mie coefficients ----------
          CALL calc_coeff_Mie(theta_table_ep, sza_table_ep, phi_table_ep, &
                           		    nvza_mie, nsza_mie, nraa_mie, &
                          		       Ret_View_angle(iPix,jPix), &
                            		           Ret_SolZen(iPix,jPix), &
                        			     RelAzAng(iPix,jPix) )

          STATUS = UVAI_Miecloud_perpixel( Ret_Lat(iPix,jPix), &
                                           Ret_Lon(iPix,jPix), &
                                        Ret_SolZen(iPix,jPix), &
                                    Ret_View_angle(iPix,jPix), &
                                          RelAzAng(iPix,jPix), &
                                          TerPress(iPix,jPix), &
                                      SnowIce_frac(iPix,jPix), &
                                 UVSurfaceAlbedo(iPix,jPix,:), &
                                              AllNormRad(1:2), & ! 354, 388
                                            NUV_AI(iPix,jPix), &
                                           NUV_COD(iPix,jPix), &
                                       NUV_CldFrac(iPix,jPix), &
                                       NUV_Residue(iPix,jPix), &
                                  NUV_Reflectivity(iPix,jPix,:) )

      ELSE
        SkipProcess = 0
        CYCLE
      ENDIF ! AllNormRad > 0 loop



      ! Cloud-free Aerosol Retrievals
      ! Compute SSA and ALH only if ref > 0 and AOD > 0 (cloud free reflectance provided by DB/DT).
      IF( NormRad(1) .GT. 0 .AND. NormRad(2) .GT. 0 .AND. &
              Ret_Tau_LandOcean(iPix,jPix,1) .GT. 0 .AND. &
              Ret_Tau_LandOcean(iPix,jPix,2) .GT. 0 .AND. &
	      Ret_Tau_LandOcean(iPix,jPix,4) .GT. 0.2) THEN  ! AOD(550) > 0.2, Ignore lower AODs

 	  Dstsmk_Index(iPix,jPix) = &
	  -10.0 * log10(Ret_ref_LandOcean(iPix,jPix,2)/Ret_ref_LandOcean(iPix,jPix,9))
             STATUS = OCI_NUV_Aer_Process(Month, uvWavel, uvWaves(:), &
                                           Ret_Lat(iPix,jPix), &
                                           Ret_Lon(iPix,jPix), &
                                        Ret_SolZen(iPix,jPix), &
                                    Ret_View_angle(iPix,jPix), &
                                          RelAzAng(iPix,jPix), &
                                          TerPress(iPix,jPix), &
                                       SurfaceType(iPix,jPix), &
                                 UVSurfaceAlbedo(iPix,jPix,:), &
                                                 NormRad(1:2), & ! 354, 388
				       AerosolType(iPix,jPix), &
                                            NUV_AI(iPix,jPix), &
                                      Dstsmk_Index(iPix,jPix), &					    
                             Ret_Tau_LandOcean(iPix,jPix,1:2), & ! 354, 388
                               Ret_Small_weighting(iPix,jPix), &
                                         NUV_SSA(iPix,jPix,:), &
                                           NUV_ALH(iPix,jPix) )

      ENDIF ! NormRad > 0 Aer-loop

      ! Above-cloud retrievals
      IF( AllNormRad(1) .GT. 0 .AND. AllNormRad(2) .GT. 0 .AND. &
             Ret_Tau_LandOcean(iPix,jPix,1) .LE. 0 .AND. &
             Ret_Tau_LandOcean(iPix,jPix,2) .LE. 0 .AND. &
	     Cloud_Frac_LandOcean(iPix,jPix) .GE. 0.9 .AND. &
	     airscocm .GT. 0) THEN 

              STATUS = OCI_NUV_ACAer_Process(Year, Month, Day, &
	  		        uvWavel, uvWaves(:), airscocm, &
                                           Ret_Lat(iPix,jPix), &
                                           Ret_Lon(iPix,jPix), &
                                        Ret_SolZen(iPix,jPix), &
                                    Ret_View_angle(iPix,jPix), &
                                          RelAzAng(iPix,jPix), &
                                          TerPress(iPix,jPix), &
                                       SurfaceType(iPix,jPix), &
                                 UVSurfaceAlbedo(iPix,jPix,:), &
                                              AllNormRad(1:2), & ! 354, 388
                                            NUV_AI(iPix,jPix), &
         			NUV_Reflectivity(iPix,jPix,:), &
                        NUV_FinalAlgorithmFlagsACA(iPix,jPix), &
                    		       NUV_ACAOD(iPix,jPix,:), &
                      		    NUV_AerCorrCOD(iPix,jPix), &
                         ApparentCloudOpticalDepth(iPix,jPix), &
                   InputAerosolSingleScattAlbACA(iPix,jPix,:), &
                       NUV_UncertaintyACAODToSSA(iPix,jPix,:), &
                         NUV_UncertaintyCODToSSA(iPix,jPix,:) )

      ENDIF ! NormRad > 0 ACAer-loop

      !
    ENDDO   ! iPix = 1, XDim-1
!
ENDDO    ! jPix = 1, YDim-1
! End Pixels and Scan line loop



! Now write NUV output file.
   IF (cfg%input_l1file /= 'NULL') THEN
     index1 = index(cfg%input_l1file, 'PACE_OCI_SIM')
     l1b_date_time_str = cfg%input_l1file(index1+13:index1+13+14)
   ELSE IF (cfg%proxy_l1file /= 'NULL') THEN
     index1 = index(cfg%proxy_l1file, 'TROP-in-Viirs_V4.1_A')
     l1b_date_time_str = cfg%proxy_l1file(index1+20:index1+20+11)
   ENDIF
       go to 1000
 NUV_outfile = trim('../L2_data/NUV_') // trim(l1b_date_time_str) // trim('.h5')   

    CALL write_NUV_output(NUV_outfile, RelAzAng, TerPress, Dstsmk_Index, &
	  Ret_Lat,Ret_Lon,Ret_SolZen,Ret_View_angle,&
          Ret_View_phi,Ret_solar_phi,Ret_Xtrack,Ret_Lines,& 
          Ret_Small_weighting,Land_sea_flag,Ret_ref_LandOcean,Ret_Tau_LandOcean,&
	  Ret_average_Omega_Ocean_UV,Ret_Index_Height,Cloud_Frac_LandOcean,&
	  NUV_AI, NUV_COD, NUV_CldFrac, NUV_SSA, NUV_ALH, AerosolType, &
	  NUV_ACAOD, NUV_AerCorrCOD, UVIR_AllRef, &
	  ApparentCloudOpticalDepth, InputAerosolSingleScattAlbACA, &
	  NUV_FinalAlgorithmFlagsACA, NUV_UncertaintyACAODToSSA, NUV_UncertaintyCODToSSA)

1000  continue

END SUBROUTINE NUV_package_main
!!
!==============================================================================
!==============================================================================
!!

SUBROUTINE make_aggres_allref(Ret_Xtrack, Ret_Lines, Ibox_Ret, &
                         UVtoSWIR_wavelengths, UVtoSWIR_Reflectances, &
                         UVIR_AllRef, Dstsmk_Index)

IMPLICIT NONE

INTEGER(KIND=4),               INTENT(IN)  :: Ret_Xtrack, Ret_Lines, Ibox_Ret
REAL(KIND=4),DIMENSION(:),     INTENT(IN)  :: UVtoSWIR_wavelengths
REAL(KIND=4),DIMENSION(:,:,:), INTENT(IN)  :: UVtoSWIR_Reflectances
REAL(KIND=4),DIMENSION(:,:,:), INTENT(OUT) :: UVIR_AllRef
REAL(KIND=4),DIMENSION(:,:),   INTENT(OUT) :: Dstsmk_Index

INTEGER(KIND=4)  :: XDim, YDim, num
INTEGER(KIND=4)  :: i, j, k, k1, x1, x2, y1, y2
REAL(KIND=4)     :: tmparr(Ibox_Ret, Ibox_Ret)

num = 1            
XDim = size(UVtoSWIR_Reflectances, 2)
YDim = size(UVtoSWIR_Reflectances, 3)

!PRINT *, 'AggRes subroutine = ', XDim, YDim, Ibox_Ret, Ret_Xtrack, Ret_Lines, UVtoSWIR_wavelengths((/2, 3, 14/))
!PRINT *, 'AggRes subroutine = ', XDim/Ibox_Ret, YDim/Ibox_Ret, Ret_Xtrack, Ret_Lines

DO i = 1, Ret_Xtrack
DO j = 1, Ret_Lines

x1 = i
y1 = j
IF (x1 /= 1) x1 = ((i-1) * Ibox_Ret) + 1
IF (y1 /= 1) y1 = ((j-1) * Ibox_Ret) + 1
x2 = x1+Ibox_Ret-1
y2 = y1+Ibox_Ret-1

!Print *, i, j, x1, x2, y1, y2, XDim, YDim

DO k = 1, 3 
  IF (k == 1) k1 = 2 	! 354
  IF (k == 2) k1 = 3	! 388
  IF (k == 3) k1 = 14	! 2250 
  tmparr(:,:) = UVtoSWIR_Reflectances(k1, x1:x2, y1:y2)
  num = (max(1,count(tmparr>0))) 
  !
  IF (num .GT. 1) THEN 
    UVIR_AllRef(i, j, k) = sum(tmparr, tmparr>0)/num
  ELSE
    UVIR_AllRef(i, j, k) = -999.
  ENDIF
  !    
  !Print *, ' '
  !Print *, k, i, j, num
  !Print *, 'Array = ', tmparr
  !Print *, 'Average = ',UVIR_AllRef(i, j, k)  
END DO ! k-loop

IF (UVIR_AllRef(i,j,2) .GT. 0 .AND. UVIR_AllRef(i,j,3) .GT. 0) THEN 
 !Dstsmk_Index(i,j) = -10.0 * LOG10(UVIR_AllRef(i,j,2)/UVIR_AllRef(i,j,3))
ELSE 
 Dstsmk_index(i,j) = -999.
ENDIF

END DO ! j-loop

END DO ! i-loop


END SUBROUTINE make_aggres_allref

!!
!==============================================================================
!==============================================================================
!!

SUBROUTINE write_NUV_output(NUV_outfile, RelAzAng, TerPress, Dstsmk_Index, &
	  Ret_Lat,Ret_Lon,Ret_SolZen,Ret_View_angle,&
          Ret_View_phi,Ret_solar_phi,Ret_Xtrack,Ret_Lines,& 
          Ret_Small_weighting,Land_sea_flag,Ret_ref_LandOcean,Ret_Tau_LandOcean,&
	  Ret_average_Omega_Ocean_UV,Ret_Index_Height,Cloud_Frac_LandOcean,&
	  NUV_AI, NUV_COD, NUV_CldFrac, NUV_SSA, NUV_ALH, AerosolType, &
	  NUV_ACAOD, NUV_AerCorrCOD, UVIR_AllRef, &
	  ApparentCloudOpticalDepth, InputAerosolSingleScattAlbACA, &
	  NUV_FinalAlgorithmFlagsACA, NUV_UncertaintyACAODToSSA, NUV_UncertaintyCODToSSA)


  IMPLICIT NONE

  Include 'output_Variables.inc'

  CHARACTER(LEN=255), INTENT(IN) :: NUV_outfile
  REAL(KIND=4), DIMENSION(:,:),   INTENT(IN) :: RelAzAng, TerPress, Dstsmk_Index
  REAL(KIND=4), DIMENSION(:,:,:), INTENT(IN) :: UVIR_AllRef
  INTEGER(KIND=2), DIMENSION(:,:),INTENT(IN) :: AerosolType

  REAL (KIND=4),   DIMENSION(:,:,:), INTENT(IN) :: InputAerosolSingleScattAlbACA
  REAL (KIND=4),   DIMENSION(:,:),   INTENT(IN) :: ApparentCloudOpticalDepth
!  INTEGER(KIND=2), DIMENSION(:,:),   INTENT(IN) :: NUV_FinalAlgorithmFlagsACA
!  REAL (KIND=4),   DIMENSION(:,:,:), INTENT(IN) :: NUV_UncertaintyACAODToSSA
!  REAL (KIND=4),   DIMENSION(:,:,:), INTENT(IN) :: NUV_UncertaintyCODToSSA
  
  INTEGER(HID_T) :: nuv_file_id
  INTEGER(HID_T) :: plist2_geo
  INTEGER(HID_T) :: plist3_wav2, plist3_wav3, plist3_wav5, plist3_wav9
  INTEGER(HID_T) :: group_geo, group_sci
  INTEGER :: hdferr1, error
  INTEGER(HSIZE_T) :: cXT, cLT, wav2, wav3, wav5, wav9

  REAL(KIND=4), ALLOCATABLE :: tmp2darr(:,:), tmp3darr(:,:,:)
!
  cLT = Ret_Lines  ! chunk along Track 404
  cXT = Ret_Xtrack    ! chunk Xtrack 400
  wav9 = 9
  wav5 = 5
  wav3 = 3
  wav2 = 2
  

  CALL h5open_f(hdferr1)
  CALL h5fcreate_f(TRIM(NUV_outfile), H5F_ACC_TRUNC_F, nuv_file_id, hdferr1)
!
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist2_geo, error)
  CALL h5pset_chunk_f(  plist2_geo, 2, (/cXT, cLT/), error)
  CALL h5pset_deflate_f(plist2_geo, 6, error) 
!
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist3_wav9, error)
  CALL h5pset_chunk_f(  plist3_wav9, 3, (/cXT, cLT, wav9/), error)
  CALL h5pset_deflate_f(plist3_wav9, 6, error) 
!
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist3_wav5, error)
  CALL h5pset_chunk_f(  plist3_wav5, 3, (/cXT, cLT, wav5/), error)
  CALL h5pset_deflate_f(plist3_wav5, 6, error) 
!
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist3_wav3, error)
  CALL h5pset_chunk_f(  plist3_wav3, 3, (/cXT, cLT, wav3/), error)
  CALL h5pset_deflate_f(plist3_wav3, 6, error) 
!
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist3_wav2, error)
  CALL h5pset_chunk_f(  plist3_wav2, 3, (/cXT, cLT, wav2/), error)
  CALL h5pset_deflate_f(plist3_wav2, 6, error) 

  
! Geolocation data
  CALL H5gcreate_f(nuv_file_id, "/GeolocationData", group_geo, hdferr1)
  ALLOCATE(tmp2darr(cXT,cLT))
  tmp2darr(:,:) = Ret_Lat(1:cXT, 1:cLT)
  CALL H5write_data(group_geo, "Latitude",            tmp2darr, dcpl = plist2_geo)
  tmp2darr(:,:) = Ret_Lon(1:cXT, 1:cLT)
  CALL H5write_data(group_geo, "Longitude",           tmp2darr, dcpl = plist2_geo)
  tmp2darr(:,:) = Ret_SolZen(1:cXT, 1:cLT)
  CALL H5write_data(group_geo, "SolarZenithAngle",    tmp2darr, dcpl = plist2_geo)
  tmp2darr(:,:) = Ret_solar_phi(1:cXT, 1:cLT)
  CALL H5write_data(group_geo, "SolarAzimuthAngle",   tmp2darr, dcpl = plist2_geo)
  tmp2darr(:,:) = Ret_View_angle(1:cXT, 1:cLT)
  CALL H5write_data(group_geo, "ViewingZenithAngle",  tmp2darr, dcpl = plist2_geo)
  tmp2darr(:,:) = Ret_View_phi(1:cXT, 1:cLT)
  CALL H5write_data(group_geo, "ViewingAzimuthAngle", tmp2darr, dcpl = plist2_geo)
  CALL H5write_data(group_geo, "RelativeAzimuthAngle", RelAzAng, dcpl = plist2_geo)
  CALL H5write_data(group_geo, "TerrainHeight", TerPress, dcpl = plist2_geo)
  CALL h5gclose_f( group_geo, hdferr1)
  

! Geophysical data
  CALL H5gcreate_f(nuv_file_id, "/ScienceData", group_sci, hdferr1)
  tmp2darr(:,:) = NUV_AI(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "NUV_AerosolIndex", tmp2darr, dcpl = plist2_geo)

  tmp2darr(:,:) = NUV_COD(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "NUV_CloudOpticalDepth", tmp2darr, dcpl = plist2_geo)

  tmp2darr(:,:) = NUV_CldFrac(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "NUV_RadiativeCloudFraction", tmp2darr, dcpl = plist2_geo)

  tmp2darr(:,:) = NUV_ALH(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "NUV_AerosolLayerHeight", tmp2darr, dcpl = plist2_geo)

  tmp2darr(:,:) = AerosolType(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "NUV_DSIAerosolType", tmp2darr, dcpl = plist2_geo)

  ALLOCATE(tmp3darr(cXT,cLT,wav5))
  tmp3darr(:,:,:) = NUV_SSA(1:cXT, 1:cLT, 1:wav5)  
  CALL H5write_data(group_sci, "NUV_AerosolSingleScattAlbedo", tmp3darr, dcpl = plist3_wav5)
  DEALLOCATE(tmp3darr)

  ALLOCATE(tmp3darr(cXT,cLT,wav3))
  tmp3darr(:,:,:) = NUV_ACAOD(1:cXT, 1:cLT, 1:wav3)  
  CALL H5write_data(group_sci, "NUV_AerosolOpticalDepthOverCloud", tmp3darr, dcpl = plist3_wav3)
  DEALLOCATE(tmp3darr)

!  PRINT *, ' WRITE H5 ==> NUV_AerCorrCOD(161, 1:2)', NUV_AerCorrCOD(161, 1:2)
  tmp2darr(:,:) = NUV_AerCorrCOD(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "NUV_AerosolCorrCloudOpticalDepth", tmp2darr, dcpl = plist2_geo)

  !PRINT *, 'Shape of Ret_Ref_LandOcean = ', Shape(Ret_Ref_LandOcean), cXT,cLT,wav9
  ALLOCATE(tmp3darr(cXT,cLT,wav9))
  tmp3darr(:,:,:) = Ret_Ref_LandOcean(1:cXT, 1:cLT, 1:wav9)
  CALL H5write_data(group_sci, "Mean_Reflectance",  tmp3darr, dcpl = plist3_wav9)
  DEALLOCATE(tmp3darr)
  
  !PRINT *, 'Shape of Ret_Tau_LandOcean = ', Shape(Ret_Tau_LandOcean), cXT,cLT,wav9
  ALLOCATE(tmp3darr(cXT,cLT,wav9))
  tmp3darr(:,:,:) = Ret_Tau_LandOcean(1:cXT, 1:cLT, 1:wav9) 
  CALL H5write_data(group_sci, "Aerosol_Optical_Depth", tmp3darr, dcpl = plist3_wav9)
  DEALLOCATE(tmp3darr)
  
  ALLOCATE(tmp3darr(cXT,cLT,wav3))
  tmp3darr(:,:,:) = UVIR_AllRef(1:cXT, 1:cLT, 1:wav3) 
  CALL H5write_data(group_sci, "UVIR_AllRef", tmp3darr, dcpl = plist3_wav3)
  DEALLOCATE(tmp3darr)
  
  tmp2darr(:,:) = Ret_Small_weighting(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "Fine-mode_Fraction_Ocean", tmp2darr, dcpl = plist2_geo)

  tmp2darr(:,:) = Dstsmk_Index(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "NUV_DustSmokeIndex", tmp2darr, dcpl = plist2_geo)

  tmp2darr(:,:) = ApparentCloudOpticalDepth(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "NUV_ApparentCloudOpticalDepth", tmp2darr, dcpl = plist2_geo)

  tmp2darr(:,:) = NUV_FinalAlgorithmFlagsACA(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "NUV_FinalAlgorithmFlagsACA", tmp2darr, dcpl = plist2_geo)

  ALLOCATE(tmp3darr(cXT,cLT,wav3))
  tmp3darr(:,:,:) = InputAerosolSingleScattAlbACA(1:cXT, 1:cLT, 1:wav3) 
  CALL H5write_data(group_sci, "NUV_InputAerosolSingleScattAlbACA", tmp3darr, dcpl = plist3_wav3)
  DEALLOCATE(tmp3darr)

  ALLOCATE(tmp3darr(cXT,cLT,wav2))
  tmp3darr(:,:,:) = NUV_UncertaintyCODToSSA(1:cXT, 1:cLT, 1:wav2) 
  CALL H5write_data(group_sci, "NUV_UncertaintyCODToSSA", tmp3darr, dcpl = plist3_wav2)
  DEALLOCATE(tmp3darr)

  ALLOCATE(tmp3darr(cXT,cLT,wav2))
  tmp3darr(:,:,:) = NUV_UncertaintyACAODToSSA(1:cXT, 1:cLT, 1:wav2) 
  CALL H5write_data(group_sci, "NUV_UncertaintyACAODToSSA", tmp3darr, dcpl = plist3_wav2)
  DEALLOCATE(tmp3darr)

  ALLOCATE(tmp3darr(cXT,cLT,wav3))
  tmp3darr(:,:,:) = Ret_average_Omega_Ocean_UV(1:cXT, 1:cLT, 1:wav3) 
  CALL H5write_data(group_sci, "DT_AerosolSingleScattAlbedo", tmp3darr, dcpl = plist3_wav3)
  DEALLOCATE(tmp3darr)

  tmp2darr(:,:) = Ret_Index_Height(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "DT_AerosolLayerHeight", tmp2darr, dcpl = plist2_geo)

  tmp2darr(:,:) = Cloud_Frac_LandOcean(1:cXT, 1:cLT)
  CALL H5write_data(group_sci, "Aerosol_Cld_Fraction_Land_Ocean", tmp2darr, dcpl = plist2_geo)

  CALL h5gclose_f( group_sci, hdferr1)

  CALL h5fclose_f(nuv_file_id, hdferr1)
  CALL h5close_f(hdferr1)

!  PRINT *, ' NUV Output file : ', NUV_outfile
     
     
END SUBROUTINE write_NUV_output

!!
!==============================================================================
!==============================================================================
!!

End Module NUV_PACKAGE_MODULES


