MODULE OCI_UVAI_MieCloudModule        

USE OCIUAAER_Config_Module
USE OCIUAAER_L1BModule
USE GetLUT_Miecloud_AI_H5module 
USE GetLUT_LER_AI_H5module
USE LookupTableModule                  
USE InterpolationModule
USE MyConstants

IMPLICIT NONE

 CHARACTER(LEN=256) :: lut_fn="xxxxxxxxxxxxxx"

CONTAINS

!!
!!
!==============================================================================
!==============================================================================
!!

SUBROUTINE OCI_UVAI_Miecloud(cfg, l1b_nXTrack, l1b_nLines, Month, &
                        Latitude, Longitude, SolarZenithAngle, ViewingZenithAngle, &
                        SolarAzimuthAngle, ViewingAzimuthAngle, TerrainHeight, &
			UVtoSWIR_Reflectances, UVAI)

 USE OCI_UV1km_DataModule


! -- Mie-AI related LUTs modules ---
 USE GetSurfAlbLUT_H5module
 USE GetLUT_LER_AI_H5module
 USE Get_TerrainPressure_H5module
 USE Get_SnowIce_module
 USE HDF5

 IMPLICIT NONE

 Include '../../Main/common_l1b_var.inc'

! Loop Counters
 INTEGER(KIND=4)           :: iPix, jPix, j
 INTEGER(KIND=4)           :: iStart = 1

 INTEGER(KIND=4)           :: STATUS, version
 REAL(KIND=4),DIMENSION(:,:),  INTENT(OUT):: UVAI
 REAL(KIND=4),DIMENSION(:,:,:),INTENT(IN) :: UVtoSWIR_Reflectances
 TYPE(ociuaaer_config_type),   INTENT(IN) :: cfg
 INTEGER(KIND=4),              INTENT(IN) :: Month
!******************************************************************************
 INTEGER(KIND=4) :: XDim, YDim

 REAL(KIND=4), DIMENSION(7) :: sza_table
 REAL(KIND=4), DIMENSION(11):: phi_table
 REAL(KIND=4), DIMENSION(14):: theta_table

 DATA theta_table/0.,12.,18.,26.,32.,36.,40.,46.,50.,54.,56.,60., 66., 72./
 DATA sza_table/0.,20.,40.,60.,66., 72.,80./
 DATA phi_table/0.,30.,60.,90.,120.,150.,160., 165.,170.,175.,180./
 !
 REAL(KIND=4), DIMENSION(14):: theta_table_ep
 REAL(KIND=4), DIMENSION(9) :: sza_table_ep
 REAL(KIND=4), DIMENSION(11):: phi_table_ep
!
 REAL(KIND=4), PARAMETER :: DTOR = (PI/180.)
 REAL(KIND=4) :: cossza

 DATA theta_table_ep /0.,12.,18.,26.,32.,36.,40.,46.,50.,54.,56.,60., 66., 72./
 DATA sza_table_ep/0.,20.,40.,60.,66., 72.,80., 84., 88./
 DATA phi_table_ep /0.,30.,60.,90.,120.,150.,160.,165.,170.,175.,180./


XDim = l1b_nXTrack
YDim = l1b_nLines

STATUS = allocate_UV1km_Data(XDim, YDim)
  IF (STATUS < 0) THEN
     PRINT *,"Error : Error allocating and initalizing UV1km data"
     CALL EXIT(1)
  ENDIF 

! Read HE4 AI Mie-Atmosphere look-up tables
! GetLUT_Miecloud_AI_H5module 
! nplev_mie, nsalb_mie, ncod_mie, nsza_mie, nvza_mie, nraa_mie
! rad354, rad388, rad_lin354_ai, rad_lin388_ai
  CALL ReadLUTAIparams(cfg%uv_ai_mielut)


! Read HE4 AI LER-Atmosphere look-up tables
! GetLUT_LER_AI_H5module 
! nsalb_ler, nplev_ler, nsza_ler, SURFALBSET_ler, PRESSURESET_ler
! rad354_ler, rad388_ler, rad_lin354_ai_ler, rad_lin388_ai_ler
  CALL ReadLUT_LER_AIparams(cfg%uv_ai_lerlut)


! Read OCI Surface AlbeDO HE4 look-up tables (Version 5 : ocean-corrected LER)
  CALL ReadSurfAlbLUTparams(cfg%uv_surfalbfile,SRFLER354, SRFLER388, GRDLAT, GRDLON)


! Read Terrain pressureH5 database [3600 longitude x 1800 latitude] bins
  CALL ReadTerrainPressparams(cfg%uv_surfprsfile,tpres, terrain_Longitude, terrain_Latitude)


! Read Snow_Ice database [360 longitude x 180 latitude x 12 month] bins
  CALL snowice_Reader(cfg%uv_snowicefile,snowice, swice_Longitude, swice_Latitude)


RelativeAzimuthAngle = SolarAzimuthAngle

! Start Pixels and Scan line loop
! Loop over Scan lines
DO jPix = 1, YDim 

! Loop over Pixels
    DO iPix = 1, XDim 

      CALL compute_relazimuth_ang(SolarAzimuthAngle(iPix,jPix), &
                                ViewingAzimuthAngle(iPix,jPix), &
                               RelativeAzimuthAngle(iPix,jPix))


      STATUS = Get_TerrainPressure(Latitude(iPix,jPix), &
                                  Longitude(iPix,jPix), &
                            TerrainPressure(iPix,jPix))
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting Terrain Pressure"
     	  CALL EXIT(1)
  	ENDIF 


      STATUS = Get_snowice_fraction(Month, Latitude(iPix,jPix), &
                                          Longitude(iPix,jPix), &
                                   SnowIce_fraction(iPix,jPix))
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting Snow/Ice Fraction"
     	  CALL EXIT(1)
  	ENDIF 
                                   
                                   
      STATUS = Get_UVSurfaceAlbedo(Month, Latitude(iPix,jPix), &
                                         Longitude(iPix,jPix), &
                         SurfaceAlbedo_Oceancorr(:,iPix,jPix))
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting Ocean corrected UV Surface albedo"
     	  CALL EXIT(1)
  	ENDIF 

      cossza = COS(DTOR * SolarZenithAngle(iPix,jPix))
      NormRadiances(1, iPix,jPix) = (UVtoSWIR_Reflectances(2,iPix,jPix)*cossza)/PI
      NormRadiances(2, iPix,jPix) = (UVtoSWIR_Reflectances(3,iPix,jPix)*cossza)/PI

      IF( (NormRadiances(1,iPix,jPix) .GT. 0) .AND. &
          (NormRadiances(2,iPix,jPix) .GT. 0) ) THEN 

          CALL calc_coeff(theta_table,sza_table,phi_table,ntheta,nsza,nphi, &
                          ViewingZenithAngle(iPix,jPix), &
                            SolarZenithAngle(iPix,jPix), &
                        RelativeAzimuthAngle(iPix,jPix))

          CALL calc_coeff_2dim(theta_table,sza_table,ntheta,nsza, &
                          ViewingZenithAngle(iPix,jPix), &
                            SolarZenithAngle(iPix,jPix) )

          ! --- Mie coefficients ----------
          CALL calc_coeff_Mie(theta_table_ep, sza_table_ep, phi_table_ep, &
                           nvza_mie, nsza_mie, nraa_mie, &
                          ViewingZenithAngle(iPix,jPix), &
                            SolarZenithAngle(iPix,jPix), &
                        RelativeAzimuthAngle(iPix,jPix))

          STATUS = UVAI_Miecloud_perpixel(Latitude(iPix,jPix), &
                                         Longitude(iPix,jPix), &
                                  SolarZenithAngle(iPix,jPix), &
                                ViewingZenithAngle(iPix,jPix), &
                              RelativeAzimuthAngle(iPix,jPix), &
                                   TerrainPressure(iPix,jPix), &
                                  SnowIce_fraction(iPix,jPix), &
                         SurfaceAlbedo_Oceancorr(:,iPix,jPix), &
                                   NormRadiances(:,iPix,jPix), & ! 354, 388
                                    UVAerosolIndex(iPix,jPix), &
                                 CloudOpticalDepth(iPix,jPix), &
                                     CloudFraction(iPix,jPix), &
                                           Residue(iPix,jPix), &
                                    Reflectivity(:,iPix,jPix) )


! 
!
      ENDIF ! Rad > 0 loop
      !
    ENDDO   ! iPix = 1, XDim-1
!
ENDDO    ! jPix = 1, YDim-1
! End Pixels and Scan line loop


!! Deallocate memory used for Terrain pressure variables
!  DEALLOCATE(tpres, terrain_Longitude, terrain_Latitude,  STAT=STATUS)

!! Deallocate snowice, swice_latitude, swice_longitude  
!  DEALLOCATE(snowice, swice_Longitude, swice_Latitude,  STAT=STATUS)
 
!! Deallocate memory used for OMI Surface Albedo HE4 LUT variables
  DEALLOCATE(SRFLER354, SRFLER388, GRDLAT, GRDLON,  STAT=STATUS)

!! Deallocate memory used for HE4 Mie-AI LUT variables
  DEALLOCATE(rad_lin354_ai, rad_lin388_ai, rad354, rad388, SURFALBSET, CODSET, STAT=STATUS)

!! Deallocate memory used for HE4 AI LER LUT variables
  DEALLOCATE(rad354_ler, rad388_ler, rad_lin354_ai_ler, rad_lin388_ai_ler, &
             SURFALBSET_ler,PRESSURESET_ler, STAT=STATUS)

END SUBROUTINE
!!
!!
!!
!!
!!
!==============================================================================
!==============================================================================
!!
FUNCTION UVAI_miecloud_perpixel(plat, plon, sun_za, sat_za, phi, pterrp, & 
                                swice, salb, normrad_obs, ainuv, &
                                codval, fcod10, ler_ai, reflectset) RESULT(STATUS)

! TITLE: Calulate the Aerosol Parameters
! NAME: L2OMUVAI_isocloud

  !USE GetLUT_module
  USE GetLUT_LER_AI_H5module
  USE LookupTableModule                 
  USE InterpolationModule

  IMPLICIT NONE  

! Function Inputs
  REAL(KIND=4),                INTENT(IN)  :: plat, plon, pterrp
  REAL(KIND=4),                INTENT(IN)  :: sun_za, sat_za, phi
  REAL(KIND=4),                INTENT(IN)  :: swice
  REAL(KIND=4),   DIMENSION(2),INTENT(IN)  :: salb
  REAL(KIND=4),   DIMENSION(2),INTENT(IN)  :: normrad_obs
  
! Function Outputs
  REAL(KIND=4),                INTENT(OUT) :: ainuv 
  REAL(KIND=4),                INTENT(OUT) :: codval
  REAL(KIND=4),                INTENT(OUT) :: fcod10
  REAL(KIND=4),                INTENT(OUT) :: ler_ai
  REAL(KIND=4), DIMENSION(2),  INTENT(OUT) :: reflectset
  
!  SnowIce Flag
  INTEGER(KIND=2)              :: SnowIceFlag 
  REAL(KIND=4)    :: swicefrac
  INTEGER(KIND=4) :: iplv, isalb, icod
!
 REAL(KIND=4)                                  :: rad_out, trm_out
 REAL(KIND=4),     DIMENSION(nsza,ntheta,nphi) :: tmprad
!
  REAL(KIND=4), DIMENSION(nsalb_mie,ncod_mie) :: radiance354_p10, radiance354_p08,radiance354_p06, &
                                                 radiance388_p10, radiance388_p08,radiance388_p06
  REAL(KIND=4), DIMENSION(ncod_mie) :: toasetsalb_388a, toasetsalb_388b, &
                                       toasetsalb_354a, toasetsalb_354b
!
  REAL(KIND=4) ::   toasetsalb_388a_ler, toasetsalb_388b_ler, &
                    toasetsalb_354a_ler, toasetsalb_354b_ler
!
  REAL(KIND=4), DIMENSION(ncod_mie)  :: intpradiance354, intpradiance388
  REAL(KIND=4), DIMENSION(ncod_mie)  :: calrad_wav1, tmpn_wav1, adjrad_wav1, fixrad_wav1
  REAL(KIND=4), DIMENSION(ncod_mie)  :: calrad_wav2, tmpn_wav2, adjrad_wav2, fixrad_wav2
  REAL(KIND=4)                       :: calrad354, calrad388
  REAL(KIND=4) :: ramcor_wav1, ramcor_wav2
!
  REAL(KIND=4), DIMENSION(nsalb_mie,ncod_mie)   :: rad_out_ai
  REAL(KIND=4) :: fractionalSalb, salb1, salb2
  REAL(KIND=4) :: fractionalSalb1, fractionalSalb2, salb11, salb12, salb21, salb22
  INTEGER(KIND=4)  :: jsalb1, jsalb2
  INTEGER(KIND=4)  :: jsalb11, jsalb12, jsalb21, jsalb22

  REAL(KIND=4) :: best354toa, best388toa
  REAL(KIND=4) :: tmpcod, cf, dif, cor_ref
!
  REAL (KIND=4), DIMENSION(2) :: out_oceanler !! Ocean LERs at 354 and 388 nm
  REAL (KIND=4), DIMENSION(2) :: newsalb

! -- LER-AI related variables ---------------------------------------------------------
  REAL(KIND=4), DIMENSION(nsalb_ler) :: radiance354_p10_ler, radiance354_p08_ler, &
                                        radiance354_p06_ler, radiance354_p02_ler, &
                                        radiance388_p10_ler, radiance388_p08_ler, &
                                        radiance388_p06_ler, radiance388_p02_ler

  REAL(KIND=4), DIMENSION(nsalb_ler) :: intpradiance354_sler, intpradiance388_sler
  REAL(KIND=4)                       :: intpradiance388_ler, intpradiance354_ler
  REAL(KIND=4), DIMENSION(nsalb_ler) :: rad_out_ai_ler
  REAL(KIND=4)                       :: ler354toa, ler388toa

  REAL(KIND=4), DIMENSION(nsalb_ler) :: calrad_wav1_ler, tmpn_wav1_ler, adjrad_wav1_ler, fixrad_wav1_ler
  REAL(KIND=4), DIMENSION(nsalb_ler) :: calrad_wav2_ler, tmpn_wav2_ler, adjrad_wav2_ler, fixrad_wav2_ler

! Pressure height  (make local to L2OMAERUVModule)
  REAL(KIND=4) :: pressure_table(2), pressure_table2(2), pressure_table3(2)
  REAL(KIND=4) :: logFacPressure

  DATA pressure_table  /1013.25, 800.0/
  DATA pressure_table2 /800.00,  600.0/
  DATA pressure_table3 /600.00,  250.0/

! For OCI Error messaging system
  INTEGER(KIND=4) :: STATUS, i

  STATUS = -1
  
! Determine the Terrain Pressure Fraction
   IF (pterrp > 800.0) THEN
      logFacPressure = (LOG(pressure_table(1)) - LOG(pterrp))           &
                      /(LOG(pressure_table(1)) - LOG(pressure_table(2)))
   ELSE IF (pterrp <= 800.0 .AND. pterrp > 600.0) THEN
      logFacPressure = (LOG(pressure_table2(1)) - LOG(pterrp))           &
                      /(LOG(pressure_table2(1)) - LOG(pressure_table2(2)))
   ELSE
      logFacPressure = (LOG(pressure_table3(1)) - LOG(pterrp))           &
                      /(LOG(pressure_table3(1)) - LOG(pressure_table3(2)))
   ENDIF
!
!   Lagrange interpolation method 
!
  DO iplv = 1, nplev_mie
      STATUS = InterpRadiance_ai(iplv,1,nsalb_mie,1,ncod_mie, rad_lin354_ai, rad_out_ai)
      SELECT CASE (iplv)
         CASE(1)
           radiance354_p10(:,:) = rad_out_ai
         CASE(2)
           radiance354_p08(:,:) = rad_out_ai
         CASE(3)
           radiance354_p06(:,:) = rad_out_ai
      END SELECT
      !
      STATUS = InterpRadiance_ai(iplv,1,nsalb_mie,1,ncod_mie, rad_lin388_ai, rad_out_ai)
      SELECT CASE (iplv)
         CASE(1)
           radiance388_p10(:,:) = rad_out_ai
         CASE(2)
           radiance388_p08(:,:) = rad_out_ai
         CASE(3)
           radiance388_p06(:,:) = rad_out_ai
      END SELECT
  ENDDO
!
! -- Do interpolation on angles for LER LUTs --------------
!
  DO iplv = 1, nplev_ler
     STATUS = InterpRadiance_ai_ler(iplv,1,nsalb_ler, rad_lin354_ai_ler, rad_out_ai_ler)
     SELECT CASE (iplv)
         CASE(1)
           radiance354_p10_ler(:) = rad_out_ai_ler
         CASE(2)
           radiance354_p08_ler(:) = rad_out_ai_ler
         CASE(3)
           radiance354_p06_ler(:) = rad_out_ai_ler
         CASE(4)
           radiance354_p02_ler(:) = rad_out_ai_ler
     END SELECT
!
     STATUS = InterpRadiance_ai_ler(iplv,1,nsalb_ler, rad_lin388_ai_ler, rad_out_ai_ler)
     SELECT CASE (iplv)
         CASE(1)
           radiance388_p10_ler(:) = rad_out_ai_ler
         CASE(2)
           radiance388_p08_ler(:) = rad_out_ai_ler
         CASE(3)
           radiance388_p06_ler(:) = rad_out_ai_ler
         CASE(4)
           radiance388_p02_ler(:) = rad_out_ai_ler
     END SELECT
!
  ENDDO


!    Interpolate radiances at 354 and 388 nm with the terrain pressure level
  DO isalb = 1,nsalb_ler
     IF (pterrp > 800.0) THEN
        STATUS = Interp1D(radiance354_p10_ler(isalb), &
                          radiance354_p08_ler(isalb), &
                          logFacPressure, intpradiance354_sler(isalb))
        STATUS = Interp1D(radiance388_p10_ler(isalb), &
                          radiance388_p08_ler(isalb), &
                          logFacPressure, intpradiance388_sler(isalb))
     ELSE IF (pterrp <= 800.0 .AND. pterrp > 600.0) THEN
        STATUS = Interp1D(radiance354_p08_ler(isalb), &
                          radiance354_p06_ler(isalb), &
                         logFacPressure, intpradiance354_sler(isalb))
        STATUS = Interp1D(radiance388_p08_ler(isalb), &
                          radiance388_p06_ler(isalb), &
                         logFacPressure, intpradiance388_sler(isalb))
     ELSE
        STATUS = Interp1D(radiance354_p06_ler(isalb), &
                          radiance354_p02_ler(isalb), &
                         logFacPressure, intpradiance354_sler(isalb))
        STATUS = Interp1D(radiance388_p06_ler(isalb), &
                          radiance388_p02_ler(isalb), &
                         logFacPressure, intpradiance388_sler(isalb))
     ENDIF
  ENDDO ! isalb-loop


! Compute Aerosol Indices ainuv and Lambert Equivalent Reflectivities

    STATUS = FindTableEntry(normrad_obs(2), intpradiance388_sler, nsalb_ler, &
                          salb1,salb2,jsalb1,jsalb2,fractionalSalb)
    STATUS = Interp1D(SURFALBSET_ler(jsalb1), &
                      SURFALBSET_ler(jsalb2), &
                      fractionalSalb, reflectset(2))
!
    STATUS = FindTableEntry(normrad_obs(1), intpradiance354_sler, nsalb_ler, &
                          salb1,salb2,jsalb1,jsalb2,fractionalSalb)
    STATUS = Interp1D(SURFALBSET_ler(jsalb1), &
                      SURFALBSET_ler(jsalb2), &
                      fractionalSalb, reflectset(1))

! --- Compute TOA radiances at 354 and 388 nm using the corrected 388 reflectivity ---------
    STATUS = FindTableEntry(reflectset(2), SURFALBSET_ler, nsalb_ler, &
                          salb1,salb2,jsalb1,jsalb2,fractionalSalb)
    STATUS = Interp1D(intpradiance388_sler(jsalb1), &
                      intpradiance388_sler(jsalb2), &
                      fractionalSalb, ler388toa)
    STATUS = Interp1D(intpradiance354_sler(jsalb1), &
                      intpradiance354_sler(jsalb2), &
                      fractionalSalb, ler354toa)
!
    ler_ai =  -100. *( alog10(normrad_obs(1)/normrad_obs(2)) - alog10(ler354toa/ler388toa) )

   newsalb(1) = salb(1)
   newsalb(2) = salb(2)
  STATUS = FindTableEntry(newsalb(1),SURFALBSET, nsalb_mie, &
                        salb11,salb12,jsalb11,jsalb12,fractionalSalb1)
  STATUS = FindTableEntry(newsalb(2),SURFALBSET, nsalb_mie, &
                        salb21,salb22,jsalb21,jsalb22,fractionalSalb2)

!
!    Interpolate radiances at 354 and 388 nm with the surface albedo and terrain pressure level
   DO icod = 1, ncod_mie
      if (salb(2) .lt. 0.20) then
         if (icod .eq. 1) then ! Use spectral dependent surface albedos for COD = 0.0 case.
             IF (pterrp > 800.0) THEN
                toasetsalb_354a(icod) = radiance354_p10(jsalb11,icod)*(1.0 - fractionalSalb1) &
			              + radiance354_p10(jsalb12,icod)* fractionalSalb1
                toasetsalb_354b(icod) = radiance354_p08(jsalb11,icod)*(1.0 - fractionalSalb1) &
				      + radiance354_p08(jsalb12,icod)* fractionalSalb1
                toasetsalb_388a(icod) = radiance388_p10(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance388_p10(jsalb22,icod)* fractionalSalb2
                toasetsalb_388b(icod) = radiance388_p08(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance388_p08(jsalb22,icod)* fractionalSalb2
             ELSE
                toasetsalb_354a(icod) = radiance354_p08(jsalb11,icod)*(1.0 - fractionalSalb1) &
				      + radiance354_p08(jsalb12,icod)* fractionalSalb1
                toasetsalb_354b(icod) = radiance354_p06(jsalb11,icod)*(1.0 - fractionalSalb1) &
				      + radiance354_p06(jsalb12,icod)* fractionalSalb1
                toasetsalb_388a(icod) = radiance388_p08(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance388_p08(jsalb22,icod)* fractionalSalb2
                toasetsalb_388b(icod) = radiance388_p06(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance388_p06(jsalb22,icod)* fractionalSalb2
              ENDIF                
         ELSE  ! Assume the same suface albedo at 388 nm
             IF (pterrp > 800.0) THEN
                toasetsalb_354a(icod) = radiance354_p10(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance354_p10(jsalb22,icod)* fractionalSalb2
                toasetsalb_354b(icod) = radiance354_p08(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance354_p08(jsalb22,icod)* fractionalSalb2
                toasetsalb_388a(icod) = radiance388_p10(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance388_p10(jsalb22,icod)* fractionalSalb2
                toasetsalb_388b(icod) = radiance388_p08(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance388_p08(jsalb22,icod)* fractionalSalb2
             ELSE
                toasetsalb_354a(icod) = radiance354_p08(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance354_p08(jsalb22,icod)* fractionalSalb2
                toasetsalb_354b(icod) = radiance354_p06(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance354_p06(jsalb22,icod)* fractionalSalb2
                toasetsalb_388a(icod) = radiance388_p08(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance388_p08(jsalb22,icod)* fractionalSalb2
                toasetsalb_388b(icod) = radiance388_p06(jsalb21,icod)*(1.0 - fractionalSalb2) &
				      + radiance388_p06(jsalb22,icod)* fractionalSalb2
             ENDIF
         ENDIF
      ELSE ! Assume the same suface albeDO at 388 nm
         IF (pterrp > 800.0) THEN
            toasetsalb_354a(icod) = radiance354_p10(jsalb21,icod)
            toasetsalb_354b(icod) = radiance354_p08(jsalb21,icod)
            toasetsalb_388a(icod) = radiance388_p10(jsalb21,icod)
            toasetsalb_388b(icod) = radiance388_p08(jsalb21,icod)
         ELSE
            toasetsalb_354a(icod) = radiance354_p08(jsalb21,icod)
            toasetsalb_354b(icod) = radiance354_p06(jsalb21,icod)
            toasetsalb_388a(icod) = radiance388_p08(jsalb21,icod)
            toasetsalb_388b(icod) = radiance388_p06(jsalb21,icod)
         ENDIF   
      ENDIF
!
      STATUS = Interp1D(toasetsalb_354a(icod), &
                        toasetsalb_354b(icod), &
                        logFacPressure, intpradiance354(icod))
      STATUS = Interp1D(toasetsalb_388a(icod), &
                        toasetsalb_388b(icod), &
                        logFacPressure, intpradiance388(icod))
   ENDDO
!

IF ( salb(2) .GT. -0.05 .and. salb(1) .GT. -0.05) THEN

! --- Compute the cloud fraction with COD = 10 (fcod10) at 388 nm --> fc = (Im - Is)/(Ic - Is) -----
    fcod10 = (normrad_obs(2) - intpradiance388(1)) / (intpradiance388(3) - intpradiance388(1))

! --- Compute Aerosol Indices ainuv and Lambert Equivalent Reflectivities
    STATUS = FindTableEntry(normrad_obs(2),intpradiance388, ncod_mie, &
                          salb1,salb2,jsalb1,jsalb2,fractionalSalb)
    STATUS = Interp1D(CODSET(jsalb1), &
                      CODSET(jsalb2), &
                      fractionalSalb, codval) 
!
! -- Use the LER model over 100% snow/ice regions (Greenland, Antarctica)
    IF (swice .ge. 100.0 .OR. pterrp < 600.0) THEN
       ainuv = ler_ai
       fcod10 = 0.0
       codval =  -9999.
    ELSE IF (fcod10 .lt. 1.0) THEN
       !=====================================================================
       !  TOA_radiance = Is * (1.0 - fc) + Icod10 * fc
       !=====================================================================
       calrad388 = intpradiance388(1)* (1.0 - fcod10) + intpradiance388(3)*fcod10
       calrad354 = intpradiance354(1)* (1.0 - fcod10) + intpradiance354(3)*fcod10
       ainuv =  -100. *( alog10(normrad_obs(1)/normrad_obs(2)) - alog10(calrad354/calrad388) )

! -- Use weighted AI values from the LER model and Mie-LER model with snow/ice fractions --
       IF (swice .GE. 1 .AND. swice .LT. 100) THEN
          swicefrac = swice/100.0
          ainuv = ler_ai*swicefrac + ainuv*(1.0 - swicefrac)
       ENDIF

    ELSE ! fcod10 above 1.0
       fcod10 = 1.0
       ! --- Compute TOA radiances at 354 and 388 nm using the Cloud Optical Depth at 388 nm ---------
       tmpcod = codval
       IF (tmpcod .gt. 100.0) THEN
          tmpcod = 100.0
          codval = tmpcod
       ENDIF
       !
       STATUS = FindTableEntry(tmpcod, CODSET, ncod_mie,salb1,salb2,jsalb1,jsalb2, &
                                                             fractionalSalb)
       STATUS = Interp1D(intpradiance388(jsalb1), &
                         intpradiance388(jsalb2), &
                         fractionalSalb, best388toa)
       STATUS = Interp1D(intpradiance354(jsalb1), &
                         intpradiance354(jsalb2), &
                         fractionalSalb, best354toa)
                         
! --- Compute the aerosol index ---------------------------------------------------
       ainuv =  -100. *( alog10(normrad_obs(1)/normrad_obs(2)) - alog10(best354toa/best388toa) )

! -- Use weighted AI values from the LER model and COD-LER model with snow/ice fractions --
       IF (swice .GE. 1 .AND. swice .LT. 100) THEN
          swicefrac = swice/100.0
          ainuv = ler_ai*swicefrac + ainuv*(1.0 - swicefrac)
       ENDIF
    ENDIF       
!
END IF  !IF SALB GT -0.05 
!

END FUNCTION UVAI_miecloud_perpixel
!!
!!
!!
!==============================================================================
!==============================================================================
!!
FUNCTION InterpRadiance_ai(iwave,nw01,nw02,ntau1,ntau2,inRad_ai, outRad_ai) RESULT(STATUS)

  USE GetLUT_Miecloud_AI_H5module
  USE LookupTableModule

  IMPLICIT NONE

  INTEGER(KIND=4),                INTENT(IN)    :: iwave
  INTEGER(KIND=4),                INTENT(IN)    :: nw01,ntau1
  INTEGER(KIND=4),                INTENT(IN)    :: nw02,ntau2
!
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(IN)    :: inRad_ai
  REAL(KIND=4), DIMENSION(:,:), INTENT(OUT)   :: outRad_ai

  INTEGER(KIND=4)    :: iw0, jw0, itau
  INTEGER(KIND=4)    :: ierr, polyInterp
  REAL(KIND=4)       :: dy
  INTEGER            :: i,k,m,n
  INTEGER(KIND=4)    :: isza,ithe,iphi

  INTEGER(KIND=4)        :: STATUS
  INTEGER(KIND=4) :: idxgral(64), idxaux, model

!  STATUS = AER_S_SUCCESS

 n = 0
 DO k=1,4
    ithe = indscn+k-1
    DO m=1,4
       iphi = indphi+m-1
       DO i=1,4
          isza = indsol+i-1
          n = n+1
!         3 -> the number of actual pressure levels (1013, 800 and 600 hPa) used on the calculations
          idxgral(n) = iwave+ (isza-1)*3*nsalb_mie*ncod_mie + &
                              (ithe-1)*3*nsalb_mie*ncod_mie*nsza_mie + &
                              (iphi-1)*3*nsalb_mie*ncod_mie*nsza_mie*nvza_mie 
       ENDDO
    ENDDO
 ENDDO

   DO iw0 = nw01,nw02
      jw0 = iw0 - nw01 + 1
      DO itau = ntau1,ntau2
         idxaux = (iw0-1)*3 + (itau-1)*3*nsalb_mie
         outRad_ai(jw0,itau) = SUM( inRad_ai(idxgral + idxaux)*cofs )
      ENDDO  ! DO itau = 1,ncod
   ENDDO  ! DO iw0  = 1,nsalb

 RETURN

 END FUNCTION InterpRadiance_ai
!!
!!
!!
!==============================================================================
!==============================================================================
!!
FUNCTION InterpRadiance_ai_ler(iwave,nw01,nw02,inRad_ai, outRad_ai) RESULT(STATUS)

  USE GetLUT_Miecloud_AI_H5module
  USE GetLUT_LER_AI_H5module
  USE LookupTableModule

  IMPLICIT NONE

  INTEGER(KIND=4),                INTENT(IN)    :: iwave
  INTEGER(KIND=4),                INTENT(IN)    :: nw01,nw02
!
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(IN)    :: inRad_ai
  REAL(KIND=4), DIMENSION(:), INTENT(OUT)   :: outRad_ai

  INTEGER(KIND=4)    :: iw0, jw0, itau
  INTEGER(KIND=4)    :: ierr, polyInterp
  INTEGER(KIND=4)    :: isza,ithe,iphi

  REAL(KIND=4)       :: dy
  INTEGER :: i,k,m,n

  INTEGER(KIND=4)        :: STATUS
  INTEGER(KIND=4) :: idxgral(64), idxaux, model

!  STATUS = AER_S_SUCCESS

 n = 0
 DO k=1,4
    ithe = indscn+k-1
    DO m=1,4
       iphi = indphi+m-1
       DO i=1,4
          isza = indsol+i-1
          n = n+1
!         4 -> the number of actual pressure levels (1013, 800, 600 and 250 hPa) used on the calculations
          idxgral(n) = iwave+ (isza-1)*4*nsalb_ler + &
                              (ithe-1)*4*nsalb_ler*nsza_mie + &
                              (iphi-1)*4*nsalb_ler*nsza_mie*nvza_mie
       ENDDO
    ENDDO
 ENDDO

   DO iw0 = nw01,nw02
      jw0 = iw0 - nw01 + 1
         idxaux = (iw0-1)*4 ! + (itau-1)*2*nsalb
         outRad_ai(jw0) = SUM( inRad_ai(idxgral + idxaux)*cofs )
   ENDDO  ! DO iw0  = 1,nsalb

 RETURN

 END FUNCTION InterpRadiance_ai_ler
!!
!!
!!
!==============================================================================
!==============================================================================
!!
SUBROUTINE calc_coeff_Mie(thenod,szanod,phinod,nthe,nsza_in,nphi_in,xtheta,xsza,xphi_in)

  USE LookupTableModule 

   IMPLICIT     NONE

   INTEGER(KIND=4) :: nthe,nsza_in,nphi_in
   REAL (KIND=4)   :: szanod(nsza_in),thenod(nthe),phinod(nphi_in)
   REAL (KIND=4)   :: xtheta,xsza,xphi_in

   REAL (KIND=4)   :: densol(4,nsza_in-3),denscn(4,nthe-3),denphi(4,nphi_in-3)
   REAL (KIND=4)   :: csza(4),ctheta(4),cphi(4)
   REAL (KIND=4)   :: xdenom, xnom
   INTEGER :: i,j,k,l,m, indmax

!     -- calculate values needed in table interpolation
!	  -- theta	  
      DO j = 1, nthe-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN
		xdenom = xdenom * (thenod(k) - thenod(i))
	    ENDIF
          END DO
          denscn(k - j + 1, j) = xdenom
        END DO
      END DO
!	  -- xsza
      DO j = 1, nsza_in-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN
		xdenom = xdenom * (szanod(k) - szanod(i))
	    ENDIF
          END DO
          densol(k - j + 1, j) = xdenom
        END DO
      END DO
!	  -- phi
      DO j = 1, nphi_in-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN 
		xdenom = xdenom * (phinod(k) - phinod(i))
	    ENDIF
          END DO
          denphi(k - j + 1, j) = xdenom
        END DO
      END DO
	  
!     -- set up index offset value for xsza, theta, and phi
!     -- theta --
	  CALL locate(thenod,nthe,xtheta,indscn)
      indscn=indscn-1
!     -- check range
      IF (indscn < 1) THEN
		indscn=1
	  ENDIF
      IF (indscn > (nthe-3)) THEN
		indscn=nthe-3
	  ENDIF

!	  -- sza --
	  CALL locate(szanod,nsza_in,xsza,indsol)
      indsol=indsol-1
!     -- check for range
      IF (indsol < 1) THEN 
		indsol=1
	  ENDIF
      IF (indsol > (nsza_in-3)) THEN
		indsol=nsza_in-3
	  ENDIF

!     -- phi --
	  CALL locate(phinod,nphi_in,xphi_in,indphi)
      indphi=indphi-1
!     -- check range
      IF (indphi < 1) THEN 
		indphi=1
	  ENDIF
      IF (indphi > (nphi_in-3)) THEN 
		indphi=nphi_in-3
	  ENDIF

!     -- computes Lagrange coefficients for sza
      indmax=indsol+3
      j=1
      DO k=indsol,indmax
         xnom=1.0
         DO i=indsol,indmax
            IF(i /= k) THEN 
	       xnom=(xsza-szanod(i))*xnom
	    ENDIF
         END DO
         csza(j)=xnom/densol(j,indsol)
         j=j+1
      END DO
!     -- computes Lagrange coefficients for theta
      j=1
      indmax=indscn+3
      DO k=indscn,indmax
         xnom=1.0
         DO i=indscn,indmax
            IF (i /= k) THEN
		xnom=(xtheta-thenod(i))*xnom
	    ENDIF
         END DO
         ctheta(j)=xnom/denscn(j,indscn)
         j=j+1
      END DO
!     -- computes Lagrange coefficients for phi
      j=1
      indmax=indphi+3
      DO k=indphi,indmax
         xnom=1.0
         DO i=indphi,indmax
            IF (i /= k) THEN
	     xnom=(xphi_in-phinod(i))*xnom
            ENDIF
         END DO
         cphi(j)=xnom/denphi(j,indphi)
         j=j+1
      END DO

!     -- store products of coefficients
      l=1
      DO i=1,4
        DO k=1,4
	 DO m=1,4
           cofs(l)=ctheta(i)*csza(m)*cphi(k)
           l=l+1
	 END DO
	END DO
      END DO 
      
END SUBROUTINE calc_coeff_Mie
!!
!!
!!
!==============================================================================
!==============================================================================
!!
 SUBROUTINE calc_coeff_2dim_Mie(thenod,szanod,nthe,nsza_in,xtheta,xsza)

  USE LookupTableModule 

	IMPLICIT     NONE

	INTEGER(KIND=4) :: nthe,nsza_in
	REAL (KIND=4) :: szanod(nsza_in),thenod(nthe)
	REAL (KIND=4) :: xtheta,xsza

	REAL (KIND=4) :: densol(4,nsza_in-3),denscn(4,nthe-3)
	REAL (KIND=4) :: csza(4),ctheta(4)
	REAL (KIND=4) :: xdenom, xnom
	INTEGER :: i,j,k,l,m, indmax

!     -- calculate values needed in table interpolation
!	  -- theta	  
      DO j = 1, nthe-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN
		xdenom = xdenom * (thenod(k) - thenod(i))
	    ENDIF
          END DO
          denscn(k - j + 1, j) = xdenom
        END DO
      END DO
!	  -- xsza
      DO j = 1, nsza_in-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN
		xdenom = xdenom * (szanod(k) - szanod(i))
	    ENDIF
          END DO
          densol(k - j + 1, j) = xdenom
        END DO
      END DO
	  
!     -- set up index offset value for xsza, theta, and phi
!     -- theta --
	  CALL locate(thenod,nthe,xtheta,indscn_2d)
      indscn_2d=indscn_2d-1
!     -- check range
      IF (indscn_2d < 1) THEN
		indscn_2d=1
	  ENDIF
      IF (indscn_2d > (nthe-3)) THEN
		indscn_2d=nthe-3
	  ENDIF

!	  -- sza --
	  CALL locate(szanod,nsza_in,xsza,indsol_2d)
      indsol_2d=indsol_2d-1
!     -- check for range
      IF (indsol_2d < 1) THEN 
		indsol_2d=1
	  ENDIF
      IF (indsol_2d > (nsza_in-3)) THEN
		indsol_2d=nsza_in-3
	  ENDIF

!     -- computes Lagrange coefficients for sza
      indmax=indsol_2d+3
      j=1
      DO k=indsol_2d,indmax
         xnom=1.0
         DO i=indsol_2d,indmax
            IF(i /= k) THEN 
	       xnom=(xsza-szanod(i))*xnom
	    ENDIF
         END DO
         csza(j)=xnom/densol(j,indsol_2d)
         j=j+1
      END DO
!     -- computes Lagrange coefficients for theta
      j=1
      indmax=indscn_2d+3
      DO k=indscn_2d,indmax
         xnom=1.0
         DO i=indscn_2d,indmax
            IF (i /= k) THEN
		xnom=(xtheta-thenod(i))*xnom
	    ENDIF
         END DO
         ctheta(j)=xnom/denscn(j,indscn_2d)
         j=j+1
      END DO

!     -- store products of coefficients
      l=1
      DO i=1,4
        DO m=1,4
           cofs_2d(l)=ctheta(i)*csza(m)	!*cphi(k)
           l=l+1
	    END DO
      END DO 

END SUBROUTINE calc_coeff_2dim_Mie
!!
!!
!!
!==============================================================================
!==============================================================================
!!
SUBROUTINE calc_coeff(thenod,szanod,phinod,nthe,nsza_in,nphi_in,xtheta,xsza,xphi_in)

  USE LookupTableModule 

    IMPLICIT     NONE

    INTEGER(KIND=4):: nthe,nsza_in,nphi_in
    REAL(KIND=4)   :: szanod(nsza_in),thenod(nthe),phinod(nphi_in)
    REAL(KIND=4)   :: xtheta,xsza,xphi_in
    REAL(KIND=4)   :: densol(4,nsza_in-3),denscn(4,nthe-3),denphi(4,nphi_in-3)
    REAL(KIND=4)   :: csza(4),ctheta(4),cphi(4)
    REAL(KIND=4)   :: xdenom, xnom
    INTEGER        :: i,j,k,l,m, indmax

!     -- calculate values needed in table interpolation
!	  -- theta	  
      DO j = 1, nthe-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN
		xdenom = xdenom * (thenod(k) - thenod(i))
	    ENDIF
          END DO
          denscn(k - j + 1, j) = xdenom
        END DO
      END DO
!	  -- xsza
      DO j = 1, nsza_in-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN
		xdenom = xdenom * (szanod(k) - szanod(i))
	    ENDIF
          END DO
          densol(k - j + 1, j) = xdenom
        END DO
      END DO
!	  -- phi
      DO j = 1, nphi_in-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN 
		xdenom = xdenom * (phinod(k) - phinod(i))
	    ENDIF
          END DO
          denphi(k - j + 1, j) = xdenom
        END DO
      END DO
	  
!     -- set up index offset value for xsza, theta, and phi
!     -- theta --
      CALL locate(thenod,nthe,xtheta,indscn)
      indscn=indscn-1
!     -- check range
      IF (indscn < 1) THEN
	  indscn=1
      ENDIF
      IF (indscn > (nthe-3)) THEN
	  indscn=nthe-3
      ENDIF

!	  -- sza --
      CALL locate(szanod,nsza_in,xsza,indsol)
      indsol=indsol-1
!     -- check for range
      IF (indsol < 1) THEN 
		indsol=1
      ENDIF
      IF (indsol > (nsza_in-3)) THEN
		indsol=nsza_in-3
      ENDIF

!     -- phi --
	  CALL locate(phinod,nphi_in,xphi_in,indphi)
      indphi=indphi-1
!     -- check range
      IF (indphi < 1) THEN 
		indphi=1
      ENDIF
      IF (indphi > (nphi_in-3)) THEN 
		indphi=nphi_in-3
      ENDIF

!     -- computes Lagrange coefficients for sza
      indmax=indsol+3
      j=1
      DO k=indsol,indmax
         xnom=1.0
         DO i=indsol,indmax
            IF(i /= k) THEN 
	       xnom=(xsza-szanod(i))*xnom
            ENDIF
         END DO
         csza(j)=xnom/densol(j,indsol)
         j=j+1
      END DO
!     -- computes Lagrange coefficients for theta
      j=1
      indmax=indscn+3
      DO k=indscn,indmax
         xnom=1.0
         DO i=indscn,indmax
            IF (i /= k) THEN
		xnom=(xtheta-thenod(i))*xnom
	    ENDIF
         END DO
         ctheta(j)=xnom/denscn(j,indscn)
         j=j+1
      END DO
!     -- computes Lagrange coefficients for phi
      j=1
      indmax=indphi+3
      DO k=indphi,indmax
         xnom=1.0
         DO i=indphi,indmax
            IF (i /= k) THEN
		xnom=(xphi_in-phinod(i))*xnom
	    ENDIF
         END DO
         cphi(j)=xnom/denphi(j,indphi)
         j=j+1
      END DO

!     -- store products of coefficients
      l=1
      DO i=1,4
        DO k=1,4
	    DO m=1,4
             cofs(l)=ctheta(i)*csza(m)*cphi(k)
             l=l+1
	   END DO
	END DO
      END DO 
	  
END SUBROUTINE calc_coeff
!!
!!
!!
!==============================================================================
!==============================================================================
!!
 SUBROUTINE calc_coeff_2dim(thenod,szanod,nthe,nsza_in,xtheta,xsza)

 USE LookupTableModule 

    IMPLICIT     NONE

    INTEGER(KIND=4) :: nthe,nsza_in
    REAL(KIND=4) :: szanod(nsza_in),thenod(nthe)
    REAL(KIND=4) :: xtheta,xsza
    REAL(KIND=4) :: densol(4,nsza_in-3),denscn(4,nthe-3)
    REAL(KIND=4) :: csza(4),ctheta(4)
    REAL(KIND=4) :: xdenom, xnom
    INTEGER      :: i,j,k,l,m, indmax

!     -- calculate values needed in table interpolation
!	  -- theta	  
      DO j = 1, nthe-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN
		xdenom = xdenom * (thenod(k) - thenod(i))
	    ENDIF
          END DO
          denscn(k - j + 1, j) = xdenom
        END DO
      END DO
!	  -- xsza
      DO j = 1, nsza_in-3
        DO k = j, j + 3
          xdenom = 1.0
          DO i = j, j + 3
            IF (i /= k) THEN
		xdenom = xdenom * (szanod(k) - szanod(i))
	    ENDIF
          END DO
          densol(k - j + 1, j) = xdenom
        END DO
      END DO
	  
!     -- set up index offset value for xsza, theta, and phi
!     -- theta --
	  CALL locate(thenod,nthe,xtheta,indscn_2d)
      indscn_2d=indscn_2d-1
!     -- check range
      IF (indscn_2d < 1) THEN
		indscn_2d=1
	  ENDIF
      IF (indscn_2d > (nthe-3)) THEN
		indscn_2d=nthe-3
	  ENDIF

!	  -- sza --
	  CALL locate(szanod,nsza_in,xsza,indsol_2d)
      indsol_2d=indsol_2d-1
!     -- check for range
      IF (indsol_2d < 1) THEN 
		indsol_2d=1
	  ENDIF
      IF (indsol_2d > (nsza_in-3)) THEN
		indsol_2d=nsza_in-3
	  ENDIF

!     -- computes Lagrange coefficients for sza
      indmax=indsol_2d+3
      j=1
      DO k=indsol_2d,indmax
         xnom=1.0
         DO i=indsol_2d,indmax
            IF (i /= k) THEN 
	       xnom=(xsza-szanod(i))*xnom
	    ENDIF
         END DO
         csza(j)=xnom/densol(j,indsol_2d)
         j=j+1
      END DO
!     -- computes Lagrange coefficients for theta
      j=1
      indmax=indscn_2d+3
      DO k=indscn_2d,indmax
         xnom=1.0
         DO i=indscn_2d,indmax
            IF (i /= k) THEN
		xnom=(xtheta-thenod(i))*xnom
	    ENDIF
         END DO
         ctheta(j)=xnom/denscn(j,indscn_2d)
         j=j+1
      END DO

!     -- store products of coefficients
      l=1
      DO i=1,4
        DO m=1,4
           cofs_2d(l)=ctheta(i)*csza(m)	!*cphi(k)
           l=l+1
	    END DO
      END DO 
      
END SUBROUTINE calc_coeff_2dim
!!
!!
!!
!==============================================================================
!==============================================================================
!!
SUBROUTINE locate(xx,n,x,j)

 IMPLICIT     NONE

 INTEGER(KIND=4) :: jl,jm,ju,j,n
 REAL (KIND=4) :: x,xx(n)

! Initialize lower and upper limits. 		
   jl=0 
   ju=n+1 

! If we are not yet done, compute a midpoint, and replace either the lower limit
! or the upper limit, as appropriate.
   DO WHILE (ju-jl.gt.1)
     jm=(ju+jl)/2
     IF ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) THEN
        jl=jm 
     ELSE
        ju=jm 
     ENDIF
   ENDDO   

! Set the output and return.
     IF (x == xx(1)) THEN
        j=1
     ELSE IF (x == xx(n)) THEN
        j=n-1
     ELSE
        j=jl
     ENDIF

END SUBROUTINE locate
!!
!!
!!
!==============================================================================
!==============================================================================
!!
SUBROUTINE compute_relazimuth_ang(SolAzAng,ViewAzAng,RelAzAng) !RESULT(STATUS)

IMPLICIT NONE

  REAL(KIND=4), INTENT(IN)  :: SolAzAng, ViewAzAng
  REAL(KIND=4), INTENT(OUT) :: RelAzAng
  REAL(KIND=4)              :: deg2radconv
  INTEGER(KIND=4)           :: STATUS

 STATUS = 1

! Calculate the Relative Azimuth Angle
  deg2radconv = 4.0*ATAN(1.0)/180.0

  IF((SolAzAng .GE. -180.0) .AND.  (SolAzAng .LE.  180.0) .AND. &
    (ViewAzAng .GE. -180.0) .AND. (ViewAzAng .LE.  180.0)) THEN

     RelAzAng = ABS(SolAzAng - ViewAzAng)
     IF(RelAzAng .GT. 180.0) RelAzAng = 360.0 - RelAzAng
     RelAzAng = 180.0 - RelAzAng
  ELSE
     !PRINT *, 'Azimuth Angles are out-of-bounds', SolAzAng, ViewAzAng
     RelAzAng = -9999.
  ENDIF

END SUBROUTINE compute_relazimuth_ang

!==============================================================================
!==============================================================================

FUNCTION OCEAN_FRESNEL_CORRECTION(gpQF,sun_za,sat_za,phi,salb) RESULT(status)

IMPLICIT NONE

  INTEGER(KIND=2), INTENT(IN)  :: gpQF
  REAL(KIND=4),    INTENT(IN)  :: sun_za, sat_za, phi
  REAL(KIND=4),    INTENT(INOUT),  DIMENSION(2)   :: salb
!
  REAL, DIMENSION(2)   :: newsalb, out_oceanler
  INTEGER              :: status

  status = -1
! --- Add the ocean Fresnel LER to each surface albedo at 354 and 388 nm, respectively ------------
  if (gpQF .eq. 17) then
      out_oceanler = AI_oceanler(sun_za,sat_za, phi)
      if (salb(1) > -0.05 .and. out_oceanler(1) > -0.05) newsalb(1) = salb(1) + out_oceanler(1) ! 354 nm
      if (salb(2) > -0.05 .and. out_oceanler(2) > -0.05) newsalb(2) = salb(2) + out_oceanler(2) ! 388 nm
  else
     newsalb(1) = salb(1)
     newsalb(2) = salb(2)
  endif

salb(1:2) = newsalb(1:2)
status = 1

RETURN

END FUNCTION OCEAN_FRESNEL_CORRECTION


!==============================================================================
!==============================================================================
!
!=======  Find ocean LER with given geometry =============================
FUNCTION AI_oceanler(sza_in, vza_in, raa_in) RESULT( out_ler )

USE Get_OceanLUT_H5module

IMPLICIT NONE
!
    REAL( KIND=4 ),                INTENT(IN)  :: sza_in, vza_in, raa_in
    REAL( KIND=4 )                             :: szain, vzain, raain
    REAL( KIND=4 ), DIMENSION(2)						   :: out_ler
    REAL( KIND=4 )  :: frac_sza, frac_vza, frac_raa
    REAL( KIND=4 )  :: slope, bias
    REAL( KIND=4 ), DIMENSION(2) :: interp_oceanler
    INTEGER( KIND=4 ) :: it, i,j,k, wave_in
    REAL( KIND=4 ), DIMENSION(16) :: sza_list=(/1.5, 6.0, 12.0,18.0,&
                                     24.0, 30.0, 36.0, 42.0,48.0, 54.0, &
                                     60.0, 66.0, 72.0, 78.0, 84.0, 88.5/)
    REAL( KIND=4 ), DIMENSION(16) :: vza_list=(/1.5, 6.0, 12.0,18.0,&
                                     24.0, 30.0, 36.0, 42.0,48.0, 54.0, &
                                     60.0, 66.0, 72.0, 78.0, 84.0, 88.5/)
    REAL( KIND=4 ), DIMENSION(16) :: raa_list=(/0.0, 12.0, 24.0, &
                               36.0, 48.0, 60.0, 72.0, 84.0, 96.0, 108.0,&
                               120.0, 132.0, 144.0, 156.0, 168.0, 180.0/)
      !-- Input angles -----------
      szain = sza_in
      vzain = vza_in
      raain = raa_in
      if (sza_in <= 1.50) szain = 1.51
      if (sza_in >= 88.50) szain = 88.49
      if (vza_in <= 1.50) vzain = 1.51
      if (vza_in >= 88.50) vzain = 88.49
      if (raa_in <= 0.0) raain = 0.01
      if (raa_in >= 180.0) raain = 179.9
!
      DO it = 1, 2
        wave_in = it
        !-- Determine interpolation indices ----
        i = 1
        DO
          if (sza_list(i) > szain) EXIT
          i = i +1
        END DO
        frac_sza = (sza_list(i) - szain)/( sza_list(i) - sza_list(i-1) )
!
        j = 1
        DO
          if (vza_list(j) > vzain) EXIT
          j = j +1
        END DO
        frac_vza = (vza_list(j) - vzain)/( vza_list(j) - vza_list(j-1) )
!
        k = 1
        DO
          if (raa_list(k) > raain) EXIT
          k = k +1
        END DO
        frac_raa = (raa_list(k) - raain)/( raa_list(k) - raa_list(k-1) )
!
        interp_oceanler(it) =  (1.0 - frac_sza)*(1.0 - frac_vza)* &
                               (1.0 - frac_raa)*oceanler(wave_in, i,j,k) &
        + (frac_sza)*(1.0 - frac_vza)*(1.0 - frac_raa)*oceanler(wave_in, i-1,j,k) &
        + (frac_sza)*(frac_vza)*(1.0 - frac_raa)*oceanler(wave_in, i-1,j-1,k) &
        + (frac_sza)*(frac_vza)*(frac_raa)*oceanler(wave_in, i-1,j-1,k-1) &
        + (frac_sza)*(1.0 - frac_vza)*(frac_raa)*oceanler(wave_in, i-1,j,k-1) &
        + (1.0 - frac_sza)*(frac_vza)*(1.0 - frac_raa)*oceanler(wave_in, i,j-1,k) &
        + (1.0 - frac_sza)*(frac_vza)*(frac_raa)*oceanler(wave_in, i,j-1,k-1) &
        + (1.0 - frac_sza)*(1.0 - frac_vza)*(frac_raa)*oceanler(wave_in, i,j,k-1)

      END DO
      
     !--- Compute LER at 354.0 nm (interpolated) and 388.0 nm (extrapolated) ----
      slope = (interp_oceanler(1) -interp_oceanler(2))/(380.0 - 340.0)
      bias = interp_oceanler(2) - slope*340.0
!
      out_ler(1) = slope*354.0 + bias  ! 354 nm LER
      out_ler(2) = slope*388.0 + bias  ! 388 nm LER
      IF (out_ler(1) < 0.0 .OR. out_ler(1) > 0.5) out_ler(1) = -999.0
      IF (out_ler(2) < 0.0 .OR. out_ler(2) > 0.5) out_ler(2) = -999.0
!
END FUNCTION AI_oceanler

!!
!==============================================================================
!==============================================================================
!!
!!
END Module OCI_UVAI_MieCloudModule
