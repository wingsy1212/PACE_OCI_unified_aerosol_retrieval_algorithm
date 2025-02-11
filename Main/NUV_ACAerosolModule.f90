MODULE NUV_ACAerosolModule       

  USE GetLUT_H5module_nc4
  USE LookupTableModule_nc4 
  USE InterpolationModule
  USE NUV_AerosolModule
  USE Get_omacaLUT7dim_H5module_nc4
  USE Get_ssaclimLUT2_H5module_nc4

 IMPLICIT NONE

 PUBLIC  :: OCI_NUV_ACAer_Process

 CONTAINS


FUNCTION OCI_NUV_ACAer_Process(yr, mon, dy, nWavel,  &
                    wavelen, coval, plat, plon, sun_za, sat_za, phi, pterrp, gpQF, &
                    salb, rad_obs, ainuv, reflect,  & 
                    AlgQF, aodtau_aac, codtau_aac_wl3, appcodtau_aac_wl3, &
                    inputssa_aca, uncaodtau_aac, unccodtau_aac, &
		    aodtau_aac_zhgt, codtau_aac_zhgt) RESULT(status)

  IMPLICIT NONE  

!==============================================================================

  INTEGER(KIND=4),             INTENT(IN)  :: yr, mon, dy
  INTEGER(KIND=4),             INTENT(IN)  :: nWavel 
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: wavelen(nWavel)
  REAL(KIND=4),                INTENT(IN)  :: coval
  REAL(KIND=4),                INTENT(IN)  :: plat, plon
  REAL(KIND=4),                INTENT(IN)  :: sun_za, sat_za, phi
  REAL(KIND=4),                INTENT(IN)  :: pterrp
  INTEGER(KIND=2),             INTENT(IN)  :: gpQF
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: salb
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: rad_obs(nWavel)
  REAL(KIND=4),   	       INTENT(IN)  :: ainuv
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: reflect(nWavel)

  REAL(KIND=4),                INTENT(OUT) :: AlgQF
  REAL(KIND=4),                INTENT(OUT) :: codtau_aac_wl3
  REAL(KIND=4),                INTENT(OUT) :: appcodtau_aac_wl3
  REAL(KIND=4), DIMENSION(3),  INTENT(OUT) :: aodtau_aac     !3-wavelengths (354, 388, 500)
  REAL(KIND=4), DIMENSION(3),  INTENT(OUT) :: inputssa_aca   !3-wavelengths (354, 388, 500)
  REAL(KIND=4), DIMENSION(2),  INTENT(OUT) :: uncaodtau_aac
  REAL(KIND=4), DIMENSION(2),  INTENT(OUT) :: unccodtau_aac
  REAL(KIND=4), DIMENSION(3,5),INTENT(OUT) :: aodtau_aac_zhgt
  REAL(KIND=4), DIMENSION(5),  INTENT(OUT) :: codtau_aac_zhgt

!==============================================================================
!  aerosol layer height
  REAL(KIND=4)                 :: zaer, zaer2
  INTEGER(KIND=2)              :: atype_aerac
  REAL(KIND=4)                 :: inssacval, ssacval354, ssacval388, ssacval500, inssacval_tmp

!  aerosol indices
  REAL(KIND=4)    	       :: aivis, ainuvLER
  REAL(KIND=4)                 :: ainuvTmp
  REAL(KIND=4), DIMENSION(4)   :: coef_nuvai
  DATA coef_nuvai/14.993834,-48.632181,-22.634173,-7.7537689/

! non-absorbing AOD at 354 and 388 nm for the Angstrom Exponent calculation
  REAL(KIND=4)                 :: nabs_aod354, nabs_aod388

!  SnowIce Flag
  INTEGER(KIND=2)              :: SnowIceFlag 

!  un-needed reflectivity
  INTEGER(KIND=4)              :: surf
  INTEGER(KIND=4)              :: indexmd(7)
  REAL(KIND=4)                 :: reflect_tmp !, tau_ret(nWavel,nzae),taufinal(nWavel)
  REAL(KIND=4)                 :: tau_nabs(nWavel,nzae)
  REAL(KIND=4)                 :: caluv, reflect471, tmpfrac
  REAL(KIND=4)                 :: ssa500_bottom_thshold
  REAL(KIND=4)        	       :: deg2radconv, sun_ga, sca
  REAL(KIND=4)                 :: SZATmp, VZATmp, PHITmp

!==============================================================================
! Aerosol Height indices, model indices
!==============================================================================
  INTEGER(KIND=4) :: izae, iw0sel

!==============================================================================
! Wavelength indexes (make local to L2OMACAModule)
!==============================================================================
  INTEGER(KIND=4) :: iwl(3),iwaveT,iwavel
  REAL(KIND=4) :: fracWave

!==============================================================================
!  Bounding values in table for input parameters
!==============================================================================
  REAL(KIND=4)    :: theta1,theta2,sza1,sza2,phi1,phi2,w1,w2
  REAL(KIND=4)    :: fractionalTheta,fractionalSZA,fractionalPHI

!==============================================================================
!  Indices for bounding values in table for input parameters
!==============================================================================
  INTEGER(KIND=4) :: jtheta, jtheta2, jsza, jsza2, jphi, jphi2, jw1

!==============================================================================
! Radiance and Transmittance from lookup tables
!==============================================================================
  REAL(KIND=4), DIMENSION(nzae,nw0sel,ntau) :: rad_out,trm_out

!==============================================================================
! Radiance, Transmittance, and TOA from Aerosol Above Clouds lookup tables
!==============================================================================

  REAL(KIND=4), DIMENSION(naod_aac,ncod_aac) :: toa388_aac, toa354_aac, ratrad
  REAL(KIND=4), DIMENSION(ncod_aac,naod_aac) :: ratrad2, toa388_aac2, fint_uvaimie_aac2
  REAL(KIND=4)                               :: ratradval
  REAL(KIND=4), DIMENSION(ncod_aac)          :: toa388_aac2_array
!
  REAL(KIND=4)                   :: aodtau_aac_wl3
  REAL(KIND=4)                   :: codtau_aac_wl3_tmp, codtau_aac_tmp
  REAL(KIND=4)                   :: aodtaumax, codtaumax

! Uncertainty array
  REAL(KIND=4), DIMENSION(3)     :: aodtau_aac_tmp
!
  REAL(KIND=4)                   :: tmplon, inzhgt
  INTEGER(KIND=4)                :: ilon, ilat
  INTEGER(KIND=4)                :: status, ierr

  INTEGER(KIND=4)                :: nzhgt = 5
  INTEGER(KIND=4)                :: izhgt
  INTEGER(KIND=2)                :: atypeTmp
  
  REAL(KIND=4), DIMENSION(5)     :: zhgt_nodes
  DATA zhgt_nodes/3.0, 6.0, 9.0, 12.0, 15.0/
  
!==============================================================================
! Radiance rations for both model and observed Radiances
!==============================================================================
  REAL(KIND=4) :: ratio_mod, ratio_obs

!==============================================================================
! Surface flux  (make local to L2OMACAModule)
!==============================================================================
  REAL(KIND=4) :: sflux_lut
  DATA sflux_lut/3.1415926/

!==============================================================================
! Pressure height  (make local to L2OMACAModule)
!==============================================================================
  REAL(KIND=4) :: pressure_table(2)
  REAL(KIND=4) :: logFacPressure
  DATA pressure_table/1013.25,800./

!==============================================================================
! Unpacking gpQF fields
!==============================================================================
   INTEGER(KIND=4), DIMENSION(6) :: gp_flags
   INTEGER(KIND=4)               :: ocean, groundqf

!==============================================================================
! Local counter
!==============================================================================
  INTEGER(KIND=4) :: i, iSSA

!==============================================================================
!  Function calls 
!==============================================================================
  status = 1
  
!==============================================================================
! Determine the Terrain Pressure Fraction
!==============================================================================
  logFacPressure = (LOG(pressure_table(1)) - LOG(pterrp))           &
                  /(LOG(pressure_table(1)) - LOG(pressure_table(2)))

!==============================================================================
! Find the two values that bracket sat_za in theta_table and their indices
!==============================================================================
  status=FindTableEntry(sat_za,theta_table,ntheta,theta1,theta2,jtheta,jtheta2,&
                                                               fractionalTheta)

!==============================================================================
! Find the two values that bracket sun_za in sza_table and their indices 
!============================================================================== 
  status=FindTableEntry(sun_za,sza_table,nsza,sza1,sza2,jsza,jsza2, &
                                                                  fractionalSZA)

!==============================================================================
! Find the two values that bracket phi in phi_table and their indices 
!============================================================================== 
  status=FindTableEntry(phi,phi_table,nphi,phi1,phi2,jphi,jphi2,fractionalPHI) 

!==============================================================================
! Get surface type
!==============================================================================
  IF (gpQF .EQ. 17) THEN
     ocean = 0
  ELSE 
     ocean = 1
  ENDIF


!=============================================================================
! Proceed with retrieval procedure if sun glint angle is above 40 degree
!=============================================================================

!==============================================================================
!Proceed with retrieval procedure if retrieved reflectivity meets the criteria 
!indicated by the cloud contamination thresholds
!==============================================================================

inzhgt = -9999.
inssacval = -9999.
ssacval388 = -9999.
aivis = -9999.

!=============================
!Calculate SunGlintAngle here.
!=============================

deg2radconv = 4.0*ATAN(1.0)/180.0
SZATmp = sun_za*deg2radconv
VZATmp = sat_za*deg2radconv
PHITmp = phi*deg2radconv

!SunGlintAngle
sun_ga = ACOS(COS(SZATmp)*COS(VZATmp) + SIN(SZATmp)*SIN(VZATmp)*COS(PHITmp))/deg2radconv


!===============
!ScatteringAngle
!===============

sca = ACOS(-COS(SZATmp)*COS(VZATmp) + SIN(SZATmp)*SIN(VZATmp)*COS(PHITmp))/deg2radconv


!=================================
aodtau_aac_wl3 = -9999.
codtau_aac_wl3 = -9999.
appcodtau_aac_wl3 = -9999.
!=================================


!500-nm AOD/COD Max LUT values...
  aodtaumax = 6.0
  codtaumax = 70.0


 aivis =  1.0
 status = CalcAerosolModel(gpQF,coval,salb(2),ainuv,aivis, atype_aerac, indexmd, &
                           mon, plat, plon, reflect(2))

! ---Reflectivity, Sun glint angle, and surface albedo thresholds---!
IF( ( (gpQF.EQ.17 .AND. reflect(2) .GT. 0.20 .AND. reflect(2) .LE. 0.30 .AND. sun_ga .GT. 20).OR. &
      (gpQF.EQ.17 .AND. reflect(2) .GT. 0.30 .AND. sun_ga .GE. 0 .AND. phi .GT.20.).OR. &
      (gpQF.NE.17 .AND. reflect(2) .GT. 0.20)) .AND. &
       salb(2) .GT. 0 .AND. salb(2) .LE. 0.20 .AND. &
       pterrp .GT. 800.0 .AND. ((atype_aerac .LE. 2).OR.(atype_aerac .EQ. 4)).AND.sun_za.LE.80) THEN

! Surface albedo (388 nm) threshold of 0.20 is expected to avoid above-cloud 
! aerosol inversion over snow and sea-ice surface.

  ! ---Calculation of UV-AI for a given range of reflectivity (0.20-0.25)---!
  IF(reflect(2) .GE. 0.20 .AND. reflect(2) .LE. 0.25) THEN
      ainuvTmp = coef_nuvai(1) + &
                 coef_nuvai(2)*reflect(2) + &
		 coef_nuvai(3)*reflect(2)**2.0 + &
		 coef_nuvai(4)*reflect(2)**3.0
  ENDIF
  
  !  -- Aerosol Above Clouds Identification---!
   IF( (reflect(2) .GE. 0.20 .AND. reflect(2) .LE. 0.25 .AND. ainuv .GE. ainuvTmp).OR. &
       (reflect(2) .GT. 0.25 .AND. ainuv.GT.-100.))THEN

! Extract aerosol layer height of pixel from climatological database (CALIOP OR MODEL)
!==============================================================================
   status = GetAerosolLayer_CALIOPHeight(mon, plat, plon, zaer2)
   IF(status /= 1) THEN
      WRITE( *,'(A)') "Error getting CALIOP Aerosol Layer Height"
      RETURN
   ENDIF

   IF (zaer2 .LT. 0) THEN 
   status = GetAerosolLayerHeight(mon, plat, plon, zaer2)
     IF(status /= 1) THEN
       WRITE( *,'(A)') "Error getting MODEL Aerosol Layer Height"
       RETURN
     ENDIF
   ENDIF
   
!==============================================================================
! Extract single scattering albedo (388 nm) of pixel from climatological database
!==============================================================================
  status = GetOMI_SSAlbedoClimValue(yr, mon, dy, plat, plon, &
                 atype_aerac, ssacval354,ssacval388,ssacval500)

  IF(status /= 1) THEN
     WRITE( *,'(A)') "Error getting Single Scattering Albedo value"
     RETURN
  ENDIF

! Stitching ssacval354, ssacval388, ssacval500 into a single array
  inputssa_aca(1) = ssacval354
  inputssa_aca(2) = ssacval388
  inputssa_aca(3) = ssacval500



! ==================================================================
! First, ACAOD and COD retrievals as a function of ALH
! ==================================================================
!
! Only for carbonaceous smoke aerosol type...
IF(atypeTmp .EQ. 1 .AND. ainuv .GE. 0.8)THEN

  DO izhgt = 1, nzhgt   !ALH Loop...

  inzhgt = zhgt_nodes(izhgt)

  CALL Interpol_aac_LUTparams(sun_za,sat_za,phi, ocean, atypeTmp, pterrp, inzhgt, inssacval, salb(1), salb(2),  &
                                    fint_rad388_aac, fint_uvaimie_aac)

   ! -- Transpose of 'fint_uvaimie_aac' and 'fint_rad388_aac'
   fint_uvaimie_aac2 = transpose(fint_uvaimie_aac)         ! [8 cod x 7 aod]
   toa388_aac2 = transpose(fint_rad388_aac)
   ratradval = rad_obs(1)/rad_obs(2)  ! ratio (354/388)

       !=================================================
       ! -- Retrieve AOD and COD @ 500nm over clouds-----
       IF(ainuv .GE. 0.8)THEN
          status = nearuv_aac2(ainuv, ainuv, fint_uvaimie_aac2, rad_obs(2), toa388_aac2, &
                              aodtau_aac_wl3, codtau_aac_wl3)
       ENDIF
       !=================================================


       IF(aodtau_aac_wl3 .GT. 0.0 .AND. aodtau_aac_wl3 .LE. aodtaumax .AND. &
          codtau_aac_wl3 .GT. 0.0 .AND. codtau_aac_wl3 .LE. codtaumax)THEN
         status = GetKRIndex(atypeTmp, aodtau_aac_wl3, inssacval, aodtau_aac_tmp)
         codtau_aac_tmp = codtau_aac_wl3
       ELSE
         aodtau_aac_tmp(:) = -9999.
         codtau_aac_tmp = -9999.
       ENDIF

       aodtau_aac_zhgt(:,izhgt) = aodtau_aac_tmp(:)
       codtau_aac_zhgt(izhgt) = codtau_aac_tmp

  ENDDO          !DO izhgt = 1, nzhgt loop ends here.
  
ENDIF            !Aerosol type selection IF condition.



! ============================================================================
! Assigning aerosol layer height
      if (zaer2 .lt. 3.0) then
         inzhgt = 3.0   ! temporary default value
      else if (zaer2 .gt. 6.0) then
         inzhgt = 6.0   ! temporary default value
      else
         inzhgt = zaer2 ! CALIOP aerosol layer height value from climatology
      endif

      if (ssacval388 .lt. 0.1 .or. ssacval388 .gt. 1.0) then
          inssacval = 0.90 ! temporary default SSA value at 388 nm 
      else
          inssacval = ssacval388
      endif

! --------------------------------------------------

        CALL Interpol_aac_LUTparams(sun_za,sat_za,phi,ocean,atype_aerac,pterrp, inzhgt, inssacval, salb(1), salb(2),  &
                                    fint_rad388_aac, fint_uvaimie_aac)
      
       ! -- Compute the LER 388 for the Aerosol Above Clouds  -----
           fint_uvaimie_aac2 = transpose(fint_uvaimie_aac)         ! [8 cod x 7 aod]
           toa388_aac2 = transpose(fint_rad388_aac)
           ratradval = rad_obs(1)/rad_obs(2)  ! ratio (354/388)


       !***********************************************************
       ! -- First, Retrieve Apparent COD @ 500nm over clouds -----
       !
       !IF(ainuv.GE.-100.)THEN
       !   status = nearuv_cld(ratradval, fint_uvaimie_aac2, rad_obs(2), toa388_aac2, &
       !                       aodtau_aac_wl3, appcodtau_aac_wl3)
       !ENDIF
       !***********************************************************
       

       !***********************************************************
       ! Cubic Spline Interpolation for ApparentCloudOpticalDepth

       IF(ainuv .GE. -100.)THEN
        !toa388_aac2(:,1) contains radiances of clouds with ACAOD=0.0
         toa388_aac2_array = toa388_aac2(:,1)
         CALL SPLINE(toa388_aac2_array,codtbl_aac,ncod_aac,rad_obs(2),appcodtau_aac_wl3)
       ENDIF
       !***********************************************************
       

       !***********************************************************
       ! -- Retrieve AOD and COD @ 500nm over clouds -----
       !
       IF(ainuv.GE.0.8)THEN
          status = nearuv_aac2(ainuv,ainuv,fint_uvaimie_aac2, rad_obs(2), toa388_aac2, &
                              aodtau_aac_wl3, codtau_aac_wl3)
	  !PRINT *, 'AC-AEROSOL MODULE = ', aodtau_aac_wl3, codtau_aac_wl3		      
       ENDIF
       !***********************************************************


       ! -- Retrieve only COD @ 500nm over clouds -----

         IF(aodtau_aac_wl3.GT.0.0 .AND. aodtau_aac_wl3.LE. aodtaumax .AND. &
            codtau_aac_wl3.GT.0.0 .AND. codtau_aac_wl3.LE. codtaumax)THEN
         status = GetKRIndex(atype_aerac,aodtau_aac_wl3,inssacval,aodtau_aac)
         ELSE
         aodtau_aac(:) = -9999.
         codtau_aac_wl3 = -9999.
         ENDIF

        IF(status /= 1) THEN
          WRITE( *,'(A)' ) "Error finding AOD and COD @388nm"
          RETURN
        ELSE
          !WRITE( *,'(A)' ) "Success finding AOD and COD @388nm"

          !----------------------------------------------------------------------------------------------------------------------------------       
          ! Modification History: 'AlgQF' value changed to '0' if 'aodtau_aac' and 'codtau_aac' are being retrieved within their valid limits

          IF( reflect(2).GE.0.25 .AND. ainuv.GE.1.3 .AND. &
	     aodtau_aac_wl3 .GT. 0 .AND. aodtau_aac_wl3 .LE. aodtaumax .AND. &
	     codtau_aac_wl3 .GT. 0 .AND. codtau_aac_wl3 .LE. codtaumax) THEN
            AlgQF = 0
          ENDIF
          
          ! AlgQF = 1
          IF( reflect(2) .GE. 0.20 .AND. reflect(2) .LT. 0.25 .AND. &
	      aodtau_aac_wl3 .GT. 0 .AND. aodtau_aac_wl3 .LE. aodtaumax .AND. &
	      codtau_aac_wl3 .GT. 0 .AND. codtau_aac_wl3 .LE. codtaumax) THEN
            AlgQF = 1
          ENDIF

          ! AlgQF = 2
          !ainuv (mie) threshold changed to 0.8 and 1.3 (ainuv (ler) threshold 0.5 and 1.0)
          IF( reflect(2) .GE. 0.25 .AND. ainuv .GE. 0.8 .AND. &
	      ainuv .LT. 1.3 .AND. aodtau_aac_wl3 .GT. 0 .AND. &
	      aodtau_aac_wl3 .LE. aodtaumax .AND. codtau_aac_wl3 .GT. 0 .AND. &
	      codtau_aac_wl3 .LE. codtaumax) THEN
            AlgQF = 2
          ENDIF

          !AlgQF = 3
          IF((sun_za .GT. 55 .AND. sca .LE. 100 .AND. ainuv .LE. 2).OR. & 
             (sat_za .GT. 55 .AND. sca .LE. 100 .AND. ainuv .LE. 2).OR. & 
             (sun_za .GT. 60 .AND. sca .LE. 130 .AND. ainuv .LE. 2).AND. &
             aodtau_aac_wl3 .GT. 0 .AND. aodtau_aac_wl3 .LE. aodtaumax.AND. &
             codtau_aac_wl3 .GT. 0 .AND. codtau_aac_wl3 .LE. codtaumax) THEN
            AlgQF = 3
          ENDIF

        ENDIF


        !Uncertainty Calculations...
        IF (aodtau_aac_wl3 .GT. 0) THEN

        !-------------
        DO iSSA = 1, 2
        !-------------

        IF (iSSA .EQ. 1) inssacval_tmp = inssacval - 0.03
        IF (iSSA .EQ. 2) inssacval_tmp = inssacval + 0.03

        CALL Interpol_aac_LUTparams(sun_za,sat_za,phi,ocean, atype_aerac, pterrp, inzhgt, inssacval_tmp, &
	                            salb(1), salb(2), fint_rad388_aac, fint_uvaimie_aac)


        ! -- Compute the LER 388 for the Aerosol Above Clouds  -----
           fint_uvaimie_aac2 = transpose(fint_uvaimie_aac)         ! [8 cod x 7 aod]
           toa388_aac2 = transpose(fint_rad388_aac)
           ratradval = rad_obs(1)/rad_obs(2)  ! ratio (354/388)


       !***********************************************************
       ! -- Retrieve AOD and COD @ 500nm over clouds -----
       !
       IF(ainuv .GE. 0.8)THEN
          status = nearuv_aac2(ainuv,ainuv,fint_uvaimie_aac2, rad_obs(2), toa388_aac2, &
                              aodtau_aac_wl3, codtau_aac_wl3_tmp)
       ENDIF
       !***********************************************************

        IF(aodtau_aac_wl3 .GT. 0.0 .AND. aodtau_aac_wl3 .LE. aodtaumax.AND. &
           codtau_aac_wl3_tmp .GT. 0.0 .AND. codtau_aac_wl3_tmp .LE. codtaumax)THEN
         status = GetKRIndex(atype_aerac,aodtau_aac_wl3,inssacval_tmp,aodtau_aac_tmp)

         uncaodtau_aac(iSSA) = (aodtau_aac_tmp(2) - aodtau_aac(2))/aodtau_aac(2)*100.
         unccodtau_aac(iSSA) = (codtau_aac_wl3_tmp - codtau_aac_wl3)/codtau_aac_wl3*100.
        ENDIF


        ENDDO   ! DO iSSA = 1, 2

        ENDIF   ! IF (aodtau_aac_wl3.GT.0) THEN


        !
        DEALLOCATE(fint_rad388_aac,fint_uvaimie_aac)
!
   END IF ! (AI .GE. 1.0)
!
ENDIF

!
END FUNCTION OCI_NUV_ACAer_Process




 FUNCTION GetKRIndex(atype_aerac,tau_wl3,inssa,tau) RESULT(status)
!==============================================================================
!
! TITLE:
!     Calulate the Aerosol Parameters at other wavelengths and get K-refractive
!     index
!
! NAME:
!     GetKrefractIndex
!
! INPUTS:
!     none
!
! OUTPUTS:
!     tau           Optical Depth of ground pixel
!     ssa           Single Scattering Albedo of ground pixel
!     ima           K-refractive index of ground pixel
!
! RETURNS:
!     1       Successful completion
!
! FUNCTIONS_CALLED:
!     findex
!     OMI_SMF_set*
!
! HISTORY:
!     14-Dec-2004 EKK Initial version
!     04-Aug-2005 JAW Modified to use 500nm tau and ssa as initial data
!
!==============================================================================

!==============================================================================
! This function interpolates Aerosol Thickness & Aerosol Single Scattering 
! Albedo at 388nm based on the values @500nm and also finds the K refractive
! index   --based on toms .
!==============================================================================

  IMPLICIT NONE  

!==============================================================================
!  Optical Depth (tau), Single Scatt. Albedo (ssa) & Krefractive Index (ima)
!==============================================================================
  REAL(KIND=4), INTENT(IN)                  :: tau_wl3,inssa
  REAL(KIND=4), INTENT(INOUT), DIMENSION(3) :: tau
  INTEGER(KIND=2), INTENT(IN)               :: atype_aerac

  INTEGER (KIND=4) :: status
  REAL(KIND=4), DIMENSION(5) :: coef_tau_wl2, coef_ssa_wl2, coef_ima_wl2, coef_ima500, &
                                coef_tau_wl1, coef_ssa_wl1, coef_ima_wl1, coef_ssa_wl2_aac_smk, coef_ssa_wl2_aac_dst
  REAL(KIND=4) :: conv_factor_wl2, conv_factor_wl1 
  
  REAL(KIND=4) :: ssa_wl3

  status = 1


 IF (atype_aerac .EQ. 1) THEN ! This is smoke
     !New set of coefficients for smoke by Hiren Jethva (USRA/GESTAR) modified on July 15, 2015
     coef_tau_wl2 = (/ 1258.4245,-5607.2688,9358.9007,-6926.2923,1917.7405 /)
     coef_ssa_wl2 = (/ 73.642888,-342.32462,596.30233,-457.93537,131.31539 /)
     coef_ima500  = (/ -19.763742,89.817572,-151.72530,113.23520,-31.563801 /)
     coef_tau_wl1 = (/ 1974.6877,-8799.8327,14684.693,-10865.865,3008.0110 /) ! valid when wl1 = 354 nm
     coef_ssa_wl1 = (/ 202.65598,-906.60985,1520.0294,-1129.3267,314.25161 /) ! valid when wl1 = 354 nm
     coef_ima_wl1 = (/ 1.5250317,0.58923326,-10.798525,14.132119,-5.4479988 /) ! valid when wl1 = 354 nm

     coef_ssa_wl2_aac_smk = (/-25.404624, 123.45677, -219.03895, 172.84429, -50.858004/)
  ELSEIF (atype_aerac .EQ. 2) THEN ! This is dust

     ! -- Coefficients for the spheroid DUST models from Hiren Jethva on July 24, 2015 --
     coef_tau_wl2 = (/ -9.6183858,  46.836601, -76.176117, 54.991750, -14.869036/)
     coef_ssa_wl2 = (/ -133.31842,  563.82931, -887.35332, 618.13487, -160.29174/)
     coef_ima500  = (/ 4.4228746,  -18.338290, 28.790138,  -20.226356,5.351622/)
     coef_ima_wl2 = (/ 9.6276955, -39.844269, 62.476901, -43.861191, 11.600838/)
     coef_tau_wl1 = (/-15.111942, 70.725592, -114.66243, 82.544919, -22.262197/) ! valid when wl1 = 354 nm
     coef_ssa_wl1 = (/-63.528719, 247.86385, -350.32612, 211.88245, -44.890627/) ! valid when wl1 = 354 nm
     coef_ima_wl1 = (/13.584530, -56.263953, 88.276371, -62.001400, 16.404418 /)! valid when wl1 = 354 nm

     coef_ssa_wl2_aac_dst = (/10.715972, -49.210211, 88.437284, -68.556941, 19.613582/)
  ELSE ! This would be industrial
     coef_tau_wl2 = (/ 0.79012994,1.3184211,-0.50053125,0.0,0.0 /)
     coef_ssa_wl2 = (/ 0.20486168,0.74138000,0.053688153,0.0,0.0 /)
     coef_ima500  = (/ 0.097054681,0.17799631,-0.55774476,0.28269098,0.0 /)
     coef_tau_wl1 = (/ 0.33975817,2.6208016,-1.09151010,0.0,0.0 /) ! valid when wl1 = 354 nm
     coef_ssa_wl1 = (/ 0.21182544,0.77122350,0.016930148,0.0,0.0 /)! valid when wl1 = 354 nm
  ENDIF
  

  !-----------------------------------------------------------------
  !Modification history: July 15, 2015 By Hiren Jethva (USRA/GESTAR)
  ! First convert the Clim or Daily value of SSA388 to SSA500 
     
  IF (atype_aerac .EQ. 1) THEN 
  ssa_wl3=coef_ssa_wl2_aac_smk(1)+coef_ssa_wl2_aac_smk(2)*inssa+coef_ssa_wl2_aac_smk(3)*inssa**2+ &
                  coef_ssa_wl2_aac_smk(4)*inssa**3+coef_ssa_wl2_aac_smk(5)*inssa**4
  ENDIF

  IF (atype_aerac .EQ. 2) THEN
  ssa_wl3=coef_ssa_wl2_aac_dst(1)+coef_ssa_wl2_aac_dst(2)*inssa+coef_ssa_wl2_aac_dst(3)*inssa**2+ &
                  coef_ssa_wl2_aac_dst(4)*inssa**3+coef_ssa_wl2_aac_dst(5)*inssa**4
  ENDIF

  !-----------------------------------------------------------------


  ! Conversion of tau500 to tau388 and tau342

  conv_factor_wl2=coef_tau_wl2(1)+coef_tau_wl2(2)*ssa_wl3+coef_tau_wl2(3)*ssa_wl3**2+ &
              coef_tau_wl2(4)*ssa_wl3**3+coef_tau_wl2(5)*ssa_wl3**4
  tau(2) = tau_wl3*conv_factor_wl2

  conv_factor_wl1=coef_tau_wl1(1)+coef_tau_wl1(2)*ssa_wl3+coef_tau_wl1(3)*ssa_wl3**2+ &
              coef_tau_wl1(4)*ssa_wl3**3+coef_tau_wl1(5)*ssa_wl3**4
  tau(1) = tau_wl3*conv_factor_wl1

  tau(3)=tau_wl3

!  !Comment out the rest of the wavelength conversion scheme: July 15, 2015
!  ! Conversion of ssa500 to ssa388 and ssa342
!
!  ssa(2) = coef_ssa_wl2(1)+coef_ssa_wl2(2)*ssa(3)+coef_ssa_wl2(3)*ssa(3)**2+ &
!           coef_ssa_wl2(4)*ssa(3)**3+coef_ssa_wl2(5)*ssa(3)**4
!
!  ssa(1) = coef_ssa_wl1(1)+coef_ssa_wl1(2)*ssa(3)+coef_ssa_wl1(3)*ssa(3)**2+ &
!           coef_ssa_wl1(4)*ssa(3)**3+coef_ssa_wl1(5)*ssa(3)**4
!
!  ! Conversion of ssa500 to imaginary refractive index at 500nm (ima500)
!
!  ima(3) = coef_ima500(1)+coef_ima500(2)*ssa(3)+coef_ima500(3)*ssa(3)**2+ &
!           coef_ima500(4)*ssa(3)**3+coef_ima500(5)*ssa(3)**4
!
!  IF (atype_aerac .EQ. 2) THEN
!
!     ! Conversion of ssa500 to imaginary refractive index at 388nm and 342nm
!     ! Use fit only for dust
!
!      ima(2) = coef_ima_wl2(1)+coef_ima_wl2(2)*ssa(3)+coef_ima_wl2(3)*ssa(3)**2+&
!               coef_ima_wl2(4)*ssa(3)**3+coef_ima_wl2(5)*ssa(3)**4
!
!      ima(1) = coef_ima_wl1(1)+coef_ima_wl1(2)*ssa(3)+coef_ima_wl1(3)*ssa(3)**2+&
!               coef_ima_wl1(4)*ssa(3)**3+coef_ima_wl1(5)*ssa(3)**4
! 
!  ELSE IF (atype_aerac .EQ. 1) THEN ! for smoke
!      ima(2) = ima(3)
!      !-- New imaginary at 354 nm from Hiren Jethva, Hampton University --------------------------------------------------
!      !ima(1) = coef_ima_wl1(1)+coef_ima_wl1(2)*ssa(3)+coef_ima_wl1(3)*ssa(3)**2+coef_ima_wl1(4)*ssa(3)**3+coef_ima_wl1(5)*ssa(3)**4 ! valid when wl1 = 354 nm
!  ELSE ! for industrial:
!      ima(2) = ima(3)
!      ima(1) = ima(3)
!  ENDIF
!
!  IF (ABS(ima(3)) .LE. 1e-4) THEN 
!     ima(3) = 0.0
!  ENDIF
!  IF (ABS(ima(2)) .LE. 1e-4) THEN 
!     ima(2) = 0.0
!  ENDIF
!

END FUNCTION GetKRIndex
!==============================================================================




 FUNCTION NEARUV_cld(muvai,uvai_aac,rad2_obs,rad2, &
                     aodtau,codtau) RESULT(status)
!==============================================================================
! TITLE:
!     Produce COD over clouds
! NAME:
!     NEARUV_cldaac
!
!      INPUTS:
!      muvai          real        Calculated UVAI without clouds 
!      uvai_aac       real        Calculated UVAI with clouds  
!      rad2           real        Calculated radiance for aerosol models at UV2
!      rad2_obs       real        Observed radiance at UV2
!
!      OUTPUTS 
!      aodtau         real        Aerosol optical thickness above clouds at UV2 (absorbing)
!      codtau         real        cloud optical thickness at UV2 (absorbing)
!
!==============================================================================

  IMPLICIT NONE

  REAL(KIND=4), INTENT(IN)                                :: muvai
  REAL(KIND=4), INTENT(IN),  DIMENSION(ncod_aac,naod_aac) :: uvai_aac, rad2
  REAL(KIND=4), INTENT(IN)                                :: rad2_obs
  REAL(KIND=4), INTENT(OUT)                               :: aodtau, codtau
!
  REAL(KIND=4)                                      :: ratio_obs, rad2wi
  REAL(KIND=4), DIMENSION(ncod_aac)                 :: rad2_mod
!  
  REAL(KIND=4), DIMENSION(naod_aac)                 :: w0g2, tauw2
  REAL(KIND=4), DIMENSION(naod_aac)                 :: rad2w
  REAL(KIND=4), DIMENSION(naod_aac)                 :: q_tmpset

!==============================================================================
! Nodal points of cloud optical depth values *** tau ***
!==============================================================================
 REAL(KIND=4), DIMENSION(10) :: codtbl_aac
 DATA codtbl_aac/0.000,2.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0/

!==============================================================================
! Tau values and indexes from tau table
!==============================================================================
  REAL(KIND=4)  :: q1, q2, q_tmp, tau
  INTEGER       :: itau, itau2, jw0, jw02
  INTEGER 	:: iw0sel, ngw02

  INTEGER(KIND=4)        :: status, ierr

  status = -1

!==============================================================================
! Initialize variables
!==============================================================================
  ngw02=0
  aodtau = -9999.
  codtau = -9999.
  ratio_obs = muvai

!  First index of naod_aac is for AOD equals zero--cloud only LUT
  DO iw0sel = 1,1

!==============================================================================
! Select the calculated UVAI_aac for this cod profile 
!==============================================================================
     rad2_mod(:) = rad2(:,iw0sel)

!==============================================================================
! Find inital tau in the table of 354 nm and 388 nm ratios
!==============================================================================
     status = FindTableEntry(rad2_obs,rad2_mod,ncod_aac,q2,q1,itau2,itau,q_tmp)

     IF((itau .EQ. itau2) .AND. (itau .EQ. 1)) itau2 = itau2+1
     IF((itau .EQ. itau2) .AND. (itau .NE. 1)) itau2 = itau2-1
     
     !===========================================================================!
     IF(itau.GE.1.AND.itau.LE.8.AND.itau2.GE.1.AND.itau2.LE.8.AND.q_tmp.LE.1)THEN
   
     q_tmp = (q1-rad2_obs)/(q1 - q2)
     q_tmpset(iw0sel) = q_tmp

!==============================================================================
! Calculate the interpolated AOD tau value
!==============================================================================
     codtau = codtbl_aac(itau) - q_tmp * (codtbl_aac(itau) - codtbl_aac(itau2))

     IF(codtau.LE.0.OR.codtau.GT.50.OR.q_tmp.GT.1)THEN
     codtau = -9999.
     ENDIF

!==============================================================================
! If initial tau and radiance are not zero add then to the new model grid
! of radiances for each Aerosol model
!
!   Check q_tmp for the AOD interpolation ONLY. No extrapolation !!! 
!==============================================================================

    ENDIF ! IF condition on ITAU and ITAU2
  
  ENDDO ! DO iw0sel = 1,nw0sel     

  status = 1

  RETURN

 END FUNCTION NEARUV_cld





 FUNCTION NEARUV_aac2(ainuv,muvai,uvai_aac,rad2_obs, rad2, &
                     aodtau,codtau) RESULT(status)
!==============================================================================
! TITLE:
!     Produce AOD and COD over clouds
! NAME:
!     NEARUV_aac
!
!      INPUTS:
!      muvai          real        Calculated UVAI without clouds 
!      uvai_aac       real        Calculated UVAI with clouds  
!      rad2           real        Calculated radiance for aerosol models at UV2
!      rad2_obs       real        Observed radiance at UV2
!
!      OUTPUTS 
!      aodtau         real        Aerosol optical thickness above clouds at UV2 (absorbing)
!      codtau         real        cloud optical thickness at UV2 (absorbing)
!==============================================================================


  IMPLICIT NONE

  REAL(KIND=4), INTENT(IN)                                :: muvai
  REAL(KIND=4), INTENT(IN),  DIMENSION(ncod_aac,naod_aac) :: uvai_aac, rad2
  REAL(KIND=4), INTENT(IN)                                :: rad2_obs
  REAL(KIND=4), INTENT(OUT)                               :: aodtau, codtau
!
  REAL(KIND=4)                                      :: ratio_obs
  REAL(KIND=4), DIMENSION(naod_aac)                 :: ratio_mod
  REAL(KIND=4)                                      :: rad2wi
! 
  REAL(KIND=4), DIMENSION(ncod_aac)                 :: w0g2, tauw2, tauw2_tmp
  REAL(KIND=4), DIMENSION(ncod_aac)                 :: rad2w, rad2w_tmp
  REAL(KIND=4), DIMENSION(ncod_aac)                 :: q_tmpset, codtbl_aac_tmp

  REAL(KIND=4)                                      :: aodtaumax, codtaumax


!==============================================================================
! Nodal points of cloud optical depth values *** tau ***
!==============================================================================
 REAL(KIND=4), DIMENSION(7) :: aodtbl_aac
 REAL(KIND=4), DIMENSION(10) :: codtbl_aac

 DATA aodtbl_aac/0.000,0.1,0.5,1.0,2.5,4.0,6.0/
 DATA codtbl_aac/0.000,2.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0/

!==============================================================================
! Tau values and indexes from tau table
!==============================================================================
  REAL(KIND=4)  :: q1, q2, q_tmp, tau, ainuv
  INTEGER 	:: iw0sel, ngw02, counter
  INTEGER       :: itau, itau2, jw0, jw02
  INTEGER(KIND=4)        :: status, ierr

  status = -1

  !500-nm AOD/COD Max LUT values...
  aodtaumax = 6.0
  codtaumax = 70.0


!==============================================================================
! Initialize variables
!==============================================================================
  ngw02=0
  aodtau = -9999.
  codtau = -9999.
  ratio_obs = muvai

!==============================================================================
! Loop over all AOD nodal points
!==============================================================================

  DO iw0sel = 1,ncod_aac

     ratio_mod(:) = uvai_aac(iw0sel,:)

!==============================================================================
! Find inital tau in the table of UV-AI (Mie) 
!==============================================================================
     status = FindTableEntry(ratio_obs,ratio_mod,naod_aac,q2,q1,itau2,itau,q_tmp)
     IF((itau .EQ. itau2) .AND. (itau .EQ. 1)) itau2 = itau2+1
     IF((itau .EQ. itau2) .AND. (itau .NE. 1)) itau2 = itau2-1
     
     !===========================================================!
     IF (itau .GE. 1 .AND. itau .LE. 7 .AND. itau2 .GE. 1 .AND. itau2 .LE. 7) THEN  
     q_tmp = (q1-ratio_obs)/(q1 - q2)
     q_tmpset(iw0sel) = q_tmp

!==============================================================================
! Calculate the interpolated AOD tau value
!==============================================================================
     tau = aodtbl_aac(itau) - q_tmp * (aodtbl_aac(itau) - aodtbl_aac(itau2))

!==============================================================================
! Interpolate on tau to find model radiance for each COD Model at AOD tau
!==============================================================================
     rad2wi = rad2(iw0sel,itau) + (rad2(iw0sel,itau2) - rad2(iw0sel,itau)) &
                                * (tau - aodtbl_aac(itau))              &
                                / (aodtbl_aac(itau2) - aodtbl_aac(itau))

!==============================================================================
! If initial tau and radiance are not zero add then to the new model grid
! of radiances for each Aerosol model
!
!   Check q_tmp for the AOD interpolation ONLY. No extrapolation !!! 
!==============================================================================
       IF ((tau .GT. 0) .AND. (tau .LT. 7) .AND. (rad2wi .GT. 0)) THEN !  388.0 nm
          ngw02 = ngw02 + 1
          w0g2(ngw02) = codtbl_aac(iw0sel)
          tauw2(ngw02) = tau
          rad2w(ngw02) = rad2wi    
       ENDIF ! IF((tau .GT. 0) .AND. (rad2wi .GT. 0)) THEN

  
    ENDIF ! IF condition on ITAU and ITAU2

  ENDDO ! DO iw0sel = 1,nw0sel     


!==============================================================================
! Interpolate on observed radiance to find the Aerosol model that fits best
!==============================================================================
  IF (ngw02 .GT. 1) THEN
     status = FindTableEntry(rad2_obs,rad2w(1:ngw02),ngw02,q1,q2,jw0,jw02,q_tmp)
     IF ((jw0 .EQ. jw02) .AND. (jw0 .EQ. 1)) jw02 = jw02+1
     IF ((jw0 .EQ. jw02) .AND. (jw02 .EQ. ngw02)) jw0 = jw02-1 
     q_tmp = (rad2_obs - q1) / (q2 - q1)


      !***New SPLINE way of performing interpolation***!

      !Checking IF observed reflectance @ 860 nm is within LUT range...
      IF (rad2_obs.GE.rad2w(1) .AND. rad2_obs.LE.rad2w(ncod_aac)) THEN

      !Preparing arrays of interpolated reflectance and COD nodes with valid positive values...
      counter = 1

      rad2w_tmp(:) = 0.0
      tauw2_tmp(:) = 0.0
      codtbl_aac_tmp(:) = 0.0


      DO iw0sel = 1,ncod_aac
      IF (rad2w(iw0sel).GT.0.AND.rad2w(iw0sel).LE.1) THEN
      rad2w_tmp(counter) = rad2w(iw0sel)
      tauw2_tmp(counter) = tauw2(iw0sel)
      codtbl_aac_tmp(counter) = codtbl_aac(iw0sel)
      counter = counter + 1
      ENDIF
      ENDDO


      !AOD_TAU
      IF (counter-1.GE.5) THEN
      CALL SPLINE(rad2w_tmp,tauw2_tmp,counter-1,rad2_obs,aodtau)

      !COD_TAU
      CALL SPLINE(rad2w_tmp,codtbl_aac_tmp,counter-1,rad2_obs,codtau)
      ELSE

      aodtau = tauw2(jw0) + q_tmp * (tauw2(jw02) - tauw2(jw0))
      codtau =  w0g2(jw0) + q_tmp *  (w0g2(jw02) -  w0g2(jw0))

      !  -- Checking interpolation indices to see out of bounds --
      IF (q_tmpset(jw0) .GT. 1.0 .OR. q_tmpset(jw02) .GT. 1.0 .OR. q_tmp .GT. 1.0 .OR. q_tmp .LT. 0.0) THEN
         aodtau = -9999.
         codtau = -9999.
      ENDIF

      ENDIF     !IF(counter-1.GE.5)THEN


     !  -- Checking the retrieved aodtau and codtau to see out of bounds --
     IF (aodtau.LE.0 .OR. aodtau.GT.aodtaumax .OR. codtau.LE.0 .OR. codtau.GT.codtaumax) THEN
        aodtau = -9999.
        codtau = -9999.
     ENDIF

     ENDIF ! IF(rad2_obs.GE.rad2w(1) .AND. rad2_obs.LE.rad2w(ncod_aac))THEN

     ELSE
     aodtau = -9999.
     codtau = -9999.

  ENDIF ! IF(ngw02.GT.1)THEN

  status = 1

  RETURN

 END FUNCTION NEARUV_aac2


!***************************************
SUBROUTINE SPLINE(XI, FI, N_NODES, X, F)
!***************************************

  IMPLICIT NONE
  INTEGER                                       :: N, M, I, K
  REAL(KIND=4), INTENT(IN)                      :: X
  INTEGER(KIND=4), INTENT(IN)                   :: N_NODES
  REAL(KIND=4)                                  :: DX, H, ALPHA, BETA, GAMMA, ETA

  REAL(KIND=4), DIMENSION (10), INTENT(IN)       :: XI, FI
  REAL(KIND=4), DIMENSION (10)                   :: P2
  REAL(KIND=4), INTENT(OUT)                     :: F


  N = N_NODES -1
  M = 2

  CALL CUBIC_SPLINE(N, XI, FI, P2)


! Find the interval that x resides
    K = 1
    DX = X-XI(1)
    DO WHILE (DX .GE. 0 .AND. K .LE. N_NODES-1)
      K = K + 1
      DX = X-XI(K)
    END DO
    K = K - 1

    !Initialize...
    F = 0.0
   IF (K .GE. 1 .AND. K .LE. N_NODES-1) THEN
! Find the value of function f(x)
    DX = XI(K+1) - XI(K)
    
    ALPHA = P2(K+1)/(6*DX)
    BETA = -P2(K)/(6*DX)
    GAMMA = FI(K+1)/DX - DX*P2(K+1)/6
    
    ETA = DX*P2(K)/6 - FI(K)/DX
    F = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) &
       +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) &
       +GAMMA*(X-XI(K))+ETA*(X-XI(K+1))
    ENDIF
    
END SUBROUTINE SPLINE


SUBROUTINE CUBIC_SPLINE (N, XI, FI, P2)
!
! Function to carry out the cubic-spline approximation
! with the second-order derivatives returned.
!
  INTEGER :: I
  INTEGER, INTENT (IN) :: N
  REAL, INTENT (IN), DIMENSION (N+1):: XI, FI
  REAL, INTENT (OUT), DIMENSION (N+1):: P2
  REAL, DIMENSION (N):: G, H
  REAL, DIMENSION (N-1):: D, B, C
!
! Assign the intervals and function differences
!
  DO I = 1, N
    H(I) = XI(I+1) - XI(I)
    G(I) = FI(I+1) - FI(I)
  END DO
!
! Evaluate the coefficient matrix elements   
  DO I = 1, N-1
    D(I) = 2*(H(I+1)+H(I)) 
!  Fix to Zero divide   
        if(H(I) .eq.0 )H(I) = .00001 
        if(H(I+1) .eq.0 )H(I+1) = .00001 
    
    B(I) = 6*(G(I+1)/H(I+1)-G(I)/H(I))
    C(I) = H(I+1)
  END DO
!
! Obtain the second-order derivatives
!
  CALL TRIDIAGONAL_LINEAR_EQ (N-1, D, C, C, B, G)
  P2(1) = 0
  P2(N+1) = 0
  DO I = 2, N
    P2(I) = G(I-1)
  END DO
END SUBROUTINE CUBIC_SPLINE
!
SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
!
! Functione to solve the tridiagonal linear equation set.
!
  INTEGER, INTENT (IN) :: L
  INTEGER :: I
  REAL, INTENT (IN), DIMENSION (L):: D, E, C, B
  REAL, INTENT (OUT), DIMENSION (L):: Z
  REAL, DIMENSION (L):: Y, W
  REAL, DIMENSION (L-1):: V, T
!
! Evaluate the elements in the LU decomposition
!
  W(1) = D(1)
  V(1)  = C(1)
  T(1)  = E(1)/W(1)
  DO I = 2, L - 1
    W(I) = D(I)-V(I-1)*T(I-1)
    V(I) = C(I)
    T(I) = E(I)/W(I)
  END DO
  W(L) = D(L)-V(L-1)*T(L-1)
!
! Forward substitution to obtain y
!
  Y(1) = B(1)/W(1)
  DO I = 2, L
    Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I)
  END DO
!
! Backward substitution to obtain z
  Z(L) = Y(L)
  DO I = L-1, 1, -1
    Z(I) = Y(I) - T(I)*Z(I+1)
  END DO
END SUBROUTINE TRIDIAGONAL_LINEAR_EQ

 
END MODULE NUV_ACAerosolModule 

