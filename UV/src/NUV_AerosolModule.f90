MODULE NUV_AerosolModule

 IMPLICIT NONE
           
  PUBLIC :: SetAerosolModel
  PUBLIC :: Get_3DModelRad_using_fmf
  PUBLIC :: OCI_NUV_Aer_Process
  PUBLIC :: InterpTransmittance
  PUBLIC :: InterpRadiance
  
  
 CONTAINS

!==============================================================================
!==============================================================================

 FUNCTION OCI_NUV_Aer_Process(Month, nWavel, wavelen, plat, plon, sun_za, sat_za, phi, pterrp, gpQF, &
                    salb, rad_obs, atype, uvai, aod, fmf, retssa, rethgt) RESULT(status)

  USE LookupTableModule_nc4 
  USE InterpolationModule
  USE Nearuv_alhssa_Module
  USE regpolymonial_predict_ssa
  

  IMPLICIT NONE  

  INTEGER(KIND=4),             INTENT(IN)  :: Month
  INTEGER(KIND=4),             INTENT(IN)  :: nWavel
  INTEGER(KIND=2),             INTENT(IN)  :: gpQF
  REAL(KIND=4),                INTENT(IN)  :: plat, plon, sun_za, sat_za, phi, pterrp
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: wavelen(nWavel)
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  ::    salb(nWavel)
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: rad_obs(nWavel)
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: AOD(9)		      
  REAL(KIND=4),                INTENT(IN)  :: uvai, fmf	
  REAL(KIND=4),                INTENT(OUT) :: retssa(5), rethgt      
  INTEGER(KIND=2),             INTENT(OUT) :: AType
  
  REAL(KIND=4), DIMENSION(nzae,nw0sel)     :: rad2d_354mod, rad2d_388mod
  REAL(KIND=4), DIMENSION(nzae,nw0sel,ntau):: rad3d_354mod, rad3d_388mod
  
  REAL(KIND=4)    :: FMF1, fmode_retssa(5), cmode_retssa(5), EAE, tmp_fmf
  REAL(KIND=4), PARAMETER :: uvai_thresh=0.8, dsi_thresh = 0
  REAL(KIND=4), PARAMETER :: fmf_thresh1 = 0.30, fmf_thresh2=0.70
  REAL(KIND=4), PARAMETER :: EAE_thresh1 = 0.40, EAE_thresh2=1.40
  INTEGER(KIND=4) :: status, indexmd(7)
  INTEGER(KIND=2) :: AType1

STATUS = -1

!=======================================================================================================
! IF Valid AOD exists Over Oceans ==> Begin the retrieval process using mix of radiances (fmode + cmode).
!=======================================================================================================
IF ( AOD(1) .GT. 0 .AND. AOD(2) .GT. 0 .AND. &
        FMF .GE. 0 .AND.    FMF .LE. 1 .AND. gpQF .EQ. 17) THEN 

 if (fmf .le. fmf_thresh1) then 
   Atype1 = 2
   tmp_fmf = 0.0
 endif
 
 if (fmf .ge. fmf_thresh2) then 
   Atype1 = 1
   tmp_fmf = 1.0
 endif
 
 if (fmf .gt. fmf_thresh1 .AND. fmf .lt. fmf_thresh2) then 
   Atype1 = 12
   tmp_fmf = (fmf-fmf_thresh1)/(fmf_thresh2-fmf_thresh1)
 endif

!Atype1 = 2
!tmp_fmf = 0.0 
status = Get_3DModelRad_using_fmf(nWavel, wavelen, sun_za, sat_za, phi, pterrp, &
                                  salb, tmp_fmf, radiance_toa)
rad3d_354mod = radiance_toa(1,:,:,:)
rad3d_388mod = radiance_toa(2,:,:,:)
status = accountforAOD_to_get_2Drad(atype1, rad3d_354mod, rad3d_388mod, &
                         AOD(1:2), tmp_fmf, rad2d_354mod, rad2d_388mod)
status = nearuv_alhssa(atype1,rad2d_354mod,rad2d_388mod,rad_obs(1),rad_obs(2), &
                             tmp_fmf, retssa, rethgt)

IF (retssa(2) .LT. 0) THEN 
status = nearuv_alhssa_alhloop(atype1,rad2d_354mod,rad2d_388mod,rad_obs(1),rad_obs(2), &
                             tmp_fmf, retssa, rethgt)
ENDIF


ENDIF


!=======================================================================================================
! IF Valid AOD exists Over Land ==> Begin the retrieval process by selecting atype based on (EAE 488-670).
!=======================================================================================================
IF ( AOD(1) .GT. 0 .AND. AOD(2) .GT. 0 .AND. gpQF .NE. 17) THEN 

 EAE = -1.0 * alog(AOD(3)/AOD(5)) / alog(488.0/670.0)
 if (uvai .le. uvai_thresh) Atype1 = 3
 if (uvai .gt. uvai_thresh .AND. &
     EAE .LE. EAE_thresh1) Atype1 = 2
 if (uvai .gt. uvai_thresh .AND. &
     EAE .GE. EAE_thresh2) Atype1 = 1
 if (uvai .gt. uvai_thresh .AND. &
     EAE .GT. EAE_thresh1  .AND. & 
     EAE .LT. EAE_thresh2) Atype1 = 12
  
  tmp_fmf = (EAE-EAE_thresh1)/(EAE_thresh2 - EAE_thresh1)
!  Atype1 = 2
!  tmp_fmf = 0.0
  if (Atype1 .EQ. 12) then 
     status = Get_3DModelRad_using_fmf(nWavel, wavelen, sun_za, sat_za, phi, pterrp, &
                                  salb, tmp_fmf, radiance_toa)
  else
     status = Get_3DModelRad_using_atyp(nWavel, wavelen, sun_za, sat_za, phi, pterrp, &
                                  salb, atype1, radiance_toa) 				    
  endif
  

rad3d_354mod = radiance_toa(1,:,:,:)
rad3d_388mod = radiance_toa(2,:,:,:)
status = accountforAOD_to_get_2Drad(atype1, rad3d_354mod, rad3d_388mod, &
                         AOD(1:2), tmp_fmf, rad2d_354mod, rad2d_388mod)
status = nearuv_alhssa(atype1,rad2d_354mod,rad2d_388mod,rad_obs(1),rad_obs(2), &
                             tmp_fmf, retssa, rethgt)
IF (retssa(2) .LT. 0) THEN 
status = nearuv_alhssa_alhloop(atype1,rad2d_354mod,rad2d_388mod,rad_obs(1),rad_obs(2), &
                             tmp_fmf, retssa, rethgt)
ENDIF


ENDIF


!=======================================================================================================
! IF UVAI is less than threshold, the aerosols are assumed at surface level (ALH=0) and retrieve SSA.
!=======================================================================================================
 IF (uvai .le. uvai_thresh .AND. retssa(2) .LT. 0) THEN
 status = nearuv_getssa_surfacealh(rad2d_354mod,rad2d_388mod,rad_obs(1),rad_obs(2), &
                                   retssa, rethgt)
 ENDIF


!=======================================================================================================
! Use polynomial coefficients derived for regional seasonal absorption models (Kayetha et al 2022, AMT).
!=======================================================================================================
 atype = Atype1
 if (retssa(2) .GT. 0.00 .AND. retssa(2) .LE. 1.00) then
    if (atype1 .NE. 12) then 
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, retssa(2), 354.0, retssa(1))
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, retssa(2), 480.0, retssa(3))
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, retssa(2), 550.0, retssa(4))
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, retssa(2), 670.0, retssa(5))
    endif
    
    if (atype1 .EQ. 12) then 
      atype = 1
      fmode_retssa(2) = retssa(2)
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, fmode_retssa(2), 354.0, fmode_retssa(1))
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, fmode_retssa(2), 480.0, fmode_retssa(3))
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, fmode_retssa(2), 550.0, fmode_retssa(4))
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, fmode_retssa(2), 670.0, fmode_retssa(5))

      atype = 2
      cmode_retssa(2) = retssa(2)
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, cmode_retssa(2), 354.0, cmode_retssa(1))
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, cmode_retssa(2), 480.0, cmode_retssa(3))
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, cmode_retssa(2), 550.0, cmode_retssa(4))
       CALL predict_wavssa_from_ssa388(plat, plon, month, atype, cmode_retssa(2), 670.0, cmode_retssa(5))

      retssa = (fmode_retssa * tmp_fmf) + (cmode_retssa * (1.0-tmp_fmf))
      atype = 12
    endif
 endif

status = 1

Return

END FUNCTION OCI_NUV_Aer_Process
		   

!==============================================================================
!==============================================================================

FUNCTION Get_3DModelRad_using_fmf(nWavel, wavelen, sun_za, sat_za, phi, pterrp, &
                                  salb, fmf, wv_mod3Drad) RESULT(status)


  USE GetLUT_H5module_nc4
  USE LookupTableModule_nc4 
  USE InterpolationModule

  IMPLICIT NONE  

! Function Inputs & output
  INTEGER(KIND=4),             INTENT(IN)  :: nWavel 
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: wavelen(nWavel)
  REAL(KIND=4),                INTENT(IN)  :: sun_za, sat_za, phi
  REAL(KIND=4),                INTENT(IN)  :: pterrp
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: salb(nWavel)
  REAL(KIND=4),    	       INTENT(IN)  :: fmf 	
!
  REAL(KIND=4),DIMENSION(nWavel,nzae,nw0sel,ntau),INTENT(OUT):: wv_mod3Drad
  INTEGER         :: status

! Profile indices for Aerosol model
  INTEGER(KIND=2) :: AType1, AType2
  INTEGER(KIND=4) :: fmode_indx(7), cmode_indx(7)

!==============================================================================
! Wavelength indexes (make local to L2OMAERUVModule)
!==============================================================================
  INTEGER(KIND=4) :: iwl(3),iwaveT,iwavel
  REAL(KIND=4) :: fracWave

!  Bounding values in table for input parameters
  REAL(KIND=4)    :: theta1,theta2,sza1,sza2,phi1,phi2,w1,w2
  REAL(KIND=4)    :: fractionalTheta,fractionalSZA,fractionalPHI
  INTEGER(KIND=4) :: jtheta, jtheta2, jsza, jsza2, jphi, jphi2, jw1

! Radiance and Transmittance from lookup tables.
  REAL(KIND=4), DIMENSION(nzae,nw0sel,ntau) :: rad_out,trm_out

! Surface flux  
  REAL(KIND=4) :: sflux_lut
  DATA sflux_lut/3.1415926/

! Pressure height  (make local to L2OMAERUVModule)
  REAL(KIND=4) :: pressure_table(2), logFacPressure
  DATA pressure_table/1013.,600./

  INTEGER(KIND=4) :: i

! Determine the Terrain Pressure Fraction and geometry indices.
  logFacPressure = (LOG(pressure_table(1)) - LOG(pterrp))           &
                  /(LOG(pressure_table(1)) - LOG(pressure_table(2)))

  status=FindTableEntry(sat_za,theta_table,ntheta,theta1,theta2, &
                                 jtheta,jtheta2,fractionalTheta)
  status=FindTableEntry(sun_za,sza_table,nsza,sza1,sza2, &
                                 jsza,jsza2,fractionalSZA)
  status=FindTableEntry(phi,phi_table,nphi,phi1,phi2, &
                                 jphi,jphi2,fractionalPHI) 

! Assume SMOKE and get SSA nodes.
AType1 = 1
status = SetAerosolModel(AType1,fmode_indx)

! Assume DUST and get SSA nodes.
AType2 = 2
status = SetAerosolModel(AType2,cmode_indx)


! Begin loop over selected wavelengths.
  DO iwavel = 1,nWavel
     iwaveT = iwavel
     iwl(iwavel) = iwaveT

!  Calculate Radiance at sea-level for Ground Pixel (using theta, sza & phi)
     status = InterpRadiance(10,iwavel,sat_za,jtheta,theta_table,sun_za,jsza, &
                           sza_table,phi,jphi,phi_table,1,nzae,fmode_indx(1), &
                                                fmode_indx(7),1,ntau,rad_out)
     radiance_geo_p10_fmode(iwaveT,:,:,:) = rad_out
     status = InterpRadiance(10,iwavel,sat_za,jtheta,theta_table,sun_za,jsza, &
                           sza_table,phi,jphi,phi_table,1,nzae,cmode_indx(1), &
                                               cmode_indx(7),1,ntau,rad_out)
     radiance_geo_p10_cmode(iwaveT,:,:,:) = rad_out
!
     status = InterpTransmittance(10,iwavel,sat_za,jtheta,theta_table,sun_za, &
                                         jsza,sza_table,1,nzae,fmode_indx(1), &
                                               fmode_indx(7),1,ntau,trm_out)
     transmittance_geo_p10_fmode(iwaveT,:,:,:) = trm_out
     status = InterpTransmittance(10,iwavel,sat_za,jtheta,theta_table,sun_za, &
                                         jsza,sza_table,1,nzae,cmode_indx(1), &
                                               cmode_indx(7),1,ntau,trm_out)
     transmittance_geo_p10_cmode(iwaveT,:,:,:) = trm_out
!
     status = InterpRadiance(6,iwavel,sat_za,jtheta,theta_table,sun_za,jsza, &
                          sza_table,phi,jphi,phi_table,1,nzae,fmode_indx(1), &
                                              fmode_indx(7),1,ntau,rad_out)
     radiance_geo_p06_fmode(iwaveT,:,:,:) = rad_out
     status = InterpRadiance(6,iwavel,sat_za,jtheta,theta_table,sun_za,jsza, &
                          sza_table,phi,jphi,phi_table,1,nzae,cmode_indx(1), &
                                              cmode_indx(7),1,ntau,rad_out)
     radiance_geo_p06_cmode(iwaveT,:,:,:) = rad_out
!
     status = InterpTransmittance(6,iwavel,sat_za,jtheta,theta_table,sun_za, &
                                        jsza,sza_table,1,nzae,fmode_indx(1), &
                                              fmode_indx(7),1,ntau,trm_out)
     transmittance_geo_p06_fmode(iwaveT,:,:,:) = trm_out
     status = InterpTransmittance(6,iwavel,sat_za,jtheta,theta_table,sun_za, &
                                        jsza,sza_table,1,nzae,cmode_indx(1), &
                                              cmode_indx(7),1,ntau,trm_out)
     transmittance_geo_p06_cmode(iwaveT,:,:,:) = trm_out
!
     radiance_geo_fmode(iwaveT,:,:,:) = radiance_geo_p10_fmode(iwaveT,:,:,:)  &
                                      -(radiance_geo_p10_fmode(iwaveT,:,:,:)  & 
                                      - radiance_geo_p06_fmode(iwaveT,:,:,:)) &
                                      * logFacPressure
     radiance_geo_cmode(iwaveT,:,:,:) = radiance_geo_p10_cmode(iwaveT,:,:,:)  &
                                      -(radiance_geo_p10_cmode(iwaveT,:,:,:)  & 
                                      - radiance_geo_p06_cmode(iwaveT,:,:,:)) &
                                      * logFacPressure
!
     transmittance_geo_fmode(iwaveT,:,:,:) = transmittance_geo_p10_fmode(iwaveT,:,:,:)  &
                                           -(transmittance_geo_p10_fmode(iwaveT,:,:,:)  &
                                           - transmittance_geo_p06_fmode(iwaveT,:,:,:)) &
                                           * logFacPressure
     transmittance_geo_cmode(iwaveT,:,:,:) = transmittance_geo_p10_cmode(iwaveT,:,:,:)  &
                                           -(transmittance_geo_p10_cmode(iwaveT,:,:,:)  &
                                           - transmittance_geo_p06_cmode(iwaveT,:,:,:)) &
                                           * logFacPressure
!
! --  Solution : use selected model indices (.fmode_indx.) correctly --
    sphalb_fmode(iwaveT,:,:,:) = sbp10(iwaveT,:,fmode_indx(1):fmode_indx(7),:)  &
                               -(sbp10(iwaveT,:,fmode_indx(1):fmode_indx(7),:)  &
                               -  sbp6(iwaveT,:,fmode_indx(1):fmode_indx(7),:)) &
                               * logFacPressure
    sphalb_cmode(iwaveT,:,:,:) = sbp10(iwaveT,:,cmode_indx(1):cmode_indx(7),:)  &
                               -(sbp10(iwaveT,:,cmode_indx(1):cmode_indx(7),:)  &
                               -  sbp6(iwaveT,:,cmode_indx(1):cmode_indx(7),:)) &
                               * logFacPressure
!
!  Calculate Radiance at the top of the atmosphere
     radiance_toa_fmode(iwaveT,:,:,:) = radiance_geo_fmode(iwaveT,:,:,:) + salb(iwavel) & 
                                 * transmittance_geo_fmode(iwaveT,:,:,:) &
                         /(1.0-salb(iwavel) * sphalb_fmode(iwaveT,:,:,:))
     radiance_toa_cmode(iwaveT,:,:,:) = radiance_geo_cmode(iwaveT,:,:,:) + salb(iwavel) & 
                                 * transmittance_geo_cmode(iwaveT,:,:,:) &
                         /(1.0-salb(iwavel) * sphalb_cmode(iwaveT,:,:,:))
     wv_mod3Drad(iwaveT,:,:,:) = ( radiance_toa_fmode(iwaveT,:,:,:) * fmf ) + &
                                 ( radiance_toa_cmode(iwaveT,:,:,:) * (1.-fmf) )
  ENDDO ! DO iwavel = 1,nWavel
! End of Wavelength Loop

wv_mod3Drad(:,:,:,:) = wv_mod3Drad(:,:,:,:)/sflux_lut
!
status = 1

END FUNCTION Get_3DModelRad_using_fmf

!
!==============================================================================
!==============================================================================
!

FUNCTION Get_3DModelRad_using_atyp(nWavel, wavelen, sun_za, sat_za, phi, pterrp, &
                                  salb, atyp, wv_mod3Drad) RESULT(status)


  USE GetLUT_H5module_nc4
  USE LookupTableModule_nc4 
  USE InterpolationModule
  USE Nearuv_alhssa_Module

  IMPLICIT NONE  

! Function Inputs & output
  INTEGER(KIND=4),             INTENT(IN)  :: nWavel 
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: wavelen(nWavel)
  REAL(KIND=4),                INTENT(IN)  :: sun_za, sat_za, phi
  REAL(KIND=4),                INTENT(IN)  :: pterrp
  REAL(KIND=4),   DIMENSION(:),INTENT(IN)  :: salb(nWavel)
!
  REAL(KIND=4),DIMENSION(nWavel,nzae,nw0sel,ntau),INTENT(OUT):: wv_mod3Drad
  INTEGER              :: status

! Profile indices for Aerosol model
  INTEGER(KIND=2), INTENT(IN) :: ATyp
  INTEGER(KIND=4) :: model_indx(7)

!==============================================================================
! Wavelength indexes (make local to L2OMAERUVModule)
!==============================================================================
  INTEGER(KIND=4) :: iwl(3),iwaveT,iwavel
  REAL(KIND=4) :: fracWave

!  Bounding values in table for input parameters
  REAL(KIND=4)    :: theta1,theta2,sza1,sza2,phi1,phi2,w1,w2
  REAL(KIND=4)    :: fractionalTheta,fractionalSZA,fractionalPHI
  INTEGER(KIND=4) :: jtheta, jtheta2, jsza, jsza2, jphi, jphi2, jw1

! Radiance and Transmittance from lookup tables.
  REAL(KIND=4), DIMENSION(nzae,nw0sel,ntau) :: rad_out,trm_out

! Surface flux  
  REAL(KIND=4) :: sflux_lut
  DATA sflux_lut/3.1415926/

! Pressure height  (make local to L2OMAERUVModule)
  REAL(KIND=4) :: pressure_table(2), logFacPressure
  DATA pressure_table/1013.,600./

  INTEGER(KIND=4) :: i

! Determine the Terrain Pressure Fraction and geometry indices.
  logFacPressure = (LOG(pressure_table(1)) - LOG(pterrp))           &
                  /(LOG(pressure_table(1)) - LOG(pressure_table(2)))

  status=FindTableEntry(sat_za,theta_table,ntheta,theta1,theta2, &
                                 jtheta,jtheta2,fractionalTheta)
  status=FindTableEntry(sun_za,sza_table,nsza,sza1,sza2, &
                                 jsza,jsza2,fractionalSZA)
  status=FindTableEntry(phi,phi_table,nphi,phi1,phi2, &
                                 jphi,jphi2,fractionalPHI) 

! Use selected atyp
status = SetAerosolModel(ATyp,model_indx)



! Begin loop over selected wavelengths.
  DO iwavel = 1,nWavel
     iwaveT = iwavel
     iwl(iwavel) = iwaveT

!  Calculate Radiance at sea-level for Ground Pixel (using theta, sza & phi)
     status = InterpRadiance(10,iwavel,sat_za,jtheta,theta_table,sun_za,jsza, &
                           sza_table,phi,jphi,phi_table,1,nzae,model_indx(1), &
                                                model_indx(7),1,ntau,rad_out)
     radiance_geo_p10_fmode(iwaveT,:,:,:) = rad_out
     status = InterpTransmittance(10,iwavel,sat_za,jtheta,theta_table,sun_za, &
                                         jsza,sza_table,1,nzae,model_indx(1), &
                                               model_indx(7),1,ntau,trm_out)
     transmittance_geo_p10_fmode(iwaveT,:,:,:) = trm_out
!
     status = InterpRadiance(6,iwavel,sat_za,jtheta,theta_table,sun_za,jsza, &
                          sza_table,phi,jphi,phi_table,1,nzae,model_indx(1), &
                                              model_indx(7),1,ntau,rad_out)
     radiance_geo_p06_fmode(iwaveT,:,:,:) = rad_out
     status = InterpTransmittance(6,iwavel,sat_za,jtheta,theta_table,sun_za, &
                                        jsza,sza_table,1,nzae,model_indx(1), &
                                              model_indx(7),1,ntau,trm_out)
     transmittance_geo_p06_fmode(iwaveT,:,:,:) = trm_out
!
     radiance_geo_fmode(iwaveT,:,:,:) = radiance_geo_p10_fmode(iwaveT,:,:,:)  &
                                      -(radiance_geo_p10_fmode(iwaveT,:,:,:)  & 
                                      - radiance_geo_p06_fmode(iwaveT,:,:,:)) &
                                      * logFacPressure
     transmittance_geo_fmode(iwaveT,:,:,:) = transmittance_geo_p10_fmode(iwaveT,:,:,:)  &
                                           -(transmittance_geo_p10_fmode(iwaveT,:,:,:)  &
                                           - transmittance_geo_p06_fmode(iwaveT,:,:,:)) &
                                           * logFacPressure
!
! --  Solution : use selected model indices (.model_indx.) correctly --
    sphalb_fmode(iwaveT,:,:,:) = sbp10(iwaveT,:,model_indx(1):model_indx(7),:)  &
                               -(sbp10(iwaveT,:,model_indx(1):model_indx(7),:)  &
                               -  sbp6(iwaveT,:,model_indx(1):model_indx(7),:)) &
                               * logFacPressure
!
!  Calculate Radiance at the top of the atmosphere
     wv_mod3Drad(iwaveT,:,:,:) = radiance_geo_fmode(iwaveT,:,:,:) + salb(iwavel) & 
                                 * transmittance_geo_fmode(iwaveT,:,:,:) &
                         /(1.0-salb(iwavel) * sphalb_fmode(iwaveT,:,:,:))
  ENDDO ! DO iwavel = 1,nWavel
! End of Wavelength Loop

wv_mod3Drad(:,:,:,:) = wv_mod3Drad(:,:,:,:)/sflux_lut
!
status = 1

END FUNCTION Get_3DModelRad_using_atyp

!
!==============================================================================
!==============================================================================
!


 FUNCTION InterpRadiance(file,iwave,theta,ThetaIdx,theta_arr,sza,SZAIdx,sza_arr, &
                         phi,PHIIdx,phi_arr,nzae1,nzae2, &
                         nw01,nw02,ntau1,ntau2,Rad) RESULT(status)
!==============================================================================
!
!	This subroutine interpolates radiance values once the Lookup tables
!	have been reduced in size and converted to vectorial form in order to
!	speed up the program by using more efficiently the indices
!
!	Written by Marcos Andrade based on f77 TOMS code (2007)
!
!==============================================================================
  USE LookupTableModule_nc4

  IMPLICIT NONE

  REAL(KIND=4),                   INTENT(IN)    :: theta,sza,phi
  INTEGER(KIND=4),                INTENT(IN)    :: ThetaIdx,SZAIdx,PHIIdx
  INTEGER(KIND=4),                INTENT(IN)    :: file, iwave
  INTEGER(KIND=4),                INTENT(IN)    :: nzae1,nw01,ntau1
  INTEGER(KIND=4),                INTENT(IN)    :: nzae2,nw02,ntau2
    REAL(KIND=4), DIMENSION(:,:,:), INTENT(OUT)   :: Rad
    
  REAL(KIND=4), DIMENSION(7)  :: sza_arr
  REAL(KIND=4), DIMENSION(11) :: phi_arr
  REAL(KIND=4), DIMENSION(14) :: theta_arr

  INTEGER(KIND=4)    :: izae, iw0, jw0, itau
  INTEGER(KIND=4)    :: ierr, polyInterp
  REAL(KIND=4)       :: dy
 
  INTEGER :: i,k,m,n
  INTEGER (KIND=4) :: isza,ithe,iphi, status
  INTEGER (KIND=4) :: idxgral(64), idxaux, model
  
  status = -1

 n = 0
 do k=1,4
    ithe = indscn+k-1
    do m=1,4
       iphi = indphi+m-1
       do i=1,4
          isza = indsol+i-1
          n = n+1
           idxgral(n) = iwave+ (isza-1)*2*nzae*nw0model*ntau + &
                        (iphi-1)*2*nzae*nw0model*ntau*nsza + &
                        (ithe-1)*2*nzae*nw0model*ntau*nsza*nphi
       enddo
    enddo
 enddo

DO izae = nzae1,nzae2
   DO iw0 = nw01,nw02
      jw0 = iw0 - nw01 + 1
      DO itau = ntau1,ntau2
 
         idxaux = (izae-1)*2 + (iw0-1)*2*nzae + (itau-1)*2*nzae*nw0model 
         IF (file == 10) THEN
            Rad(izae,jw0,itau) = SUM( radlin_p10(idxgral + idxaux)*cofs )
         ELSE
            Rad(izae,jw0,itau) = SUM( radlin_p06(idxgral + idxaux)*cofs )
         END IF
      ENDDO  ! DO itau = 1,ntau
   ENDDO  ! DO iw0  = 1,nw0
ENDDO  ! DO izae = 1,nzae

 STATUS = 1
 
 RETURN

 END FUNCTION InterpRadiance
!
!==============================================================================
!==============================================================================
!
 FUNCTION InterpTransmittance(file,iwave,theta,ThetaIdx,theta_arr,sza,SZAIdx, &
                              sza_arr,nzae1,nzae2, &
                              nw01,nw02,ntau1,ntau2,Trans) RESULT(status)
!==============================================================================
!
!	This subroutine interpolates transmittance values once the lookup tables
!	have been reduced in size and converted to vectorial form in order to
!	speed up the program by using more efficiently the indices
!
!	Written by Marcos Andrade based on f77 TOMS code (2007)
!
!==============================================================================
  USE LookupTableModule_nc4

  IMPLICIT NONE

  REAL(KIND=4),                   INTENT(IN)    :: theta,sza
  INTEGER(KIND=4),                INTENT(IN)    :: ThetaIdx,SZAIdx
  INTEGER(KIND=4),                INTENT(IN)    :: file, iwave
  INTEGER(KIND=4),                INTENT(IN)    :: nzae1,nw01,ntau1
  INTEGER(KIND=4),                INTENT(IN)    :: nzae2,nw02,ntau2
  REAL(KIND=4), DIMENSION(:,:,:), INTENT(OUT)   :: Trans

  REAL(KIND=4), DIMENSION(7)  :: sza_arr
  REAL(KIND=4), DIMENSION(14) :: theta_arr

  INTEGER(KIND=4)    :: izae, iw0, jw0, itau
  INTEGER(KIND=4)    :: ierr, polyInterp
  REAL(KIND=4)       :: dy
 
  INTEGER :: i,k,m,n
  INTEGER (KIND=4) :: isza,ithe,iphi,status
  INTEGER (KIND=4) :: idxgral(16), idxaux
  
  status = -1

  n = 0
  DO k=1,4
     ithe = indscn_2d+k-1
     DO i=1,4
        isza = indsol_2d+i-1
        n = n+1
        idxgral(n) = iwave+ (isza-1)*2*nzae*nw0model*ntau + &
                    (ithe-1)*2*nzae*nw0model*ntau*nsza
     ENDDO
  ENDDO

  DO izae = nzae1,nzae2
     DO iw0 = nw01,nw02
        jw0 = iw0 - nw01 + 1
        DO itau = ntau1,ntau2
 
           idxaux = (izae-1)*2 + (iw0-1)*2*nzae + (itau-1)*2*nzae*nw0model 
           IF (file == 10) THEN
              Trans(izae,jw0,itau) = SUM( tralin_p10(idxgral + idxaux)*cofs_2d )
           ELSE
              Trans(izae,jw0,itau) = SUM( tralin_p06(idxgral + idxaux)*cofs_2d )
           END IF

        ENDDO  ! DO itau = 1,ntau
     ENDDO  ! DO iw0  = 1,nw0
  ENDDO  ! DO izae = 1,nzae

 STATUS = 1
 RETURN

 END FUNCTION InterpTransmittance
!
!==============================================================================
!==============================================================================
!  
 
  FUNCTION CalcAerosolModel(gpQF,coval,salb,AeroIndexNUV,AeroIndexVIS, AeroType,indexmd, &
                           mon,xlat,xlon,reflect2) RESULT(status)

  USE Nearuv_alhssa_Module
  
  IMPLICIT NONE  

!==============================================================================
! Function Inputs
!==============================================================================
!  Aerosol Index (UV, VIS)
  REAL(KIND=4),    INTENT(IN)  :: COval 
  REAL(KIND=4),    INTENT(IN)  :: AeroIndexNUV , AeroIndexVIS
  REAL(KIND=4),                    INTENT(IN)  :: xlat
  REAL(KIND=4),                    INTENT(IN)  :: xlon
  REAL(KIND=4),                    INTENT(IN)  :: reflect2
  INTEGER(KIND=4),                 INTENT(IN)  :: mon
  REAL(KIND=4),                    INTENT(IN)  :: salb
  INTEGER(KIND=2),                 INTENT(IN)  :: gpQF

!==============================================================================
! Function Outputs
!==============================================================================
  INTEGER(KIND=2),               INTENT(OUT) :: AeroType
  INTEGER(KIND=4),               INTENT(OUT) :: indexmd(7)

!==============================================================================
! Aerosol index for near UV and Visible wavelengths
!==============================================================================
  REAL(KIND=4) :: AeroIndexVIS0, AeroIndexNUV0(1)
  REAL(KIND=4) :: delR, recoval 
  REAL(KIND=4) :: slope1, slope2, gain1, gain2

! CO thesholds for N and S. hemispheres 
  REAL(KIND=4) :: cothres1, cothres2
  
!==============================================================================
! Latitude and Longitude indexes in the surface category catalog
!==============================================================================
  REAL(KIND=4)    :: tmplon
  INTEGER(KIND=4) :: ilat,ilon, status, ierr

  status = 1

!==============================================================================
! Minimun value for the Aerosol Index
!============================================================================== 
   IF (gpQF .EQ. 16 .OR. gpQF .EQ. 7) THEN
      AeroIndexVIS0 = VisAIThreshold(1)
   ELSEIF (gpQF .EQ. 17) THEN
      AeroIndexVIS0 = VisAIThreshold(2)
   ELSE
      AeroIndexVIS0 = VisAIThreshold(3)
   ENDIF
   AeroIndexNUV0 = UVAIThreshold

!============================================================================
! New aerosol type selection scheme with the AIRS CO and UV Aerosol Index 
!============================================================================
recoval = coval!/1.0E38

! -- Use the interpolated CO threshold in between 10 N and 10 S latitudes ----
IF (xlat .gt. 10.0) then
    cothres1 = 2.0
else if (xlat .lt. -10.0) then
    cothres1 = 1.6
else
    slope1 = (2.0 - 1.6)/(10 + 10)
    gain1 = 2.0 - slope1*(10)
    cothres1 = gain1 + slope1* xlat
endif

!PRINT *, coval, recoval, xlat, cothres1

! -- Interpolated CO threshold between 10N and 10S for high CO values only ---------------
IF ( (AeroIndexNUV .GE. 0.8) .AND. (recoval .GT. cothres1) .AND. (xlat .LT. -10.0) ) THEN
     AeroType = 1
     status = SetAerosolModel(AeroType,indexmd)
     IF(status /= 1) THEN
        WRITE( *,'(A)' ) "Error setting AerosolType"
        RETURN
     ENDIF
ELSE IF ( (AeroIndexNUV .GE. 0.8) .AND. (recoval .GT. cothres1) .AND. (xlat .GT. 10.0) ) THEN
     AeroType = 1
     status = SetAerosolModel(AeroType,indexmd)
     IF(status /= 1) THEN
        WRITE( *,'(A)' ) "Error setting AerosolType"
        RETURN
     ENDIF

!!  -- Use interpolated CO thrsholds between 10N and 10S latitudes -------------------------------------------
ELSE IF ( (AeroIndexNUV .GE. 0.8) .AND. (recoval .GT. cothres1) .AND. (xlat .GE. -10.0 .AND. xlat .LE. 10.0) ) THEN
     AeroType = 1
     status = SetAerosolModel(AeroType,indexmd)
     IF(status /= 1) THEN
        WRITE( *,'(A)' ) "Error setting AerosolType"
        RETURN
     ENDIF
!! SULFATE model is used when high CO level & Low AI < 0.5  is found -------------------------------------------
ELSE IF ( (AeroIndexNUV .LT. 0.8) .AND. (recoval .GT. cothres1) .AND. (xlat .LT. -10.0) ) THEN
     AeroType = 3
     status = SetAerosolModel(AeroType,indexmd)
     IF(status /= 1) THEN
        WRITE( *,'(A)' ) "Error setting AerosolType"
        RETURN
     ENDIF
ELSE IF ( (AeroIndexNUV .LT. 0.8) .AND. (recoval .GT. cothres1) .AND. (xlat .GT. 10.0) ) THEN
     AeroType = 3
     status = SetAerosolModel(AeroType,indexmd)
     IF(status /= 1) THEN
        WRITE( *,'(A)' ) "Error setting AerosolType"
        RETURN
     ENDIF
ELSE IF ( (AeroIndexNUV .LT. 0.8) .AND. (recoval .GT. cothres1) .AND. (xlat .GE. -10.0 .AND. xlat .LE. 10.0) ) THEN
     AeroType = 3
     status = SetAerosolModel(AeroType,indexmd)
     IF(status /= 1) THEN
        WRITE( *,'(A)' ) "Error setting AerosolType"
        RETURN
     ENDIF
! -----------------------------------------------------------------------------------------------
ELSE IF ( AeroIndexNUV .GE. 0.8) THEN
     AeroType = 2
     status = SetAerosolModel(AeroType,indexmd)
    IF(status /= 1) THEN
       WRITE( *,'(A)' ) "Error setting AerosolType"
       RETURN
    ENDIF
ELSE IF ( (gpQF .EQ. 16 .OR. gpQF .EQ. 7) ) THEN
     AeroType = 2
     status = SetAerosolModel(AeroType,indexmd)
    IF(status /= 1) THEN
       WRITE( *,'(A)' ) "Error setting AerosolType"
       RETURN
    ENDIF
ELSE
  AeroType = 3
  status = SetAerosolModel(AeroType,indexmd)
  IF(status /= 1) THEN
     WRITE( *,'(A)' ) "Error setting AerosolType"
     RETURN
  ENDIF
ENDIF
!
status = 1
   
RETURN

END FUNCTION CalcAerosolModel

!
!==============================================================================
!==============================================================================
!


 
 FUNCTION ChoseAerosolModel(mon,xlat,xlon,ref_sfc,ainuv,AeroType,indexmd, &
                            reflect2, salb) RESULT(status)
!
! TITLE:
!     Chose the Aerosol model based on geographical considerations
!
! NAME:
!     ChoseAerosolModel
!
! INPUTS:
!     mon           The month of the observation
!     xlat          Latitude of the Ground Pixel
!     xlon          Longitude of the Ground Pixel
!     ref_sfc       Surface Category array
!     ainuv         The NUV Aerosol index
!
! OUTPUTS:
!     AeroType      The Aerosol Type (smoke, dust, industrial)
!     indexmd       The calcualted aerosol model

  USE Nearuv_alhssa_Module
  
  IMPLICIT NONE

! Function Inputs
  INTEGER(KIND=4),                 INTENT(IN)  :: mon
  REAL(KIND=4),                    INTENT(IN)  :: xlat
  REAL(KIND=4),                    INTENT(IN)  :: xlon
  INTEGER(KIND=4), DIMENSION(:,:), INTENT(IN)  :: ref_sfc
  REAL(KIND=4),                    INTENT(IN)  :: ainuv
  REAL(KIND=4),                    INTENT(IN)  :: reflect2
  REAL(KIND=4),                    INTENT(IN)  :: salb

! local variable
  REAL(KIND=4)                                 :: wfact, wai

! Function Outputs
  INTEGER(KIND=2), INTENT(OUT)                 :: AeroType
  INTEGER(KIND=4), INTENT(OUT)                 :: indexmd(7)

! Latitude and Longitude indexes in the surface category catalog
  REAL(KIND=4)    :: tmplon
  INTEGER(KIND=4) :: ilat,ilon
  INTEGER(KIND=4) :: surf

! Boundaries for dust region east of the Sahara desert
  REAL(KIND=4)                :: d1lat1
  REAL(KIND=4), PARAMETER     :: d1lat2 =  55.0, d1lat2sah = 30.0 
  REAL(KIND=4), PARAMETER     :: d1lon1 = -85.0, d1lon2    = 30.0

  REAL(KIND=4), DIMENSION(12) :: d1lat1clm
  DATA d1lat1clm/13.,12.,11.,5.,5.,5.,5.,5.,5.,5.,12.,13./

  INTEGER(KIND=4)           :: status

  status = -1
  tmplon=0.0d00

!==============================================================================
!  If NUV Aerosol Index is less than 0.0 the set AerosolType to Industrial
!==============================================================================
  IF(ainuv .LT. 0.4) THEN
     AeroType = 3
     status = SetAerosolModel(AeroType,indexmd)
     IF(status /= 1) THEN
        WRITE( *,'(A)' ) "Error setting AerosolType - IND"
        RETURN
     ENDIF
  ENDIF

!==============================================================================
!  For all NH latitudes greater than 55N set AerosolType to Smoke
!==============================================================================
     IF(xlat .GT. 55.0) THEN
        AeroType = 1
        status = SetAerosolModel(AeroType,indexmd)
        IF(status /= 1) THEN
           WRITE( *,'(A)' ) "Error setting AerosolType - SMK"
           RETURN
        ENDIF
     ENDIF

!==============================================================================
!  For the NH Pacific Ocean between 25N and 55N and Japan set AerosolType to Dust
!==============================================================================
     IF(((xlat .LE.  55.0) .AND. (xlat .GT.  25.0)) .AND. &
        ((xlon .GE. 45.0) .AND. (xlon .LE. 180.0))) THEN
        AeroType = 2
        status = SetAerosolModel(AeroType,indexmd)
        IF(status /= 1) THEN
           WRITE( *,'(A)' ) "Error setting AerosolType - DST"
           RETURN
        ENDIF
     ENDIF

     IF(((xlat .LE.   55.0) .AND. (xlat .GT.   25.0)) .AND. &
        ((xlon .GE. -180.0) .AND. (xlon .LE. -115.0))) THEN
        AeroType = 2
        status = SetAerosolModel(AeroType,indexmd)
        IF(status /= 1) THEN
           WRITE( *,'(A)' ) "Error setting AerosolType - DST"
           RETURN
        ENDIF
     ENDIF

!============================================================
! Use a weighted UV AI to separate Smoke from Dust
!============================================================
    wfact = ((salb-reflect2)/reflect2)
    wai  = ainuv*wfact
!
    IF(((xlat .LE.  15.0) .AND. (xlat .GT.  -10.0)) .AND. &
        ((xlon .GE. -60.0) .AND. (xlon .LE. 40.0))) THEN
        wfact = ((salb-reflect2)/reflect2)
        wai  = ainuv*wfact
        IF (wai .LT. -1.3) THEN
            AeroType = 1
            status = SetAerosolModel(AeroType,indexmd)
            IF(status /= 1) THEN
               WRITE( *,'(A)' ) "Error setting AerosolType - SMK"
               RETURN
            ENDIF
        ELSE  
            AeroType = 2
            status = SetAerosolModel(AeroType,indexmd)
            IF(status /= 1) THEN
               WRITE( *,'(A)' ) "Error setting AerosolType - DST"
               RETURN
            ENDIF
        END IF
    ENDIF

!==============================================================================
!  For all SH latitudes south of -40S set AerosolType to Dust
!==============================================================================
     IF(xlat .LT. -40.0) THEN
        AeroType = 2
        status = SetAerosolModel(AeroType,indexmd)
        IF(status /= 1) THEN
           WRITE( *,'(A)' ) "Error setting AerosolType - DST"
           RETURN
        ENDIF
     ENDIF

!==============================================================================
!  Get the surface category out of the surface category reference catalog
!==============================================================================

  tmplon = xlon
  IF(xlon .LT. 0) tmplon=xlon+360.
  ilon = INT(tmplon) + 1
  ilat = INT(90.-xlat) + 1
  if (ilat .gt. 180) ilat = 180
  if (ilon .gt. 360) ilon = 360
  surf = ref_sfc(ilat,ilon) 

  IF (surf .EQ. 17) THEN

!==============================================================================
!  For all SH oceans between the Equator and 40S set AerosolType to Smoke
!==============================================================================
     IF((xlat .LT. 0.0) .AND. (xlat .GE. -40.0)) THEN
        AeroType = 1
        status = SetAerosolModel(AeroType,indexmd)
        IF(status /= 1) THEN
           WRITE( *,'(A)' ) "Error setting AerosolType - SMK"
           RETURN
        ENDIF
     ENDIF

!==============================================================================
!  For Northern Hemisphere Indian Ocean, Pacific Ocean up to 25N, 
!   and Gulf of Mexico set AerosolType to Smoke                    
!==============================================================================
     IF(((xlat .LE. 25.0) .AND. (xlat .GE.   0.0)) .AND. &
        ((xlon .GE. 80.0) .AND. (xlon .LE. 180.0))) THEN
        AeroType = 1
        status = SetAerosolModel(AeroType,indexmd)
        IF(status /= 1) THEN
           WRITE( *,'(A)' ) "Error setting AerosolType - SMK"
           RETURN
        ENDIF
     ENDIF
 
     IF(((xlat .LE.   25.0) .AND. (xlat .GE.   0.0)) .AND. &
        ((xlon .GE. -180.0) .AND. (xlon .LE. -85.0))) THEN
        AeroType = 1
        status = SetAerosolModel(AeroType,indexmd)
        IF(status /= 1) THEN
           WRITE( *,'(A)' ) "Error setting AerosolType - SMK"
           RETURN
        ENDIF
     ENDIF

!==============================================================================
!  For the Mediterranean and bodies of water in Eurasia and Arabian Peninsula 
!  set AerosolType to Dust
!==============================================================================
     IF(((xlat .LE. 55.0) .AND. (xlat .GT. 10.0)) .AND. &
        ((xlon .GE. 30.0) .AND. (xlon .LT. 80.0))) THEN
        AeroType = 2
        status = SetAerosolModel(AeroType,indexmd)
        IF(status /= 1) THEN
           WRITE( *,'(A)' ) "Error setting AerosolType - DST"
           RETURN
        ENDIF
     ENDIF
 
!==============================================================================
!  For Atlantic Ocean East of the Sahara desert set AerosolType to Dust
!   (Note: Southern boundary of this region changes with the month.) 
!==============================================================================
     d1lat1=d1lat1clm(mon)

     IF(((xlat .GE. d1lat1) .AND. (xlat .LE. d1lat2)) .AND. &
        ((xlon .GE. d1lon1) .AND. (xlon .LE. d1lon2))) THEN
        AeroType = 2
        status = SetAerosolModel(AeroType,indexmd)
        IF(status /= 1) THEN
           WRITE( *,'(A)' ) "Error setting AerosolType - DST"
           RETURN
        ENDIF
     ENDIF

  ENDIF !IF(surf .EQ. 17) THEN

!==============================================================================
!  For Atlantic Ocean East of the Sahara desert set AerosolType to Dust
!   (Note: Southern boundary of this region changes with the month.) 
!==============================================================================
  d1lat1=d1lat1clm(mon)

  IF(((xlat .GE. d1lat1) .AND. (xlat .LE. d1lat2sah)) .AND. &
     ((xlon .GE. d1lon1) .AND. (xlon .LE. d1lon2   ))) THEN
     AeroType = 2
     status = SetAerosolModel(AeroType,indexmd)
     IF(status /= 1) THEN
        WRITE( *,'(A)' ) "Error setting AerosolType - DST"
        RETURN
     ENDIF
  ENDIF

!==============================================================================
!  For Arid and semi-Arid regions set AerosolType to Dust
!==============================================================================
  IF((surf .EQ. 16) .or. (surf .EQ. 7)) THEN
     AeroType = 2
     status = SetAerosolModel(AeroType,indexmd)
     IF(status /= 1) THEN
        WRITE( *,'(A)' ) "Error setting AerosolType - DST"
        RETURN
     ENDIF
  ENDIF


  AeroType = 1 
  status = SetAerosolModel(AeroType,indexmd)
  IF(status /= 1) THEN
     WRITE( *,'(A)' ) "Error setting AerosolType - SMK"
     RETURN
  ENDIF

  status = 1

  RETURN

END FUNCTION ChoseAerosolModel

!
!==============================================================================
!==============================================================================
!
 FUNCTION SetAerosolModel(atype,indexmd) RESULT(status)

! TITLE:
!     Select the Aerosol model based on whether the aerosol is smoke, dust or 
!     industial pollutants.
! NAME:
!     SetAerosolModel
! INPUTS:
!     atype         The Aerosol type
!

  USE Nearuv_alhssa_Module
  
  IMPLICIT NONE  

! Function Inputs
  INTEGER(KIND=2), INTENT(IN)  :: atype

! Function Outputs
  INTEGER(KIND=4), INTENT(OUT) :: indexmd(7)

! Profile indices for Aerosol model
  INTEGER(KIND=4)              :: im, nw0sel, status

   status = -1
   nw0sel = SIZE(indexmd)

!==============================================================================
!   Choose carbonaceous models
!==============================================================================
   IF(atype .EQ. 1)THEN
      DO im = 1,nw0sel
            indexmd(im) = indexcs(im)
              w0sel(im) = w0model(indexmd(im)) 
           w0388sel(im) = w0388model(indexmd(im)) 
           Kr500sel(im) = Kr500model(indexmd(im)) 
           Kr388sel(im) = Kr388model(indexmd(im)) 
        tauRatiosel(im) = tauRatiomodel(indexmd(im)) 
      ENDDO ! DO im = 1,nw0sel

   ENDIF

!==============================================================================
!   Choose dust models
!==============================================================================
   IF(atype .EQ. 2)THEN
     DO im = 1,nw0sel
            indexmd(im) = indexds(im)
              w0sel(im) = w0model(indexmd(im)) 
           w0388sel(im) = w0388model(indexmd(im)) 
           Kr500sel(im) = Kr500model(indexmd(im)) 
           Kr388sel(im) = Kr388model(indexmd(im)) 
        tauRatiosel(im) = tauRatiomodel(indexmd(im)) 
     ENDDO ! DO im = 1,nw0sel

   ENDIF

!==============================================================================
!    Choose industrial pollutant
!==============================================================================
   IF(atype .EQ. 3)THEN
      DO im = 1,nw0sel
            indexmd(im) = indexind(im)
              w0sel(im) = w0model(indexmd(im)) 
           w0388sel(im) = w0388model(indexmd(im)) 
           Kr500sel(im) = Kr500model(indexmd(im)) 
           Kr388sel(im) = Kr388model(indexmd(im)) 
        tauRatiosel(im) = tauRatiomodel(indexmd(im)) 
      ENDDO ! DO im = 1,nw0sel

   ENDIF

   status = 1
   
   RETURN

END FUNCTION SetAerosolModel

!
!==============================================================================
!==============================================================================
!


END MODULE NUV_AerosolModule        
