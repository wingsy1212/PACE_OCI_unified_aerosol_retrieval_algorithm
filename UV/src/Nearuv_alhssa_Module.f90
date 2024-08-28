MODULE Nearuv_alhssa_Module

!==============================================================================
!
! FILENAME: Nearuv_alhssa_Module.f90
!
! AUTHOR: Vinay Kayetha (SSAI/NASA-GSFC), May 2022.
!
! Functions in this module:
!	accountforAOD_to_get_2Drad
! 	retreive_alhssa_hgtsol
! 	nearuv_alhssa
!	nearuv_getssa_surfacealh
!
!==============================================================================
 IMPLICIT NONE


!==============================================================================
! Nodes of w0 at 354, 388, 500
!==============================================================================
  REAL(KIND=4), DIMENSION(7)  :: w0sel, w0388sel, w0354sel
  
  REAL(KIND=4), DIMENSION(21), PARAMETER :: w0354model = &
                (/0.7498, 0.8074, 0.8472, 0.8806, 0.9321, 0.9622, 1.0000, & 
                  0.7806, 0.8082, 0.8549, 0.8879, 0.9435, 0.9696, 1.0000, & 
                  0.8002, 0.8261, 0.8540, 0.8841, 0.9172, 0.9543, 1.0000/)  

  REAL(KIND=4), DIMENSION(21), PARAMETER :: w0388model = &
  		(/0.7792, 0.8377, 0.8760, 0.9053, 0.9488, 0.9723, 1.0000, &
                  0.7806, 0.8082, 0.8549, 0.8879, 0.9435, 0.9696, 1.0000, &
                  0.8186, 0.8425, 0.8680, 0.8954, 0.9254, 0.9591, 1.0000/)

  REAL(KIND=4), DIMENSION(21), PARAMETER :: w0model = &
             (/0.8627, 0.9105, 0.9364, 0.9543, 0.9780, 0.9862, 1.0000, & ! Dust
               0.8265, 0.8486, 0.8785, 0.9117, 0.9603, 0.9789, 1.0000, & ! Smoke
               0.8531, 0.8723, 0.8926, 0.9182, 0.9422, 0.9689, 1.0000/)  ! Sulfate

                  
!==============================================================================
! Nodes of KrefractiveIndex at 354, 388, 500
!==============================================================================
  REAL(KIND=4), DIMENSION(7)  :: Kr500sel, Kr388sel, Kr354sel

  REAL(KIND=4), DIMENSION(21), PARAMETER :: Kr354model = &                  
  		 (/0.02303, 0.01279, 0.00832, 0.00561, 0.00256, 0.00128, 0.00000, &
                   0.05760, 0.04800, 0.03600, 0.02400, 0.01200, 0.00600, 0.00000, &
                   0.01200, 0.01000, 0.00800, 0.00600, 0.00400, 0.00200, 0.00000/)

  REAL(KIND=4), DIMENSION(21), PARAMETER :: Kr388model = &
  		(/0.01662, 0.00923, 0.00600, 0.00405, 0.00185, 0.00092, 0.00000, &
                  0.04800, 0.04000, 0.03000, 0.02000, 0.01000, 0.00500, 0.00000, &
                  0.01200, 0.01000, 0.00800, 0.00600, 0.00400, 0.00200, 0.00000/)

  REAL(KIND=4), DIMENSION(21), PARAMETER :: Kr500model = &
  	         (/0.00720, 0.00400, 0.00260, 0.00176, 0.00080, 0.00040, 0.0000, &
                   0.02880, 0.02400, 0.01800, 0.01200, 0.00600, 0.00300, 0.0000, &
                   0.01850, 0.01550, 0.01250, 0.00900, 0.00600, 0.00300, 0.0000/)


!==============================================================================
! Ratio of optical depth *** tauRatiomodel ***
!==============================================================================
  REAL(KIND=4)  :: tauRatiosel(7), tauRatiomodel(21)
  
  REAL(KIND=4), DIMENSION(7), PARAMETER :: tau_table = (/0.00, 0.10, 0.50, 1.00, 2.50, 4.00, 6.00/)
  
  REAL(KIND=4), DIMENSION(21), PARAMETER :: tauRatio388_500 = &
                     (/1.166, 1.165, 1.165, 1.165, 1.165, 1.165, 1.165, & !Dust
                       1.501, 1.517, 1.538, 1.560, 1.483, 1.492, 1.502, & !Smoke
                       1.507, 1.576, 1.582, 1.588, 1.595, 1.601, 1.608/)   !Sulphate

  REAL(KIND=4), DIMENSION(21), PARAMETER :: tauRatio354_500 = &		       
  		     (/1.233, 1.234, 1.234, 1.234, 1.234, 1.234, 1.234, & 
                       1.702, 1.725, 1.756, 1.788, 1.664, 1.677, 1.690, &
                       1.813, 1.822, 1.831, 1.840, 1.850, 1.859, 1.870/) 


!==============================================================================
! Height nodal points
!==============================================================================
  REAL(KIND=4), DIMENSION(5), PARAMETER :: hgt_table = &
  				(/0.000, 1.500, 3.000, 6.000, 10.000/)
                    
!==============================================================================
! Aerosol types: indexds -- dust (ainuv: +, aivis: +) 
!                indexcs -- carbonaeous(smoke) (ainuv: +, aivis: -) 
!               indexind -- industrial pollutants (ainuv: -, aivis: -) 
!==============================================================================
  INTEGER(KIND=4), DIMENSION(7), PARAMETER :: indexds = (/1,2,3,4,5,6,7/)
  INTEGER(KIND=4), DIMENSION(7), PARAMETER :: indexcs = (/8,9,10,11,12,13,14/)
  INTEGER(KIND=4), DIMENSION(7), PARAMETER :: indexind = (/15,16,17,18,19,20,21/)

  REAL(KIND=4), PARAMETER :: UVAIThreshold = 1.0
  REAL(KIND=4), DIMENSION(3), PARAMETER :: VisAIThreshold = (/2.5,0.5,1.5/)


 CONTAINS

!==============================================================================
!==============================================================================

FUNCTION accountforAOD_to_get_2Drad(atype,rad3d_wv1,rad3d_wv2,UVAOD,&
                                    fmf,rad2d_wv1,rad2d_wv2) RESULT(status)

USE InterpolationModule

IMPLICIT NONE

  INTEGER(KIND=2),                INTENT(IN)  :: AType
  REAL(KIND=4), DIMENSION(2),     INTENT(IN)  :: UVAOD
  REAL(KIND=4),                   INTENT(IN)  :: fmf
  REAL(KIND=4), DIMENSION(5,7,7), INTENT(IN)  :: rad3d_wv1,rad3d_wv2
  REAL(KIND=4), DIMENSION(5,7)  , INTENT(OUT) :: rad2d_wv1,rad2d_wv2
!
  INTEGER                    :: status, i, j
  REAL(KIND=4), DIMENSION(7) :: rad354_allaodnodes, rad388_allaodnodes
  REAL(KIND=4), DIMENSION(7) :: tau354_table, tau388_table
  REAL(KIND=4), DIMENSION(7) :: tmpRatio354_500, tmpRatio388_500
  REAL(KIND=4)               :: a1, a2, frac, tmp1, tmp2
  INTEGER(KIND=4)            :: ja1, ja2


IF (atype .EQ. 1) THEN 
  tau354_table = tau_table * tauRatio354_500(8:14)
  tau388_table = tau_table * tauRatio388_500(8:14)
ELSE IF (atype .EQ. 2) THEN 
  tau354_table = tau_table * tauRatio354_500(1:7)
  tau388_table = tau_table * tauRatio388_500(1:7)
ELSE IF (atype .EQ. 3) THEN 
  tau354_table = tau_table * tauRatio354_500(15:21)
  tau388_table = tau_table * tauRatio388_500(15:21)
ELSE IF (atype .EQ. 12) THEN 
  tmpRatio354_500 = (tauRatio354_500(1:7)*(1.-fmf)) + (tauRatio354_500(8:14)*fmf)
  tmpRatio388_500 = (tauRatio388_500(1:7)*(1.-fmf)) + (tauRatio388_500(8:14)*fmf)
  tau354_table = tau_table * tmpRatio354_500
  tau388_table = tau_table * tmpRatio388_500
ENDIF

rad2d_wv1(:,:) = 0.0
rad2d_wv2(:,:) = 0.0

DO i = 1, 5 ! hgt-loop 
  DO j = 1, 7 ! ssa-loop
  rad354_allaodnodes = rad3d_wv1(i,j,:)
  rad388_allaodnodes = rad3d_wv2(i,j,:)

    IF ( UVAOD(1) .GE. minval(tau354_table) .AND. & 
         UVAOD(1) .LE. maxval(tau354_table) ) THEN
    status = FindTableEntry(UVAOD(1),tau354_table,7,a1,a2,ja1,ja2,frac)
    status = Interp1D(rad354_allaodnodes(ja1), rad354_allaodnodes(ja2), frac, tmp1)
    rad2d_wv1(i,j) = tmp1
    ENDIF

    IF ( UVAOD(2) .GE. minval(tau388_table) .AND. & 
         UVAOD(2) .LE. maxval(tau388_table) ) THEN
    status = FindTableEntry(UVAOD(2),tau388_table,7,a1,a2,ja1,ja2,frac)
    status = Interp1D(rad388_allaodnodes(ja1), rad388_allaodnodes(ja2), frac, tmp2)
    rad2d_wv2(i,j) = tmp2
    ENDIF

  ENDDO ! ssa-loop

ENDDO ! hgt-loop


RETURN

END FUNCTION accountforAOD_to_get_2Drad

!==============================================================================
!==============================================================================

FUNCTION retrieve_alhssa_hgtsol(aodvshgt,ssavshgt,aod388,retssa,rethgt)RESULT(status)

USE InterpolationModule

IMPLICIT NONE

  REAL(KIND=4), DIMENSION(5), INTENT(IN)  :: aodvshgt, ssavshgt
  REAL(KIND=4),               INTENT(IN)  :: aod388
  REAL(KIND=4),               INTENT(OUT) :: retssa, rethgt
!
  REAL, DIMENSION(5)   :: tmpaodvshgt
  INTEGER              :: status
  REAL(KIND=4)         :: a1, a2, fracaod, minaod
  INTEGER(KIND=4)      :: ja1, ja2, nhgt
!
  rethgt = -9999.
  retssa = -9999.
  nhgt = 5

  IF ( aod388 .GE. minval(aodvshgt) .AND. &
       aod388 .LE. maxval(aodvshgt) .AND. minval(aodvshgt) .GT. 0 ) THEN
    tmpaodvshgt = aodvshgt
    ! Since AOD decreases with height use (1-frac) 
    ! IMPORTANT NOTE: FindTableEntry changes the actual data if it is NOT ascending.
    status = FindTableEntry(aod388,tmpaodvshgt,nhgt,a1,a2,ja1,ja2,fracaod)
    status = Interp1D(hgt_table(ja1), hgt_table(ja2), (1.-fracaod), rethgt)
    status = Interp1D(ssavshgt(ja1), ssavshgt(ja2), (1.-fracaod), retssa)
    !
  ENDIF

status = 1

RETURN

END FUNCTION retrieve_alhssa_hgtsol

!==============================================================================
!==============================================================================

 FUNCTION NEARUV_ALHSSA(AType,rad1,rad2,rad1_obs,rad2_obs, &
                        fmf, retssa2,rethgt2) RESULT(status)

!-------------------------------------------------------------------------------!
! TITLE      : NEARUV_ALHSSA							!
!										!
! AUTHOR     : Vinay Kayetha (SSAI/NASA-GSFC), May 2022.			!
!										!
! DESCRIPTION: Uses AOD to retrieve ALH and SSA388 simutaneously.		!
!										!
! INPUTS     :	W0 	Selected aerosol model ssa nodes			!
!		RAD 	Model Norm.Radiances at UV1, UV2 (354 and 388 nm)	!
!			2D array as function of (nhgt, nssa)			!
!		RAD_OBS	Measured Norm.Radiances at UV1, UV2 (354 and 388 nm)	!
!										!
! OUTPUTS    : 	RETSSA	Retrieved SSA at 388 nm					!
!		RETHGT	Retrieved HGT.						!
!-------------------------------------------------------------------------------!

  USE InterpolationModule

  IMPLICIT NONE

  INTEGER(KIND=2),              INTENT(IN)  :: AType
  REAL(KIND=4), INTENT(IN),  DIMENSION(:,:) :: rad1, rad2
  REAL(KIND=4), INTENT(IN)                  :: rad1_obs, rad2_obs
  REAL(KIND=4), INTENT(IN)                  :: fmf
!
  REAL(KIND=4), INTENT(OUT)                 :: retssa2(5), rethgt2
!
  REAL(KIND=4)                              :: ratio_obs
  REAL(KIND=4), DIMENSION(7)                :: hgtloop, radw2_hgt
  REAL(KIND=4), DIMENSION(5)                :: ratio_mod1
  REAL(KIND=4), DIMENSION(7)                :: w0388, w0388_f, w0388_c
!
  REAL(KIND=4)    :: q2, q1, q_tmp, ssa500_m2, ssa388_m2, ssa354_m2, hgt, rad2wi
  INTEGER(KIND=4) :: idx2, idx1, issa, newssa_id, xx
!
  INTEGER :: STATUS
  REAL(KIND=4) :: tmphgt(5), tmpssa(7), tmpssa388(7), tmpssa354(7)
  REAL(KIND=4), DIMENSION(5) :: coef_ssa_wl2, coef_ssa_wl1


 IF (atype .EQ. 1 ) THEN ! This is smoke
	w0388(1:7) = w0388model(8:14)
  ELSEIF (atype .EQ. 2 ) THEN ! This is dust
        w0388(1:7) = w0388model(1:7)
  ELSEIF (atype .EQ. 3 ) THEN  ! This is urban/industrial
        w0388(1:7) = w0388model(15:21)
  ELSEIF (atype .EQ. 12 ) THEN  ! Mixture/Smoke + Dust
        w0388_c(1:7) = w0388model(1:7)
        w0388_f(1:7) = w0388model(8:14)
        w0388(1:7) = (w0388_f(1:7) * fmf) + (w0388_c(1:7) * (1.0-fmf))
  ELSE 
     PRINT *, 'Invalid Aerosol Type or Mixing'
     CALL EXIT(1)
  ENDIF

! Calculate the ratio of radiances at 354 nm and 388 nm
  ratio_obs = rad1_obs/rad2_obs
  hgtloop(:) = -9999.
  radw2_hgt(:) = -9999.
  newssa_id = 0
  ssa388_m2 = -9999. 

! Print *, 'NUV_alh ATYPE = ',atype
! Print *, 'SSA Loop = ', atype, w0388(1:7), w0388model(1:21)

! METHOD-2 start with SSA-loop
! Begin loop over all SSA nodes
  DO issa = 1, 7
    
    !Print *, 'SSA loop ', issa, ratio_obs, minval(ratio_mod1), maxval(ratio_mod1), ratio_mod1(:)
    
     ratio_mod1(:) = rad1(:,issa)/rad2(:,issa)
     IF ( ratio_obs .GE. minval(ratio_mod1) .AND. &
          ratio_obs .LE. maxval(ratio_mod1) ) THEN
  
       ! ratio_mod and rad2 decreases with increasing height.
       ! CAUTION: FindTableEntry - changes the input vector to ascending order. 
       ! For descending order, use (1-frac) with the function Interp1D.
       status = FindTableEntry(ratio_obs,ratio_mod1,5,q1,q2,idx1,idx2,q_tmp)
       status = Interp1D(hgt_table(idx1), hgt_table(idx2), (1.-q_tmp), hgt)
       !PRINT *, q2, q1, idx2, idx1, (1.-q_tmp), hgt
       status = FindTableEntry(hgt,hgt_table,5,q1,q2,idx1,idx2,q_tmp)
       status = Interp1D(rad2(idx1,issa), rad2(idx2,issa), q_tmp, rad2wi)

       IF (hgt .GE. minval(hgt_table) .AND. hgt .LE. maxval(hgt_table) .AND. rad2wi .GT. 0) THEN
         newssa_id = newssa_id + 1
         tmpssa388(newssa_id) = w0388(issa)
         hgtloop(newssa_id) = hgt
         radw2_hgt(newssa_id) = rad2wi
       ENDIF
     ENDIF

  ENDDO ! DO issa = 1, 7
! End loop over SSA nodes     

  IF(newssa_id.GT.1 .AND. rad2_obs .GE. minval(radw2_hgt(1:newssa_id)) .AND. &
                          rad2_obs .LE. maxval(radw2_hgt(1:newssa_id)) )THEN
     status = FindTableEntry(rad2_obs,radw2_hgt(1:newssa_id),newssa_id,q1,q2,idx1,idx2,q_tmp)
     IF ( (q_tmp .GE. 0 .OR. q_tmp .LE. 1) ) THEN
       status = Interp1D(tmpssa388(idx1), tmpssa388(idx2), q_tmp, ssa388_m2)
       status = Interp1D(hgtloop(idx1), hgtloop(idx2), q_tmp, rethgt2)
       retssa2(2) = ssa388_m2
     ELSE 
       retssa2(:) = -9999.
       rethgt2 = -9999.
     ENDIF
  ELSE
      retssa2(:) = -9999.
      rethgt2 = -9999.
  ENDIF ! IF(newhgt_id.GT.1)THEN

  STATUS = 1

  RETURN

 END FUNCTION NEARUV_ALHSSA
 
!==============================================================================
!==============================================================================

 FUNCTION nearuv_getssa_surfacealh(rad1,rad2,rad1_obs,rad2_obs, &
                        retssa2,rethgt2) RESULT(status)

  USE InterpolationModule

  IMPLICIT NONE

  REAL(KIND=4), INTENT(IN),  DIMENSION(:,:) :: rad1, rad2
  REAL(KIND=4), INTENT(IN)                  :: rad1_obs, rad2_obs
!
  REAL(KIND=4), INTENT(OUT)                 :: retssa2(5), rethgt2
!
  REAL(KIND=4)                              :: ratio_obs, ssa388_m2
  REAL(KIND=4), DIMENSION(7)                :: rad2_mod, w0388

!
  REAL(KIND=4)    :: q2, q1, q_tmp
  INTEGER(KIND=4) :: idx2, idx1
!
  INTEGER :: STATUS

! This will always be Urban
  w0388 = w0388model(15:21)
  
! Retrieve SSA assuming ALH=0
  ratio_obs = rad1_obs/rad2_obs
  rad2_mod(:) = rad2(1,:)
  
  retssa2(2) = -9999.
  rethgt2 = -9999.
     
  IF ( rad2_obs .GE. minval(rad2_mod) .AND. &
       rad2_obs .LE. maxval(rad2_mod) ) THEN
  
    status = FindTableEntry(rad2_obs,rad2_mod,7,q1,q2,idx1,idx2,q_tmp)
    status = Interp1D(w0388(idx1), w0388(idx2), q_tmp, ssa388_m2)
 
    IF ( ssa388_m2 .GE. minval(w0388) .AND. &
         ssa388_m2 .LE. maxval(w0388) ) THEN
       retssa2(2) = ssa388_m2
       rethgt2 = 0.
    ENDIF
    
  ENDIF
  
  STATUS = 1

  RETURN

 END FUNCTION nearuv_getssa_surfacealh
 
!==============================================================================
!==============================================================================


 FUNCTION NEARUV_ALHSSA_ALHLOOP(AType,rad1,rad2,rad1_obs,rad2_obs, &
                        fmf, retssa2,rethgt2) RESULT(status)

!-------------------------------------------------------------------------------!
! TITLE      : NEARUV_ALHSSA							!
!										!
! AUTHOR     : Vinay Kayetha (SSAI/NASA-GSFC), Dec 2023.			!
!										!
! DESCRIPTION: Uses AOD to retrieve ALH and SSA388 simutaneously.		!
!										!
! INPUTS     :	W0 	Selected aerosol model ssa nodes			!
!		RAD 	Model Norm.Radiances at UV1, UV2 (354 and 388 nm)	!
!			2D array as function of (nhgt, nssa)			!
!		RAD_OBS	Measured Norm.Radiances at UV1, UV2 (354 and 388 nm)	!
!										!
! OUTPUTS    : 	RETSSA	Retrieved SSA at 388 nm					!
!		RETHGT	Retrieved HGT.						!
!-------------------------------------------------------------------------------!

  USE InterpolationModule

  IMPLICIT NONE

  INTEGER(KIND=2),              INTENT(IN)  :: AType
  REAL(KIND=4), INTENT(IN),  DIMENSION(:,:) :: rad1, rad2
  REAL(KIND=4), INTENT(IN)                  :: rad1_obs, rad2_obs
  REAL(KIND=4), INTENT(IN)                  :: fmf
!
  REAL(KIND=4), INTENT(OUT)                 :: retssa2(5), rethgt2
!
  REAL(KIND=4)                              :: ratio_obs
  REAL(KIND=4), DIMENSION(7)                :: hgtloop, radw2_hgt
  REAL(KIND=4), DIMENSION(7)                :: ratio_mod1
  REAL(KIND=4), DIMENSION(7)                :: w0500, w0388, w0354
!
  REAL(KIND=4)    :: q2, q1, q_tmp, ssa500_m2, ssa388_m2, ssa354_m2, tmp_w0, rad2wi
  INTEGER(KIND=4) :: idx2, idx1, issa, newssa_id, xx, ihgt
!
  INTEGER :: STATUS
  REAL(KIND=4) :: tmphgt(5), tmpssa(7), tmpssa388(7), tmpssa354(7)
  REAL(KIND=4), DIMENSION(5) :: coef_ssa_wl2, coef_ssa_wl1

!  Print *, 'NUV_alh ATYPE = ', atype
!  Print *, 'ALHLoop subroutine w0388model = ', w0388model
      
 IF (atype .EQ. 1 ) THEN ! This is smoke
     DO xx = 1, 7 
      w0388(xx) = w0388model(xx + 7)
     ENDDO
  ELSEIF (atype .EQ. 2 ) THEN ! This is dust
     DO xx = 1, 7 
      w0388(xx) = w0388model(xx + 0)
     ENDDO
  ELSEIF (atype .EQ. 3 ) THEN  ! This would be industrial
     DO xx = 1, 7 
      w0388(xx) = w0388model(xx + 14)
     ENDDO 
  ELSEIF (atype .EQ. 12 ) THEN 
     DO xx = 1, 7 
     !Print*, 'ALH Loop = ', xx, fmf, w0388model(xx + 7), w0388model(xx + 0)
      w0388(xx) = (w0388model(xx + 7) * fmf) + (w0388model(xx + 0) * (1.-fmf))
     ENDDO
  ELSE 
     PRINT *, 'Invalid Aerosol Type or Mixing'
     CALL EXIT(1)
  ENDIF

! Calculate the ratio of radiances at 354 nm and 388 nm
  ratio_obs = rad1_obs/rad2_obs
  hgtloop(:) = -9999.
  radw2_hgt(:) = -9999.
  newssa_id = 0
  ssa388_m2 = -9999. 

! METHOD-1 start with HGT-loop
! Begin loop over all HGT nodes
  DO ihgt = 1, 5
     
     ratio_mod1(:) = rad1(ihgt, :)/rad2(ihgt, :)
     IF ( ratio_obs .GE. minval(ratio_mod1) .AND. &
          ratio_obs .LE. maxval(ratio_mod1) ) THEN
  
       ! ratio_mod and rad2 decreases with increasing height.
       ! CAUTION: FindTableEntry - changes the input vector to ascending order. 
       ! For descending order, use (1-frac) with the function Interp1D.
             
       status = FindTableEntry(ratio_obs,ratio_mod1,7,q1,q2,idx1,idx2,q_tmp)
!Print *, 'HGT loop ', ihgt, ratio_obs, minval(ratio_mod1), maxval(ratio_mod1), ratio_mod1(:), idx1, idx2, q1, q2 
       status = Interp1D(w0388(idx1), w0388(idx2), q_tmp, tmp_w0)
       
       status = FindTableEntry(tmp_w0, w0388, 7, q1,q2,idx1,idx2,q_tmp)
       status = Interp1D(rad2(ihgt,idx1), rad2(ihgt,idx2), q_tmp, rad2wi)

       IF (tmp_w0 .GE. minval(w0388) .AND. tmp_w0 .LE. maxval(w0388) .AND. rad2wi .GT. 0) THEN
         newssa_id = newssa_id + 1
         tmpssa388(newssa_id) = tmp_w0
         hgtloop(newssa_id) = hgt_table(ihgt)
         radw2_hgt(newssa_id) = rad2wi
       ENDIF
     ENDIF

  ENDDO ! DO ihgt = 1, 5
! End loop over HGT nodes   


  IF(newssa_id.GT.1 .AND. rad2_obs .GE. minval(radw2_hgt(1:newssa_id)) .AND. &
                          rad2_obs .LE. maxval(radw2_hgt(1:newssa_id)) )THEN
     status = FindTableEntry(rad2_obs,radw2_hgt(1:newssa_id),newssa_id,q1,q2,idx1,idx2,q_tmp)
     IF ( (q_tmp .GE. 0 .OR. q_tmp .LE. 1) ) THEN
       status = Interp1D(tmpssa388(idx1), tmpssa388(idx2), q_tmp, ssa388_m2)
       status = Interp1D(hgtloop(idx1), hgtloop(idx2), q_tmp, rethgt2)
       retssa2(2) = ssa388_m2
     ELSE 
       retssa2(:) = -9999.
       rethgt2 = -9999.
     ENDIF
  ELSE
      retssa2(:) = -9999.
      rethgt2 = -9999.
  ENDIF ! IF(newhgt_id.GT.1)THEN

  STATUS = 1

  RETURN

 END FUNCTION NEARUV_ALHSSA_ALHLOOP
 
!==============================================================================
!==============================================================================

END MODULE Nearuv_alhssa_Module

