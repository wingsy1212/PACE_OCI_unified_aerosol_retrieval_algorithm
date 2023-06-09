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
!	predict_ssa_frompoly
!	nearuv_getssa_surfacealh
!
!==============================================================================
 IMPLICIT NONE

 CONTAINS

!==============================================================================
!==============================================================================

FUNCTION accountforAOD_to_get_2Drad(atype,indexmd,rad3d_wv1,rad3d_wv2,UVAOD,&
                                    fmf,rad2d_wv1,rad2d_wv2) RESULT(status)

USE InterpolationModule

IMPLICIT NONE

  INTEGER(KIND=2),                INTENT(IN)  :: AType
  INTEGER(KIND=4),                INTENT(IN)  :: indexmd(7)
  REAL(KIND=4), DIMENSION(2),     INTENT(IN)  :: UVAOD
  REAL(KIND=4),                   INTENT(IN)  :: fmf
  REAL(KIND=4), DIMENSION(5,7,7), INTENT(IN)  :: rad3d_wv1,rad3d_wv2
  REAL(KIND=4), DIMENSION(5,7)  , INTENT(OUT) :: rad2d_wv1,rad2d_wv2
!
  INTEGER                    :: status, i, j
  REAL(KIND=4), DIMENSION(7) :: rad354_allaodnodes, rad388_allaodnodes
  REAL(KIND=4), DIMENSION(7) :: tau354_table, tau388_table, tau_table
  REAL(KIND=4), DIMENSION(7) :: tmpRatio354_500, tmpRatio388_500
  REAL(KIND=4), DIMENSION(21):: tauRatio354_500, tauRatio388_500
  REAL(KIND=4)               :: a1, a2, frac, tmp1, tmp2
  INTEGER(KIND=4)            :: ja1, ja2

  DATA tau_table/0.00, 0.10, 0.50, 1.00, 2.50, 4.00, 6.00/
  DATA tauRatio388_500/1.166,1.166,1.166,1.166,1.166,1.166,1.166, & !Dust
                       1.501,1.517,1.538,1.560,1.483,1.492,1.502, & !Smoke
                       1.507,1.576,1.582,1.588,1.595,1.601,1.608/   !Sulphate
  DATA tauRatio354_500/1.233,1.233,1.233,1.233,1.233,1.233,1.233, & 
                       1.702,1.725,1.756,1.788,1.663,1.677,1.690, &
                       1.813,1.822,1.831,1.840,1.850,1.859,1.870/ 

IF (atype .GE. 1 .AND. atype .LE. 3 .AND. fmf .LT. 0) THEN 
  tau354_table = tau_table * tauRatio354_500(indexmd(1:7))
  tau388_table = tau_table * tauRatio388_500(indexmd(1:7))
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

    IF ( UVAOD(1) .GE. minval(tau_table) .AND. & 
         UVAOD(1) .LE. maxval(tau_table) ) THEN
    status = FindTableEntry(UVAOD(1),tau_table,7,a1,a2,ja1,ja2,frac)
    status = Interp1D(rad354_allaodnodes(ja1), rad354_allaodnodes(ja2), frac, tmp1)
    rad2d_wv1(i,j) = tmp1
    ENDIF

    IF ( UVAOD(2) .GE. minval(tau_table) .AND. & 
         UVAOD(2) .LE. maxval(tau_table) ) THEN
    status = FindTableEntry(UVAOD(2),tau_table,7,a1,a2,ja1,ja2,frac)
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
  REAL(KIND=4), DIMENSION(5) :: hgt_table
  DATA hgt_table/0.000, 1.500, 3.000, 6.000, 10.000/
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

 FUNCTION NEARUV_ALHSSA(AType,w0,rad1,rad2,rad1_obs,rad2_obs, &
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
  REAL(KIND=4), INTENT(IN),  DIMENSION(:)   :: w0
  REAL(KIND=4), INTENT(IN),  DIMENSION(:,:) :: rad1, rad2
  REAL(KIND=4), INTENT(IN)                  :: rad1_obs, rad2_obs
  REAL(KIND=4), INTENT(IN)                  :: fmf
!
  REAL(KIND=4), INTENT(OUT)                 :: retssa2(5), rethgt2
!
  REAL(KIND=4)                              :: ratio_obs
  REAL(KIND=4), DIMENSION(7)                :: hgtloop, radw2_hgt
  REAL(KIND=4), DIMENSION(5)                :: ratio_mod1
  REAL(KIND=4), DIMENSION(7)                :: w0500, w0388, w0354
  REAL(KIND=4), DIMENSION(21)               :: w0model, w0388model, w0354model
!
  REAL(KIND=4)    :: q2, q1, q_tmp, ssa500_m2, ssa388_m2, ssa354_m2, hgt, rad2wi
  INTEGER(KIND=4) :: idx2, idx1, issa, newssa_id, xx
!
  INTEGER :: STATUS
  REAL(KIND=4) :: tmphgt(5), tmpssa(7), tmpssa388(7), tmpssa354(7)
  REAL(KIND=4), DIMENSION(5) :: coef_ssa_wl2, coef_ssa_wl1

! Height nodal points
  REAL(KIND=4), DIMENSION(5) :: hgt_table
  DATA hgt_table/0.000, 1.500, 3.000, 6.000, 10.000/

  DATA w0354model/0.7498,0.8076,0.8467,0.8803,0.9358,0.9583,1.000, &
                  0.7175,0.7502,0.7984,0.8499,0.9227,0.9548,1.000, &
                  0.8000,0.8267,0.8528,0.8850,0.9167,0.9543,1.000/
  DATA w0388model/0.7792,0.8380,0.8753,0.9054,0.9516,0.9692,1.000, &
                  0.7795,0.8113,0.8509,0.8899,0.9442,0.9679,1.000, &
                  0.8185,0.8432,0.8669,0.8963,0.9249,0.9592,1.000/
  DATA w0model/0.8627, 0.9105, 0.9364, 0.9543, 0.9780, 0.9862, 1.0000, & ! Dust
               0.8265, 0.8486, 0.8785, 0.9117, 0.9603, 0.9789, 1.0000, & ! Smoke
               0.8531, 0.8723, 0.8926, 0.9182, 0.9422, 0.9689, 1.0000/   ! Sulfate

! Print *, 'NUV_alh ATYPE = ',atype
 
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
      w0388(xx) = (w0388model(xx + 7) * fmf) + (w0388model(xx + 0) * (1.-fmf))
     ENDDO
  ELSE 
     PRINT *, 'Invalid Aerosol Type or Mixing'
     CALL EXIT(1)
  ENDIF

! Calculate the ratio of radiances at 354 nm and 388 nm
  ratio_obs = rad1_obs/rad2_obs
  hgtloop(:) = -999.
  radw2_hgt(:) = -999.
  newssa_id = 0
  ssa388_m2 = -999. 


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
       retssa2(:) = -999.
       rethgt2 = -999.
     ENDIF
  ELSE
      retssa2(:) = -999.
      rethgt2 = -999.
  ENDIF ! IF(newhgt_id.GT.1)THEN

  STATUS = 1

  RETURN

 END FUNCTION NEARUV_ALHSSA
 
!==============================================================================
!==============================================================================

SUBROUTINE predict_ssa_frompoly(month, atype, retssa2)

IMPLICIT NONE

  INTEGER(KIND=4),      INTENT(IN)     :: Month
  INTEGER(KIND=2),      INTENT(IN)     :: AType
  REAL(KIND=4),         INTENT(INOUT)  :: retssa2(5)
!
  REAL(KIND=4), DIMENSION(3) :: coeffs
  REAL(KIND=4) :: ssa388, wav1, ssawav, offset, adjssawav
  
  IF (atype .EQ. 1 ) THEN ! This is smoke
	if (month .GE. 3 .AND. month .LE. 5) then ! MAM
	  coeffs = (/0.651702, 0.001097, -1.0542E-06/)
	elseif (month .GE. 6 .AND. month .LE. 8) then ! JJA
	  coeffs = (/0.626015, 0.001222, -1.2471E-06/)
	elseif (month .GE. 9 .AND. month .LE. 11) then ! SON
	  coeffs = (/0.643045, 0.001135, -1.1513E-06/)
	else ! DJF
	  coeffs = (/0.720132, 0.000731, -7.0638E-07/)
	endif

  ELSEIF (atype .EQ. 2 ) THEN ! This is dust
	if (month .GE. 3 .AND. month .LE. 5) then ! MAM
	  coeffs = (/0.591968, 0.001211, -9.7449E-07/)
	elseif (month .GE. 6 .AND. month .LE. 8) then ! JJA
	  coeffs = (/0.580360, 0.001254, -1.0101E-06/)
	elseif (month .GE. 9 .AND. month .LE. 11) then ! SON
	  coeffs = (/0.675821, 0.000924, -7.3562E-07/)
	else ! DJF
	  coeffs = (/0.776112, 0.000550, -3.8337E-07/)
	endif

  ELSEIF (atype .EQ. 3 ) THEN  ! This would be industrial
	if (month .GE. 3 .AND. month .LE. 5) then ! MAM
	  coeffs = (/0.794680, 0.000605, -7.0880E-07/)
	elseif (month .GE. 6 .AND. month .LE. 8) then ! JJA
	  coeffs = (/0.789785, 0.000612, -7.3643E-07/)
	elseif (month .GE. 9 .AND. month .LE. 11) then ! SON
	  coeffs = (/0.728115, 0.000744, -8.1387E-07/)
	else ! DJF
	  coeffs = (/0.806732, 0.000370, -4.1067E-07/)
	endif
  ELSE
  
  RETURN
  ENDIF

ssa388 = coeffs(1) + (coeffs(2) * 388.0) + (coeffs(3) * 388.0**2.0)
offset = retssa2(2) - ssa388

wav1 = 354.0
ssawav = coeffs(1) + (coeffs(2) * wav1) + (coeffs(3) * wav1**2.0)
retssa2(1) = ssawav + offset
if (retssa2(1) .gt. 1.0) retssa2(1) = 0.9999

wav1 = 480.0
ssawav = coeffs(1) + (coeffs(2) * wav1) + (coeffs(3) * wav1**2.0)
retssa2(3) = ssawav + offset
if (retssa2(3) .gt. 1.0) retssa2(3) = 0.9999

wav1 = 550.0
ssawav = coeffs(1) + (coeffs(2) * wav1) + (coeffs(3) * wav1**2.0)
retssa2(4) = ssawav + offset
if (retssa2(4) .gt. 1.0) retssa2(4) = 0.9999

wav1 = 670.0
ssawav = coeffs(1) + (coeffs(2) * wav1) + (coeffs(3) * wav1**2.0)
retssa2(5) = ssawav + offset
if (retssa2(5) .gt. 1.0) retssa2(5) = 0.9999

END SUBROUTINE predict_ssa_frompoly

!==============================================================================
!==============================================================================

 FUNCTION nearuv_getssa_surfacealh(w0388,rad1,rad2,rad1_obs,rad2_obs, &
                        retssa2,rethgt2) RESULT(status)

  USE InterpolationModule

  IMPLICIT NONE

  REAL(KIND=4), INTENT(IN)     :: w0388(7)
  REAL(KIND=4), INTENT(IN),  DIMENSION(:,:) :: rad1, rad2
  REAL(KIND=4), INTENT(IN)                  :: rad1_obs, rad2_obs
!
  REAL(KIND=4), INTENT(OUT)                 :: retssa2(5), rethgt2
!
  REAL(KIND=4)                              :: ratio_obs, ssa388_m2
  REAL(KIND=4), DIMENSION(7)                :: rad2_mod

!
  REAL(KIND=4)    :: q2, q1, q_tmp
  INTEGER(KIND=4) :: idx2, idx1
!
  INTEGER :: STATUS


! Retrieve SSA assuming ALH=0
  ratio_obs = rad1_obs/rad2_obs
  rad2_mod(:) = rad2(1,:)
  
  retssa2(2) = -999.
  rethgt2 = -999.
     
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

END MODULE Nearuv_alhssa_Module

