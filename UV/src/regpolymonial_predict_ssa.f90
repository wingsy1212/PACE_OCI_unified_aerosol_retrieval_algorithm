MODULE regpolymonial_predict_ssa

!==============================================================================
!
! FILENAME: regpolymonial_predict_ssa.f90
!
! AUTHOR: Vinay Kayetha (SSAI/NASA-GSFC), Nov 2023.
!
! Functions in this module:
!	predict_wavssa_from_ssa388
!
!==============================================================================

 INTEGER(KIND=4), PARAMETER              :: n_seas = 4
 INTEGER(KIND=4), PARAMETER              :: n_coef = 3 

 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: default_smk 
 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: default_dst
 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: default_urb 

 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: reg1_smk  ! Sahel
 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: reg2_smk  ! South_Africa
 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: reg3_smk  ! Middle_East
 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: reg4_smk  ! E_China
 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: reg5_smk  ! N_India
 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: reg6_smk  ! Australia

 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: reg1_dst  ! Middle_East
 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: reg2_dst  ! E_China
 REAL(KIND=4), DIMENSION(n_seas, n_coef) :: reg3_dst  ! N_India

! Default models     
 DATA default_smk(1, 1:3) /0.6454, 0.0011, -1.00E-06/
 DATA default_smk(2, 1:3) /0.7261, 0.0008, -7.00E-07/  ! Verified 
 DATA default_smk(3, 1:3) /0.6438, 0.0012, -1.00E-06/
 DATA default_smk(4, 1:3) /0.6669, 0.0011, -1.00E-06/

 DATA default_dst(1, 1:3) /0.7575, 0.0006, -4.00E-07/
 DATA default_dst(2, 1:3) /0.6686, 0.0009, -7.00E-07/
 DATA default_dst(3, 1:3) /0.6110, 0.0011, -8.00E-07/
 DATA default_dst(4, 1:3) /0.6837, 0.0009, -7.00E-07/

 DATA default_urb(1, 1:3) /0.8067, 0.0004, -4.11E-07/
 DATA default_urb(2, 1:3) /0.7947, 0.0006, -7.09E-07/
 DATA default_urb(3, 1:3) /0.7898, 0.0006, -7.36E-07/
 DATA default_urb(4, 1:3) /0.7281, 0.0007, -8.14E-07/

! Regional Smoke Models
 ! Sah & Sahel [0, -30, 40, 25]
 DATA reg1_smk(1, 1:3) /0.7980, 0.0003, -3.00E-07/
 DATA reg1_smk(2, 1:3) /0.7261, 0.0008, -7.00E-07/
 DATA reg1_smk(3, 1:3) /0.6952, 0.0009, -1.00E-06/
 DATA reg1_smk(4, 1:3) /0.6669, 0.0011, -1.00E-06/

 ! ! S.Africa [-40, 0, 0, 100]
 DATA reg2_smk(1, 1:3) /0.6454, 0.0011, -1.00E-06/
 DATA reg2_smk(2, 1:3) /0.7261, 0.0008, -7.00E-07/
 DATA reg2_smk(3, 1:3) /0.5151, 0.0015, -2.00E-06/
 DATA reg2_smk(4, 1:3) /0.4490, 0.0018, -2.00E-06/

 ! Middle_East [0, 25, 40, 60]
 DATA reg3_smk(1, 1:3) /0.8241, 0.0003, -3.00E-07/
 DATA reg3_smk(2, 1:3) /0.5893, 0.0014, -1.00E-06/
 DATA reg3_smk(3, 1:3) /0.6868, 0.0010, -1.00E-06/
 DATA reg3_smk(4, 1:3) /0.5930, 0.0013, -1.00E-06/

 ! E_China [35, 70, 50, 150]
 DATA reg4_smk(1, 1:3) /0.5989, 0.0012, -1.00E-06/
 DATA reg4_smk(2, 1:3) /0.6491, 0.0011, -1.00E-06/
 DATA reg4_smk(3, 1:3) /0.6561, 0.0012, -1.00E-06/
 DATA reg4_smk(4, 1:3) /0.5408, 0.0015, -1.00E-06/

 ! N_India[0, 60, 35, 100]
 DATA reg5_smk(1, 1:3) /0.5355, 0.0016, -2.00E-06/
 DATA reg5_smk(2, 1:3) /0.5306, 0.0015, -1.00E-06/
 DATA reg5_smk(3, 1:3) /0.4400, 0.0020, -2.00E-06/
 DATA reg5_smk(4, 1:3) /0.5621, 0.0015, -2.00E-06/

 ! Australia [-90, 100, 0, 180]
 DATA reg6_smk(1, 1:3) /0.6454, 0.0011, -1.00E-06/
 DATA reg6_smk(2, 1:3) /0.7261, 0.0008, -7.00E-07/
 DATA reg6_smk(3, 1:3) /0.6438, 0.0012, -1.00E-06/
 DATA reg6_smk(4, 1:3) /0.7703, 0.0005, -6.00E-07/


! Regional Dust Models
 ! Middle_East [0, 25, 40, 60]
 DATA reg1_dst(1, 1:3) /0.6594, 0.0011, -9.00E-07/
 DATA reg1_dst(2, 1:3) /0.4818, 0.0017, -1.00E-06/
 DATA reg1_dst(3, 1:3) /0.5782, 0.0013, -1.00E-06/
 DATA reg1_dst(4, 1:3) /0.6286, 0.0012, -9.00E-07/

 ! E_China [35, 70, 50, 150]
 DATA reg2_dst(1, 1:3) /0.7575, 0.0006, -4.00E-07/
 DATA reg2_dst(2, 1:3) /0.5977, 0.0012, -1.00E-06/
 DATA reg2_dst(3, 1:3) /0.6016, 0.0011, -9.00E-07/
 DATA reg2_dst(4, 1:3) /0.6837, 0.0009, -7.00E-07/

 ! N_India[0, 60, 35, 100]
 DATA reg3_dst(1, 1:3) /0.7575, 0.0006, -4.00E-07/
 DATA reg3_dst(2, 1:3) /0.4832, 0.0016, -1.00E-06/
 DATA reg3_dst(3, 1:3) /0.5072, 0.0016, -1.00E-06/
 DATA reg3_dst(4, 1:3) /0.6837, 0.0009, -7.00E-07/


CONTAINS

!=====================================================================
!=====================================================================

SUBROUTINE predict_wavssa_from_ssa388(plat, plon, mon, atype, ssa388, inwave, wavssa)
 
IMPLICIT NONE

 INTEGER(KIND=2),   INTENT(IN)    :: atype
 INTEGER(KIND=4),   INTENT(IN)    :: mon
 REAL(KIND=4),      INTENT(IN)    :: plat, plon
 REAL(KIND=4),      INTENT(IN)    :: ssa388, inwave
 REAL(KIND=4),      INTENT(OUT)   :: wavssa

 INTEGER(KIND=4):: seas_id
 REAL(KIND=4)   :: tmp_wavssa, tmp_ssa388, offset
 REAL(KIND=4)   :: polycoeff(3), seascoeff(n_seas, n_coef)

  IF (mon .EQ. 12 .OR. mon .EQ.  1 .OR. mon .EQ.  2) seas_id = 1
  IF (mon .EQ.  3 .OR. mon .EQ.  4 .OR. mon .EQ.  5) seas_id = 2
  IF (mon .EQ.  6 .OR. mon .EQ.  7 .OR. mon .EQ.  8) seas_id = 3
  IF (mon .EQ.  9 .OR. mon .EQ. 10 .OR. mon .EQ. 11) seas_id = 4
  
  
  SELECT CASE (atype)
      CASE(1) ! Smoke
	IF ( plat .GT. 0.00 .AND. plat .LE. 40.00 .AND. &
	     plon .GT. -30.00 .AND. plon .LE. 25.00 ) THEN  
	     seascoeff(:,:) = reg1_smk
	ELSE IF ( plat .GT. -40.00 .AND. plat .LE. 0.00 .AND. &
	          plon .GT. 0.00 .AND. plon .LE. 100.00 ) THEN 
	     seascoeff(:,:) = reg2_smk		  
	ELSE IF ( plat .GT. 0.00 .AND. plat .LE. 40.00 .AND. &
	          plon .GT. 25.00 .AND. plon .LE. 60.00 ) THEN 
	     seascoeff(:,:) = reg3_smk		  
	ELSE IF ( plat .GT. 35.00 .AND. plat .LE. 50.00 .AND. &
	          plon .GT. 70.00 .AND. plon .LE. 150.00 ) THEN 
	     seascoeff(:,:) = reg4_smk		  
	ELSE IF ( plat .GT.  0.00 .AND. plat .LE. 35.00 .AND. &
	          plon .GT. 60.00 .AND. plon .LE. 100.00 ) THEN  
	     seascoeff(:,:) = reg5_smk
	ELSE IF ( plat .GT. -90.00 .AND. plat .LE. 0.00 .AND. &
	          plon .GT. 110.00 .AND. plon .LE. 180.00 ) THEN 
	     seascoeff(:,:) = reg6_smk
	ELSE
	     seascoeff(:,:) = default_smk	
	ENDIF	

      CASE(2, 4) ! Dust
	IF ( plat .GT. 0.00 .AND. plat .LE. 40.00 .AND. &
	          plon .GT. 25.00 .AND. plon .LE. 60.00 ) THEN ! Middle_East 
	     seascoeff(:,:) = reg1_dst		  
	ELSE IF ( plat .GT. 35.00 .AND. plat .LE. 50.00 .AND. &
	          plon .GT. 70.00 .AND. plon .LE. 150.00 ) THEN ! E_China 
	     seascoeff(:,:) = reg2_dst		  
	ELSE IF ( plat .GT.  0.00 .AND. plat .LE. 35.00 .AND. &
	          plon .GT. 60.00 .AND. plon .LE. 100.00 ) THEN  ! N_India
	     seascoeff(:,:) = reg3_dst
	ELSE
	     seascoeff(:,:) = default_dst	
	ENDIF	
      
      CASE(3) ! Urban	
	     seascoeff(:,:) = default_urb               
  END SELECT
 
 polycoeff(1:3) = seascoeff(seas_id, :)   
 tmp_wavssa = polycoeff(1) + inwave*polycoeff(2) + inwave**2 *polycoeff(3)
 tmp_ssa388 = polycoeff(1) +  388.0*polycoeff(2) +  388.0**2 *polycoeff(3)
 offset = ssa388 - tmp_ssa388
 wavssa = tmp_wavssa + offset
 if (wavssa .GT. 1.0) wavssa=0.999999
 
END SUBROUTINE predict_wavssa_from_ssa388

!=====================================================================
!=====================================================================

END MODULE regpolymonial_predict_ssa
