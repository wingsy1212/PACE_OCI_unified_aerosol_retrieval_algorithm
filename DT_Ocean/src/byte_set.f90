!****************************************************************
      SUBROUTINE BYTE_SET(QA_V,Bit_SP,QA_Byte)
      IMPLICIT NONE
      SAVE
!-------------------------------------------------------------------
!
! !F77
!
! !DESCRIPTION:
!                 This subroutine is to set QA bit given nth byte
!                 Note: the QA byte setting starts from the rightmost
!                       (not the leftmost) bit
!
!        1st word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
!        2nd word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
!        3rd word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
!        4th word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
!        5th word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
!
! !INPUT PARAMETERS:
!
!        QA_V       QA parameter value
!        Bit_SP     Bit starting position of Ith QA parameter
!                    (see MODIS atmosphere QA plan)
!
! !OUTPUT PARAMETERS:
!
!        QA_Byte    Byte set for quality control
 
      Intrinsic ibset,ibclr
      INTEGER QA_V,Bit_SP,Byte_Temp
      BYTE QA_Byte
      Byte_Temp=QA_Byte

      IF(QA_V.EQ.0) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)

      ELSE IF(QA_V.EQ.1) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)

      ELSE IF(QA_V.EQ.2) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)

      ELSE IF(QA_V.EQ.3) THEN

        Byte_temp = ibset(Byte_Temp,Bit_SP)
        Byte_temp = ibset(Byte_Temp,Bit_SP+1)

      ELSE IF(QA_V.EQ.4) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)

      ELSE IF(QA_V.EQ.5) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)

      ELSE IF(QA_V.EQ.6) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)

      ELSE IF(QA_V.EQ.7) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)

      ELSE IF(QA_V.EQ.8) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibclr(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.9) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibclr(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.10) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibclr(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.11) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibclr(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.12) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.13) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.14) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.15) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ENDIF

      QA_Byte=Byte_Temp

      RETURN
      END


