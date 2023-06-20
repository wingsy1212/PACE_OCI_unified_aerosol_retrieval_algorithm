
 module Linear_interpolation
      implicit none

contains
!*********************************************************************

      SUBROUTINE INTERP(M,XBAR,X,Y,Y1)

!---------------------------------------------------------------------
! !F90      
!
! !DESCRIPTION: This subroutine is a general purpose routine and
!              interpolates linearly.  Value of y1 is interpolated
!             for x1.no extrapolation is allowed.
!
! !INPUT PARAMETERS:I,M,X1,X,Y
!
! !OUTPUT PARAMETERS:Y1,LOPT
! 
!----------------------------------------------------------------------

      IMPLICIT NONE
       Save
      INTEGER IL,LL,LOPT,M
      REAL PINTEN,PPHI,SINTEN,SPHI
      REAL  X(M),Y(M),Y1,XBAR,Diff

      Y1=0.0
      LOPT=0
      LL=M-1
        DO 230 IL=1,LL
!        Extrapolation on lower bound
       IF(XBAR .LE.X(1))THEN
          PPHI=X(1)
          SPHI=X(2)
          PINTEN=Y(1)
          SINTEN=Y(2)
           Diff=(SPHI-PPHI)
           if(Diff .Le.  0.0) Diff=1
          Y1=PINTEN+((SINTEN-PINTEN)*((XBAR-PPHI)/Diff))
               LOPT=1 
           RETURN

!        Extrapolation on upper bound
       ELSEIF(XBAR .GE.X(LL+1)) THEN
           PPHI=X(LL)
           SPHI=X(LL+1)
         PINTEN=Y(LL)
         SINTEN=Y(LL+1)
         Diff=(SPHI-PPHI)
           if(Diff .Le. 0.0) Diff=1
          Y1=PINTEN+((SINTEN-PINTEN)*((XBAR-PPHI)/Diff))
              LOPT=1 
           RETURN

!      interpolation
       ELSEIF (XBAR .GE.X(IL) .AND.XBAR.LE. X(IL+1)) THEN
         PPHI=X(IL)
         SPHI=X(IL+1)
         PINTEN=Y(IL)
         SINTEN=Y(IL+1)
         Diff=(SPHI-PPHI)
           if(Diff .Le. 0.0) Diff=1
          Y1=PINTEN+((SINTEN-PINTEN)*((XBAR-PPHI)/Diff))
          LOPT=1  
           RETURN 
         ENDIF 
  230    CONTINUE 
           end subroutine  INTERP
            end module  Linear_interpolation
            
            
            