      module compute_center_Glint_Angle
       implicit none

contains   
! ---------------------------------------------------------------------
! COMPUTE_GLINTANGLE    
       
          SUBROUTINE  COMPUTE_GLINTANGLE(MTHET0,MTHET,MDPHI,&
                 GLINT_ANGLE,QA_Flag_Ocean) 

       IMPLICIT NONE
       SAVE

       include 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'

       REAL GLINT_ANGLE, MTHET0,MTHET,MDPHI
       INTEGER QA_Flag_Ocean(12)

       GLINT_ANGLE = 0.0
        IF(MTHET0.GT.0.0.AND.MTHET.GT.0.0.AND.MDPHI.GT.0.0) THEN
         GLINT_ANGLE=(COS(MTHET0*DTR))*(COS(MTHET*DTR))&
                  +((SIN(MTHET0*DTR))*(SIN(MTHET*DTR))&
                  *( COS(MDPHI*DTR)))
         GLINT_ANGLE = (ACOS(GLINT_ANGLE))*RTD
        ENDIF

       IF(GLINT_ANGLE .GE.GLINT_THRESHOLD) then
         QA_Flag_Ocean(9)=1
       ELSE
         QA_Flag_Ocean(9)=0
       ENDIF 
       RETURN
       end  subroutine COMPUTE_GLINTANGLE
           end module  compute_center_Glint_Angle
           
           