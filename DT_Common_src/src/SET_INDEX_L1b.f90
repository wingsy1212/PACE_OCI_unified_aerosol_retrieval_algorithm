     
     module SET_INDEX_L1b
 
      implicit none

contains   
! ---------------------------------------------------------------------
! SET_INDEX

     
      SUBROUTINE SET_INDEX(START_500,END_500,START_250,END_250,&
                         START_1KM,END_1KM,IDATA)
      IMPLICIT  NONE
      INCLUDE 'read_Sat_MODIS.inc'
      INTEGER START_500,END_500,START_250,END_250,IDATA,START_1KM
      INTEGER END_1KM
      IF( IDATA .EQ.1) THEN
        START_1KM=1
        END_1KM=Iline
        START_500=1
        END_500=Iline*2
        START_250=1
        END_250=Iline*4
      ELSE
        START_500=START_500+(Iline*2)
        END_500=END_500+(Iline*2)
        START_250=START_250+(Iline*4)
        END_250=END_250+(Iline*4)
        START_1KM=START_1KM+Iline
        END_1KM=END_1KM+Iline
      ENDIF 
       end  subroutine SET_INDEX
           end module  SET_INDEX_L1b
           
