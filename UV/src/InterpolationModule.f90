MODULE InterpolationModule
!==============================================================================
!     
! FILENAME:
!     InterpolationModule.f90
!
! DESCRIPTION:
!     This module contains the functions that are used for interpolation.
!
! AUTHORS:
!     Ellyne Kinney / Science Systems and Applications, Inc.
!
! HISTORY:
!     18-Feb-2005 EKK Initial version
!==============================================================================
 IMPLICIT NONE

 CONTAINS

 FUNCTION FindTableEntry(Value,TableValues,nTableEntries,BoundVal1,BoundVal2, &
                        IndexVal1,IndexVal2,FractionalVal)RESULT(status)
!==============================================================================
!
! TITLE:
!     Finds Bounding values for specified number in the lookup table 
!
! NAME:
!     FindTableEntry
!
! INPUTS:
!     Value                     Input value
!     TableValues               Array of values in lookup table
!     nTableEntries             Number of values in array
!
! OUTPUTS:
!     BoundVal1, BoundVal2      The Bounding value of 'Value' in 'TableValues'
!     IndexVal1, IndexVal2      Indices of 'BoundVal1' and 'BoundVal2'
!     FractionalVal             The fractional distance of Value between
!                                BoundVal1 & BoundVal2
! HISTORY:
!     22-Feb-2005 EKK  Changed to F90 module. Removed Implict variables.
!==============================================================================

  IMPLICIT NONE

  REAL(KIND=4),    INTENT(IN)               :: Value
  REAL(KIND=4),                DIMENSION(:) :: TableValues
  INTEGER(KIND=4)                           :: nTableEntries
  REAL(KIND=4),    INTENT(OUT)              :: BoundVal1, BoundVal2
  INTEGER(KIND=4), INTENT(OUT)              :: IndexVal1, IndexVal2
  REAL(KIND=4),    INTENT(OUT)              :: FractionalVal
  INTEGER(KIND=4)                           :: DesFlag

  INTEGER(KIND=4)                           :: iTable
  REAL(KIND=4)   , DIMENSION(nTableEntries) :: tmpArray
  REAL(KIND=4)                              :: tmpValue

  INTEGER(KIND=4)        :: STATUS

  DesFlag=0
  tmpArray=0.0d00
  tmpValue=0.0d00

  IF(TableValues(1) .GT. TableValues(nTableEntries))THEN
    DesFlag=1
    DO iTable = 1,nTableEntries
      tmpArray(iTable) = TableValues(nTableEntries+1-iTable)
    ENDDO
    DO iTable = 1,nTableEntries
      TableValues(iTable) = tmpArray(iTable)
    ENDDO
  ENDIF

  IF((Value .GE. TableValues(1)) .AND. &              ! == Interpolate ==
     (Value .LE. TableValues(nTableEntries))) THEN

     DO iTable  = 1,nTableEntries-1
       IF((Value .GE. TableValues(iTable)) .AND. &
          (Value .LE. TableValues(iTable+1))) EXIT
     ENDDO

     BoundVal1 = TableValues(iTable)
     BoundVal2 = TableValues(iTable+1)
     IndexVal1 = iTable
     IndexVal2 = iTable + 1

  ELSE ! IF((Value .GE. TableValues(1)) .AND...       ! == Extrapolate ==

     IF(Value .LT. TableValues(1))THEN
        BoundVal1 = TableValues(1)
        BoundVal2 = TableValues(2)
        IndexVal1 = 1
        IndexVal2 = 1
     ELSE ! IF(Value .LT. TableValues(1))THEN  
        IF(Value .GT. TableValues(nTableEntries))THEN
           BoundVal1 = TableValues(nTableEntries-1)
           BoundVal2 = TableValues(nTableEntries)
           IndexVal1 = nTableEntries-1
           IndexVal2 = nTableEntries
        ENDIF ! IF(Value .GT. TableValues(nTableEntries))THEN
     ENDIF ! IF(Value .LT. TableValues(1))THEN

  ENDIF ! IF((Value .GE. TableValues(1)) .AND...

  FractionalVal = (Value - BoundVal1)/(BoundVal2 - BoundVal1)
  IF(DesFlag.eq.0) RETURN

  tmpValue  = IndexVal1
  IndexVal1 = nTableEntries+1-IndexVal2
  IndexVal2 = nTableEntries+1-tmpValue

  tmpValue  = BoundVal1
  BoundVal1 = BoundVal2
  BoundVal2 = tmpValue

  RETURN

 END FUNCTION FindTableEntry

 FUNCTION Interp1D(BoundVal1, BoundVal2, FractionalVal, Value) RESULT(status)
!==============================================================================
!
! TITLE:
!     Interpolate over a 1D array
! NAME:
!     Interp1D
! INPUTS:
!     BoundVal1, BoundVal2
!     FractionalVal
! OUTPUTS:
!     Value
! HISTORY:
!     02-Mar-2005 EKK  Initial Version.
!
!==============================================================================

  IMPLICIT NONE

  REAL(KIND=4), INTENT(IN)  :: BoundVal1, BoundVal2
  REAL(KIND=4), INTENT(IN)  :: FractionalVal
  REAL(KIND=4), INTENT(OUT) :: Value

  INTEGER(KIND=4)           :: STATUS

  IF(BoundVal1 .LT. BoundVal2) THEN            ! Ascending
     Value = BoundVal1 + FractionalVal * (BoundVal2 - BoundVal1)
  ELSE ! IF(BoundVal1 .LT. BoundVal2) THEN     ! Descending
     Value = BoundVal1 - FractionalVal * (BoundVal1 - BoundVal2)
  ENDIF ! IF(BoundVal1 .LT. BoundVal2) THEN

  RETURN

 END FUNCTION Interp1D


END MODULE InterpolationModule
