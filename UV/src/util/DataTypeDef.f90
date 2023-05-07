MODULE DataTypeDef
!///////////////////////////////////////////////////////////////////////////////
! $Id: DataTypeDef.f90 155 2017-05-23 17:11:03Z jyli $
!-------------------------------------------------------------------------------
! D E S C R I P T I O N :
!    Define common data types that will be used for all my Fortran programs.
!
!
!
! A U TH O R: 
!    Jason Li (SSAI)
!///////////////////////////////////////////////////////////////////////////////


! Initially, I  borrowed the idea from "Numerical Recipes in Fortran 90". In 
! 2017, I decided to use ISO_FORTRAN_ENV module instead. It provides a bunch of
! scalar default-integer named constants. Similar idea, but now as a part of 
! the language standard since Fortran 2003. Note that kind value can be negative
! if a target platform does not support the particular kind. 
!
! INT8, INT16, INT32, INT64:
!   Kind type parameters to specify an INTEGER type with a storage size of 8, 16
!   32 and 64 bits. 
!
! REAL32, REAL64, REAL128:
!   Kind type parameters to specify a REAL type with a storage size of 32, 64, 
!   and 128 bits. 
!

use, intrinsic :: ISO_FORTRAN_ENV, only: I1B => INT8,  &
                                         I2B => INT16, &
                                         I4B => INT32, &
                                         I8B => INT64, &
                                         SP => REAL32, &
                                         DP => REAL64, &
                                         QP => REAL128

IMPLICIT NONE
SAVE

!
! *** Integer data types:
!

!INTEGER,  PARAMETER :: I8B = SELECTED_INT_KIND(18) ! (+/-) 9223372036854775807
!INTEGER,  PARAMETER :: I4B = SELECTED_INT_KIND(9)  ! (+/-)          2147483647
!INTEGER,  PARAMETER :: I2B = SELECTED_INT_KIND(4)  ! (+/-)               32767
!INTEGER,  PARAMETER :: I1B = SELECTED_INT_KIND(2)  ! (+/-)                 127

!INTEGER,  PARAMETER :: I8B = INT64
!INTEGER,  PARAMETER :: I4B = INT32
!INTEGER,  PARAMETER :: I2B = INT16
!INTEGER,  PARAMETER :: I1B = INT8



!
! *** Floating point data types: 
!

! SP = SELECTED_REAL_KIND(6) : (+/-) 3.4028235E+38
!    = SELECTED_REAL_KIND(6,37) 
!    = KIND(1.0)
!INTEGER,  PARAMETER :: SP = SELECTED_REAL_KIND(6,37)
!INTEGER,  PARAMETER :: SP = REAL32

! DP = SELECTED_REAL_KIND(15): (+/-) 1.7976931348623167E+308
!    = SELECTED_REAL_KIND(15,307)
!    = KIND(1.0D0)  
!INTEGER,  PARAMETER :: DP = SELECTED_REAL_KIND(15,307)
!INTEGER,  PARAMETER :: DP = REAL64


!
! *** Complex data type:
!

INTEGER,  PARAMETER :: SPC = KIND((1.0,1.0))
INTEGER,  PARAMETER :: DPC = KIND((1.0D0,1.0D0))


END MODULE DataTypeDef
