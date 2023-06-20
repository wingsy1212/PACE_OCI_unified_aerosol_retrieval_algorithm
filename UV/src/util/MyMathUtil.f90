MODULE MyMathUtil
!///////////////////////////////////////////////////////////////////////////////
! $Id: MyMathUtil.f90 155 2017-05-23 17:11:03Z jyli $
!-------------------------------------------------------------------------------
! D E S C R I P T I O N:
!    A collection of useful math utilities.
!
!    function radians_sp_scalar(deg) result (rad)
!    function radians_sp_1d(deg) result (rad)
!    function interpol(x, fx, xout, fillvalue) result (fxout)
!    function interpol2(x, fx, xout) result (fxout)
!    function radians_sp_scalar(deg) result (rad)
!    function radians_sp_1d(deg) result (rad)
!    function locate(xx,x)
!    function interpol(x, fx, xout, fillvalue) result (fxout)
!    function interpol2(x, fx, xout) result (fxout)
!
! A U T H O R: 
!    Jason Li (SSAI)
!///////////////////////////////////////////////////////////////////////////////


USE DataTypeDef
USE MyConstants, ONLY: PI

interface radians
   module procedure radians_sp_scalar
   module procedure radians_sp_1d
end interface

contains


!...............................................................................

function radians_sp_scalar(deg) result (rad)

implicit none

! real(SP), parameter, save :: DTOR = PI / 180.0
real(SP), parameter :: DTOR = PI / 180.0

real(SP), intent(in) :: deg
real(SP) :: rad

rad = deg * DTOR

end function radians_sp_scalar

!...............................................................................

function radians_sp_1d(deg) result (rad)

implicit none

! real(SP), parameter, save :: DTOR = PI / 180.0
real(SP), parameter :: DTOR = PI / 180.0

real(SP), intent(in), dimension(:) :: deg
real(SP), dimension(size(deg)) :: rad

rad = deg * DTOR

end function radians_sp_1d

!...............................................................................

function locate(xx,x)

IMPLICIT NONE
REAL(SP), DIMENSION(:), INTENT(IN) :: xx
REAL(SP), INTENT(IN) :: x
INTEGER(I4B) :: locate
INTEGER(I4B) :: n,jl,jm,ju
LOGICAL :: ascnd
n=size(xx)
ascnd = (xx(n) >= xx(1))
jl=0
ju=n+1
do
   if (ju-jl <= 1) exit
   jm=(ju+jl)/2
   if (ascnd .eqv. (x >= xx(jm))) then
      jl=jm
   else
      ju=jm
   end if
end do
if (x == xx(1)) then
    locate=1
else if (x == xx(n)) then
    locate=n-1
else
    locate=jl
end if
end function locate

!...............................................................................

function interpol(x, fx, xout, fillvalue) result (fxout)

USE MyConstants, only: FILLVALUE_SP
implicit none

real(SP), dimension(:), intent(in) :: x, fx
real(SP), intent(in) :: xout
real(SP), intent(in), optional :: fillvalue
real(SP) :: fxout

integer(I4B) :: ilo, ihi, nx
real(SP) :: dx, frac, tolerance, fill

nx = size(x)

if(present(fillvalue)) then
   fill = fillvalue
else
   fill = FILLVALUE_SP
endif
fxout = fill

ilo = locate(x, xout)
if(ilo > 0 .and. ilo < nx) then
   tolerance = epsilon(xout)
   ihi = ilo + 1
   dx = x(ihi) - x(ilo)
   if(abs(dx) > tolerance) then ! non-zero
      if(abs(fx(ilo)-fill)>tolerance .and. abs(fx(ihi)-fill)>tolerance) then
         frac =  (xout - x(ilo)) / dx
         fxout = (1.0-frac)*fx(ilo) + frac*fx(ihi)
      endif
   else
      fxout = fx(ilo)
   endif
endif    

end function interpol

!...............................................................................

function interpol2(x, fx, xout) result (fxout)

! This is interpol() with extrapolation.

implicit none

real(SP), dimension(:), intent(in) :: x, fx
real(SP), intent(in) :: xout
real(SP) :: fxout

integer(I4B) :: ilo, ihi, nx
real(SP) :: slope


nx = size(x)

ilo = locate(x, xout)

if(ilo <= 0) then
   slope = ( fx(1)-fx(2) ) / ( x(1)-x(2) ) 
   fxout = slope * (xout-x(1)) + fx(1)
else if(ilo >= nx) then
   slope = ( fx(nx)-fx(nx-1) ) / ( x(nx)-x(nx-1) ) 
   fxout = slope * (xout-x(nx)) + fx(nx)
else
   fxout = interpol(x, fx, xout) 
endif    

end function interpol2

END MODULE

