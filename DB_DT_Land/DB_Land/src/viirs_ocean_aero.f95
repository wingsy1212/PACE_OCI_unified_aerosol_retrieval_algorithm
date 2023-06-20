module viirs_ocean_aero

  implicit none
  private

                              
  public  ::  unload_viirs_ocean_aerosol_luts
  public  ::  calc_median
  public  ::  viirs_ocean_cell_output
  public  ::  run_ocean_retrieval
  public  ::  check_status

! -- ocean aerosol variables
  real, dimension(:,:,:), allocatable         ::  plut
  real, dimension(:), allocatable             ::  plut_sum
  integer, parameter                  ::  r8 = selected_real_kind(p=13)

  type  ::  viirs_ocean_input
    real                                ::  lat
    real                                ::  lon
    real                                ::  sza
    real                                ::  vza
    real                                ::  raa
    real                                ::  ws

    real                                ::  m03_refl
    real                                ::  m04_refl
    real                                ::  m05_refl
    real                                ::  m07_refl
    real                                ::  m08_refl
    real                                ::  m10_refl
    real                                ::  m11_refl
  end type viirs_ocean_input

  type  ::  viirs_ocean_cell_output
    real, dimension(:,:,:), allocatable       ::  aot
    real, dimension(:,:), allocatable         ::  aot550
    real, dimension(:,:), allocatable         ::  ae
    real, dimension(:,:), allocatable         ::  fmf
    real, dimension(:,:), allocatable         ::  aot550_stdv
    real, dimension(:,:,:), allocatable       ::  error
    integer, dimension(:,:), allocatable      ::  npixels
    real, dimension(:,:), allocatable         ::  ss
    integer, dimension(:,:), allocatable      ::  model_flag
    integer, dimension(:,:), allocatable      ::  alg_flag
    real, dimension(:,:), allocatable         ::  opt_err
  end type viirs_ocean_cell_output
    
  type  ::  viirs_ocean_lut

    integer                                   ::  nsza
    integer                                   ::  nvza
    integer                                   ::  nraa
    integer                                   ::  naot
    integer                                   ::  nfmf
    integer                                   ::  nwspd
    integer                                   ::  nchl
    integer                                   ::  nwave

    real, dimension(:,:), allocatable         ::  aot488
    real, dimension(:,:), allocatable         ::  aot551
    real, dimension(:,:), allocatable         ::  aot675
    real, dimension(:,:), allocatable         ::  aot865
    real, dimension(:,:), allocatable         ::  aot1240
    real, dimension(:,:), allocatable         ::  aot1600
    real, dimension(:,:), allocatable         ::  aot2250
    real, dimension(:,:), allocatable         ::  ae

    real, dimension(:), allocatable           ::  sza
    real, dimension(:), allocatable           ::  vza
    real, dimension(:), allocatable           ::  raa
    real, dimension(:), allocatable           ::  aot
    real, dimension(:), allocatable           ::  fmf
    real, dimension(:), allocatable           ::  wspd
    real, dimension(:), allocatable           ::  chl
    real, dimension(:), allocatable           ::  wave
    real, dimension(:), allocatable           ::  coeff

    real, dimension(:,:,:,:,:,:,:), allocatable ::  lut_m03
    real, dimension(:,:,:,:,:,:,:), allocatable ::  lut_m04
    real, dimension(:,:,:,:,:,:,:), allocatable ::  lut_m05
    real, dimension(:,:,:,:,:,:,:), allocatable ::  lut_m07
    real, dimension(:,:,:,:,:,:), allocatable   ::  lut_m08
    real, dimension(:,:,:,:,:,:), allocatable   ::  lut_m10
    real, dimension(:,:,:,:,:,:), allocatable   ::  lut_m11

  end type viirs_ocean_lut

  type(viirs_ocean_lut)   ::  dust_lut
  type(viirs_ocean_lut)   ::  fine_lut
  type(viirs_ocean_lut)   ::  mari_lut
  type(viirs_ocean_lut)   ::  mix_lut


  contains

! -- perform aerosol retrieval over water
! -- status = 0 -- no errors
! -- status = -1 -- retrieval failed.
! -- status = 1 -- water too turbid or glinty
!
! -- algorithm flag values (alg) 0 = regular, 1 = IR (turbid) retrieval
  type(viirs_ocean_cell_output) function run_ocean_retrieval(lat, lon, sza, vza, raa,  &
  &         r488, r550, r675, r865, r1240, r1600, r2250, ws, ls, bpmask, status, turbid_res, &
  &         bathy, chl, debug, debug_data,elev,cell_resolution,nrt,platform,cxscan,cscan) result(out)

    implicit none
    character(len=*)                        ::  platform
    character(len=255)                      ::  str
    real, dimension(:,:), intent(in)        ::  lat         ! geolocation
    real, dimension(:,:), intent(in)        ::  lon
    real, dimension(:,:), intent(in)        ::  sza
    real, dimension(:,:), intent(in)        ::  vza         ! geometry
    real, dimension(:,:), intent(in)        ::  raa
    real, dimension(:,:), intent(in)        ::  r488
    real, dimension(:,:), intent(in)        ::  r550
    real, dimension(:,:), intent(in)        ::  r675
    real, dimension(:,:), intent(in)        ::  r865        ! TOA reflectivity
    real, dimension(:,:), intent(in)        ::  r1240
    real, dimension(:,:), intent(in)        ::  r1600
    real, dimension(:,:), intent(in)        ::  r2250
    real, dimension(:,:), intent(in)        ::  ws          ! wind speed
    real, dimension(:,:), intent(in)        ::  chl       ! Climatological Chl
    integer, dimension(:,:), intent(in)     ::  ls          ! land/sea mask
    integer, dimension(:,:), intent(in)     ::  bpmask      ! bad pixel mask
    integer, dimension(:,:), intent(in)     ::  bathy       ! bathymetry
    integer, intent(out)                    ::  status
    type(viirs_ocean_cell_output), dimension(:), intent(out), optional ::  debug_data
    logical, intent(in), optional           ::  debug
    real, dimension(:,:), intent(inout)     ::  turbid_res
    real, dimension(:,:), intent(in)        ::  elev
    real(kind=r8)             				      ::  opt_err 
    integer, intent(in)                     ::  cxscan,cscan


    integer, dimension(2)                   ::  dims, cdims
    
    real, dimension(7)                      ::  tref           ! input TOA reflectivity
    real                                    ::  taot550
    real, dimension(7)                      ::  taot           ! output ocean AOT
    real(kind=r8), dimension(7)             ::  terr           ! output ocean pixel error
    real                                    ::  tfmf
    real                                    ::  tae
    real                                    ::  tss, vza_limit,sza_limit
    integer                                 ::  talg,cell_resolution
    real                                    ::  tres          ! temp for turbid residual test

    integer                                 ::  i, j, k
    integer                                 ::  i1, i2, j1, j2
    integer                                 ::  imod, nrt,oaod_count
    
    real, dimension(:,:), allocatable       ::  aot550, topt_err
    real, dimension(:,:,:), allocatable     ::  aot
    real, dimension(:,:), allocatable       ::  ae
    real, dimension(:,:), allocatable       ::  fmf
    real, dimension(:,:), allocatable       ::  ss
    real, dimension(:,:,:), allocatable     ::  err 
    integer, dimension(:,:), allocatable    ::  alg
    real, dimension(4)                      ::  test_ss

    real                                    ::  glint
    
    type(viirs_ocean_cell_output)           ::  tmp_out
    type(viirs_ocean_cell_output)           ::  dust1, fine1, mari1, mix1

    type(viirs_ocean_lut)                   ::  lut
    integer                                 ::  min_model
    
    
    dims = (/size(lat,1), size(lat,2)/)   
!     cdims = ceiling(dims / float(cell_resolution))
!     cdims = floor(dims / float(cell_resolution))
    cdims =(/cxscan,cscan/)

!    print *, 'dims: ', dims
!    print *, 'cdims: ', cdims
        
!   -- first, check that we actually have water pixels to retrieve over

!     -- no ocean pixels, initialize output and return.
      status = init_viirs_ocean_cell_output(out, dims(1), dims(2), cdims(1), cdims(2), 7)
      call check_status(status,'ERROR: Failed to allocate final aggregated output: ')      

      out%aot(:,:,:)        = -999.0
      out%aot550(:,:)       = -999.0
      out%ae(:,:)           = -999.0
      out%fmf(:,:)          = -999.0
      out%aot550_stdv(:,:)  = -999.0
      out%error(:,:,:)      = -999.0
      out%npixels(:,:)      = 0
      out%ss(:,:)           = -999.0
      out%model_flag(:,:)   = -999
      out%alg_flag(:,:)     = -999
      out%opt_err(:,:)      = -999.0


    status = 0
    return

  end function run_ocean_retrieval



  subroutine unload_viirs_ocean_aerosol_luts(status)
    implicit none
    
    integer, intent(inout)  ::  status
  
    status = 0
    return
    
  end subroutine unload_viirs_ocean_aerosol_luts

! -- perform a binary search of array xx for j such that x lies between xx(j) and xx(j+1).
! -- see Numerical Recipes in Fortran, Second Edition, p.110
! -- returns  values:
! --     0 if xx(1) > x and size(xx) if xx(size(xx)) < x
! --    -1, 1 if x < xx(1) or x > xx(size(xx)) respectively.
! --    -2 if xx is not monotonic.
  integer function locate(xx,x,status) result(j)
    implicit none
    
    real, dimension(:), intent(in)    ::  xx
    real, intent(in)                  ::  x
    integer, intent(inout)            ::  status
    
    integer                           ::  jl, jm, ju
    integer                           ::  i, n
    
    n = size(xx)
    status = 0

!   -- check for monoticity in xx    
    do i = 1, n-1
      if (.NOT. (xx(1) <= xx(2) .AND. xx(i) <= xx(i+1))) then
        print *, "ERROR: xx in locate() not monotonic."
        status = -2
        return
      end if
    end do
      
!   -- start binary search.
    jl = 0
    ju = n + 1
    do
      if (ju-jl > 1) then
        jm = (ju + jl) / 2
        if ((xx(n) >= xx(1)) .EQV. (x >= xx(jm))) then
          jl = jm
        else
          ju = jm
        end if
      else
        exit
      endif
    end do

!   -- check endpoint equality, otherwise use jl from above.
    if (x == xx(1)) then 
      j = 1
    else if (x == xx(n)) then
      j = n-1
    else
      j = jl
    end if

!   -- set status, indicate success or appropriate failure condition.
    status = 0
    if (j >= n) status = 1
    if (j < 1) status = -1
    
    return
     
  end function locate

! -- calculate the glint reflectance in terms of I/F


      
  integer function linfit(x, y, r) result(status)
    implicit none
    
    real, dimension(:), intent(in)    ::  x
    real, dimension(:), intent(in)    ::  y
    real, dimension(2), intent(inout) ::  r
    
    real                              ::  sx, sy
    real                              ::  sxx, syy, sxy

    integer                           ::  i, n

    n   = size(x)
    sx  = 0.0
    sy  = 0.0
    sxy = 0.0
    sxx = 0.0
    do i = 1, n
      sx = sx + x(i)
      sy = sy + y(i)
      sxx = sxx + (x(i) * x(i))
      sxy = sxy + (x(i) * y(i))
      syy = syy + (y(i) * y(i))
    end do
    
    r(2) = ((n*sxy) - (sx*sy))/((n*sxx)-(sx*sx))
    r(1) = (sy/n)-(r(2)*sx/n)
    
    status = 0
    return
    
  end function linfit


! -- calculate the median of the 2D in_data. Assumes -999.0 as fill value.
	real function calc_median(in_data, status) result(medval)
	  implicit none

	  real, dimension(:,:), intent(in)  ::  in_data
	  integer, intent(inout)            ::  status

	  real, dimension(:), allocatable   ::  median_tmp
	  
	  integer                           ::  defcnt, icnt
	  integer                           ::  i, j
	  
	  defcnt = count(in_data > -900.0)
	  
	  if (defcnt > 15000) then
	    print *, "ERROR: Input array is too big. Must be less than 15000 valid elements."
	    status = -1
	    return
	  end if
	  if (defcnt == 0) then
	    medval = -999.0
	    status = 0
	    return
	  end if
	  
!   -- if there's only a single valid value, just return it.
	  if (defcnt == 1) then
	    medval = minval(in_data, in_data > -900.0)
	    status = 0
	    return
	  end if
	  
!   -- if there's numerous values but they're all the same, return the value.
    if (minval(in_data, in_data > -900.0) == maxval(in_data, in_data > -900.0)) then
      medval = minval(in_data, in_data > -900.0)
      status = 0
      return
    end if
	  
	  allocate(median_tmp(defcnt), stat=status)
	  call check_status(status,"ERROR: Failed to allocate array for median calculation: "   ) 

!   -- copy valid data from 2D in_data to 1D median_tmp array for processing.
	  icnt = 1
	  do j = 1, size(in_data,2)
	    do i = 1, size(in_data,1)
	      if (in_data(i,j) > -900.0) then
	         median_tmp(icnt) = in_data(i,j)
	        icnt = icnt + 1
	      end if
	    end do
	  end do
	  
	  call Sort(median_tmp, defcnt)
	  medval = median(median_tmp, defcnt)

!   -- clean up and return
	  deallocate(median_tmp, stat=status)
	  call check_status(status, "ERROR: Failed to deallocate array for median calculation: "  ) 
	  
	  return
	  
	end function calc_median

! -- calculate the standard deviation of the 2D in_data. Assumes -999.0 as fill value.
	real function calc_stdv(in_data, status) result(stdv)
	  implicit none
	  
	  real, dimension(:,:), intent(in)  ::  in_data
	  integer, intent(inout)            ::  status
	  	  
	  real                              ::  avg
	  real                              ::  var
	  real                              ::  diff
	  integer                           ::  defcnt
	  integer                           ::  icnt
	  integer                           ::  i, j
	  
	  defcnt = count(in_data > -900.0)
	  
	  if (defcnt > 15000) then
	    print *, "ERROR: Input array is too big. Must be less than 15000 valid elements."
	    status = -1
	    return
	  end if
	  if (defcnt <= 1) then
	    stdv = -999.0
	    status = 0
	    return
	  end if

	  avg = sum(in_data, mask=in_data > -900) / real(defcnt)  
    
    var = 0.0
	  do j = 1, size(in_data,2)
	    do i = 1, size(in_data,1)
	      if (in_data(i,j) > -900.0) then
	        diff = in_data(i,j) - avg
	        var  = var + (diff * diff)
	      end if
	    end do
	  end do
	  var = var / real(defcnt-1)
	  stdv = sqrt(var)

	  return

	end function calc_stdv

	integer function init_viirs_ocean_cell_output(t, i, j, x, y, z) result(status)
	  implicit none

	  type(viirs_ocean_cell_output), intent(inout)  ::  t
	  integer, intent(in)                           ::  i     ! pixel resolution
	  integer, intent(in)                           ::  j
	  integer, intent(in)                           ::  x     ! cell resolution
	  integer, intent(in)                           ::  y
	  integer, intent(in)                           ::  z     ! n bands

	  allocate(               &
    &  t%aot(x,y,z),        &
    &  t%aot550(x,y),       &
    &  t%ae(x,y),           &
    &  t%fmf(x,y),          &
    &  t%aot550_stdv(x,y),  &
    &  t%error(i,j,z),      &
    &  t%npixels(x,y),      &
    &  t%ss(x,y),           &
    &  t%model_flag(x,y),   &
    &  t%alg_flag(x,y),     &
    &  t%opt_err(x,y),      &
    &  stat=status)
    call check_status(status, "ERROR: Failed to allocate arrays for ocean output products: "  ) 

    return


	end function init_viirs_ocean_cell_output

	integer function release_viirs_ocean_cell_output(t) result(status)
	  implicit none

	  type(viirs_ocean_cell_output), intent(inout)  ::  t

	  deallocate(                &
	  &  t%aot,                  &
	  &  t%aot550,               &
	  &  t%ae,                   &
	  &  t%fmf,                  &
	  &  t%aot550_stdv,          &
	  &  t%error,                &
	  &  t%npixels,              &
	  &  t%ss,                   &
	  &  t%model_flag,           &
	  &  t%alg_flag,             &
	  &  t%opt_err,              &
	  &  stat=status)
	  call check_status(status, "ERROR: Failed to deallocate VIIRS ocean cell output arrays: "  ) 

	  return
	end function release_viirs_ocean_cell_output

  
  subroutine check_status(status, str)
    integer, intent (in) :: status
    character(len=*), intent (in)    :: str
    
    if(status /= 0) then 
      print *, str, status
    end if
  end subroutine check_status

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse


SUBROUTINE M22INV (A, AINV, OK_FLAG)

  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(2,2), INTENT(IN)  :: A
  DOUBLE PRECISION, DIMENSION(2,2), INTENT(OUT) :: AINV
  LOGICAL, INTENT(OUT) :: OK_FLAG

  DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
  DOUBLE PRECISION :: DET
  DOUBLE PRECISION, DIMENSION(2,2) :: COFACTOR


  DET =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

  IF (ABS(DET) .LE. EPS) THEN
	 AINV = 0.0D0
	 OK_FLAG = .FALSE.
	 RETURN
  END IF

  COFACTOR(1,1) = +A(2,2)
  COFACTOR(1,2) = -A(2,1)
  COFACTOR(2,1) = -A(1,2)
  COFACTOR(2,2) = +A(1,1)

  AINV = TRANSPOSE(COFACTOR) / DET

  OK_FLAG = .TRUE.

  RETURN

END SUBROUTINE M22INV

integer function calc_minmax_ocean(refl412, mm_ratio) result(status)
  implicit none
  
  real, dimension(:,:), intent(in)    ::  refl412
  real, dimension(:,:), intent(inout) ::  mm_ratio

  real        ::  xmin, xmax
  integer     ::  i, j, i1, i2, j1, j2
  
  status = -1
  
  do j = 1, size(refl412,2)
    do i = 1, size(refl412,1)
      if (refl412(i,j) < -900.0) continue
      
      i1 = max(i-1,1)
      i2 = min(i+1,size(refl412,1))
      j1 = max(j-1,1)
      j2 = min(j+1,size(refl412,2))
      xmin = -999.0 
      xmax = -999.0
      if (count(refl412(i1:i2,j1:j2) > -900.0) > 0) then
        xmin = minval(refl412(i1:i2,j1:j2), refl412(i1:i2,j1:j2) > -900.0)
        xmax = maxval(refl412(i1:i2,j1:j2), refl412(i1:i2,j1:j2) > -900.0)
        mm_ratio(i,j) = xmax/xmin
      end if    
      
    end do
  end do
  
  status = 0
  return

end function calc_minmax_ocean


! --------------------------------------------------------------------
! REAL FUNCTION  Median() :
!    This function receives an array X of N entries, copies its value
! to a local array Temp(), sorts Temp() and computes the median.
!    The returned value is of REAL type.
! --------------------------------------------------------------------

   REAL FUNCTION  Median(X, N)
      IMPLICIT  NONE
      REAL, DIMENSION(1:), INTENT(IN)    :: X
      INTEGER, INTENT(IN)                :: N
      REAL, DIMENSION(1:N)            :: Temp
      INTEGER                            :: i

      DO i = 1, N                       ! make a copy
         Temp(i) = X(i)
      END DO
!       CALL  Sort(Temp, N)               ! sort the copy
      IF (MOD(N,2) == 0) THEN           ! compute the median
         Median = (Temp(N/2) + Temp(N/2+1)) / 2.0
      ELSE
         Median = Temp(N/2+1)
      END IF
   END FUNCTION  Median


! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      Real, DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      real                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)               ! assume the first is the min
      Location = Start                  ! record its position
      DO i = Start+1, End               ! start with next elements
         IF (x(i) < Minimum) THEN       !   if x(i) less than the min?
            Minimum  = x(i)             !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      real, INTENT(INOUT) :: a, b
      real                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      real, DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1                  ! except for the last
         Location = FindMinimum(x, i, Size)     ! find min from this to last
         CALL  Swap(x(i), x(Location))  ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort
end module viirs_ocean_aero
