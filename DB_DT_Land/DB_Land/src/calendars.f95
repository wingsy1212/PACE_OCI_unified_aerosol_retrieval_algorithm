module calendars
  
  implicit none

! Restrict access to module components. 
! All are private unless otherwise stated.
  private 
  
  public  ::  datetime, gdatetime, operator(-), operator(+), operator(<), &
  &   operator(>), operator(/=), operator(==), operator(<=), operator(>=)
  
  public  ::  day_of_week_from_fixed, timezone_from_longitude,universal_from_local, &
  &   local_from_universal, standard_from_universal, universal_from_standard
  
  public  ::  gregorian_leap_year, fixed_from_gregorian, time_from_gregorian, &
  &   gregorian_from_time, gregorian_year_start, gregorian_year_end,          &
  &   gregorian_year_range, gregorian_year_from_fixed, gregorian_from_fixed,  &
  &   gregorian_date_difference, gregorian_from_doy, doy_from_gregorian,      &
  &   doy_from_fixed, fixed_from_doy, tai93_from_fixed, fixed_from_tai93,     &
  &   season_from_doy, seconds_from_fixed
  

! Define custom data types
! Fixed date and time
  type :: datetime
    integer   ::  day         ! in days
    integer   ::  time        ! in milliseconds
  end type datetime
  
! Gregorian date and time
  type :: gdatetime
    integer   ::  year
    integer   ::  month
    integer   ::  mday
    integer   ::  hour
    integer   ::  min
    integer   ::  sec
    integer   ::  msec
    integer   ::  tzone       ! in hours, tzone = localtime - utc, range=[0,23]
  end type gdatetime

! Define calendar epochs
  integer, private, parameter   :: GREG_EPOCH   = 1
  integer, private, parameter   :: FIXED_EPOCH  = 0

! Gregorian constants
  integer, parameter  ::  JAN  = 1
  integer, parameter  ::  FEB  = 2
  integer, parameter  ::  MAR  = 3
  integer, parameter  ::  APR  = 4
  integer, parameter  ::  MAY  = 5
  integer, parameter  ::  JUN  = 6
  integer, parameter  ::  JUL  = 7
  integer, parameter  ::  AUG  = 8
  integer, parameter  ::  SEP  = 9
  integer, parameter  ::  OCT  = 10
  integer, parameter  ::  NOV  = 11
  integer, parameter  ::  DEC  = 12
  
  integer, parameter  ::  SUN = 0
  integer, parameter  ::  MON = SUN + 1
  integer, parameter  ::  TUE = SUN + 2
  integer, parameter  ::  WED = SUN + 3
  integer, parameter  ::  THU = SUN + 4
  integer, parameter  ::  FRI = SUN + 5
  integer, parameter  ::  SAT = SUN + 6
    
!---------------------------------------------------------------------------------------------------    
  interface operator (-)
    module procedure subtract
  end interface
  
  interface operator (+)
    module procedure add
  end interface

  interface operator (>)
    module procedure greater_than
  end interface
  
  interface operator (>=)
    module procedure greater_than_equal
  end interface
  
  interface operator (<)
    module procedure less_than
  end interface
  
  interface operator (<=)
    module procedure less_than_equal
  end interface
  
  contains

!---------------------------------------------------------------------------------------------------
! FIXED FUNCTIONS
!---------------------------------------------------------------------------------------------------
    
  type (datetime) function subtract(dt1, dt2)
    type (datetime), intent(in) :: dt1, dt2
    
    subtract = datetime(dt1%day - dt2%day, dt1%time - dt2%time)
    if (subtract%time >= 86400000) then
      subtract%day = subtract%day + int(subtract%time/86400000)
      subtract%time = mod(subtract%time,86400000)
    end if
    if (subtract%time < 0) then
      subtract%day = subtract%day + int((subtract%time)/86400000) - 1
      subtract%time = 86400000 + mod(subtract%time,86400000)
    end if
    
    if (subtract%time == 86400000) then
      subtract%time = 0
      subtract%day = subtract%day + 1
    end if
    
    return
    
  end function subtract
  
  type (datetime) function add(dt1, dt2)
    type (datetime), intent(in) :: dt1, dt2

    add = datetime(dt1%day + dt2%day, dt1%time + dt2%time)
    if (add%time >= 86400000) then
      add%day = add%day + int(add%time/86400000)
      add%time = mod(add%time,86400000)
    end if
    if (add%time < 0) then
      add%day = add%day + int(add%time/86400000) - 1
      add%time = 86400000 + mod(add%time, 86400000)
    end if

    return
    
  end function add

  logical function greater_than(dt1, dt2)
    type (datetime), intent(in) :: dt1, dt2
    
    greater_than = .false.
    if (dt1%day > dt2%day) then
      greater_than = .true.
    else if (dt1%day == dt2%day) then
      if (dt1%time > dt2%time) then
        greater_than = .true.
      end if
    end if
    
    return
  
  end function greater_than
  
  logical function greater_than_equal(dt1, dt2)
    type (datetime), intent(in) :: dt1, dt2
    
    greater_than_equal = .false.
    if (dt1 > dt2) then
      greater_than_equal = .true.
    else if (dt1%day == dt2%day .AND. dt1%time == dt2%time) then
      greater_than_equal = .true.
    end if
    
    return
  
  end function greater_than_equal
  
  logical function less_than(dt1, dt2)
    type (datetime), intent(in) :: dt1, dt2
    
    less_than = .false.
    if (dt1%day < dt2%day) then
      less_than = .true.
    else if (dt1%day == dt2%day) then
      if (dt1%time < dt2%time) then
        less_than = .true.
      end if
    end if
    
    return
  
  end function less_than
  
  logical function less_than_equal(dt1, dt2)
    type (datetime), intent(in) :: dt1, dt2
    
    less_than_equal = .false.
    if (dt1 < dt2) then 
      less_than_equal = .true.
    else if (dt1%day == dt2%day .AND. dt1%time == dt2%time) then
      less_than_equal = .true.
    end if
  
    return
  end function less_than_equal
  
  integer function rd(day)
    implicit none
    integer, intent(in) ::  day
    
    rd = day - FIXED_EPOCH
  
  end function rd

  integer function day_of_week_from_fixed(dt1)
    implicit none
    
    type(datetime), intent(in)  ::  dt1
    
    day_of_week_from_fixed = mod(dt1%day - FIXED_EPOCH - SUN, 7)
    
    return
  end function day_of_week_from_fixed
  
! Return the value of dt1 in seconds
  real function seconds_from_fixed(dt1)
    implicit none
    
    type(datetime), intent(in)  ::  dt1
    
    seconds_from_fixed = dt1%day * 86400 + dt1%time/1000.0
    
    return
  end function seconds_from_fixed
	
!---------------------------------------------------------------------------------------------------
! TIME FUNCTIONS
!---------------------------------------------------------------------------------------------------
	
! returns local mean time offset from UTC in milliseconds
	integer function timezone_from_longitude(lon)
	  implicit none
	  
	  real, intent(in)  :: lon
	  
	  integer, parameter :: MSEC_PER_LON = 240000
	  
	  timezone_from_longitude =  lon * MSEC_PER_LON
	  
	  return
	
	end function timezone_from_longitude
	
	type(datetime) function universal_from_local(dt1, lon)
	  implicit none
	  
	  type(datetime), intent(in)  ::  dt1
	  real, intent(in)            ::  lon
	  	  
	  universal_from_local = dt1 - datetime(0, timezone_from_longitude(lon))
	 
	  return
	  
	end function universal_from_local
	
	type(datetime) function local_from_universal(dt1, lon)
	  implicit none
	  
	  type(datetime), intent(in)  ::  dt1
	  real, intent(in)            ::  lon
	  
	  type(datetime)  ::  dt2
	  
	  local_from_universal = dt1 + datetime(0, timezone_from_longitude(lon))
	  
	  return
	
	end function local_from_universal
	
! tz = UTC offset where standard_time = universal_time + tz (in milliseconds)
	type(datetime) function standard_from_universal(dt1, tz)
	  type(datetime), intent(in)  ::  dt1
	  integer, intent(in)         ::  tz
	  
	  standard_from_universal = dt1 + datetime(0, tz)
	  return
	  
	end function standard_from_universal
	
	type(datetime) function universal_from_standard(dt1, tz)
	  type(datetime), intent(in)  ::  dt1
	  integer, intent(in)         ::  tz
	  
	  universal_from_standard = dt1 - datetime(0, tz)
	  return
	  
	end function universal_from_standard
	
!---------------------------------------------------------------------------------------------------
! GREGORIAN FUNCTIONS
!---------------------------------------------------------------------------------------------------

  logical function gregorian_leap_year(yr)
    implicit none
    integer, intent(in)   ::  yr
    
    if ((mod(yr,4) == 0 .AND. mod(yr,100) /= 0) .OR. (mod(yr,400) == 0)) then
       gregorian_leap_year = .true.
    else
       gregorian_leap_year = .false.
    end if
    
    return
    
  end function gregorian_leap_year
    
    
  type(datetime) function fixed_from_gregorian(gdt1)
    implicit none
  
    type(gdatetime), intent(in)   ::  gdt1
    type(datetime)                ::  dt1
    integer :: offset
      
      if (gdt1%month <= 2) then
        offset = 0
      else if (gregorian_leap_year(gdt1%year)) then
        offset = -1
      else
        offset = -2
      end if
      
      fixed_from_gregorian%day = GREG_EPOCH - 1 + 365*(gdt1%year - 1) + floor((gdt1%year - 1)/4.0) - &
      &       floor((gdt1%year-1)/100.0) + floor((gdt1%year - 1)/400.0) +                             &
      &       floor((1.0/12.0)*(367*gdt1%month - 362)) + offset + gdt1%mday
      fixed_from_gregorian%time = time_from_gregorian(gdt1)
      
      return
      
  end function fixed_from_gregorian
  
  integer function time_from_gregorian(gdt1)
    implicit none
    
    type(gdatetime), intent(in) ::  gdt1
    
    time_from_gregorian = gdt1%hour*3600000 + gdt1%min*60000 + gdt1%sec*1000 + gdt1%msec
  
    return
  end function time_from_gregorian
  
  type(gdatetime) function gregorian_from_time(dt1)
    implicit none
    
    type(datetime), intent(in)  ::  dt1
    type(gdatetime)             :: gdt1
    
    gdt1%hour = dt1%time / 3600000
    gdt1%min = (dt1%time - gdt1%hour*3600000) / 60000
    gdt1%sec = (dt1%time - gdt1%hour*3600000 - gdt1%min*60000) / 1000
    gdt1%msec = dt1%time - gdt1%hour*3600000 - gdt1%min*60000 - gdt1%sec*1000
    
    gregorian_from_time = gdt1
    
    return
  end function gregorian_from_time
  
  type(datetime) function gregorian_year_start(yr)
    implicit none
    integer, intent(in) ::  yr
    type(gdatetime)     ::  gdt1
    
    gdt1 = gdatetime(yr, JAN, 1, 0, 0, 0, 0, 0)
    gregorian_year_start = fixed_from_gregorian(gdt1)
    
    return
  end function gregorian_year_start
    
  type(datetime) function gregorian_year_end(yr)
    implicit none
    integer, intent(in) ::  yr
    type(gdatetime)     ::  gdt1
    
    gdt1 = gdatetime(yr, DEC, 31, 23, 59, 59, 86400000, 0)
    gregorian_year_end = fixed_from_gregorian(gdt1)
    
    return
    
  end function gregorian_year_end
  
  subroutine gregorian_year_range(yr, yrange)
    integer, intent(in)   :: yr
    integer, dimension(:), intent(inout)  :: yrange
    type(datetime)  :: dt1, dt2
    
    dt1 = gregorian_year_start(yr)
    dt2 = gregorian_year_end(yr)
    yrange = (/dt1%day, dt2%day/)
    
    return
    
  end subroutine gregorian_year_range
  
  integer function gregorian_year_from_fixed(dt1)
    implicit none
      
    type (datetime), intent(in) :: dt1
      
    integer   ::  d0, d1, d2, d3
    integer   ::  n1, n4, n100, n400
    integer   ::  year
      
    d0 = dt1%day - GREG_EPOCH
    n400 = d0/146097
    d1  = mod(d0, 146097)
    n100 = d1/36524
    d2 = mod(d1, 36524)
    n4 = d2/1461
    d3 = mod(d2, 1461)
    n1 = d3/365
    year = 400*n400 +  100*n100 + 4*n4 + n1
      
    if (n100 == 4 .OR. n1 == 4) then 
      gregorian_year_from_fixed = year
    else
      gregorian_year_from_fixed = year + 1
    endif
      
    return
  end function gregorian_year_from_fixed
  
  type(gdatetime) function gregorian_from_fixed(dt1)
    implicit none
    
    type(datetime), intent(in)  ::  dt1
    
    integer   ::  yr, pdays, correction, mo, mdy
    type(gdatetime) ::  gdt1
    type(datetime)  ::  dt2
    
    yr = gregorian_year_from_fixed(dt1)
    dt2 = gregorian_year_start(yr)
    pdays = dt1%day - dt2%day
    
    gdt1 = gdatetime(yr, MAR, 1, 0, 0, 0, 0, 0)
    dt2 = fixed_from_gregorian(gdt1)
    if (dt1%day < dt2%day) then
      correction = 0
    else if (gregorian_leap_year(yr)) then
      correction = 1
    else
      correction = 2
    end if
    
    mo = floor((1.0/367.0) * (12 * (pdays+correction) + 373))
    
    gdt1 = gdatetime(yr, mo, 1, 0, 0, 0, 0, 0)
    dt2 = fixed_from_gregorian(gdt1)
    mdy = 1 + dt1%day - dt2%day
    
    gdt1 = gregorian_from_time(dt1)
    
    gregorian_from_fixed = gdatetime(yr, mo, mdy, gdt1%hour, gdt1%min, gdt1%sec, &
   &  gdt1%msec, 0)
    
    return
  end function gregorian_from_fixed
  
  type(datetime) function gregorian_date_difference(gdt1, gdt2)
    implicit none
    
    type(gdatetime), intent(in) :: gdt1, gdt2
    
    gregorian_date_difference = fixed_from_gregorian(gdt2) - fixed_from_gregorian(gdt1)
  
  end function gregorian_date_difference
  
  type(gdatetime) function gregorian_from_doy(yr, doy)
    implicit none
    
    integer, intent(in) ::  yr, doy
    type(datetime)      ::  dt1
    type(gdatetime)     ::  gdt1
    
    dt1 = fixed_from_doy(yr, doy)
    gregorian_from_doy = gregorian_from_fixed(dt1)
    
    return
    
  end function gregorian_from_doy
  
  integer function doy_from_gregorian(gdt1)
    implicit none
    
    type(gdatetime), intent(in) ::  gdt1
    type(gdatetime)             ::  gdt2
    type(datetime)              ::  dt1
    
    gdt2 = gdatetime(gdt1%year-1, DEC, 31, 0, 0, 0, 0, 0)
    dt1 = gregorian_date_difference(gdt2, gdt1)
    doy_from_gregorian = dt1%day
    
    return
    
  end function doy_from_gregorian
  
  integer function doy_from_fixed(dt1)
    implicit none
    
    type(datetime), intent(in) ::  dt1
    
    integer   ::  d0, d1, d2, d3
    integer    ::  n1, n100
    
    d0 = dt1%day - GREG_EPOCH
    d1 = mod(d0,146097)
    n100 = d1/36524
    d2 = mod(d1, 36524)
    d3 = mod(d2, 1461)
    n1 = d3/365
  
    if (n1 /= 4 .AND. n100 /= 4) then
      doy_from_fixed = mod(d3,365) + 1
    else
      doy_from_fixed = 366
    end if
  
    return
  
  end function doy_from_fixed
  
  type(datetime) function fixed_from_doy(yr, doy)
    implicit none
    
    integer, intent(in) ::  doy, yr
    type(gdatetime)     ::  gdt1
    type(datetime)      ::  dt1
    
    gdt1 = gdatetime(yr-1, DEC, 31, 0, 0, 0, 0, 0)
    fixed_from_doy = fixed_from_gregorian(gdt1) + datetime(doy,0)
    
    return
  
  end function fixed_from_doy
  
! Compare given year and day of year to beginning and ending of each season.
! Return the season: 1 = Winter, 2 = Spring, 3 = Summer, 4 = Fall
  integer function season_from_doy(yr, doy)
    implicit none
    
    integer, intent(in)   ::  yr
    integer, intent(in)   ::  doy
    
    type(datetime)  ::  dt1, sstart, sstop
    type(gdatetime) ::  gdt1
    
    dt1 = fixed_from_doy(yr,doy)
    gdt1 = gregorian_from_fixed(dt1)
    
    select case (gdt1%month)
     case(1,2,12)                   ! Winter
      season_from_doy = 1
     
     case(3:5)                      ! Spring
      season_from_doy = 2
     
     case(6:8)                      ! Summer
      season_from_doy = 3
     
     case(9:11)                     ! Fall
      season_from_doy = 4
     
     case default                   ! Uh oh
      print *, "ERROR: Invalid month specified: ", gdt1%month
      season_from_doy = -1
    
    end select
    
!   Winter
!    sstart = fixed_from_gregorian(gdatetime(yr-1,DEC,1,0,0,0,0,0))
!    sstop = fixed_from_gregorian(gdatetime(yr,MAR,1,0,0,0,0,0))
!    if (dt1 >= sstart .AND. dt1 < sstop) then
!      season_from_doy = 1
!      return
!    end if
!    
!   Spring
!    sstart = fixed_from_gregorian(gdatetime(yr,MAR,1,0,0,0,0,0))
!    sstop = fixed_from_gregorian(gdatetime(yr,JUN,1,0,0,0,0,0))
!    if (dt1 >= sstart .AND. dt1 < sstop) then
!      season_from_doy = 2
!      return
!    end if
!    
!   Summer
!    sstart = fixed_from_gregorian(gdatetime(yr,JUN,1,0,0,0,0,0))
!    sstop = fixed_from_gregorian(gdatetime(yr,SEP,1,0,0,0,0,0))
!    if (dt1 >= sstart .AND. dt1 < sstop) then
!      season_from_doy = 3
!      return
!    end if
!    
!   Fall
!    sstart = fixed_from_gregorian(gdatetime(yr,SEP,1,0,0,0,0,0))
!    sstop = fixed_from_gregorian(gdatetime(yr,NOV,1,0,0,0,0,0))
!    if (dt1 >= sstart .AND. dt1 < sstop) then
!      season_from_doy = 4
!      return
!    end if
!    
!   Something went wrong.
!    season_from_doy = -1
    return
    
  end function season_from_doy
!---------------------------------------------------------------------------------------------------
! TAI-93 FUNCTIONS
!---------------------------------------------------------------------------------------------------
! TAI-93 is a time standard used in EOS data files.  Time is specified as seconds since
!   1993-01-01, aka January 1, 1993.

! TODO: set up and return a status code somehow.
  integer function tai93_from_fixed(dt)
    type(datetime), intent(in)  :: dt
    
    type(datetime)    ::  tai_epoch, dt1
    type(gdatetime)   ::  gdt1
     
    tai93_from_fixed = -999
    gdt1 = gdatetime(2043, 1, 1, 0, 0, 0, 0, 0)
    dt1 = fixed_from_gregorian(gdt1)
    if (dt >= dt1) then
      print *, "ERROR: Specified datetime out of range for TAI93 values (max=2043Jan1)"
      return
    end if
    
    gdt1 = gdatetime(1993, 1, 1, 0, 0, 0, 0, 0)
    tai_epoch = fixed_from_gregorian(gdt1)
    dt1 = dt - tai_epoch
    
    tai93_from_fixed = dt1%day * 86400 + int(dt1%time/1000)
    
    return
  end function tai93_from_fixed
  
  type(datetime) function fixed_from_tai93(tai)
    integer, intent(in) :: tai
    type(datetime)  ::  tai_epoch
    type(gdatetime) ::  gdt1
    integer         ::  dtai, mstai
    gdt1 = gdatetime(1993,1,1,0,0,0,0,0)
    tai_epoch = fixed_from_gregorian(gdt1)
    
    dtai = tai / 86400
    mstai = mod(tai,86400)*1000
    fixed_from_tai93 = tai_epoch + datetime(dtai, mstai)
    
    return
  end function fixed_from_tai93
  
end module calendars
