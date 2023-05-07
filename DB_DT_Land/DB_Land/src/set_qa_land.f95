subroutine set_qa_land(cxscan,nc,sdv,aot550,nedge,noob,smoke,high_alt_smoke_count,gzone,dstar,flag2,platform)
    implicit none
    
    integer, intent(in)         ::  cxscan
    integer(kind=2), intent(in) ::  nc
    real, intent(in)            ::  sdv
    real, intent(in)            ::  aot550
    integer(kind=2), intent(in) ::  nedge
    integer, intent(in)         ::  noob            ! 0 = majority of retrievals within table limits
    integer, intent(in)         ::  smoke           ! number of smoke pixels in cell
    integer, intent(in)         ::  high_alt_smoke_count  
    integer, intent(in)         ::  gzone           ! geographic zone index
    real, intent(in)            ::  dstar           ! Dstar parameter
    integer, intent(inout)      ::  flag2
    real                        ::  ncadj   
    character(len=*)            ::  platform
 
! debug
!    print *, 'lsflag, scat_ang ',lsflag,scat_ang
!    if (sfc_typ  >  -99) then 
!      print *, 'sfc_typ ', sfc_typ
!    endif
! debug
!   -- initialize QA=1 as default
    flag2 = 1

    ncadj = 1.0
    if (cxscan <= 126 .or. cxscan > 274) ncadj = 7.0/8.0
    if (cxscan <= 80 .or. cxscan > 320) ncadj = 6.0/8.0
    
    if (platform .eq. 'AHI' .or. platform .eq. 'GOES' ) ncadj = 2./8.0

!   -- for ConUS (zones 13, 18), India, Australia, if smoke count gt 0, return QA=3.
    if (gzone == 13 .OR. gzone == 18 .OR. gzone == 15 .OR. gzone == 12) then
      if (high_alt_smoke_count > 0) then
        if (nc > 6*ncadj) flag2 = 2
        if (nc > 19*ncadj) flag2 = 3

        return

      end if

     if (gzone == 13 .OR. gzone == 18) then
      if (smoke > 0) then
        if (nc > 13*ncadj) flag2 = 2
        if (nc > 19*ncadj) flag2 = 3
        return
      end if
     end if

     if (gzone == 15 .OR. gzone == 12) then
      if (smoke > 20) then
        if (nc > 20*ncadj) flag2 = 2
        if (nc > 30*ncadj) flag2 = 3
        return
      end if
     end if


    end if

!   -- lower pixel thresholds over Middle East region. (jlee added)
    if (gzone == 30 .or. gzone == 32 .or. gzone == 33) then
      if (nc > 20*ncadj) flag2 = 2
      if (nc >= 30*ncadj) flag2 = 3
      return
    end if

!     -- do some special processing for N.Africa (gzone 1-5)                                 
!      if ((gzone >= 1 .AND. gzone <= 5) .AND. Dstar > 1.2) then
    if (Dstar > 1.2) then
      if (nc > 24*ncadj) flag2 = 2
      if (nc >= 36*ncadj) flag2 = 3
      return
    end if    

!   -- lower pixel thresholds for high D* pixels, exclude N.Africa zones.      
    if (Dstar >= 1.06 .AND. (gzone < 1 .OR. (gzone > 5 .AND. gzone /= 26 .AND. gzone /= 27))) then
      if (nc > 24*ncadj) flag2 = 2
      if (nc >= 36*ncadj) flag2 = 3
      return
    end if      

!   if pass below checks, increase to qa=2
    if (noob == 1) go to 10           ! majority of retrievals were oob (out of bounds)
    
    if (nedge >= 3) go to 10
    if (aot550 < 1.0 .and. sdv > 0.18) go to 10 
    if (aot550 >= 1.0 .and. sdv > aot550*0.18) go to 10
    
    if (nc > 24*ncadj) flag2 = 2
   
!   if pass below checks, increase to qa=3
    if (gzone /= 14) then               ! over S. America, don't check this.
      if (noob == 1) go to 10           ! majority of retrievals were oob (out of bounds)
    end if
      
    if (nedge >= 1) go to 10
!      if (sdv > 0.15  .and.  aot550 < 0.8) go to 10
      
    if (aot550 < 1.0 .and. sdv > 0.15) go to 10
    if (aot550 >= 1.0 .and. sdv > aot550*0.15) go to 10
    
    if (nc >= 36*ncadj) flag2 = 3
10  continue

    return
    
   

  end subroutine set_qa_land
