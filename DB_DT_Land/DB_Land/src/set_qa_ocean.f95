subroutine set_qa_ocean(lat, glint, aot550, nc, aot550_stdv, n_total, retrerr, alg_flag, flag2,platform)
    implicit none
    
    real, intent(in)            ::  lat           ! latitude
    real, intent(in)            ::  glint  
    real, intent(in)            ::  aot550
    integer, intent(in)         ::  nc            ! number of pixels use for ocean retr
    real, intent(in)            ::  aot550_stdv
    integer, intent(in)         ::  n_total       ! total number of valid pixel data (includes clouds)
    real, intent(in)            ::  retrerr
    integer, intent(in)         ::  alg_flag      ! ocean retrieval algorithm flag
    integer, intent(inout)      ::  flag2
    character(len=*)                        ::  platform

!   -- initialize QA=1 as default
    flag2 = 1
!   -- ocean        
!   if pass below checks, increase to qa=3
    
    if (platform .eq. 'VIIRS') then 
      if (alg_flag >= 1) then   ! all turbid or a mix of turbid and clean pixels
        if (retrerr < 10.0 .AND. (nc/real(n_total)) >= 0.50 .AND. aot550 < 4.95 .AND. aot550_stdv < 0.5) then
          flag2 = 3
            !new flag for sunglint area in NRT retrieval
            if (alg_flag == 3) then 
              flag2 = 2 
            end if
        end if
      else
        if (retrerr < 10.0 .AND. (nc/real(n_total)) >= 0.30 .AND. aot550 < 4.95 .AND. aot550_stdv < 0.5) then
          flag2 = 3
            !new flag for sunglint area in NRT retrieval
            if (alg_flag == 3) then 
              flag2 = 2 
            end if
        end if
      end if
    else  !else of if (platform .eq. 'VIIRS') then 

!reterr calculation for ABI AHI has problem, need to be updated future    
      if (alg_flag >= 1) then   ! all turbid or a mix of turbid and clean pixels
        if ((nc/real(n_total)) >= 0.50 .AND. aot550 < 4.95 .AND. aot550_stdv < 0.5) then
          flag2 = 3
            !new flag for sunglint area in NRT retrieval
            if (alg_flag == 3) then 
              flag2 = 2 
            end if
        end if
      else
        if ((nc/real(n_total)) >= 0.30 .AND. aot550 < 4.95 .AND. aot550_stdv < 0.5) then
          flag2 = 3
            !new flag for sunglint area in NRT retrieval
            if (alg_flag == 3) then 
              flag2 = 2 
            end if
        end if
      end if
          
    
    end if  !else of if (platform .eq. 'VIIRS') then 
    
    return

  end subroutine set_qa_ocean
