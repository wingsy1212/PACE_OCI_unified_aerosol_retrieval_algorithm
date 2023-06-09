module screening
  use viirs_db_utils, only: calc_minmax_ratio,            & 
                            calc_smoke_mask,              &
                            calc_smoke_ae_mask,           &
        									  calc_pyrocb_mask,             &		
								            calc_high_alt_smoke_mask
                            
  use viirs_config, only: viirs_config_type
  use modis_surface, only:  get_LER650,                 &
                            latlon_to_index_ler,        &
                            get_geographic_zone,        &
                            get_sfc_elev_std                          

  use read_l1b, only        : viirs_db_svm

  implicit none

  private

  public  ::  cloud_screen,cloud_screen2
  
  contains  
  
  
  
  integer function cloud_screen(viirs_data,config, wv,platform,lcld_mask, ocld_mask,&
    & month,gzflg, snow_mask, snow_mask2,sfcstd, sr650, minmaxbt,dt_cld_reform) result(status)
    
    implicit none
    character(len=*) , intent(in)                   ::  platform
    type (viirs_db_svm), intent(in)                 ::  viirs_data
    type (viirs_config_type), intent(in)            ::  config
    real, dimension(:,:), intent(in)                ::  wv         ! in cm
    integer, dimension(:,:), intent(in)             ::  dt_cld_reform
   
    real, dimension(:,:), allocatable               ::  minmax
    real, dimension(:,:), allocatable, intent(out)  ::  minmaxbt
    real, dimension(:,:), allocatable, intent(out)  ::  sfcstd, sr650

    integer, dimension(:,:), allocatable, intent(out)   ::  lcld_mask, ocld_mask
    integer, dimension(:,:), allocatable, intent(out)   ::  snow_mask, snow_mask2, gzflg, month
    
    integer                               ::  i,j,ilat, ilon,min_flag
    real                                  ::  ndsi
    real                                  ::  ndvi_lower, ndvi_upper
    real                                  ::  dd    
    
    status  = 0
  
    !-----------------------------------------------------------------------------------------
    ! -- cloud filtering
    !-----------------------------------------------------------------------------------------
    ! Note: Requires reflectances that match MODIS units - (pi*I)/(F*cos(sza))

      allocate(minmax(viirs_data%xscan, viirs_data%scan),minmaxbt(viirs_data%xscan, viirs_data%scan), &
      &   lcld_mask(viirs_data%xscan, viirs_data%scan), ocld_mask(viirs_data%xscan, viirs_data%scan), stat=status)
      if (status /= 0) then 
        print *, "ERROR: Failed to allocate min/max array: ", status
        stop
      end if
      lcld_mask(:,:) = 0
      ocld_mask(:,:) = 0

      allocate(snow_mask(viirs_data%xscan, viirs_data%scan), stat=status)
      if (status /= 0) then
        print *, "ERROR: Unable to allocate snow mask array: ", status
        stop
      end if
      snow_mask(:,:) = 0

      allocate(snow_mask2(viirs_data%xscan, viirs_data%scan), month(viirs_data%xscan, viirs_data%scan),  stat=status)
      if (status /= 0) then
        print *, "ERROR: Unable to allocate snow mask array: ", status
        stop
      end if
      snow_mask2(:,:) = 0
      month(:,:) = config%month
  
      allocate(sr650(viirs_data%xscan, viirs_data%scan), gzflg(viirs_data%xscan, viirs_data%scan),  &
      &        sfcstd(viirs_data%xscan, viirs_data%scan), stat=status)
      if (status /= 0) then
        print *, "ERROR: Unable to allocate sr650 array: ", status
        stop
      end if
      gzflg(:,:) = -999
      sr650(:,:) = -999.0  
      sfcstd(:,:) = -999.0

      do j = 1, viirs_data%scan
        do i = 1, viirs_data%xscan
      
          if (viirs_data%lat(i,j) < -900.0 .OR. viirs_data%lon(i,j) < -900.0) cycle
      
          status = latlon_to_index_LER(viirs_data%lat(i,j), viirs_data%lon(i,j), ilat, ilon)
          if (status /= 0) then
            print *, "ERROR: Failed to convert lat/lon to LER index: ", status, viirs_data%lat(i,j), viirs_data%lon(i,j)
            stop
          end if
          
          sr650(i,j) = get_LER650(ilat, ilon, viirs_data%ndvi(i,j), viirs_data%sca(i,j),  &
          &             viirs_data%raa(i,j),min_flag)/100.0    ! convert output to value from 0.0 to 1.0
      
          gzflg(i,j) = get_geographic_zone(viirs_data%lat(i,j), viirs_data%lon(i,j), status)
          if (status /= 0) then
            print *, "ERROR: Failed to get geographic zone: ", viirs_data%lat(i,j), viirs_data%lon(i,j), status
            stop
          end if

          sfcstd(i,j) = get_sfc_elev_std(viirs_data%lat(i,j), viirs_data%lon(i,j), status)
          if (status /= 0) then
            print *, "ERROR: Failed to get geographic zone: ", viirs_data%lat(i,j), viirs_data%lon(i,j), status
            stop
          end if
      
        end do
      end do

      if (platform .eq. 'VIIRS') status = calc_minmax_ratio(viirs_data%m01_refl, minmax)
!       if (platform .eq. 'VIIRS') status = calc_minmax_ratio(viirs_data%m15_bt, minmaxbt)
      if (platform .eq. 'AHI' .or. platform .eq. 'GOES') status = calc_minmax_ratio(viirs_data%m03_refl, minmax)
      if (platform .eq. 'AHI' .or. platform .eq. 'GOES') status = calc_minmax_ratio(viirs_data%m15_bt, minmaxbt)

      if (status /= 0) then
        print *, "ERROR: Failed to calculate min/max ratio: ", status
        stop
      end if
      
      if (platform .eq. 'VIIRS') then 
        where(minmax > 1.2) lcld_mask = 1
        where(sr650 >= 0.08 .AND. minmax > 1.15) lcld_mask = 1
      end if
      if (platform .eq. 'AHI' .or. platform .eq. 'GOES') then
        where(minmaxbt > 1.04 .OR. minmax > 1.7) lcld_mask = 1
        where(viirs_data%m15_bt < 280. .AND. minmax > 1.2) lcld_mask = 1

        where(gzflg == 12 .AND. viirs_data%m15_bt >= 298.) lcld_mask = 0 
        where(gzflg == 12 .AND. minmaxbt > 1.03.AND. viirs_data%m15_bt < 298.) lcld_mask = 1 
        where(gzflg == 12 .AND. minmax > 1.7) lcld_mask = 1 
      end if

      if (platform .eq. 'GOES') then
        where(viirs_data%lon < -30.0 .AND.  &
        &      minmax > 1.3) lcld_mask = 1
      end if
            
      do j = 1, viirs_data%scan
        do i = 1, viirs_data%xscan
      
          if (viirs_data%lat(i,j) < -900.0 .OR. viirs_data%lon(i,j) < -900.0) cycle

    !     -- skip bow-tie deletion pixels 
          if (platform .eq. 'VIIRS') then 
            if (viirs_data%m01_refl(i,j) < -900.0) cycle
          end if
          if (platform .eq. 'AHI' .or. platform .eq. 'GOES') then 
            if (viirs_data%m03_refl(i,j) < -900.0) cycle
          end if
     
    ! Over bright surfaces, apply RR1.38/0.66 and BTD11-12 threshold !JH
    ! test addition - 6/7/2011
          if (sr650(i,j) .ge. 0.08) then

!             if (viirs_data%dstar(i,j) .lt. 1.12 .AND. viirs_data%sza(i,j) .lt. 68.0 .OR. &
!             & (viirs_data%dstar(i,j) .lt. 1.2 .AND. viirs_data%sza(i,j) .gt. 68.0)) then
!               ! Thin Cirrus Screening, and Thin Cirrus Over-Correction Rectification
!               if (platform .eq. 'VIIRS' .or. platform .eq. 'GOES') then
!                  if (viirs_data%m09_refl(i,j).gt.0.018 .and. wv(i,j).ge.0.9) lcld_mask(i,j) = 1			
!               end if
!           
!               if (gzflg(i,j) == 24) then    ! Taklimakan
!                 if (platform .eq. 'VIIRS' .or. platform .eq. 'GOES') then
!                     if (viirs_data%m09_refl(i,j).gt.0.04 .and. wv(i,j).lt.0.9 .and. &
!                 &   viirs_data%btd11(i,j).gt.-2.5 .and. viirs_data%m05_refl(i,j).lt.0.45) lcld_mask(i,j) = 1
!                 end if
!                 if (viirs_data%m15_bt(i,j) .lt. 265.0) lcld_mask(i,j) = 1	! Cloud Screening
!               else                          ! Everywhere else
!                 if (platform .eq. 'VIIRS' .or. platform .eq. 'GOES') then
!                   if (viirs_data%m09_refl(i,j).gt.0.018 .and. wv(i,j).lt.0.9 .and. &
!                 &   viirs_data%btd11(i,j).gt.-1.0 .and. viirs_data%m05_refl(i,j).lt.0.55) &
!                 & 	 lcld_mask(i,j) = 1
!                 end if
!                 if (platform .eq. 'VIIRS' .and. viirs_data%m15_bt(i,j) .lt. 270.0) lcld_mask(i,j) = 1	! Cloud Screening
!                 if (platform .ne. 'VIIRS' .and. viirs_data%m15_bt(i,j) .lt. 265.0) lcld_mask(i,j) = 1	! Cloud Screening
!               end if
!           
!               if (platform .eq. 'VIIRS' .and. sr650(i,j) .ge. 0.16) then
!                 if (viirs_data%m15_bt(i,j) .ge. 270.0 .and. viirs_data%m15_bt(i,j) .lt. 281.0 .and. &
!               &    viirs_data%btd11(i,j) .gt. -0.5) lcld_mask(i,j) = 1   ! Cloud Edge Contamination Removal
!               end if
!             endif
          end if

          if (platform .eq. 'GOES') then         
            if (viirs_data%m09_refl(i,j)/cos(viirs_data%sza(i,j)*3.14159/180.).gt.0.03)     lcld_mask(i,j) = 1 ! Cloud Screening ?
          end if
                          
    !     Over dark surfaces      
!           if (sr650(i,j).lt.0.08) then
!             if (viirs_data%dstar(i,j) .lt. 1.12 .AND. viirs_data%sza(i,j) .lt. 68.0 .OR. &
!             & (viirs_data%dstar(i,j) .lt. 1.2 .AND. viirs_data%sza(i,j) .gt. 68.0)) then  
!               if (platform .eq. 'VIIRS') then         
!                 dd = viirs_data%m05_refl(i,j)/viirs_data%m01_refl(i,j)
!                 if (viirs_data%m11_refl(i,j).gt.0.36 .and. dd.lt.0.95)     lcld_mask(i,j) = 1 ! Cloud Screening ?
!                 !Thin Cirrus Screening, and Thin Cirrus Over-Correction Rectification
!                 if (viirs_data%m09_refl(i,j).gt.0.018 .and. wv(i,j).gt.0.4) lcld_mask(i,j) = 1              
!               end if
!               if (platform .eq. 'AHI' .or. platform .eq. 'GOES') then         
!                 if (viirs_data%m11_refl(i,j).gt.0.36)     lcld_mask(i,j) = 1 ! Cloud Screening ?
!               end if
!               if (platform .eq. 'GOES') then         
!                 !if (viirs_data%m11_refl(i,j)/cos(viirs_data%sza(i,j)*3.14159/180.).gt.0.05)     lcld_mask(i,j) = 1 ! Cloud Screening ?
!                 !if (viirs_data%m11_refl(i,j)/cos(viirs_data%sza(i,j)*3.14159/180.).gt.0.09)     lcld_mask(i,j) = 1 ! Cloud Screening ?
!                 if (minmax(i,j) > 1.15)     lcld_mask(i,j) = 1 ! Cloud Screening ?
!                 if (viirs_data%sza(i,j) > 80. .and. minmax(i,j) > 1.1)     lcld_mask(i,j) = 1 ! Cloud Screening ?
!               end if
!               
!               if (platform .ne. 'PACE') then                             
!                 if (viirs_data%m15_bt(i,j) .lt. 270.0)                          lcld_mask(i,j) = 1 !Cloud Screening
!               end if
!             end if 
!           end if
      
          ! High AMF
          if (platform .eq. 'VIIRS' .or. platform .eq. 'GOES') then
           if (viirs_data%lat(i,j) > 60.0 .and. viirs_data%amf(i,j) > 7.0) then  
             if (viirs_data%m09_refl(i,j) > 0.018 .and. wv(i,j) > 0.9) lcld_mask(i,j) = 1
             if (viirs_data%m03_refl(i,j) > 0.4) lcld_mask(i,j) = 1
           end if
          end if

    ! -- only apply these cloud tests over region 5 (N. Africa) and southern Africa.
    ! -- developed using MYD021KM.A2011171.A12[25|30]*.hdf as test case.
    
    ! -- over *really* bright surfaces, except Arabian Peninsula (6 <= gzflg <= 11)
!           if (platform .eq. 'VIIRS') then    
!           if (gzflg(i,j) < 6 .OR. gzflg(i,j) > 11) then 
!            if (sr650(i,j) > 0.25) then
!              if (viirs_data%dstar(i,j) .lt. 0.85.and. wv(i,j).lt.2.3) lcld_mask(i,j) = 1 
!            end if
!           end if
!           end if             

!           if ((gzflg(i,j) == 5 .OR. gzflg(i,j) == 1 .OR. gzflg(i,j) == 26 .OR. &
!           &  gzflg(i,j) == 27) .OR. (viirs_data%lat(i,j) < 15.0 .AND. &
!           &  ((viirs_data%lon(i,j) >-20.0 .AND. viirs_data%lon(i,j) < 55.0) .AND. &
!           &  gzflg(i,j) == -999))) then
!   
!             if (sr650(i,j) < 0.08 .AND. viirs_data%dstar(i,j) .lt. 1.12) then 
!               if (viirs_data%m03_refl(i,j) > 0.15) then 
!                 if (viirs_data%m11_refl(i,j) > 0.16)  lcld_mask(i,j) = 1   ! used to be 0.16
!               end if
!             end if
!             if (sr650(i,j) > 0.08 .AND. viirs_data%dstar(i,j) .lt. 1.12) then
!               if (viirs_data%m03_refl(i,j) > 0.25) then
!                 if (viirs_data%btd11(i,j) > 3.0) lcld_mask(i,j) = 1 
!               end if
!     !          if (viirs_data%m03_refl(i,j) > 0.20 .AND. viirs_data%m15_bt(i,j) < 284.0) then
!     !            lcld_mask(i,j) = 1
!     !          end if
!             end if
!           end if
               
    ! -- try to detect snow cover and skip pixel.
    ! -- source: http://modis-snow-ice.gsfc.nasa.gov/?c=atbd&t=atbd
          ndsi = -999.0
          
          if (platform .ne. 'GOES') then
          
             if (viirs_data%m04_refl(i,j) > -900.0 .AND. viirs_data%m10_refl(i,j) > -900.0) then
               ndsi = (viirs_data%m04_refl(i,j) - viirs_data%m10_refl(i,j)) / &
               & (viirs_data%m04_refl(i,j) + viirs_data%m10_refl(i,j))
             else
               lcld_mask(i,j) = 1 
             end if

!              if (ndsi >= 0.35 .AND. viirs_data%m07_refl(i,j) > 0.11 .AND. &   !  modified by CH
!              &  viirs_data%m15_bt(i,j) < 286.0) then
!                snow_mask(i,j) = 1
!              end if
             
!              if (ndsi > -0.2 .and. viirs_data%m07_refl(i,j) > 0.11 .and. viirs_data%m15_bt(i,j) < 286.0) then
!                snow_mask2(i,j) = 1
!              endif   
             
             if (platform .eq. 'AHI') then !this was commented out on VIIRS DB, should be check on VIIRS
                if (ndsi >= 0.15 .AND. viirs_data%m07_refl(i,j) > 0.11 .AND. &   
                &  viirs_data%m15_bt(i,j) < 286.0) then
                  snow_mask(i,j) = 1
                end if
 
                ndvi_lower = (-0.5*ndsi) + 0.3
                ndvi_upper = 0.75*ndsi + 0.3
                if (ndsi > 0.0 .AND. ndsi < 0.15 .AND. viirs_data%m07_refl(i,j) > 0.11 .AND. &
                & viirs_data%m15_bt(i,j) < 286.0 .AND. viirs_data%ndvi(i,j) > ndvi_lower .AND. &
                & viirs_data%ndvi(i,j) < ndvi_upper) then
                  snow_mask(i,j) = 1
                end if 
             end if                       
          end if
 
          !merge DT cld mask to DB mask
!           if (viirs_data%lat(i,j) > 50.0 .and. viirs_data%land_mask(i,j) < 2 .and. lcld_mask(i,j)==0 &
!             & .and. dt_cld_reform(i,j)==0) then 
!             lcld_mask(i,j)  = 1
! !             print *, 'use DT cloud'
!           end if
!                 
        end do
      end do 
  
    ! -- sync up land cloud mask with snow mask.
      where (snow_mask == 1) lcld_mask = 1
!      print *, 'cloud filtering complete.'		
  
    ! -- if D* is high enough (indicates aerosols roughly), skip spatvar checks and
    ! -- other cloud filters to avoid losing the aerosol plume.
    ! -- Exclude N. Africa regions (gzone 1-5) though. D* not dependable at high SZA's, so excluded.
!       where ((viirs_data%dstar >= 1.06 .AND. viirs_data%sza < 68.0) .AND. &
!       &  (gzflg < 1 .OR. (gzflg > 5 .AND. gzflg /= 26 .AND. gzflg /= 27 .AND. gzflg /= 24))) lcld_mask = 0 
  
    ! -- reset the land cloud mask over oceans.
      where(viirs_data%land_mask >= 2) lcld_mask = 0     ! reset land mask over ocean.
  
  end function cloud_screen
  
  
  integer function cloud_screen2(viirs_data,platform,ler412,ler488,ler670, sr650, sfcstd, snow_mask, &
    & snow_mask2,ocld_mask, gzflg,lcld_mask,month,lc,lskip_mask, oskip_mask,smoke_mask, &
    & smoke_ae_mask,pyrocb_mask	,high_alt_smoke_mask,m09_noiof,m07_noiof,dt_cld_reform) result(status)  
    implicit none
    
    character(len=*) , intent(in)             :: platform
    type (viirs_db_svm), intent(in)           :: viirs_data
    real, dimension(:,:), intent(in)          :: ler412,ler488,ler670,sr650, sfcstd,m09_noiof,m07_noiof
    integer, dimension(:,:), intent(in)       :: snow_mask, snow_mask2, gzflg, month 

    integer, dimension(:,:), intent(inout)    :: lcld_mask,lc,ocld_mask,dt_cld_reform
    integer, dimension(:,:), allocatable      :: tmp_mask
    
    real, dimension(:,:), allocatable         :: minmax_ler
    integer, dimension(:,:), allocatable, intent(out)      :: smoke_mask, smoke_ae_mask
    integer, dimension(:,:), allocatable, intent(out)      :: pyrocb_mask	,high_alt_smoke_mask	
    integer, dimension(:,:), allocatable, intent(inout) :: lskip_mask, oskip_mask

    integer                                   :: i,j,i1, j1,i2,j2
    real                                      :: rat488670, rat412488,thold 

    status  = 0
    
    allocate(minmax_ler(viirs_data%xscan, viirs_data%scan), stat=status)
    if (status /= 0) then 
      print *, "ERROR: Failed to allocate min/max array: ", status
      stop
    end if     
! -- one final set of cloud checks. Down here since we have to use the LER412 value.
  do j = 1, viirs_data%scan
    do i = 1, viirs_data%xscan

      if (ler412(i,j) > 50.0) lcld_mask(i,j) = 1

      if (platform .eq. 'AHI') then
          if (ler488(i,j) > 50.0 .and. viirs_data%sza(i,j) < 70.0) lcld_mask(i,j) = 1
      end if

      if (platform .eq. 'GOES') then
          if (ler488(i,j) > 50.0 .or. ler488(i,j) < 1.0) lcld_mask(i,j) = 1
          if (sr650(i,j) < 0.1 .and. ler488(i,j) > 30.0) lcld_mask(i,j) = 1
      end if
      
! modified by CH 2/18/2020

      if (sr650(i,j) .lt. 0.1) then
         if (gzflg(i,j) == 12) then
         if (viirs_data%m11_refl(i,j) .gt. 0.05 .AND. ler412(i,j) .gt. 12.0 .AND.  &
        &    ler488(i,j)/ler670(i,j) < 0.7) lcld_mask(i,j) = 1
         else
         if (viirs_data%m11_refl(i,j) .gt. 0.05 .AND. ler412(i,j) .gt. 12.0) lcld_mask(i,j) = 1
         end if
      end if
      
      if (sr650(i,j).lt.0.08) then
        ! if (viirs_data%dstar(i,j) .lt. 1.12) then
!           if (viirs_data%m15_bt(i,j) .gt. 270 .AND. viirs_data%m15_bt(i,j) .lt. 274 .AND. &
!           & ler412(i,j) .gt. 40.0) then
!             lcld_mask(i,j) = 1
!           end if
! 
!           if (viirs_data%sza(i,j) > 72.0) then
!             if (viirs_data%m15_bt(i,j) .gt. 270 .AND. viirs_data%m15_bt(i,j) .lt. 275 .AND. &
!             & ler412(i,j) .gt. 20.0) then
!               lcld_mask(i,j) = 1
!             end if
!           end if
! 
!         end if
      end if
    
    end do
  end do
  if (platform .eq. 'VIIRS') then 
  
  !this screening should be modified for AHI and ABI later, temporarily it is only applied for VIIRS
  
   ! 28 February 2018 JLee, additional snow filter using stricter snow mask and ler412 minmax
     status = calc_minmax_ratio(ler412, minmax_ler)
     if (status /= 0) then
       print *, "ERROR: Failed to calculate min/max ratio of ler412: ", status
       stop
     end if

   !-- N. America
     where (snow_mask2 == 1 .and. gzflg == 13 .and. viirs_data%ndvi < 0.35 .and. &
     &      sfcstd <= 50 .and. minmax_ler > 1.20) lcld_mask = 1

     ! higher ndvi and minmax thresholds for spring at high latitudes
     where (viirs_data%lat > 45 .and. month >= 2 .and. month <= 6 .and. &
     &      snow_mask2 == 1 .and. gzflg == 13 .and. viirs_data%ndvi < 0.45 .and. &
     &      sfcstd <= 50 .and. minmax_ler > 1.30) lcld_mask = 1 

   !-- Rest of the world
     where (snow_mask2 == 1 .and. gzflg /= 13 .and. viirs_data%ndvi < 0.35 .and. &
     &      sfcstd <= 50 .and. minmax_ler > 1.40) lcld_mask = 1

     ! higher ndvi threshold and lower minmax for spring at high latitudes
     where (viirs_data%lat > 50 .and. month >= 2 .and. month <= 6 .and. &
     &      snow_mask2 == 1 .and. gzflg /= 13 .and. viirs_data%ndvi < 0.45 .and. &
     &      sfcstd <= 50 .and. minmax_ler > 1.30) lcld_mask = 1

   !-- Mountainous regions
     where (snow_mask2 == 1 .and. sfcstd > 50 .and. minmax_ler > 1.80) lcld_mask = 1
  end if 
! -- detect anomalous D* values in large dust storm over Taklimakan region and reset
! -- the cloud mask to 0.
  do j = 1, viirs_data%scan
    do i = 1, viirs_data%xscan
    
      if (gzflg(i,j) == 24) then
        if (ler412(i,j) > -900.0 .AND. ler488(i,j) > -900.0 .AND. ler670(i,j) > -900.0) then
          if (ler488(i,j)/ler670(i,j) < 0.65 .AND. &      !viirs_data%m16_bt(i,j) < 280.0 .AND. 
          &   ler412(i,j) > 18.0) then
            lcld_mask(i,j) = 0
          end if
        end if
      end if

    end do
  end do

!-----------------------------------------------------------------------------------------
! -- cloud adjacency
!-----------------------------------------------------------------------------------------
! -- if a pixel is next to a cloudy or snowy pixel, mark it as cloudy too.
! -- only applied to some zones. 
  allocate(tmp_mask(viirs_data%xscan, viirs_data%scan), stat=status)
  if (status /= 0) then 
    print *, "ERROR: Failed to allocate temporary mask array: ", status
    stop
  end if
  tmp_mask(:,:) = 0
  
  do j = 1, viirs_data%scan
    do i = 1, viirs_data%xscan
    
      if ((gzflg(i,j) == 13 .OR. gzflg(i,j) == 18 .OR.  gzflg(i,j) == 29 .OR. gzflg(i,j) == 31 .OR. & !N. America, Fresno, Mexico City, Yuma
      &    gzflg(i,j) == 14) .and. snow_mask2(i,j) == 0) cycle  ! S. America
      
      i1 = max(1, i-1)
      i2 = min(i+1, viirs_data%xscan)
      j1 = max(1, j-1)
      j2 = min(j+1, viirs_data%scan)
      
      if (lcld_mask(i,j) /= 1 .AND. count(lcld_mask(i1:i2,j1:j2) == 1) > 0) tmp_mask(i,j) = 1
    end do
  end do
  where (tmp_mask == 1) lcld_mask = 1
  deallocate(tmp_mask, stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to deallocate temporary mask: ", status
    stop
  end if
  
!-----------------------------------------------------------------------------------------
! -- smoke mask
!-----------------------------------------------------------------------------------------
! -- detect smoke with 412 LER and 2.25 um TOA reflectances. Reset land
! -- cloud mask to 0 if smoke found.
! -- Used later as well in QA section.
  allocate(smoke_mask(viirs_data%xscan, viirs_data%scan), 	&
  &        smoke_ae_mask(viirs_data%xscan, viirs_data%scan),&
  &				pyrocb_mask(viirs_data%xscan,viirs_data%scan),		&
  				high_alt_smoke_mask(viirs_data%xscan, viirs_data%scan), stat=status) 
  if (status /= 0) then 
    print *, "ERROR: Failed to allocate smoke mask array: ", status
    stop
  end if
  smoke_mask(:,:) = 0
  smoke_ae_mask(:,:) = 0
  pyrocb_mask(:,:) = 0
  high_alt_smoke_mask(:,:) = 0
  
!   if (platform .eq. 'VIIRS') then
!    print *, 'entering calc_smoke_mask...'		  
    status = calc_smoke_mask(ler412, ler488, ler670, viirs_data%m09_refl, &
    & viirs_data%m11_refl, viirs_data%m15_bt, gzflg, lc, sr650, smoke_mask)
    if (status /= 0) then
      print *, "ERROR: Failed to calculate the smoke flag: ", status
      stop
    end if
  
!     status = calc_pyrocb_mask(ler412, ler488, ler670, viirs_data%m09_refl, &
!         viirs_data%m11_refl, viirs_data%m14_bt, pyrocb_mask)
!       
!     if (status /= 0) then
!       print *, "ERROR: Failed to calculate smoke mask: ", status
!       stop
!     end if
!   
!     status = calc_high_alt_smoke_mask(ler412, ler488, ler670, viirs_data%m09_refl, &
!         viirs_data%m11_refl, viirs_data%m15_bt, gzflg, sr650, lc, high_alt_smoke_mask)
!     if (status /= 0) then
!       print *, "ERROR: Failed to calculate smoke mask: ", status
!       stop
!     end if
!   end if

  do j = 1, viirs_data%scan
    do i = 1, viirs_data%xscan

          if (gzflg(i,j) == 31 .OR. (gzflg(i,j) == 13 .AND. sr650(i,j)> 0.16)) then
            rat488670 = ler488(i,j) / ler670(i,j)
            rat412488 = ler412(i,j) / ler488(i,j)

            if (platform .eq. 'VIIRS') then
              if (ler488(i,j) > 16.0 .AND. (rat488670 > 0.9 .AND. rat488670 < 1.1) .AND. &
            &  (rat412488 > 0.7 .AND. rat412488 < 0.88) .AND. ler488(i,j) < 40.0 .AND. &
            &   viirs_data%m11_refl(i,j) < 0.22 .AND. viirs_data%m09_refl(i,j) > 0.0032) high_alt_smoke_mask(i,j) = 1
            end if
          end if
    end do
  end do

!  print *, 'done with smoke mask...'	

!   all_output_datasets(16)  = dataset('smoke_mask', smoke_mask)
!   all_output_datasets(17)  = dataset('high_alt_smoke_mask', high_alt_smoke_mask)
!   all_output_datasets(18)  = dataset('pyrocb_mask', pyrocb_mask)
!   all_output_datasets(19)  = dataset('sr650', sr650)
!   all_output_datasets(20)  = dataset('lc', lc)
  
 where (smoke_mask == 1)
   lcld_mask = 0
 end where		

  where (viirs_data%amf < 7.0)
    where (smoke_mask == 1)          lcld_mask = 0
    where (pyrocb_mask == 1)         lcld_mask = 0
    where (high_alt_smoke_mask == 1) lcld_mask = 0
  end where
 
  do j = 1, viirs_data%scan
    do i = 1, viirs_data%xscan
      !merge DT cld mask to DB mask
!       if (viirs_data%lat(i,j) > 50.0 .and. viirs_data%land_mask(i,j) < 2 .and. lcld_mask(i,j)==0 &
!         & .and. dt_cld_reform(i,j)==0) then 
!         lcld_mask(i,j)  = 1
!       end if

      !DT cloud mask over dark surface
       thold  = 0.08
       if (viirs_data%lat(i,j) > 50.0) then 
          thold  = thold+(viirs_data%lat(i,j)-50.)/400. 
       endif 
       if (viirs_data%land_mask(i,j) < 2 .and. sr650(i,j)< thold .and. dt_cld_reform(i,j)==0) lcld_mask(i,j)  = 1
       if (viirs_data%land_mask(i,j) < 2 .and. sr650(i,j)< thold .and. dt_cld_reform(i,j)==1) lcld_mask(i,j)  = 0
      
      !merge DB DT cloud mask
      dt_cld_reform(i,j)  = lcld_mask(i,j)
    end do
  end do 

  
! -- detect smoke over Australia and ConUS, to be used later to reset
! -- the AE to 1.8 (see modis.f90). This test isn't very discriminating and
! -- "detects" smoke along with clouds and dark vegetated regions. So look out!
!   if (platform .eq. 'VIIRS') then 
!   status = calc_smoke_ae_mask(ler412, ler488, ler670, viirs_data%m09_refl, &
!   &   viirs_data%dstar, gzflg, smoke_ae_mask)
!   if (status /= 0) then
!     print *, "ERROR: Failed to calculate the smoke flag2: ", status
!     stop
!   end if
!   end if
    
!-----------------------------------------------------------------------------------------
! -- ocean cloud filters
!-----------------------------------------------------------------------------------------  
!   if (platform .eq. 'VIIRS') then 
!     status = ocean_cloud_filter(viirs_data%lat, viirs_data%m01_refl, viirs_data%m03_refl, &
!     &                           viirs_data%m08_refl, viirs_data%m09_refl, &
!     &                           viirs_data%sza, viirs_data%land_mask, ocld_mask)
!     if (status /= 0) then
!       print *, "ERROR: Ocean cloud filter failed: ", status
!       stop
!     end if
!   end if 
! !   if (platform .eq. 'AHI') then 
! !     status = ocean_cloud_filter_ahi(viirs_data%lat, viirs_data%m03_refl,viirs_data%m09_refl, &
! !     &                           viirs_data%sza, viirs_data%land_mask, ocld_mask,platform)
! !     if (status /= 0) then
! !       print *, "ERROR: Ocean cloud filter failed: ", status
! !       stop
! !     end if
! !   end if  
! 
!   if (platform .eq. 'GOES' .or. platform .eq. 'AHI') then 
!     status = ocean_cloud_filter_goes(viirs_data%lat, viirs_data%m07_refl,viirs_data%m1120_bt, &
!     &                           viirs_data%m11_refl,viirs_data%m09_refl,ler488,viirs_data%sza,&
!     &                           viirs_data%vza, viirs_data%land_mask, ocld_mask,platform,&
!     &                           m09_noiof,m07_noiof,viirs_data%xscan, viirs_data%scan,&
!     &                           viirs_data%dstar,viirs_data%btd11,viirs_data%btd4,viirs_data%btd8)
!     if (status /= 0) then
!       print *, "ERROR: Ocean cloud filter failed: ", status
!       stop
!     end if
!   end if         
  ! -- reset ocean cloud mask over land.
    where(viirs_data%land_mask <= 1) ocld_mask = 0     
 
!-----------------------------------------------------------------------------------------
! -- combine all of our masks into one.
!-----------------------------------------------------------------------------------------
 allocate(lskip_mask(viirs_data%xscan, viirs_data%scan), &
 &        oskip_mask(viirs_data%xscan, viirs_data%scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate bad pixel masks: ", status
    stop
  end if
  lskip_mask(:,:) = 0
  oskip_mask(:,:) = 0 
  if (platform .eq. 'VIIRS') then 
    where (viirs_data%land_mask < -900 .OR. lcld_mask == 1 .OR. viirs_data%m01_refl < -900.0)
      lskip_mask = 1
    end where
    where (viirs_data%land_mask < -900 .OR. ocld_mask == 1 .OR. viirs_data%m01_refl < -900.0) 
      oskip_mask = 1
    end where
  end if   
  if (platform .eq. 'AHI' .or. platform .eq. 'GOES') then 
    where (viirs_data%land_mask < -900 .OR. lcld_mask == 1 .OR. viirs_data%m03_refl < -900.0)
      lskip_mask = 1
    end where
    where (viirs_data%land_mask < -900 .OR. ocld_mask == 1 .OR. viirs_data%m03_refl < -900.0) 
      oskip_mask = 1
    end where
  end if 
! -- skip pixels where land cover == 5, permanent snow and ice.
  where (lc == 5)
    lskip_mask = 1
  end where
! -- skip pixels with D* > 10.0 -- indicated very high cirrus clouds, not aerosols.
!   where (viirs_data%dstar > 10.0)
!     lskip_mask = 1
!   end where  
  
end function cloud_screen2
  
  
end module screening
