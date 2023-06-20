! @TODO: search() can't handle data that is at the last node of the table
! TODO: Convert all of this to F95

! NOTES: Some "3.5" are sprinkled in here. The AOT tables only currently go up to AOT 3.5, but we
! extrapolate out to a maximum of 5.0. However, we don't want to change the AE or SSA values, so
! in some cases, the AOT values were limited to 3.5 or 3.5 is used in order to maintain the AE and SSA
! calculations in here. At the end, the spectral AOT values are truncated to 3.5 if needed when returned
! to modis.f90. Again, this is done in order to maintain the current AE and SSA values.
      subroutine find_v_viirs(realbuf, tmpvg, outbuf, qa_flag, px_elev, lm, windsp, wv,platform)
      
      use landcover, only: get_landcover
      use calendars, only: gdatetime, season_from_doy, gregorian_from_doy
      use modis_surface, only: get_brdfcorr_sr,
     1                         get_aot500,
     1                         terrain_flag, 
     1                         terrain_flag_new,
     1                         get_geographic_zone,
     1                         get_background_aod, ! 9 January 2018 JLee
     1                         get_LER412, 
     1                         get_LER470, 
     1                         get_LER650, 
     1                         get_modis_LER865,
     1                         get_swir_coeffs412,
     1                         get_swir_coeffs470,
     1                         get_swir_stderr412, 
     1                         get_swir_stderr470,
     1                         get_swir_range
        
      use viirs_aerosol_luts, only:   aero_412,      
     1                                aero_470,
     1                                aero_650,
     1                                new_intep,
     1                                aero_412_abs,
     1                                aero_470_abs,
     1                                aero_412_dust,
     1                                aero_470_dust,
     1                                aero_650_dust,
     1                                aero_412_abs_dust,
     1                                aero_470_abs_dust     

      include 'newaottbl.inc'
      
      real realbuf(26), outbuf(22)
      integer lm      !land mask
c
c     Taking into account the terrain effect on Rayleigh
c
c     ---input parameters
      character(len=*)    ::   platform
      real pdiff(3), Dstar1, tmpvg(7), bt11
c     ---intermediate parameters
      integer doy, ilat, ilon, ilat6, ilon6
      real terrain_flag_new5
      real trflg, pteran, x1, x2, x3, xr470
      real aot_mod(6),r470new, xday, r470ss2 
      real sza, scat_ang_2, amf
      real  outbufvg(21)
      real  model_frac, aod_frac
      
      real            ::  toa_ndvi, tmp412, tmp470, tmp670, ndvi_thold
      real						::	px_elev, bgaod
      integer         ::  lc, gzflg, season, mod_sfc, gzflg_sav, regid
      type(gdatetime) ::  gdt1
      integer         ::  status,min_flag, dum_flag
      
      real            ::  r412_tbl, r470_tbl, r650_tbl, windsp, wv
      real            ::  r412_135, r470_135, r650_135
      
      real    ::  tau_x412, tau_x412ss, tau_x412ss2, tau_x412_91
      real    ::  tau_x412ss_995, tau_x412ss2_995, tau_x412ss_96
      real    ::  tau_x412ss_97, tau_x412ss_94, tau_x412ss_95
      real    ::  tau_x412ss_98, tau_x412ss_91, tau_smoke, tau_x412ss2_98
      
      real    ::  tau_x412_new_91, tau_x412_new_93, tau_x412_new_94
      real    ::  tau_x412_new_96, tau_x412_new_995
      real    ::  tau_x470_new_91, tau_x470_new_92, tau_x470_new_93
      real    ::  tau_x470_new_94, tau_x470_new_95, tau_x470_new_96, tau_x470_new_995

      real    ::  tau_x412_dust_91, tau_x412_dust_93, tau_x412_dust_94
      real    ::  tau_x412_dust_96, tau_x412_dust_995
      real    ::  tau_x470_dust_91, tau_x470_dust_92, tau_x470_dust_93
      real    ::  tau_x470_dust_94, tau_x470_dust_95, tau_x470_dust_96, tau_x470_dust_995

      real    ::  swir_coeffs412(3), swir_coeffs470(3)
      real    ::  swir_stderr412, swir_stderr470
      real    ::  swir_range(2)      
 
      integer ::  tau_x412_flag, tau_x412ss_flag, tau_x412ss2_flag, tau_x412_91_flag
      integer ::  tau_x412_flag_dust, tau_x470_flag_dust, tau_x470_flag_veg 
    
      real    ::  tau_x470, tau_x470ss, tau_x470ss2, tau_x470_new
      integer ::  tau_x470_flag2, tau_x470ss_flag, tau_x470_new_flag, tau_x470_new_91_flag
      real    ::  tau_x650
      integer ::  tau_x650_flag
      integer ::  tau_x412_flag2, tau_x470_flag, tau_x412_flag_91

      real    ::  refp_2100
      real    ::  r470sv, r412sv, r470sv_veg, r412sv_veg
      real    ::  tau_x412sv94, tau_x412sv96, tau_x412sv98
      real    ::  tau_x470sv94, tau_x470sv96, tau_x470sv995
      real    ::  tau_x470sv96_dust
      integer ::  tau_flagsv      
      logical ::  abs_aero_flag
      
      logical ::  sr_fail_flag, do_veg
      integer ::  brdf_flag       
            
      logical dflag, dflag2
      logical ::  debug, use_alternate_brdf
      integer ::  lprint
      integer ::  alg_flag, old_alg_flag
      
c     ---output parameters
      real xtau(3), alpha, ssa(3), tau550, sfc_typ
      integer qa_flag(4)
c     ---common parameter      
      common  /xday/ doy
			
			integer handle_lut_out_of_bounds
	    
	    lprint = -999
      debug = .false.
      
c-----------------------------------------------------------------------
c      Initialization
c-----------------------------------------------------------------------

      mm = 10     ! solar zenith
      nn = 46     ! satellite zenith
      ll = 30     ! rel. azimuth
      ma = 10     ! tau
      mw =  8     ! ssa
      tau_x470_flag_veg = -999  ! out-of-bound flag 
      min_flag  = 0
      dum_flag  = 0
c-----------------------------------------------------------------------
c      Start processing the data
c-----------------------------------------------------------------------

      use_alternate_brdf = .false.
      
c     Load the input data into local storage

      xlat        = realbuf(1)
      xlong       = realbuf(2)
      sza         = realbuf(3)
      xthet       = realbuf(4)  !VZA
      xphi        = realbuf(5)  !RAA
      ref650      = realbuf(6)
      stdv        = realbuf(13)
      toa_ndvi    = realbuf(14)
      bt11        = realbuf(15)
      Dstar1      = realbuf(16)
      pdiff(1)    = realbuf(18)
      pdiff(2)    = realbuf(19)
      pdiff(3)    = realbuf(20)
      pteran      = realbuf(21)

!     dump debug output
      if (abs(xlat-28.1) < 0.005 .AND. abs(xlong-(-16.6)) < 0.005) then
       debug = .true.
       lprint = 1
      end if
!      if (abs(xlat-12.75) > 0.05 .OR. abs(xlong-23.92) > 0.07) go to 10
      
      if (xphi.gt.179.99) go to 10
      if (xphi.lt.6.0) xphi = 6.
 
c     -- sun glint mask
      cc     = 3.14159/180.
      psi    = acos(cos(sza*cc)*cos(xthet*cc) +
     1         sin(sza*cc)*sin(xthet*cc)*cos(xphi*cc))
      glint_ang = psi/cc
c      if (abs(psi/cc).lt.35.0) go to 10

c     -- scattering angle (scat_ang)
      cc     = 3.14159/180.
      psi    = acos(cos(sza*cc)*cos(xthet*cc) -
     1         sin(sza*cc)*sin(xthet*cc)*cos(xphi*cc))
      scat_ang = 180. - psi/cc
c      if (scat_ang .gt. 175.) go to 10

c     -- air mass factor
      amf = 1.0/cos(sza*cc)+1.0/cos(xthet*cc)

      ilat = floor(xlat*10.0) + 900 + 1
      if (ilat > 1800) ilat = 1800
      if (ilat < 1)    ilat = 1
      
      ilon = floor(xlong*10.0) + 1800 + 1
      if (ilon > 3600) ilon = 3600
      if (ilon < 1)    ilon = 1

      ilat6 = floor(xlat/0.06) + 1500 + 1
      if (ilat > 3000) ilat = 3000
      if (ilat < 1)    ilat = 1

      ilon6 = floor(xlong/0.06) + 3000 + 1
      if (ilon > 6000) ilon = 6000
      if (ilon < 1)    ilon = 1      

      season = season_from_doy(2005, doy)   ! year isn't currently available in here.
      gdt1 = gregorian_from_doy(2005,doy)   ! just use 2005 for now.
            
      trflg = terrain_flag(ilon,ilat)
      terrain_flag_new5 = terrain_flag_new(ilon,ilat)
      gzflg = get_geographic_zone(xlat, xlong, status)      
      if (status /= 0) then
        print *, "ERROR: Failed to get geographic zone: ", xlat, xlong, status
        return
      end if

      gzflg_sav = gzflg

      regid = regid_2(ilon,ilat) ! 15 March 2018 JLee added

      ! 9 January 2018 JLee, read bg_aod
      bgaod = get_background_aod(xlat, xlong, season, status)
      if (status /= 0) then
        print *, "ERROR: Failed to get background aod: ", xlat, xlong, status
        return
      end if
             
      xr470 = get_LER470(ilat, ilon, toa_ndvi, scat_ang, xphi,dum_flag)/100.0
      
      do i = 1, 3
      if (pdiff(i).lt.0.0.and.pdiff(i).gt.-1.E-4) 
     1    pdiff(i) = 0.0
!      if (trflg.lt.0.0) pdiff(i) = 0.0
      enddo
      if (xr470.lt.0.) pdiff(2) = 0.
      
!     -- disable pressure corrections
!      pdiff(:) = 0.0
      
      x1 = sza
      x2 = xthet
      x3 = xphi
      refl1 = realbuf(11)          ! 412 nm
      refl3 = realbuf(9)           ! 470 nm
      refl6 = realbuf(7)           ! 650 nm      
      refl11 = realbuf(12)         ! 2250 nm
      refl1_newgc = realbuf(10)    ! 412 nm
      refl3_newgc = realbuf(8)     ! 470 nm


!     -- apply pressure correction to TOA reflectance over Morocco and West Sudan.
      if (px_elev > 750.0 .AND. (xlat > 28.0 .AND. xlat < 37.0 .AND. xlong > -12.0 .AND. 
     1    xlong < 10.0)) then 
        refl3 = realbuf(9) + pdiff(2)/1.30    ! 470 nm
      end if
      
      if (px_elev > 900.0 .AND. (xlat > 10.5 .AND. xlat < 19.5 .AND. xlong > 20.5 .AND. 
     1    xlong < 29.0)) then
        refl3 = realbuf(9) + pdiff(2)/1.50    ! 470 nm
       end if
             
      if (refl6 > 0.2) go to 10

      rr412_mod     = realbuf(22)
      rr470_mod     = realbuf(23)
      band26        = realbuf(17)   !CL_flag
      rat           = refl6 / refl1
      rat1          = rr470_mod / rr412_mod
      rat_650_470   = ref650 / rr470_mod
      rat_650_412   = ref650 / rr412_mod

			if (debug) then
        print *, 'find_v, in: ', xlat, xlong, sza, xthet, xphi, realbuf(11), realbuf(9), realbuf(7), Dstar1, lm
        print *, 'ler, 412, 470, 650: ', realbuf(22), realbuf(23), realbuf(6)
      endif
      
 5    continue
c---------------------------------------------------------------
c     Initialization
c---------------------------------------------------------------
c     --intermediate parameters
      dflag      = .false.
      dflag2     = .false.
      w0_x       = -999.
      w0_int     = -999.
      w0_int_470 = -999.
      w0_x412    = -999.
      w0_x470    = -999.
      
      tau_x412          = -999.
      tau_x412_91       = -999.
      tau_x412ss        = -999.
      tau_x412ss_995    = -999.0
      tau_x412ss2_995   = -999.0
      tau_x412ss2_98   = -999.0
      tau_x412ss_96     = -999.0
      tau_x412ss_97     = -999.0
      tau_x412ss_94     = -999.0
      tau_x412ss_95     = -999.0
      tau_x412ss_98     = -999.0
      tau_x412ss_91     = -999.0
      tau_smoke         = -999.0
      tau_x412_new_91   = -999.0
      tau_x412_new_93   = -999.0
      tau_x412_new_94   = -999.0
      tau_x412_new_96   = -999.0
      tau_x412_new_995  = -999.0
      tau_x470_new      = -999.0
      tau_x470_new_91   = -999.0
      tau_x470_new_92   = -999.0
      tau_x470_new_93   = -999.0
      tau_x470_new_94   = -999.0
      tau_x470_new_95   = -999.0
      tau_x470_new_96   = -999.0
      tau_x470_new_995  = -999.0
      
      tau_x470          = -999.
      tau_x470ss        = -999.
      tau_x470ss2       = -999.
      tau_x650          = -999.
      tau_x412_flag     = -999
      tau_x470_flag     = -999
      tau_x650_flag     = -999
      tau_x412_flag2    = -999
      tau_x470_flag2    = -999
      tau_x412_flag_91  = -999
      xxrat             = -999.
      xxrat2            = -999.
      aot               = -999.
      r412              = -999.0
      r412new           = -999. 
      r412ss            = -999.0
      r412ss2           = -999.0
      r470              = -999.0
      r470ss            = -999.0
      r470ss2           = -999.0
      r470new           = -999.
      r470_sav          = -999.0
      r650              = -999.0
      tau_x412_new      = -999.
      tau_x470_new      = -999.
      
      tau_x412ss2       = -999.
      tau_x470_new_91   = -999.
      sr_fail_flag      = .false.
      abs_aero_flag     = .false.
      mod_sfc           = -999
     
c     -- 2.2 um surface database aod parameters 
      r412sv            = -999.0
      r470sv            = -999.0
      r412sv_veg        = -999.0
      r470sv_veg        = -999.0
      tau_x412sv94      = -999.0
      tau_x412sv96      = -999.0
      tau_x412sv98      = -999.0
      tau_x470sv94      = -999.0
      tau_x470sv96      = -999.0
      tau_x470sv995     = -999.0
      tau_x470sv96_dust = -999.0
      tau_flagsv        = -999 
      swir_coeffs412(:) = -999.0
      swir_coeffs470(:) = -999.0
      swir_stderr412    = -999.0
      swir_stderr470    = -999.0
      swir_range(:)     = -999.0
      
c     -- output parameters
      tau550     = -999.
      alpha      = -999.
      sfc_typ    = -999.
      alg_flag   = 0 
      old_alg_flag   = 0 
      
      do i = 1, 3
         xtau(i) = -999.
         ssa(i) = -999.
      enddo

      do i = 1, 4
         qa_flag(i) = 0
      enddo

      do i = 1, 6
         aot_mod(i) = -999.
      enddo
      
      do i = 1,22
         outbuf(i) = -999.
      enddo
      outbufvg(:) = -999.0
      
      xday = real(doy)

c---------------------------------------------------------------
c     Screen for pixels outside reasonable ranges of reflectance
c---------------------------------------------------------------
c      if (refl1.gt.0.0.and.refl1.lt.0.09.and.
c     1    refl6.gt.0.0.and.refl6.lt.0.14) go to 11
c      if (refl1.gt.0.09.and.refl1.lt.0.50.and.
c     1    res.gt.6.0) go to 11
c      go to 10

11    continue
      
!     -- over land, skip pixels with SZA > 72.0!
      if (platform .eq. 'VIIRS' .and. sza > 84.0) go to 10
C     if (platform .eq. 'AHI' .and. sza > 74.0) go to 10
      if (platform .eq. 'AHI' .and. sza > 84.0) go to 10
      if (platform .eq. 'GOES' .and. sza > 84.0) go to 10

      if (platform .eq. 'AHI' .and. xthet > 80.0) go to 10     

c--------------------------------------------------------
c   Load Surface Reflectance
c--------------------------------------------------------
c     -- get base surf. reflc. values from surf. coeff. tables and save.
      r412_tbl = get_LER412(ilat, ilon, toa_ndvi, scat_ang, xphi,min_flag)
      r470_tbl = get_LER470(ilat, ilon, toa_ndvi, scat_ang, xphi,min_flag)
      r650_tbl = get_LER650(ilat, ilon, toa_ndvi, scat_ang, xphi,min_flag)
      r865_tbl = get_modis_LER865(ilat, ilon)  ! No 865 yet for VIIRS. Just use MODIS value.
      
      r412_135 = get_LER412(ilat, ilon, toa_ndvi, 135.0, xphi,dum_flag)
      r470_135 = get_LER470(ilat, ilon, toa_ndvi, 135.0, xphi,dum_flag)
      r650_135 = get_LER650(ilat, ilon, toa_ndvi, 135.0, xphi,dum_flag)
      r865_135 = get_modis_LER865(ilat, ilon)
      
c     -- set to surface reflectance values to default table values.
      r412 = r412_tbl
      r470 = r470_tbl
      r650 = r650_tbl
      r865 = r865_tbl

      if (r865.lt.12.0.and.glint_ang.lt.30.0) go to 10    ! sun glint mask

      if (debug) print *, 'glint, scat. ang, r865: ', glint_ang, scat_ang, r865
      if (debug) print *, 'ilat, ilon, toa_ndvi, xphi: ', ilat, ilon, toa_ndvi, xphi
      if (debug) print *, "r412, r470, r650: ", r412, r470, r650

c     -- out of scope, skip.      
c      if (r412.gt.30.0) go to 10
c      if (r470.gt.50.0) go to 10
c      if (r650.gt.60.0) go to 10
      
!     -- initialize other surface variables
      r470ss    = r470_tbl
      r470new   = r470_tbl
      
      r412ss    = r412_tbl
      r412ss2   = r412_tbl
      r412new   = r412_tbl
      
!---------------------------------------------------------------------------------------------------
! START NEW AERONET-DERIVED SURFACE REFLECTIVITY
!---------------------------------------------------------------------------------------------------
!     -- get landcover, then surface reflectivity for pixel..  
      lc = get_landcover(xlat, xlong, status)
      if (status .ne. 0) then
        print *, "ERROR: Failed to get land cover value: ", i, j, xlat, xlong, status
        return
      end if

!     -- turn off tropical Sahel transition zone (zone 27) if it's not winter.
      if (gzflg == 27 .AND. (xday >= 60 .AND. xday < 335)) gzflg = 5
      if (gzflg == 27 .AND. realbuf(22) > 8.0) gzflg = 26

      brdf_flag = -1
      tmp412=-999.0 ; tmp470=-999.0 ; tmp670=-999.0
      brdf_flag = get_brdfcorr_sr(xlat, xlong, xphi, scat_ang, xthet, amf, px_elev, gdt1%month, 
     &						toa_ndvi, stdv, gzflg, lc, bgaod, tmp412, tmp470, tmp670, use_alternate_brdf=use_alternate_brdf, 
     &            debug=debug)
      !brdf_flag == -1,-2 : error
      !brdf_flag == 0  : AERONET brdf
      !brdf_flag == 1  : BRDF from table
      if (brdf_flag == 0 .OR. brdf_flag == 1) then
        r412 = tmp412
        r470 = tmp470
        if (tmp670 > -900.0) then
          r650 = tmp670
        end if
      end if
      if (debug) then
        print *, 'surface reflc, 412, 490, 670: ', r412, r470, r650
      end if

      if (gzflg == 31 .and. px_elev >= 750.0) gzflg = 13

      if (gzflg == 18 .and. toa_ndvi < 0.3 .and. lc /= 6 .and. lc /= 4 .and.
     1    realbuf(22) < 12.0 .and. realbuf(22)/realbuf(23) < 0.8) then ! modified by CH 2/3/2017
      r412 = r412 * 1.6
      r470 = r470 * 1.3
      endif

      
!     -- insert special function to calculate surface reflectivity over N.Africa coast 
!     --  e.g. Blida, Saada -- geozone 2.  Overriding values from BRDF in modis_surface.f95.
!      if (gzflg == 2) then
!       mod_sfc = 8
!       call rsfc470(dflag2,mod_sfc,xday,xthet,xphi,scat_ang,r470)
!       if (dflag2) go to 10
!      endif

!     -- special function surface reflectivity for Hamim and all of Arabian Peninsula.
      if (gzflg >= 6 .AND. gzflg <= 11 .AND. gzflg /= 10) then 
        if (r650_135 > 32.0 .AND. (r650_135/r412_135) > 3.7) then
          
!         -- r412_135+0.5 == value from table was a bit too low here. Adjust upwards to
!         -- bring the AOT back down.
          mod_sfc = 9
          call newsfc412_arab(dflag2,mod_sfc,xday,ilat,ilon,xthet,xphi,
     &                        scat_ang,terrain_flag_new5,r412_135,r412new)
            if ((r412_135 < 7.25) .OR. r650_135 > 42.0 .OR. Dstar1 >1.03 .OR. use_alternate_brdf) then
              continue
            else
              r412 = r412new
            end if
        end if
      end if

      
!     -- use surface table value as last resort.  Skip if undefined.
      if (r412 < -900.0) r412 = r412_tbl
      if (r470 < -900.0) r470 = r470_tbl
      if (r650 < -900.0) r650 = r650_tbl

!     -- check if we have any surface information. 
!     -- if not, just skip the pixel if we have an AERONET BRDF but still 
!     -- failed (brdf_flag == -2). 
!     -- Otherwise, set sr_fail_flag and try the veg retrieval.
      if (r412 < -900.0 .OR. r470 < -900.0 .OR. r650 < -900.0) then
        if (brdf_flag == -2) then
          go to 10
        else   
          sr_fail_flag = .true.
          go to 637               ! go to vege. retrieval, no need for SR database values.
        end if
      end if

!     -- 2.2 um surface database approach first. Otherwise, no retrieval
!     results for certain conditions due to a lot of go to statement
!     below

! jlee added to test swir vs. vis surface table for 488 nm surface
! reflectance
      ! output: r470sv, 488 nm surface reflectance

      !if(xlong.gt.-20.0.and.xlong.lt.70.0.and.xlat.gt.5.0.and.xlat.lt.45.0) then
        swir_coeffs412 = get_swir_coeffs412(ilat6,ilon6)
        swir_stderr412 = get_swir_stderr412(ilat6,ilon6)
        swir_coeffs470 = get_swir_coeffs470(ilat6,ilon6)
        swir_stderr470 = get_swir_stderr470(ilat6,ilon6)
        swir_range = get_swir_range(ilat6,ilon6)
  
!        print *, swir_coeffs(:), swir_stderr, swir_range(:)

        refp_2100 = refl11*3.14159/cos(sza*cc)*100 !reflectance in percentage, not I over F
        if(refp_2100.ge.(swir_range(1)-2).and.refp_2100.le.(swir_range(2)+2)) then
          if(swir_coeffs470(1).gt.-900.0.and.swir_coeffs470(2).gt.-900.0.and.
     1       swir_coeffs470(3).gt.-900.0.and.swir_stderr470.lt.1.0) then

            r470sv = swir_coeffs470(1)
     1             + swir_coeffs470(2)*refp_2100
     1             + swir_coeffs470(3)*refp_2100*refp_2100
            r470sv_veg = r470sv
          endif
          !==========
          if(swir_coeffs412(1).gt.-900.0.and.swir_coeffs412(2).gt.-900.0.and.
     1       swir_coeffs412(3).gt.-900.0.and.swir_stderr412.lt.1.0) then

            r412sv = swir_coeffs412(1)
     1             + swir_coeffs412(2)*refp_2100
     1             + swir_coeffs412(3)*refp_2100*refp_2100
            r412sv_veg = r412sv
          endif
        endif
      !endif

      ! AOD using 490 nm band
      if (r470sv > 0.0 .and. r470sv < 24.0) then
        if (r470sv < 1.0) r470sv = 1.0
        !========
        imod = 2                                ! w0 = 0.94
        call aero_470(dflag,refl3_newgc,x1,x2,x3,mm,nn,ll,ma,
     1       imod,r470sv,tau_x470sv94,tau_flagsv,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg,tau_flagsv,tau_x470sv94)
        if (status /= 0) then
          print *,
     1    "ERROR: Failed to check/reset AOT out of bounds condition: ",status
          return
        endif
        !========
        imod = 3                                ! w0 = 0.96
        call aero_470(dflag,refl3_newgc,x1,x2,x3,mm,nn,ll,ma,
     1       imod,r470sv,tau_x470sv96,tau_flagsv,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg,tau_flagsv,tau_x470sv96)
        if (status /= 0) then
          print *,
     1    "ERROR: Failed to check/reset AOT out of bounds condition: ",status
          return
        endif
        !========
        imod = 4                                ! w0 = 0.995
        call aero_470(dflag,refl3_newgc,x1,x2,x3,mm,nn,ll,ma,
     1       imod,r470sv,tau_x470sv995,tau_flagsv,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg,tau_flagsv,tau_x470sv995)
        if (status /= 0) then
          print *,
     1    "ERROR: Failed to check/reset AOT out of bounds condition: ",status
          return
        endif
        !========
      endif

      ! AOD using 412 nm band
      if (refl1_newgc > -900.) then 
      if (r412sv > 0.0 .and. r412sv < 20.0) then
        if (r412sv < 1.0) r412sv = 1.0
        imod = 5                                ! w0 = 0.94
        call aero_412(dflag,refl1_newgc,x1,x2,x3,mm,nn,ll,ma,
     1       imod,r412sv,tau_x412sv94,tau_flagsv,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg,tau_flagsv,tau_x412sv94)
        if (status /= 0) then
          print *,
     1    "ERROR: Failed to check/reset AOT out of bounds condition: ",status
          return
        endif
        !========
        imod = 6                                ! w0 = 0.96
        call aero_412(dflag,refl1_newgc,x1,x2,x3,mm,nn,ll,ma,
     1       imod,r412sv,tau_x412sv96,tau_flagsv,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg,tau_flagsv,tau_x412sv96)
        if (status /= 0) then
          print *,
     1    "ERROR: Failed to check/reset AOT out of bounds condition: ",status
          return
        endif
        !========
        imod = 7                                ! w0 = 0.98
        call aero_412(dflag,refl1_newgc,x1,x2,x3,mm,nn,ll,ma,
     1       imod,r412sv,tau_x412sv98,tau_flagsv,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg,tau_flagsv,tau_x412sv98)
        if (status /= 0) then
          print *,
     1    "ERROR: Failed to check/reset AOT out of bounds condition: ",status
          return
        endif
      endif
      endif

      !r470sv adjustment for tau_x470sv96_dust
      ddx = 0.
      if (px_elev < 500.0) then
        if (xphi < 90.0 .and. r470sv > 0.0 .and. r470sv < 24.0) then
        ddx = xthet*1.5/65.0
        if (r412_tbl > 12.0) ddx = xthet*0.5/65.0
        if (r412_tbl > 14.0) ddx = 0.0
        r470sv = r470sv - ddx
        endif
        if (xphi >= 90.0 .and. r412_tbl > 12.0 .and. r470sv > 0.0 .and. r470sv < 24.0) then
        ddx = xthet*0.8/65.0
        r470sv = r470sv + ddx
        endif
      endif

      ! output: tau_x470sv, 470 nm AOD from default aerosol model (SSA = 0.96)
      if (r470sv > 0.0 .and. r470sv < 24.0) then
        if (r470sv < 1.0) r470sv = 1.0
        imod = 3                                ! w0 = 0.96, dust
        call aero_470_dust(dflag,refl3_newgc,x1,x2,x3,mm,nn,ll,ma,
     1       imod,r470sv,tau_x470sv96_dust,tau_flagsv,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg,tau_flagsv,tau_x470sv96_dust)
        if (status /= 0) then
          print *, 
     1    "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        endif
      endif
! end jlee added
 
      r650_sav = r650
      r470_sav = r470
      r412_sav = r412
            
c     --- DIRECT COPY/PASTE FROM SEAWIFS CODE
!-- 470 nm

      r470new = r470_sav
      
!     -- start tweaking r412ss, r412ss2, r470ss---------------------------------------
!     -- new bi-directional surface 2
      ddx = 1.8
      if (xday > 150.0 .and. xday < 258.0)  ddx = 1.0   ! Jun, Jul
 
      ddx2 = 1.0
      if (xday > 150.0 .and. xday < 258.0)  ddx2 =  1.5 ! Jun, Jul
 
      r470ss  = r470_tbl 
      r470ss2 = r470_tbl
 
      dda1 = 2.0
      dda2 = 1.0
      if (scat_ang >= 100.0 .and. scat_ang < 155.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/55.
       if (scat_ang >= 155.0)
     1   r470ss = r470_tbl+dda1
      r470ss = r470ss  - ddx - ddx2
 
      dda1 = 2.0
      dda2 = 1.0
      if (r470_tbl > 9.0)  then
       if (scat_ang >= 100.0 .and. scat_ang < 155.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/55.
       if (scat_ang >= 155.0)
     1   r470ss = r470_tbl+dda1
      r470ss = r470ss  - ddx - ddx2
      endif
 
!     Libya - Egypt 1
      if (xlat >21. .and. xlong > 12.0) then
      dda1 = 2.5
      if (r412_tbl > 12.0) dda1 = 3.0
      dda2 = 1.0
      if (r412_tbl > 10.0)  then
       if (scat_ang >= 100.0 .and. scat_ang < 155.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/55.
       if (scat_ang >= 155.0)
     1   r470ss = r470_tbl+dda1
      r470ss = r470ss  - ddx 
      endif
      endif
 
!     Libya - Egypt 2
      if (xlat >15. .and. xlong > 22.0) then
      dda1 = 2.5
      if (r412_tbl > 12.0) dda1 = 3.0
      dda2 = 1.0
      if (r412_tbl > 10.0)  then
       if (scat_ang >= 100.0 .and. scat_ang < 155.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/55.
       if (scat_ang >= 155.0)
     1   r470ss = r470_tbl+dda1
      r470ss = r470ss  - ddx 
      endif
      endif
 
!     Central Algeria
      if (xday > 243.0 .and. xday < 274.0)  then     ! Sept
      if (xlat >22.5 .and. xlat <30. .and. xlong > -4.0.and. xlong < 1.) then
      dda1 = 2.5
      if (r412_tbl > 10.)  then
       if (scat_ang >= 100.0 .and. scat_ang < 155.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/55.
       if (scat_ang >= 155.0)
     1   r470ss = r470_tbl+dda1
      r470ss = r470ss  - ddx - ddx2
      endif
      endif
      endif
 
!     N. Algeria 1
      if (xlat >31.5 .and. xlong > 4. .and. xlong < 10.) then
      dda1 = 1.5
      dda2 = 1.0
      if (r412_tbl > 9.0)  then
       if (scat_ang >= 100.0 .and. scat_ang < 155.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-100.)/55.
       if (scat_ang >= 155.0)
     1   r470ss = r470_tbl+dda1
      r470ss = r470ss  - ddx 
      endif
      endif
 
!     NE Mauritania 1
      if (xlat>22.5 .and. xlat<30.0.and. xlong>-11.0.and. xlong< -5.) then
      if (r412_tbl > 9.)  then
        r470ss = r470ss + 0.0
      endif
      endif
 
!     NE Mauritania 2
      if (xlat>22.5 .and. xlat<26.0.and. xlong>-13.0.and. xlong< -11.001) then
      if (r412_tbl > 9.)  then
        r470ss = r470ss + 0.0
      endif
      endif
 
!     NE Mauritania 3
      if (xlat>20.0 .and. xlat<30.0.and. xlong>-20.0.and. xlong< -12.0) then
      if (r412_tbl > 9.)  then
        r470ss = r470ss + 0.0
      endif
      endif
 
!     NE Mauritania 3
      if (xlat>15.0 .and. xlat<20.0 .and. xlong< -14.9) then
      if (r412_tbl > 10.5)  then
        if (xday >= 244.0 .and. xday < 274.0)  r470ss = r470_tbl - ddx - ddx2
      endif
      endif
 
 
      if (xphi > 90.0) then
        dda1 = 0.0
        dda2 = 1.4
        if (xday > 150.0 .and. xday < 258.0) dda2 = 2.0
        if (xthet >= 30.0 .and. xthet < 70.0) dda1 = (xthet-30.0) * dda2/40.
        if (xthet >= 70.0) dda1 = dda2
      r470ss = r470ss + dda1
      endif
 
      if (Dstar1 > 1.1) r470ss = r470_tbl - ddx - ddx2
      if (Dstar1 > 1.06.and.xday > 150.0 .and. xday < 258.0) r470ss = r470_tbl - ddx - ddx2
      if (Dstar1 > 1.04.and.wv>1.5.and.xday > 150.0 .and. xday < 258.0) r470ss = r470_tbl - ddx - ddx2

      if (gzflg >= 6 .and. gzflg <= 11) then
 
      ddx = 0.5 +0.8
      if (xday > 59.0 .and. xday < 151.0)   ddx = 0.5+0.5  ! Mar, Apr, May
      if (xday > 150.0 .and. xday < 258.0)  ddx =  1.0+1.5 ! Jun, Jul, Aug
      if (xday > 243.0 .and. xday < 335.0)  ddx =  0.7+1.5 ! Sep, Oct, Nov
      if (xday > 243.0 .and. xday < 335.0.and.r412_tbl >11.0)  ddx =  0.7+1.0 ! Sep, Oct, Nov
 
      ddx2 = 0.0
      if (xday > 150.0 .and. xday < 258.0.and.r412_tbl <11.0)  ddx2 =  0.5 ! Jun, Jul, Aug
 
      r470ss  = r470 - ddx
      r470ss2 = r470
 
      dda1 = r470_tbl*0.1 +0.5
      dda2 = r470_tbl*0.0
      if (scat_ang >= 90.0 .and. scat_ang < 170.0)
     1   r470ss = r470_tbl + dda1*(scat_ang-90.)/80.
      if (scat_ang >= 170.0)
     1   r470ss = r470_tbl+dda1+dda2*(scat_ang-170.)/10.
      r470ss = r470ss  - ddx - ddx2

      dda1 = r470_tbl*0.1 +1.0
      dda2 = r470_tbl*0.0
      if (r412_tbl > 11.0) dda1 = r470_tbl*0.1 +1.0
      if (r412_tbl > 12.0) dda1 = r470_tbl*0.1 +1.0
      if (xday > 120.0 .and. xday < 152.0 .and. r412_tbl > 11.0)  dda1 = dda1 - 1.5  ! May
      if (xday > 59.0 .and. xday < 244.0.and.r412_tbl > 10.3) dda1 = r470_tbl*0.1 +1.5  ! spring, summer
      if (xday > 243.0 .and. xday < 335.0.and.r412_tbl > 11.0) dda1 = r470_tbl*0.1 +1.0  ! fall

       if (r412_tbl > 9.0)  then
       scat_ang_2 = scat_ang
      if (scat_ang_2 >= 90.0 .and. scat_ang_2 < 160.0)
     1   r470ss = r470_tbl + dda1*(scat_ang_2-90.)/70.
       if (scat_ang_2 >= 160.0)
     1   r470ss = r470_tbl+dda1+dda2*(scat_ang_2-160.)/20.
      r470ss = r470ss  - ddx  - ddx2
       endif
 
      if (xday > 150.0 .and. xday < 244.0) then       ! Jun,Jul, Aug
          if (Dstar1 > 1.1)  r470ss = r470_tbl
      endif
 
      if (xday > 243.0 .and. xday < 274.0.and.Dstar1 > 1.02) then  ! Sept
          r470ss = r470_tbl-1.0
      endif
      if (xday > 120.0 .and. xday < 152.0.and.Dstar1 > 1.08) then  ! May
          r470ss = r470_tbl-1.0
      endif
 
!     Southern Arabian Pen
      if (xday > 152.0 .and. xday < 244.0) then   ! Jun, Jul, Aug
 
      ddx  = 0.0
      if (xday > 152.0 .and. xday < 182.) ddx  = 0.5    ! Jun
      if (xday > 181.0 .and. xday < 244.0) ddx  = 1.5   ! Jul, Aug
      dda1 = 1.5
      dda2 = 0.0
      if (r412_tbl > 11.0) dda1 = 1.8
      if (r412_tbl > 12.0) dda1 = 2.0
      if (xday > 181.0 .and. xday < 244.0) then ! Jul, Aug
        if (r412_tbl > 11.0) dda1 = 1.8
        if (r412_tbl > 12.0) dda1 = 2.0
      endif
 
      if (xlat < 18.0.and. xlong <= 49.0) then
      if (r412_tbl > 11.)  then
       if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1   r470ss = r470_tbl-ddx + dda1*(scat_ang-100.)/75.
       if (scat_ang >= 175.0)
     1   r470ss = r470_tbl-ddx +dda1+dda2*(scat_ang-175.)/5.
      endif
      endif
 
      if (xlat < 18.5 .and. xlong > 49.0.and. xlong <= 52.5) then
      if (r412_tbl > 11.)  then
       if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1   r470ss = r470_tbl-ddx + dda1*(scat_ang-100.)/75.
       if (scat_ang >= 175.0)
     1   r470ss = r470_tbl-ddx +dda1+dda2*(scat_ang-175.)/5.
      endif
      endif
 
      if (xlat < 23.0 .and. xlong > 52.5) then
      if (r412_tbl > 11.)  then
       if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1   r470ss = r470_tbl-ddx + dda1*(scat_ang-100.)/75.
       if (scat_ang >= 175.0)
     1   r470ss = r470_tbl-ddx +dda1+dda2*(scat_ang-175.)/5.
      endif
      endif
 
      endif      ! end of southern Arabian Pen
 
      dda2 = 1.01
      if (xday > 151.0 .and. xday < 244.0) dda2 = 1.1
      if (xphi > 90.0 .and. Dstar1 <dda2) then
        ddx = 0.0
        dda1 = 1.4
        if (xday > 151.0 .and. xday < 244.0) dda1 = 1.8
        if (xday > 151.0 .and. xday < 244.0.and.r412_tbl > 10.3.and.(rr470_mod-r470ss) > 0.0) dda1 = 2.2
        if (xday > 243.0 .and. xday < 335.0) dda1 = 1.8
 !        if (xday > 243.0 .and. xday < 335.0) dda1 = 2.5
        if (xday > 243.0 .and. xday < 335.0.and.r412_tbl > 10.3.and.(rr470_mod-r470ss) > 0.0) dda1 = 2.5
        if (xday > 243.0 .and. xday < 335.0.and.r412_tbl > 11.) dda1 = 1.8
        if (xthet >= 30.0 .and. xthet < 70.0) ddx = (xthet-30.0) * dda1/40.
        if (xthet >= 70.0) ddx = dda1
      r470ss = r470ss + ddx
      endif
 
      endif
 
      if (r470 >= 24.0)   r470 = 23.9
      if (r470 < 1.0)     r470 = 1.0
      if (r470ss >= 24.0) r470ss = 23.9
      if (r470ss < 1.0 .AND. (gzflg > 0 .AND. (gzflg <= 11 .OR. gzflg == 27))) go to 10
      
      if (lprint > 0) print *,'final r470,r470ss =', r470,r470ss
 
!-- 412 nm
 
      ddx = 0.0
      if (xday > 32.0 .and. xday < 60.0)  ddx = 0.5            ! Feb
      if (xday > 59.0 .and. xday < 152.0)   ddx = 1.0 + 0.4+0.3    ! Mar,Apr, May
      if (xday > 151.0 .and. xday < 213.0)  ddx = 1.0 + 0.8+0.5    ! Jun, Jul
      if (xday > 212.0 .and. xday < 244.0)  ddx = 0.5+0.5          ! Aug
      if (xday > 243.0 .and. xday < 335.0)  ddx = 0.5+0.5          ! Sep, Oct, Nov
      if (xday > 273.0 .and. xday < 305.0)  ddx = 0.7+0.5          ! Oct
 
      ddx1 = 0.5
      if (xday < 32.0)   ddx1 = 0.5                         ! Jan
      if (xday > 59.0 .and. xday < 152.0)   ddx1 = 0.7      ! Mar,Apr, May
      if (xday > 151.0 .and. xday < 244.0)  ddx1 = 0.5      ! Jun,Jul,Aug
      if (xday > 243.0 .and. xday < 305.0)  ddx1 = 0.5      ! Sep, Oct
      if (xday > 304.0)  ddx1 = 0.5                         ! Nov, Dec
 
      ddx2 = 0.7
      if (xday < 32.0)   ddx2 = 0.7                         ! Jan
      if (xday > 59.0 .and. xday < 152.0)   ddx2 = 0.7      ! Mar,Apr,May
      if (xday > 243.0 .and. xday < 305.0)  ddx2 = 0.7      ! Sep, Oct
      if (xday > 304.0)  ddx2 = 0.7                         ! Nov, Dec
 
      ddx3 = 0.3
      if (xday > 0.0 .and. xday < 152.0)   ddx3 = 1.0       ! Jan,Feb,Mar,Apr,May
      ddx33 = 0.4
      if (xday > 32.0 .and. xday < 60.0)    ddx33 = 0.7     ! Feb
      if (xday > 59.0 .and. xday < 152.0)   ddx33 = 1.0     ! Mar,Apr, May
      if (xday > 151.0 .and. xday < 244.0)  ddx33 = 0.7     ! Jun,Jul,Aug
      if (xday > 243.0 .and. xday < 305.0)  ddx33 = 0.4     ! Sep, Oct
      if (xday > 304.0) ddx33 = 0.6                         ! Nov, Dec
 
      r412ss2 = r412_tbl - ddx
 
      if (scat_ang >= 120.0 .and. scat_ang < 170.0) 
     1   r412ss2 = r412_tbl- ddx + ddx1*(scat_ang-120.)/50.
       if (scat_ang >= 170.0)  
     1   r412ss2 = r412_tbl- ddx+ddx1+0.0*(scat_ang-170.)/10.
 
      if (r412_tbl > 9.0)  then
       if (scat_ang >= 120.0 .and. scat_ang < 170.0) 
     1   r412ss2 = r412_tbl- ddx + ddx2*(scat_ang-120.)/50.
       if (scat_ang >= 170.0)  
     1     r412ss2 = r412_tbl- ddx+ddx2+0.0*(scat_ang-170.)/10.
      endif
 
      if (r412_tbl > 11.0)  then
       if (scat_ang >= 110.0 .and. scat_ang < 160.0) 
     1   r412ss2 = r412_tbl - ddx+ ddx3*(scat_ang-110.)/50.
       if (scat_ang >= 160.0)  
     1   r412ss2 = r412_tbl- ddx+ddx3+0.0*(scat_ang-160.)/20.
      endif
 
      if (r412_tbl > 12.0)  then
       if (scat_ang >= 110.0 .and. scat_ang < 160.0) 
     1    r412ss2 = r412_tbl - ddx+ ddx33*(scat_ang-110.)/50.
       if (scat_ang >= 160.0)  
     1    r412ss2 = r412_tbl- ddx+ddx33+0.0*(scat_ang-160.)/20.
      endif
 
      if (xphi > 90.0) then
        ddx = 0.0
        ddx33 = 1.4
        if (xday > 243.0 .and. xday < 305.0) ddx33 = 0.0
        if (r412_tbl > 12.0) then
            ddx33 = 2.0
            if (xday > 243.0 .and. xday < 305.0) ddx33 = 1.5
        endif
        if (xthet < 70.0) ddx = xthet * ddx33/70.
        if (xthet >= 70.0) ddx = ddx33
      r412ss2 = r412ss2 + ddx
      endif
            
!     -- new bi-directional surface 2
 
      ddx = 0.0
      if (xday > 0.0 .and. xday < 32.0)   ddx = 0.5 + 0.7    ! Jan
      if (xday > 31.0 .and. xday < 60.0)  ddx = 0.5 + 0.4    ! Feb
      if (xday > 59.0 .and. xday < 152.0)  ddx = 0.5 + 0.4+0.8    ! Mar, Apr, May
      if (xday > 151.0 .and. xday < 244.0) ddx = 0.5 + 0.4+1.0    ! Jun, Jul, Aug
      if (xday > 243.0 .and. xday < 305.0)  ddx = 0.5 + 0.4+1.0     ! Sep, Oct
      if (xday > 304.0 .and. xday < 335.0)  ddx = 0.5 + 0.4     ! Nov
      if (xday > 334.0)                     ddx = 0.5 + 0.4     ! Dec
      ddx2 = 1.0
      if (xday > 59.0 .and. xday < 152.0)   ddx2 = 1.0    ! Mar, Apr,May
      if (xday > 151.0 .and. xday < 244.0)  ddx2 = 1.2    ! Jun,Jul,Aug
      if (xday > 243.0 .and. xday < 305.0)  ddx2 = 0.7    ! Sep, Oct
      if (xday > 304.0 .and. xday < 335.0)  ddx2 = 1.0    ! Nov
      ddx3 = 1.5
      if (xday < 32.0)  ddx3 = 1.0    
      if (xday > 151.0 .and. xday < 244.0)  ddx3 = 1.5    ! Jun,Jul,Aug
      if (xday > 243.0 .and. xday < 305.0)  ddx3 = 0.8    ! Sep, Oct
      ddx33 = 0.0
      if (xday > 151.0 .and. xday < 244.0)  ddx33 = 0.0    ! Jun,Jul,Aug
 
      if (xday > 151.0 .and. xday < 244.0) then
        dda1 = rr412_mod - r412_tbl  
        if (dda1 <-1.5 .and.Dstar1 >1.0 )  ddx2 = 0.    ! Jun,Jul,Aug
      endif
 
      r412ss = r412_tbl - ddx
 
      if (xday > 151.0 .and. xday < 274.0)  then   
       dd = (xlong - 8.0) /4.0
       if (xlong > 12.0) dd = 1.
       if (xlong <= 8.0) dd = 0.
       dda1 = ddx2
       ddx2 = dda1*dd
      endif
 
      if (scat_ang >= 100.0 .and. scat_ang < 160.0)       
     1   r412ss = r412_tbl- ddx + ddx2*(scat_ang-100.)/60. 
      if (scat_ang >= 160.0)  
     1   r412ss = r412_tbl- ddx+ddx2+ddx33*(scat_ang-160.)/20.
 
      if (xday > 243.0 .and. xday < 305.0) then
      if (scat_ang >= 140.0 .and. scat_ang < 160.0)
     1   r412ss = r412_tbl- ddx + ddx2*(scat_ang-140.)/20.
      if (scat_ang >= 160.0)
     1   r412ss = r412_tbl- ddx+ddx2
       endif
 
      if (r412_tbl > 12.0)  then
       if (scat_ang >= 120.0 .and. scat_ang < 170.0) 
     1    r412ss = r412_tbl- ddx + ddx3*(scat_ang-120.)/50.
       if (scat_ang >= 170.0)  
     1    r412ss = r412_tbl- ddx+ddx3+0.0*(scat_ang-170.)/10.
      endif
 
      if (xday > 151.0 .and. xday < 244.0) then
      if (r412_tbl > 10.0)  r412ss = r412ss + 0.5
      endif
 
!     E. Central Algeria, long narrow desert
      dda1 = 3.2
      dda2 = 9.5
      dda3 = 100.0
      dda4 = 0.5
      dda5 = 70.
       if (xday > 181.0 .and. xday < 305.0)  then       ! summer, Sep, Oct
        dda1 = 2.5
        dda2 = 9.5
        dda3 = 100.0
        dda4 = 0.0
        dda5 = 70.
       endif
       if (xday > 151.0 .and. xday < 182.0)  dda1 = 1.5
       if (xday > 243.0 .and. xday < 274.0)  dda1 = 0.7
       if (xday > 59.0 .and. xday < 152.0)   dda1 = 1.0
 
       if (xlat>27. .and. xlat<36. .and. xlong>2.5 .and. xlong<11.5) then
       if (r412_135 > dda2)  then
        if (scat_ang >= dda3 .and. scat_ang < 170.0) 
     1    r412ss = r412_tbl- ddx + dda1*(scat_ang-dda3)/dda5
        if (scat_ang >= 170.0)  
     1    r412ss = r412_tbl- ddx+dda1+dda4*(scat_ang-170.)/10.
       endif
       endif
 
!     Libya - Egypt 1
      if (xlat >20. .and. xlong > 14.9) then
      dda1 = r412_tbl*0.08 +0.5
      if (r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 0.8
      if (xday > 59.0 .and. xday < 91.0 .and.r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 1.5   ! Mar
      if (xday > 151.0 .and. xday < 244.0 .and.r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 0.4   ! summer
      dda2 = r412_tbl*0.0
      if (r412_tbl > 9.6)  then
       if (scat_ang >= 110.0 .and. scat_ang < 160.0)  
     1   r412ss = r412_tbl- ddx + dda1*(scat_ang-110.)/50.
       if (scat_ang >= 160.0)  
     1   r412ss = r412_tbl- ddx+dda1+dda2*(scat_ang-160.)/20.
        if (xday > 151.0 .and. xday < 244.0) then
          if (xday > 151.0 .and. xday < 182.0) then
             if (r412_tbl >10.0.and.xlat>20.0 .and. xlat<24.8.and. xlong>24.7.and.xlong<30.0) go to 211
             if (r412_tbl >10.0.and.xlat>17.5 .and. xlat<=20.0.and. xlong>24.0.and.xlong<30.0) go to 211
          endif
         dda1 = r412_tbl- ddx+r412_tbl*0.08 +0.5
         if (scat_ang >= 110.0 .and. scat_ang < 160.0) r412ss = dda1+0.8*(scat_ang-110.)/50.
         if (scat_ang >= 160.0) r412ss = dda1+0.8
         if (r412_tbl > 14.0) then
          dda1 = r412_tbl- ddx+r412_tbl*0.08 +0.5
          if (scat_ang >= 110.0 .and. scat_ang < 160.0) r412ss = dda1+1.0*(scat_ang-110.)/50.
          if (scat_ang >= 160.0) r412ss = dda1+1.0
         endif
        endif
      endif
      endif
211   continue
 
!     Libya - Egypt 2
      if (xlat >15. .and. xlong > 22.0) then
      dda1 = r412_tbl*0.08 +0.5
      if (r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 0.8
      if (xday > 59.0 .and. xday < 91.0 .and.r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 1.5   ! Mar
      if (xday > 151.0 .and. xday < 244.0 .and.r412_tbl > 12.0) dda1 = r412_tbl*0.08 + 0.4   ! summer
      dda2 = r412_tbl*0.0
      if (r412_tbl > 9.6)  then
       if (scat_ang >= 110.0 .and. scat_ang < 160.0)  
     1   r412ss = r412_tbl- ddx + dda1*(scat_ang-110.)/50.
       if (scat_ang >= 160.0)  
     1   r412ss = r412_tbl- ddx+dda1+dda2*(scat_ang-160.)/20.
        if (xday > 151.0 .and. xday < 244.0) then
          if (xday > 151.0 .and. xday < 182.0) then
             if (r412_tbl >10.0.and.xlat>20.0 .and. xlat<24.8.and. xlong>24.7.and.xlong<30.0) go to 212
             if (r412_tbl >10.0.and.xlat>17.5 .and. xlat<=20.0.and. xlong>24.0.and.xlong<30.0) go to 212
          endif
         dda1 = r412_tbl- ddx+r412_tbl*0.08 +0.5
         if (scat_ang >= 110.0 .and. scat_ang < 160.0) r412ss = dda1+0.8*(scat_ang-110.)/50.
         if (scat_ang >= 160.0) r412ss = dda1+0.8
         if (r412_tbl > 14.0) then
          dda1 = r412_tbl- ddx+r412_tbl*0.08 +0.5
          if (scat_ang >= 110.0 .and. scat_ang < 160.0) r412ss = dda1+1.0*(scat_ang-110.)/50.
          if (scat_ang >= 160.0) r412ss = dda1+1.0
         endif
        endif
      endif
      endif
212   continue
 
!     N. Algeria 1
      if (xlat >31.5 .and. xlong > 4. .and. xlong < 10.) then
      dda1 = r412_tbl*0.08 + 0.5
      dda2 = r412_tbl*0.05 + 0.5
      if (xday > 151.0 .and. xday < 244.0) then
      dda1 = r412_tbl*0.08 -0.3
      dda2 = 0.0
      endif
      if (xday > 243.0 .and. xday < 274.0) then
      dda1 = r412_tbl*0.08 -0.3
      dda2 = 0.0
      endif
      if (r412_tbl > 9.4)  then
       if (scat_ang >= 120.0 .and. scat_ang < 170.0) 
     1   r412ss = r412_tbl - ddx+0.3+ dda1*(scat_ang-120.)/50.
       if (scat_ang >= 170.0)  
     1   r412ss = r412_tbl- ddx+0.3+dda1+dda2*(scat_ang-170.)/10.
      endif
      endif
33    continue
 
!     N. Algeria 2
      if (xlat >30.0 .and. xlat <32.0 .and. xlong > 4.8 .and. xlong < 7.) then
      dda1 = r412_tbl*0.08 + 0.5
      dda2 = 0.0
      if (xday > 151.0 .and. xday < 244.0) then
      dda1 = r412_tbl*0.08 -0.3
      dda2 = 0.0
      endif
      if (xday > 243.0 .and. xday < 274.0) then
      dda1 = r412_tbl*0.08 -0.3
      dda2 = 0.0
      endif
      if (r412_tbl > 9.4)  then
        if (scat_ang >= 100.0 .and. scat_ang < 160.0) 
     1    r412ss = r412_tbl- ddx+0.3 + dda1*(scat_ang-100.)/60.
      if (scat_ang >= 160.0)  
     1    r412ss = r412_tbl- ddx+0.3+dda1+dda2*(scat_ang-160.)/20.
      endif
      endif
 
!     Chad - Libya border
      if (xday > 31.0 .and. xday < 60.0)  go to 36  ! use for all months except for Feb
      if (xlat >20. .and. xlat <25. .and. xlong >15.0.and. xlong <17.5) then
      if (r412_tbl > 12.0)  then
        r412ss = r412ss + 0.5
      endif
      endif
36    continue
 
      dda2 = 0.5
      if (xday > 59.0 .and. xday < 121.0) dda2 = 1.0
!     NE Mauritania 1
      if (xlat>22.5 .and. xlat<30.0.and. xlong>-11.0.and. xlong< -5.) then
      if (r412_tbl > 9.)  then
        r412ss = r412ss + dda2
        if (xday > 151.0 .and. xday < 274.0)  r412ss = r412_tbl-1.5      ! summer
      endif
      endif
 
      dda2 = 0.5
      if (xday > 59.0 .and. xday < 121.0) dda2 = 1.0
!     NE Mauritania 2
      if (xlat>22.5 .and. xlat<26.0.and. xlong>-13.0.and. xlong< -11.001) then
      if (r412_tbl > 9.)  then
        r412ss = r412ss + dda2
        if (xday > 151.0 .and. xday < 274.0)  r412ss = r412_tbl-1.5      ! summer
      endif
      endif
 
      dda2 = 0.0
      if (xday > 59.0 .and. xday < 121.0) dda2 = 1.0
!     NE Mauritania 3
      if (xlat>=20.0 .and. xlat<29.0.and. xlong< -12.5) then
      if (r412_tbl > 10.5)  then
        r412ss = r412ss + dda2
        if (xday > 151.0 .and. xday < 274.0)  r412ss = r412_tbl-1.5      ! summer
      endif
      endif
 
      dda2 = 0.0
      if (xday > 59.0 .and. xday < 121.0) dda2 = 1.0
!     NE Mauritania 3
      if (xlat>15.0 .and. xlat<20.0 .and. xlong< -14.9) then
      if (r412_tbl > 10.5)  then
        r412ss = r412ss + dda2
        if (xday > 151.0 .and. xday < 244.0)  r412ss = r412_tbl-1.5      ! summer
        if (xday >= 244.0 .and. xday < 274.0)  r412ss = r412ss-1.5
      endif
      endif
 
!     Lake Chad
      if (xlat >10.0 .and. xlat <21.0 .and. xlong > 10.0 .and. xlong < 20.0) then
      if (r412_tbl > 12.)  then
        if (xday > 181.0 .and. xday < 244.0)  r412ss = r412_tbl         ! summer
      endif
      endif
 
!     Morocco 1
      if (xlat>30.0 .and. xlat<36.0 .and. xlong<= -7.5) then
      if (r412_tbl > 5.8)  then
        if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1    r412ss = r412_tbl + 4.0*(scat_ang-100.)/75.
      if (scat_ang >= 175.0)
     1    r412ss = r412_tbl+4.0+0.0*(scat_ang-175.)/5.
      endif
      endif
 
!     Morocco 2
      if (xlat>31.5 .and. xlat<36.0 .and. xlong<= -5.0.and. xlong> -7.5) then
      if (r412_tbl > 5.8)  then
        if (scat_ang >= 100.0 .and. scat_ang < 175.0)
     1    r412ss = r412_tbl + 4.0*(scat_ang-100.)/75.
      if (scat_ang >= 175.0)
     1    r412ss = r412_tbl+4.0+0.0*(scat_ang-175.)/5.
      endif
      endif
 
 
      if (xday > 151.0 .and. xday < 244.0)  go to 37   ! winter, spring
 
      dda2 = 0.5
      if (xday > 243.0 .and. xday < 274.0)  dda2 =  0.0     ! Sept
!     Central Algeria
      if (xlat >22.5 .and. xlat <30. .and. xlong > -4.0.and. xlong < 1.) then
      if (r412_tbl > 10.)  then
        r412ss = r412ss + dda2
      endif
      endif
 
!     Mali - Algeria border
      if (xlat>20. .and. xlat<22.501.and. xlong>-2.0.and. xlong< 2.) then
      if (r412_tbl > 10.)  then
        r412ss = r412ss + dda2-1.0
      endif
      endif
 
      if (xday > 243.0 .and. xday < 274.0) then
      if (xlat>16. .and. xlat<22.501.and. xlong>-12.5.and. xlong< 2.) then
      if (r412_tbl > 12.)  then
        r412ss = r412ss - 1.0
      endif
      endif
      endif
 
37    continue
 
      if (xday > 59.0 .and. xday < 244.0) then    ! summer
        if (xphi > 90.0) then
        ddx = 0.0
        ddx33 = 0.6
        if (r412_tbl > 12.0) ddx33 = 0.9
        if (xthet < 70.0) ddx = xthet * ddx33/70.
        if (xthet >= 70.0) ddx = ddx33
        r412ss = r412ss + ddx
        endif
 
      if (xday > 59.0 .and. xday < 152.0) go to 38
        if (r412_tbl > 9.6) then
        if (xlat >20. .and. xlong > 14.9) go to 38
        if (xlat >15. .and. xlong > 22.0) go to 38
        endif
      endif
38    continue
      
!     -- special treatment for gzone 21 to alleviate a hot spot.
      if (gzflg == 21) then
        if (xphi > 90.0) then
        ddx = 0.0
        ddx33 = 3.0
        if (xthet < 70.0) ddx = xthet * ddx33/70.
        if (xthet >= 70.0) ddx = ddx33
        r412 = r412 + ddx
        endif
      endif
      
      if (lprint > 0) print *,'scat_ang,b,f r412ss =',r412_tbl,scat_ang,r412ss
 
      if (r412ss2 > 20.0) r412ss2 = 19.9
      if (r412ss2 < 1.0) r412ss2 = 1.0      
      if (r412 > 20.0 .and. r412 < 40.0) r412 = 19.9
      if (r412 > 40.0) go to 10
      if (r412 < 1.0)  r412 = 1.0
      if (r412ss >= 20.0 .and. r412ss < 40.0) r412ss = 19.9
      if (r412ss < 1.0)  r412ss = 1.0
      
      if (r470 >= 24.0)   r470 = 23.9
      if (r470 < 1.0)     r470 = 1.0
      if (r470ss >= 24.0) r470ss = 23.9
      if (r470ss < 1.0)   r470ss = 1.0
      
      if (lprint > 0) print *,'final r412,r412ss =',r412,r412ss
!      write (888,'(3(F10.5,1X))') xlat, xlong, r412ss
      
      if (r650 >= 47.0) r650 = 46.9
      if (r650 < 1.0)   r650 = 1.0

C     -- END DIRECT COPY SEAWIFS CODE

      if (r412.gt.12.0.and.r412.lt.80.) qa_flag(4) = 2
            
c--------------------------------------------------------
c       Cloud Screening
c--------------------------------------------------------
c-- band26 values should be fill value (-999) so most of these
c-- test are non-functional. Confirmed w/ Dr. Hsu 2013-03-01 --CB.
      if (platform .eq. 'VIIRS') then 
         if (debug) print *,'--- start cloud screening ----'
         if (debug) print *, 'band26, realbuf(7), realbuf(11): ', band26, realbuf(7), realbuf(11)
         if (band26.lt.0.0.and.realbuf(11).lt.0.0.and. realbuf(7).gt.0.0) go to 620
         if (band26.gt.0.0.and.realbuf(11).lt.0.0.and. realbuf(7).gt.0.0) go to 10
         rat = ref650/ rr412_mod
         rat1 = rr470_mod / rr412_mod
         if (debug) print *, 'rat, rat1, ref650: ', rat, rat1, ref650
         if (ref650.gt.45.0.and.rat.lt.1.4) go to 10
         if (ref650.gt.56.0.and.rat.lt.1.3) go to 10
      
         if (debug) print *, 'trflg, rr412_mod, rat1, r412: ', trflg, rr412_mod, rat1, r412
         if (trflg.gt.0.0.and.rr412_mod.gt.10.0.and.rat1.gt.1.25) go to 50

         if (r412.gt.6.0) then
           if (band26.gt.0.0.and.rr412_mod.gt.6.0) go to 10
         else
           if (band26.gt.0.0.and.rr412_mod.gt.10.0.and.rat1.lt.1.25) go to 10
         endif
         if (debug) print *, '--- end cloud screening ---'
      end if      
50    continue

c--------------------------------------------------------
c   Special Case: For thick strong-absorbing dust plume, 
c                 go to 3-channels
c--------------------------------------------------------

      if (ref650 > 32.0 .and. (gzflg >= 6 .and. gzflg /= 21 .and. gzflg /= 24 
     1    .and. gzflg /= 13 .and. gzflg /= 31) .and. Dstar1 > 1.1) then
        xxrat = 0.8
      
        abs_aero_flag = .true.    ! -- set flag for use below.

        go to 620
      endif
 80   continue
      
!     -- vegetated area in zone 1 next to DMN_Maine_Soroa site in summer. Tweak.
      if (gzflg == 1 .AND. (gdt1%month >= 6 .AND. gdt1%month <= 8)) then
        if (r650_135 < 16.0) then
          if (toa_ndvi < 0.18) then
            r412 = min(r412 + 2.0, 20.0)
            r470new = min(r470new + 2.0, 24.0)
            r650 = min(r650 + 2.0, 47.0)
          else
            r412 = max(r412 - 2.0, 1.0)
            r470new = max(r470new - 2.0, 1.0)
            r650 = max(r650 - 2.0, 1.0)
          end if
        end if
      end if
     
c--------------------------------------------------------
c     Preliminary Retrieval on AOT
c--------------------------------------------------------

c
c     For Moderate AOT, Use 412 - 470 nm Pair
c
c     Retrieving 470 nm AOT
c
      if (debug) print *, '--- starting aot470 --- '
      refl = refl3
      x3 = xphi

      tau_x470_flag = -999
      tau_x470_flag2 = -999
      
!     -- retrieve using most aerosol models -- in preparation for call to get_aot500.
      tau_x470_new_91 = -999.0 ; tau_x470_new_92 = -999.0 ; tau_x470_new_93 = -999.0
      tau_x470_new_94 = -999.0 ; tau_x470_new_95 = -999.0 ; tau_x470_new_96 = -999.0 
      tau_x470_new_995 = -999.0

      imod = 3                                ! w0 = 0.96
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1             imod,r470,tau_x470,tau_x470_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag, tau_x470)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
            
!      if (gzflg < 0 .or. gzflg > 11) go to 81    ! gzflg > 11 everywhere except N. Africa/SAP.

      imod = 3                                ! w0 = 0.96
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1              imod,r470ss,tau_x470ss,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470ss)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      if (tau_x470.lt.0.0601.and.Dstar1.gt.1.1.and.rat1.gt.1.6)
     1    then
      imod = 3                                ! w0 = 0.96
      r470ss = r470ss - 1.0
      if (r470ss .lt. 1.0) r470ss = 1.0
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1              imod,r470ss,tau_x470ss,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      go to 313
      endif
 
      if (tau_x470.lt.0.5.and.Dstar1.gt.0.98.and.rat1.gt.1.46)
     1    then
      imod = 3                                ! w0 = 0.96
      r470ss = r470ss - 1.0
      if (Dstar1.gt.1.04) r470ss = r470ss - 0.5
      if (r470ss .lt. 1.0) r470ss = 1.0
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1              imod,r470ss,tau_x470ss2,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      endif
313   continue
      
      if (r470new >= 24.0) go to 81
			if (r470new .lt. 1.0) r470new = 1.0
			
      imod = 3                                ! w0 = 0.96
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1              imod,r470new,tau_x470_new,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      tau_x470_new_96 = tau_x470_new
      
      imod = 4                                ! w0=0.995
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_995,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_995)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      
      imod = 2      
      model_frac = 0.5                        ! w0=0.95
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_95,tau_x470_flag2,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_95)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 2                                ! w0=0.94
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_94,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_94)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 1          
      model_frac = 2.0/3.0                    ! w0=0.93
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_93,tau_x470_flag2,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_93)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 1
      model_frac = 1.0/3.0                    ! w0=0.92
      call aero_470(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r470new,tau_x470_new_92,tau_x470_flag2,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_92)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod =  1                               ! w0 = 0.91
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,          
     1             imod,r470new,tau_x470_new_91,tau_x470_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_91)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
 
      if (r470 < 12.0 .and. rr470_mod > 11.0) then
        rat_470_412 = rr470_mod / rr412_mod
        if (rat_470_412 > 1.85) then
          imod = 1                              ! w0 = 0.91
          call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,                 
     1                  imod,r470new,tau_x470_new_91,tau_x470_flag2,trflg,0.0,debug) 
          if (dflag) go to 10                                                 
          status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_91)
          if (status /= 0) then
            print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
            return
          end if
        endif
      endif
 
      if (xthet > 62.0) then
        rat_470_412 = rr470_mod / rr412_mod
        if (rat_470_412 > 1.7) then
          imod = 1                              ! w0 = 0.91
          call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,             
     1                  imod,r470new,tau_x470_new_91,tau_x470_flag2,trflg,0.0,debug)
          if (dflag) go to 10
          status = handle_lut_out_of_bounds(gzflg, tau_x470_flag2, tau_x470_new_91) 
          if (status /= 0) then
            print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
            return
          end if
        end if
      end if

!--------------------------------------------------------
! 470 nm aod using dust table

      tau_x470_flag_dust = -999

!     -- retrieve using most aerosol models -- in preparation for call to get_aot500.
      tau_x470_dust_91 = -999.0 ; tau_x470_dust_92 = -999.0 ; tau_x470_dust_93 = -999.0
      tau_x470_dust_94 = -999.0 ; tau_x470_dust_95 = -999.0 ; tau_x470_dust_96 = -999.0
      tau_x470_dust_995 = -999.0

      imod = 3                                ! w0 = 0.96
      call aero_470_dust(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1                   imod,r470new,tau_x470_dust_96,tau_x470_flag_dust,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag_dust, tau_x470_dust_96)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod = 4                                ! w0=0.995
      call aero_470_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1                   imod,r470new,tau_x470_dust_995,tau_x470_flag_dust,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag_dust, tau_x470_dust_995)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod = 2
      model_frac = 0.5                        ! w0=0.95
      call aero_470_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1                   imod,r470new,tau_x470_dust_95,tau_x470_flag_dust,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag_dust, tau_x470_dust_95)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod = 2                                ! w0=0.94
      call aero_470_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1                   imod,r470new,tau_x470_dust_94,tau_x470_flag_dust,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag_dust, tau_x470_dust_94)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod = 1
      model_frac = 2.0/3.0                    ! w0=0.93
      call aero_470_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1                   imod,r470new,tau_x470_dust_93,tau_x470_flag_dust,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag_dust, tau_x470_dust_93)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod = 1
      model_frac = 1.0/3.0                    ! w0=0.92
      call aero_470_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1                   imod,r470new,tau_x470_dust_92,tau_x470_flag_dust,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag_dust, tau_x470_dust_92)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod =  1                               ! w0 = 0.91
      call aero_470_dust(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1                   imod,r470new,tau_x470_dust_91,tau_x470_flag_dust,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x470_flag_dust, tau_x470_dust_91)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      if (r470 < 12.0 .and. rr470_mod > 11.0) then
        rat_470_412 = rr470_mod / rr412_mod
        if (rat_470_412 > 1.85) then
          imod = 1                              ! w0 = 0.91
          call aero_470_dust(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1                       imod,r470new,tau_x470_dust_91,tau_x470_flag_dust,trflg,0.0,debug)
          if (dflag) go to 10
          status = handle_lut_out_of_bounds(gzflg, tau_x470_flag_dust, tau_x470_dust_91)
          if (status /= 0) then
            print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
            return
          end if
        endif
      endif

      if (xthet > 62.0) then
        rat_470_412 = rr470_mod / rr412_mod
        if (rat_470_412 > 1.7) then
          imod = 1                              ! w0 = 0.91
          call aero_470_dust(dflag,refl,x1,x2,x3,mm,nn,ll,ma, 
     1                       imod,r470new,tau_x470_dust_91,tau_x470_flag_dust,trflg,0.0,debug)
          if (dflag) go to 10
          status = handle_lut_out_of_bounds(gzflg, tau_x470_flag_dust, tau_x470_dust_91)
          if (status /= 0) then
            print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
            return
          end if
        end if
      end if
! end 470 nm aod using dust table

      if (debug) print *, '--- end aot470 ---'
81    continue

C      if (lprint > 0) print *,'tau_x470,tau_x470ss,tau_x470_new=',
C     1    tau_x470,tau_x470ss,tau_x470_new
     

c--------------------------------------------------------
c          Retrieving 412 nm AOT
c--------------------------------------------------------
 
      refl = refl1
      x3 = xphi

      tau_x412_flag = -999
      tau_x412_flag2 = -999
      
      r412new = r412                         ! transitional zone
      
      imod = 5                                ! w0 = 0.94
      if (refl < -900.) go to 605 !skip 412nm retrieval
      call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,   
     1              imod,r412,tau_x412,tau_x412_flag2,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag2, tau_x412)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
!---
      if (r412 > 11.0 .and. tau_x412 < 0.4) then
        tau_x412_91 = tau_x412 * 2.
        go to 630
      else if (r412 > 12.0) then
        tau_x412_91 = tau_x412 * 2.
        go to 630
      endif

      refl = refl1
      x3 = xphi
      if (refl < -900.) go to 605 !skip 412nm retrieval
      
      tau_x412_flag_91 = -999

      imod = 4                               ! w0 = 0.91
      call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1             imod,r412,tau_x412_91,tau_x412_flag_91,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag_91, tau_x412_91)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
630   continue
      
      tau_x412_flag = -999
      tau_x412_new_91 = -999.0 ; tau_x412_new_93 = -999.0 ; tau_x412_new_94 = -999.0
      tau_x412_new_96 = -999.0 ; tau_x412_new_995 = -999.0
      
      tau_x412ss_995 = -999.0 ; tau_x412ss2_995 = -999.0 ; tau_x412ss_96 = -999.0
      tau_x412ss_97 = -999.0 ; tau_x412ss_94 = -999.0 ; tau_x412ss_95 = -999.0
      tau_x412ss_98 = -999.0
    
      if (r412new < 1.0) go to 631 
      if (r412new >= 20.0) go to 631
      if (refl < -900.) go to 631 !skip 412nm retrieval
      
      imod = 8      ! w0=0.995
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_995,tau_x412_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_995)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
!      print *
!      print *, 'aero_412, tau_x412_new_995: ', refl, sza, xthet, xphi, imod, r412new, tau_x412_new_995, tau_x412_flag, trflg, dflag
!      print *
      
      imod = 6      ! w0=0.96
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_96,tau_x412_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_96)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 5      ! w0=0.94
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_94,tau_x412_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_94)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 4      ! w0=0.93
      model_frac = 2.0/3.0
      if (dflag) print *,'calling aero_412, model_frac: imod, model_frac: ', imod, model_frac
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_93,tau_x412_flag,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_93)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
      
      imod = 4      ! w0=0.91
      call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,      
     1             imod,r412new,tau_x412_new_91,tau_x412_flag,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412_new_91)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
            
!     -- arabian peninsula
      if (gzflg >= 6 .and. gzflg <= 11) then
        tau_x412_new = tau_x412_new_94
      end if

      if (gzflg == 6) then
        if (xday >= 60.0 .and. xday < 274.0) go to 602
        tau_x470_new = tau_x412_new
      end if
      
602   continue

      if (gzflg == 8) then
        if (xday >= 182.0 .and. xday < 274.0) go to 603
        tau_x470_new = tau_x412_new
      end if
603   continue

      if (gzflg == 7) tau_x470_new = tau_x412_new

      if (gzflg == 9) tau_x470_new = tau_x412_new

      if (gzflg == 10) then
        if (xday >= 60.0 .and. xday < 274.0) go to 605
        tau_x470_new = tau_x412_new
      end if 
605   continue

!---------------------------------------------------------------
! 412 nm aod using dust table

      tau_x412_flag_dust = -999
      tau_x412_dust_91 = -999.0 ; tau_x412_dust_93 = -999.0 ; tau_x412_dust_94 = -999.0
      tau_x412_dust_96 = -999.0 ; tau_x412_dust_995 = -999.0

      if (refl < -900.) go to 631 !skip 412nm retrieval
      imod = 8      ! w0=0.995
      call aero_412_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1              imod,r412new,tau_x412_dust_995,tau_x412_flag_dust,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag_dust, tau_x412_dust_995)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod = 6      ! w0=0.96
      call aero_412_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1              imod,r412new,tau_x412_dust_96,tau_x412_flag_dust,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag_dust, tau_x412_dust_96)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod = 5      ! w0=0.94
      call aero_412_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1              imod,r412new,tau_x412_dust_94,tau_x412_flag_dust,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag_dust, tau_x412_dust_94)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod = 4      ! w0=0.93
      model_frac = 2.0/3.0
      call aero_412_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1              imod,r412new,tau_x412_dust_93,tau_x412_flag_dust,trflg,model_frac,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag_dust, tau_x412_dust_93)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if

      imod = 4      ! w0=0.91
      call aero_412_dust(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma,
     1              imod,r412new,tau_x412_dust_91,tau_x412_flag_dust,trflg,0.0,debug)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x412_flag_dust, tau_x412_dust_91)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
        return
      end if
! end 412 nm aod using dust table

631   continue

      if ((gzflg <= 5 .and. gzflg > 0) .OR. gzflg == 27) then

        tau_x412_flag = -999
!        if (r412ss > 19.89) print *, 'out, r412ss: ', r412ss,xlat,xlong

        imod = 4                                ! w0 = 0.91
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss_91,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_91)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
        
        imod = 5                                ! w0 = 0.94
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
  
        tau_x412ss_94 = tau_x412ss
        imod = 8                                ! w0=0.995
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss_995,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_995)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
        
        imod =  7     ! w0=0.98
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss_98,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_98)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
        

        imod = 6      ! w0=0.96
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss,tau_x412ss_96,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_96)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
        
        imod = 5      ! w0=0.95
        model_frac = 0.5
        call aero_412(dflag,refl,sza,xthet,xphi,mm,nn,ll,ma, 
     1      imod,r412ss,tau_x412ss_95,tau_x412_flag,trflg,model_frac,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss_95)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if

!       General surfaces,  tau_x412ss = default (0.94)
        if (xday < 60.0 .or. xday > 334.0) then                 ! Dec, Jan, Feb
          if (tau_x412ss_97 <0.6) tau_x412ss = tau_x412ss_96     
        endif
        if (xday > 243.0 .and. xday < 274.0) then               ! Sep
          if (tau_x412ss_97 <0.6) tau_x412ss = tau_x412ss_96     
        endif
        if (xday > 273.0 .and. xday < 335.0) then               ! Oct, Nov
          if (tau_x412ss < 0.5)   tau_x412ss = tau_x412ss_96
        endif
        if (xday > 59.0 .and. xday < 121.0) then               ! Mar,Apr
          if (tau_x412ss < 0.5)   tau_x412ss = tau_x412ss_96
        endif

!      Bright Surfaces        
        ddx = 10.5
        if (xday > 151.0 .and. xday < 244.0) ddx = 12.0         ! Jun, Jul, Aug
        if (xday > 243.0 .and. xday < 274.0) ddx = 12.0         ! Sep
        if (r412_tbl > ddx ) then                
          if (tau_x412ss_97 < 0.6) tau_x412ss = tau_x412ss_96 
          if (xday > 243.0 .or. xday < 60.0) then             ! Sep,Nov,Dec,Jan,Feb
            if (tau_x412ss_995 < 0.5) tau_x412ss = tau_x412ss_995 
            if (xday > 273.0 .and. xday < 305.0) then             ! Oct
              if (tau_x412ss_97 < 0.6) tau_x412ss = tau_x412ss_96 
            endif   
            if (xlat >10.0 .and. xlat <21.0 .and. xlong > 10.0 .and. xlong < 20.0)   
     1          tau_x412ss = tau_x412ss_96                        
          endif   
          if (xday > 59.0 .and. xday < 121.0) then             ! Mar,Apr
            tau_x412ss = tau_x412ss_96
            if (xlat >10.0 .and. xlat <21.0 .and. xlong > 10.0 .and. xlong < 20.0)
     1          tau_x412ss = tau_x412ss_96
            if (Dstar1 < 1.01 .and. tau_x412ss>0.6) tau_x412ss = tau_x412ss_995
            if (Dstar1 > 1.01 .and. Dstar1 < 1.04 .and. tau_x412ss>0.6) tau_x412ss = tau_x412ss_98
            endif   
          endif   

!      Deserts near Egypt and Sudan        
        if (r412_tbl > 10.5 .and. Dstar1 > 1.04) then                
           if (xlat >20.0 .and. xlat <25.0 .and. xlong > 25.0 .and. xlong < 30.0)
     1         tau_x412ss = tau_x412ss_94
        endif   

!      Deserts near Libya and Egypt         
        ddx = 1.01
        if (xday > 151.0 .and. xday < 305.0) ddx = 1.04

        dda1 = tau_x412ss

        if (xlat >20. .and. xlong > 12.9) then               
          if (r412_tbl > 8.5) then
            if (xlat <22. .and. xlong < 17.0) go to 632
              if (Dstar1 .lt. ddx) then
              dd = (xlong - 11.9) /3.0
              if (xlong > 14.9) dd = 1.
              tau_x412ss = tau_x412ss-(tau_x412ss-tau_x412ss_995)*dd
              dda2 = tau_x412ss
              dd = (xlat - 18.0) /4.0
              if (xlat > 22.0) dd = 1.
              tau_x412ss = dda1-(dda1-dda2)*dd
              endif   
          endif   
        endif   
        if (xlat >15. .and. xlat <=20. .and. xlong > 22.0) then                
          if (r412_tbl > 8.5) then
              if (Dstar1 .lt. ddx) then
              dda2 = tau_x412ss_995
              dd = (xlat - 18.0) /4.0
              if (xlat > 22.0) dd = 1.
              if (xlat < 18.0) dd = 0.
              tau_x412ss = dda1-(dda1-dda2)*dd
          endif
        endif   
        endif   
        if (xlat >19.0 .and. xlat <=20.0 .and. xlong > 19.8) then
          if (r412_tbl > 8.5) then
              if (Dstar1 .lt. ddx) then
              dda2 = tau_x412ss_995
              dd = (xlat - 18.0) /4.0
              if (xlat > 22.0) dd = 1.
              tau_x412ss = dda1-(dda1-dda2)*dd
              endif
          endif
        endif

        if (xday > 243.0 .and. xday < 274.0) then
        if (xlat >19.0 .and. xlat <=21.0 .and. xlong > 19.8.and. xlong <23.5) then
          if (r412_tbl > 8.5) then
              if (Dstar1 .lt. ddx) tau_x412ss = tau_x412ss_995
          endif
          endif
        endif   

!     NE Mauritania 3
      if (xday > 243.0 .and. xday < 274.0) then
       if (xlat>22.5 .and. xlat<30.0.and. xlong>-11.0.and. xlong< -5.) then
       if (r412_tbl > 9.0)  tau_x412ss = tau_x412ss_94
       endif
       if (xlat>22.5 .and. xlat<26.0.and. xlong>-13.0.and. xlong< -11.001) then
       if (r412_tbl > 9.0)  tau_x412ss = tau_x412ss_94
       endif
       if (xlat>=20.0 .and. xlat<29.0.and. xlong< -12.5) then
       if (r412_tbl > 10.5)  tau_x412ss = tau_x412ss_94
       endif
       if (xlat>15.0 .and. xlat<20.0 .and. xlong< -14.9) then
       if (r412_tbl > 10.5)  tau_x412ss = tau_x412ss_94
       endif
      endif

      if (xday > 151.0 .and. xday < 274.0) go to 632         ! Jun, Jul, Aug

!     NE Mauritania 1
      if (xlat>20.0 .and. xlat<30.0 .and. xlong< -12.5) then
      if (r412_tbl > 9.)  then
          if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
      endif
      endif
 
!     NE Mauritania 2
      if (xlat>26.0 .and. xlat<30.0.and. xlong>=-12.5.and. xlong< -11.001) then
      if (r412_tbl > 9.)  then
          if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
      endif
      endif

!     Morocco 1
      if (xlat>30.0 .and. xlat<36.0 .and. xlong<= -7.5) then
      if (r412_tbl > 5.8)  then
          if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
      endif
      endif   

!     Morocco 2
      if (xlat>31.5 .and. xlat<36.0 .and. xlong<= -5.0.and. xlong> -7.5) then
      if (r412_tbl > 5.8)  then
          if (Dstar1 .lt. 1.01) tau_x412ss = tau_x412ss_995
      endif
        endif   

632     continue

        go to 635
633     continue
        if (tau_x412ss_97 < 0.6) tau_x412ss = tau_x412ss_96
        if (xday > 181.0 .and. xday < 274.0) then               ! Jul, Aug, Sep
          if (tau_x412ss_995 <0.5) tau_x412ss = tau_x412ss_995
        endif
        if (xday > 273.0 .and. xday < 335.0) then               ! Oct, Nov
          if (tau_x412ss_995 <0.5) tau_x412ss = tau_x412ss_995
        endif   
        if (xday > 59.0 .and. xday < 121.0) then               ! Mar,Apr
          if (tau_x412ss_995 <0.5) tau_x412ss = tau_x412ss_995
        endif

635   continue
      endif   
       
      if (gzflg >= 6 .and. gzflg <= 11) then

        tau_x412_flag = -999
        if (r412ss2 > 20.0) print *, 'out, r412ss2: ', r412ss2
        if (r412ss2 < 1.0) print *, 'less, r412ss2: ', r412ss2

        imod = 5                                ! w0 = 0.94
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,      
     1               imod,r412ss2,tau_x412ss2,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss2)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if

        imod =  7                               ! w0=0.98
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1               imod,r412ss2,tau_x412ss2_98,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss2_98)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if

        imod =  8                               ! w0=0.995
        call aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1               imod,r412ss2,tau_x412ss2_995,tau_x412_flag,trflg,0.0,debug)
        if (dflag) go to 10
        status = handle_lut_out_of_bounds(gzflg, tau_x412_flag, tau_x412ss2_98)
        if (status /= 0) then
          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
          return
        end if
      end if
 
c  If surface is vegetated, derive the surface reflectivity at 490 nm and 670 nm
c    from 865 nm surface LER and NDVI.  Specifically exclude Sahara and Arabia via gzflg
c    flag.  Exclude areas with AERONET-based BRDF info via brdf_flag == 0.
c--------------------------------------------------------------------------------------------------
 637  continue
 
!     -- set up our NDVI threshold over which the veg code can run.
!     -- if we failed to get surface info above (sr_fail_flag == .true.),
!     -- let veg try no matter the NDVI (unless superseded by another test
!     -- later on.
      ndvi_thold = 0.1
      if (sr_fail_flag .eqv. .true.) then
        ndvi_thold = -999.0
      end if
      
      if (lc == 6) then
        ndvi_thold = 0.2
      end if

!     -- increase NDVI threshold over high elevation or barren areas
!     -- in the United States (13)
      if (gzflg == 13) then   
        if (px_elev >= 500.0 .OR. lc == 6) then
            ndvi_thold = 0.3
        end if
      endif
      
!     -- increase NDVI threshold over high elevation areas in North and South America.
!     -- this covers Central America which does not currently have a region assigned to it.
      if (gzflg < -900 .AND. (xlat < 30.0 .AND. xlong < -30.0)) then
        if (px_elev >= 500) then
          ndvi_thold = 0.3
        end if
      end if

!     -- decide if we need to run the veg retrieval.
      do_veg = .false.
      if (toa_ndvi >= ndvi_thold .AND. brdf_flag /= 0) then
        do_veg = .true.
      end if

!     -- override do_veg from above to exclude vegetated algorithm from the following
!     -- zones.
      if (gzflg > 0 .AND. gzflg <= 11) do_veg = .false.    ! N. Africa, Arabia
      if (gzflg == 15) do_veg = .false.                   ! India, Kanpur region 
      if (gzflg == 20) do_veg = .false.                   ! India, Thar desert
      if (gzflg == 19) do_veg = .false.                   ! India, Pune region 
      if (gzflg == 34) do_veg = .false.                   ! India, high elevation
      if (gzflg == 21) do_veg = .false.                   ! Asia dust transport region
      if (gzflg == 23) do_veg = .false.                   ! Iraq, Iran region
      if (gzflg == 24) do_veg = .false.                   ! Taklimakan desert
      if (gzflg == 27) do_veg = .false.                   ! Sahel transition region

!     -- some special logic over China in Winter.
      if (gzflg == 16 .AND. (xday >= 335 .OR. xday < 60)) then
        if (toa_ndvi > 0.2) then
          if (lc == 1 .OR. lc == 2 .OR. lc == 4) then
            do_veg = .true.
          end if
          if (lc == 3 .AND.  xlat < 25.0) then
	          do_veg = .true.
          end if
        end if
      end if

!     -- some special logic over India, only test NDVI regardless of
!     brdf_flag, adjustment of ndvi_thold above won't trigger veg
!     retrievals due to the brdf_flag condition, 3 January 2018 JLee
      if ((gzflg == 15 .or. gzflg == 19) .and. toa_ndvi > 0.4) then
        do_veg = .true.
      end if

!     -- special logic over N. America urban areas
      if (gzflg == 13 .and. toa_ndvi > 0.4) then 
        do_veg = .true.
      end if

!!     -- special logic over Sahel region
!      if ((gzflg == 26 .or. gzflg == 27) .and. toa_ndvi > 0.4) then
!        do_veg = .true.
!      end if    
      
			if (debug) print *, "ndvi, gzflg, brdf_flag, sr_fail_flag, do_veg: ", toa_ndvi, gzflg, brdf_flag, sr_fail_flag, do_veg

c     -- skip pixel if no surf. reflc. value and not suitable for vege. retrieval.
      if ((do_veg .eqv. .false.) .AND. (sr_fail_flag .eqv. .true.)) go to 10
      
      if (do_veg) then ! .OR. oob_samer) then
        alg_flag = 10 
        old_alg_flag = 1 
c      if (ndvi670 >= 0.1 .AND. (gzflg < -900 .OR. gzflg > 11) .AND. brdf_flag /= 0) then
c        print *, 'doing veg retrieval: ', itmp, jtmp, alg_flag
c       Get current season and tweak for swf_aero_veg input.
c        iopss = season_from_doy(yr, doy)  
c        iopss = iopss - 1
c        if (iopss == 0) iopss = 4

c... *******************************************************************
c... *******************************************************************
c... *******************************************************************
c-----------------------------------------------------------------------
c... do retrieval over vegetated surfaces (NDVI=>0.1)
c-----------------------------------------------------------------------
c...  imod=2  ! ssa_490=0.94 --> You may open this line and remove
c...                             do loop  "do 2500 imod=1,4", or use
c...  --> "call swf_aero_veg(nvalx470,nvalx650,iopss,2,sza,xthet,"
c...  ------------------------------------------------------
c        gdt1 = gregorian_from_doy(yr,doy)
c        ioprg=0   ! region index initialization @@new@@
c        idx=int((xlong-(-180.0))/0.10)+1    ! @@new@@
c        idy=int((xlat-(-90.0))/0.10)+1     ! @@new@@
c        if(idx.ge.1.and.idx.le.3600.and.idy.ge.1.and.idy.le.1800) then   ! @@new@@
c           ioprg=veg_regions(idx,idy)   ! @@new@@
c        else          ! @@new@@
c           ioprg=0    ! @@new@@
c        endif         ! @@new@@
c        imod = 2
c        select case (ioprg)
c          case (6)                        ! S. Africa
c            select case (gdt1%month)  
c              case (6:11)       ! Jun-Nov use more absorbing model (0.89).
c                imod = 1
c              case default
c                imod = 2
c            end select
c          case default
c            imod = 2
c        end select
        
c        tau_x470_flag = -999
c        tau_x650_flag = -999
        
        tau_x470_flag_veg = -999
        call find_v_veg(gdt1%month,season,realbuf,tmpvg,
     1       r412sv_veg,r470sv_veg,gzflg,outbufvg,tau_x470_flag_veg, platform)
        if (outbufvg(7) < 0.0 .and. tau_x470_flag_veg == 1) go to 805

c       translate outbufvg back to local variables.
        do i=1,3
          xtau(i) = outbufvg(i)
          ssa(i) = outbufvg(i+3)
          outbuf(i+3) = ssa(i)
        enddo
        tau550 = outbufvg(7)
        alg_flag = 11
        alpha = outbufvg(8)
        r412 = -999.0
        r470 = outbufvg(11)
        r650 = outbufvg(12)
        tau_x650_flag = outbufvg(10)
        xthet = outbufvg(13)
        scat_ang = outbufvg(14)
        sfc_typ = outbufvg(15)
        
        if (debug) print *, "veg final: aot550, aot: ", tau550, xtau(1), xtau(2), xtau(3)

        ! jlee test 20170726, use 2.2 um surface reflectance for this
        ! condition, Dr. Hsu found a high bias over bright vegetated
        ! surfaces in CONUS, Africa, etc.
        if ((outbufvg(3).gt.0.2.and.alpha.lt.-0.5).or.
     &      (outbufvg(3).gt.0.2.and.(realbuf(23)/realbuf(6)).lt.0.7).or.
     &      (outbufvg(2).gt.0.4.and.(realbuf(23)/realbuf(6)).lt.0.7)) then
           tau550 = tau_x470sv96
           alg_flag = 12
           old_alg_flag = 0
        endif

c     Return fill for 412 nm, no dark target retrieval yet.
c        xtau(1) = -999.0
c        tau_x412_flag2 = -999
c        if (lprint > 0) print *, 'after veg: ', i, j, xtau(2), xtau(4), tau550, & 
c        & alpha, tau_x470_flag, tau_x650_flag
c        
c!       -- reset AOT to zone-specific value if we ran off the LUT.
c        status = handle_lut_out_of_bounds(gzflg, tau_x470_flag, xtau(2))
c        if (status /= 0) then
c          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
c          return
c        end if
c        status = handle_lut_out_of_bounds(gzflg, tau_x650_flag, xtau(4))
c        if (status /= 0) then
c          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
c          return
c        end if
c        status = handle_lut_out_of_bounds(gzflg, tau_x470_flag, tau550)
c        if (status /= 0) then
c          print *, "ERROR: Failed to check/reset AOT out of bounds condition: ", status
c          return
c        end if
        
      
c        write (500,'(7(F10.5,1X))') xlat, xlong, r470, r650, xtau(2), xtau(4), tau550
        
c        if (lprint > 0) print *, 'after LUT check, veg: ', i, j, xlat, xlong, ndvi670, &
c        &   r470, r650, xtau(2), xtau(4), tau550, alpha
        
c!     Check consistency
        if (xtau(2) < 0.0 .or. xtau(3) < 0.0 .or. tau550 < 0.0) then
          xtau(1)       = -999.0
          xtau(2)       = -999.0
          xtau(3)       = -999.0
          ssa(1)        = -999.0
          ssa(2)        = -999.0
          ssa(3)        = -999.0
          tau550        = -999.0
          alpha         = -999.0
          r412          = -999.0
          r470          = -999.0
          r650          = -999.0
          tau_x470_flag = -999
          tau_x650_flag = -999
          alg_flag      = -999
          old_alg_flag  = -999
          outbufvg(:)   = -999.0 
        end if
        
        if (debug) print *, 'final: ', xlat, xlong, xtau(2), xtau(3), tau550, alpha
c        write(888,'(7(F12.5), I3)') xlat, xlong, xtau(2), xtau(4), r470, r650, z_ndvi, proc_flag(i,j)
c-----------------------------------------------------------------------
c... end of retrieval over vegetated surfaces
c-----------------------------------------------------------------------
c... *******************************************************************
c... *******************************************************************
c... *******************************************************************

c         r470 = get_ssr_490(yr, doy, ler650(i,j)*100.0, ler865(i,j)*100.0, status)
c        if (status < 0) then
c          print *, "ERROR: Unable to derive vegetated surface "//
c     1      "reflectivity for 490 nm: ", status
c          stop
c          cycle
c        end if
c       
c        r650 = get_ssr_670(yr, doy, ler650(i,j)*100.0, ler865(i,j)*100.0, status)
c        if (status < 0) then
c          print *, "ERROR: Unable to derive vegetated surface "//
c     1      "reflectivity for 670 nm: ", status
c          stop
c          cycle
c        end if
        go to 865
      end if
      
      alg_flag = 0
      old_alg_flag = 0

c      if (lprint > 0) print *,'tau_x412,tau_x412_new,tau_x412ss,tau_x412ss2=',  
c     1    tau_x412,tau_x412_new,tau_x412ss,tau_x412ss2
c--------------------------------------------------------
c          Retrieving 412 nm SSA (412 - 470 nm)
c--------------------------------------------------------
      if (platform .eq. 'AHI') go to 620
      if (Dstar1.gt.1.05.and.rat1.gt.1.6) go to 620
      if (tau_x470.lt.0.2.and.tau_x470.gt.0.0) go to 805
      if (rr412_mod.gt.20.0.and.tau_x470.lt.0.0) go to 620
      if (tau_x470.lt.0.0) go to 620

      refl = refl1
      x3 = xphi
      tau_x = tau_x470
      w0_x = -999.         

      if (tau_x  >=  3.5) tau_x = 3.5-0.0001     ! search2 dies if tau_x >= tau(10), 3.5
      call aero_412_abs(dflag,refl,x1,x2,x3,mm,nn,ll,
     1    r412,tau_x,w0_x)
      if (dflag) go to 10

      w0_int = w0_x

      if (w0_x.lt.0.0) go to 10
 
      w0_int_470 = w0_x +(0.976 -w0_x)*(470.-412.)/(650.-412.)
      if (w0_int_470.lt.0.0) w0_int_470 = -999.

c--------------------------------------------------------
c          Retrieving 650 nm AOT
c--------------------------------------------------------
      if (Dstar1 > 1.1 .and. xlong < -30.0) go to 620
      if (xthet.gt.60.0.and.xphi.lt.90.0.and.tau_x470.gt.0.5)
     1    go to 620
      if (xlong < -30.0) then	! mostly N.America
        if (r650.lt.30.0.and.tau_x470.lt.0.7.and.tau_x470.gt.0.0)
     1    go to 805
      else
        if (r650.lt.30.0.and.r650.gt.15.0.and.tau_x470.lt.0.7.and.tau_x470.gt.0.0)
     1    go to 805
     	end if
c      if (r650.lt.30.0.and.tau_x470.lt.0.7.and.tau_x470.gt.0.0)
c     1    go to 805
      if (r650.ge.30.0.and.tau_x470.lt.1.0.and.tau_x470.gt.0.0)
     1    go to 805
      if (scat_ang.gt.165.0.and.tau_x470.lt.1.8) go to 805

620   continue
      refl = refl6
      x3 = xphi

      tau_x650_flag = -999

      call aero_650(dflag,refl,x1,x2,x3,mm,nn,ll,ma,    
     1    r650,tau_x650,tau_x650_flag,tau_x470_flag2,   
     1    tau_x412, tau_x470,tau_x412_flag_91,trflg)
      if (dflag) go to 10
      status = handle_lut_out_of_bounds(gzflg, tau_x650_flag, tau_x650)
      if (status /= 0) then
        print *, "ERROR: Failed to check/reset AOT out of bounds 
     1            condition: ", status
        return
      end if

!--------------------------------------------------------
!          Retrieving 412 and 470 nm SSA
!--------------------------------------------------------
      if (tau_x650 < 0.0) go to 805
      if (platform .eq. 'AHI') go to 805
      
c      if (tau_x650 > 3.5) go to 805

!--- 412 nm

      refl = refl1
      x3 = xphi
      tau_x = tau_x650
      if (rat1.lt.1.8.and.tau_x470ss.gt.0.8) tau_x = tau_x470ss
      w0_x = -999.
      if (tau_x > 3.5) tau_x = 3.5
      if (tau_x  > 3.5) then
         goto 805
      endif
      if (tau_x  >=  3.5) tau_x = 3.5-0.0001 

      call aero_412_abs(dflag,refl,x1,x2,x3,mm,nn,ll,   
     1                 r412,tau_x,w0_x)

      if (dflag) go to 10
      w0_x412 = w0_x
      if (w0_x < 0.0) go to 10
 
!--- 470 nm 
 
      refl = refl3
      x3 = xphi
      tau_x = tau_x650
      w0_x = -999.
      if (tau_x  >=  3.5) tau_x = 3.5-0.0001 
! ADDED BY COREY
      if (r470  >  24.0) go to 805
! END ADDED BY COREY
      call aero_470_abs(dflag2,refl,x1,x2,x3,mm,nn,ll,  
     1                 r470,tau_x,w0_x470)
      if (dflag2) go to 10
 
      if (w0_x470 < 0.0) w0_x470 = -999.

 805  continue
!
!     -- Selecting Models
!

      call aero_mod (tau_x412,tau_x470,tau_x650,tau_x412_91,aot_mod)

      if (xlong < -30.0 .and. xlat >= 10.0 .and. Dstar1 > 1.1 .and. 
     1    realbuf(23) > 12.0 ) then
        if ((gzflg == 13 .or. gzflg == 31) .and. realbuf(22)/realbuf(6) < 0.3) go to 807
          tau_x412 = tau_x650
          tau_x470 = tau_x650
          aot = tau_x650
          go to 860
      endif
 807  continue

      if (tau_x412 > 0.0 .and. tau_x470 > 0.0) xxrat = tau_x412 / tau_x470
      if (tau_x650 > 0.0 .and. tau_x470 > 0.0) xxrat2 =  tau_x470 / tau_x650

!     -- just a hack to avoid AE changes now that AOT values can be extrapolated
!     -- past 3.5 whereas before they couldn't and would be truncated to 3.5.
      if (tau_x412 > 0.0 .AND. tau_x470 > 0.0) xxrat = min(tau_x412, 3.5)/min(tau_x470, 3.5)
      if (tau_x650 > 0.0 .and. tau_x470 > 0.0) xxrat2 = min(tau_x470, 3.5)/min(tau_x650, 3.5)
      
      if (xxrat < 0.0) then
         alpha   = -999.
         go to 806
      end if

      dd      = alog(412./490.)
      alpha   = alog(xxrat)
      alpha   = -1.*alpha/dd
      
 806  continue

      dd    = ref650 / rr470_mod
      dd1   = ref650 / rr412_mod
      dd2   = rr470_mod / rr412_mod
      sfcdd = r650 / r412
      view  = xthet

      if (r650 < 30.0 .and. r650 > 9.0 .and. sfcdd < 2.4)  go to 850
      if (r650 < 30.0 .and. r650 > 9.0 .and. sfcdd > 2.4 .and. 
     1   sfcdd < 3.2 .and. trflg > 0.)  go to 850

      if (platform .eq. 'VIIRS') aot = aot_mod(3)
      if (platform .eq. 'AHI') aot = aot_mod(6)
      
      if (px_elev > 500.0 .and. r650 > 10.0 .and. xlong < -30.0) aot = tau_x412
      if (px_elev > 500.0 .and. 1./rat_650_470 < 0.6 .and. tau_x650 > 0.6) then
      if (xlong < -30.0 .and. xlat < 10.0) go to 10
      if (xlong < -30.0 .and. xlat >= 10.0 .and. lc /= 6) go to 10
      endif

      if (tau_x412_91 < 0.0 .and. tau_x650 > 0.0 .and. r650 < 30.0) aot = aot_mod(1)
      if (tau_x412_91 < 0.0 .and. r650 >= 30.0 .and. tau_x412 > 0.0) aot = aot_mod(4)
      if (aot < 0.3 .and. aot > 0.0) go to 860
      if (w0_int < 0.92 .and. aot > 0.0) aot = aot_mod(2)
      if (r650 < 30.0 .and. tau_x412_91 > 1.2                       
     1      .and. tau_x412_91 > tau_x650 .and. tau_x650 > tau_x470  
     1      .and. w0_int >=  0.939) aot = aot_mod(1)
      if (tau_x412_91 < 0.0 .and. tau_x470 < 0.0 .and. tau_x650 > 0.0) aot = aot_mod(1)
      if (aot < 0.0 .and. tau_x650 > 0.0) aot = aot_mod(1)
      if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 > 1.35) alpha = -0.4
      if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 <= 1.35) alpha = 1.8
      if (tau_x412 < 0.0 .and. tau_x470 > 0.0 .and. dd1 <= 1.4) alpha = 1.8
      if (aot < 0.0) go to 10
      go to 860

850   continue

      aot = aot_mod(2)

      if (tau_x650 > aot .and. xlong > -30.0) aot = aot_mod(1)
      if (px_elev > 500.0 .and. r650 > 10.0 .and. xlong < -30.0) aot = tau_x412
      if (px_elev > 500.0 .and. 1./rat_650_470 < 0.6 .and. tau_x650 > 0.6) then
      if (xlong < -30.0 .and. xlat < 10.0) go to 10
      if (xlong < -30.0 .and. xlat >= 10.0 .and. lc /= 6) go to 10
      endif
      
      if (platform .eq. 'VIIRS') then
       if (aot < 0.0) go to 10
!      modified by CH 2/18/2020
       if (gzflg == 12 .and. 1/dd2 < 0.95 .and. 1/dd > 0.8 .and. aot > 1.0) go to 851

       if (gzflg /= 16 .and. tau_x650 > 1.0 .and. dd2 < 1.1 .and. r650 > 15.) go to 10

       if (tau_x412 > 0.0 .and. tau_x650 > 1.9 .and. dd2 < 1.2) go to 10
       if (tau_x412 > 0.0 .and. tau_x470 > 1.9 .and. tau_x650 > 0.4    
     1     .and. dd2 < 1.2 .and. gzflg /= 31 .and. gzflg /= 13) go to 10
       if (gzflg /= 16 .and. tau_x412 > 0.0 .and. aot > 1.0 .and. dd2 < 1.1 .and.
     1    tau_x650 > 0.4 .and. r650 > 13.0 .and. r650 <= 15.0 .and. gzflg /= 31
     1    .and. gzflg /= 13) go to 10

       if (tau_x412 > 0.0 .and. r412 > 18.0 .and. tau_x650 > 0.4   
     1     .and. dd < 1.4 .and. r650 <= 15.) go to 10

       if (gzflg /= 16 .and. tau_x412 > 1.0 .and. xxrat > 1.2 .and. r650 > 15.0
     1    .and. tau_x650 > 0.4 .and. gzflg /= 31 .and. gzflg /= 13) go to 10

       if (tau_x650 > 1.5 .and. xxrat > 1.05 .and. tau_x412 < tau_x650 .and. tau_x412 > 0.0) go to 10
       if (tau_x650 > 1.5 .and. tau_x412 < 0.0 .and. tau_x470 < 0.0 .and.  
     1   dd2 < 1.6 .and. w0_x412 > 0.96) go to 10
       if (tau_x650 > 1.2 .and. tau_x412 < 0.0 .and. xxrat2 > 1.2 .and. w0_x412 > 0.97) go to 10
      
851   continue
       if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 > 1.35)   alpha = -0.4
       if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 <= 1.35)  alpha = 1.8
       if (tau_x412 < 0.0 .and. tau_x470 > 0.0 .and. dd1 <= 1.4)   alpha = 1.8

       if (tau_x650 > 2.0 .and. tau_x412 < 0.0 .and. tau_x470 > 2.0      
     1     .and. tau_x470 < 3.0 .and. xxrat2 < 1.2 .and. xxrat2 > 1.0) 
     1    go to 10 
      
       if (tau_x650 > 2.0 .and. tau_x412 < 0.0 .and. tau_x470 > 3.0  
     1     .and. xxrat2 < 1.45 .and. xxrat2 > 1.0) go to 10

       if (tau_x412 < 0.0 .and. tau_x470 > 0.0 .and. xxrat2 > 2.)    
     1    aot = aot_mod(5)
       if (tau_x412 > 1.5 .and. tau_x470 > 0.0 .and. tau_x650 < 0.3) 
     1    aot = aot_mod(5)
       if (alpha > 1.0 .and. tau_x470 > 0.2) aot = aot_mod(5)
       if (alpha > 1.0 .and. tau_x470 <= 0.2) aot = aot_mod(5)*0.75
       if (tau_x412 < 0.0 .and. tau_x470 < 0.0 .and. dd1 <= 1.35)  aot = aot_mod(6)
      end if 

860   continue
      tau550 = aot
      !brdf_flag == 0  : AERONET brdf
      !brdf_flag == 1  : BRDF from table
      !alg_flag == 1  : AERONET brdf
      !alg_flag == 2  : BRDF from table
      !alg_flag == 3  : min SR from db
      
      if (brdf_flag > -1) then   
        old_alg_flag = 0   
        alg_flag = 1
        alg_flag = alg_flag + brdf_flag
        if (alg_flag == 1 .and. min_flag ==1) alg_flag  = 3
      endif
      if (brdf_flag < 0) then   
        old_alg_flag = 0   
        alg_flag = 2
        alg_flag = alg_flag + min_flag
      endif      
      

c     --- Additional Cloud Screening
!      if ((gzflg <= 6 .OR. gzflg == 27) .and. (gzflg > 0 .and. gzflg /= 2)) then
!      if (tau_x412ss > 0.6.and.rat1 <1.4.and.Dstar1 < 1.0.and.r412_tbl < 10.) go to 10
!      endif
      
      !     -- skip all of the stuff below, just use aot.
      if (abs_aero_flag .eqv. .true.) then
        goto 864
      end if

      ! Overwrite DB AOD (NDVI < NDVI_thold and no AERONET BRDF)
      ! Default, not using default because we do not know the impact on other regions
!      aot = tau_x470_new_96
!      if (aot < 0.7) then
!        model_frac = 1.0-aot/0.7
!        aot = tau_x470_new_96*model_frac + tau_x470_new_94*(1.0-model_frac)
!      else
!        aot = tau_x470_new_94
!      endif

      ! North America, test needed 
!      if (tau_x470sv_94 > 0.0 .and. tau_x470sv_96 > 0.0 .and. tau_x470sv_995 > 0.0) then 
!        if (regid == 1) then  
!          aot = tau_x470sv_995
!          if (aot < 0.7) then 
!            model_frac = 1.0-aot/0.7
!            aot = tau_x470sv_995*model_frac + tau_x470sv_96*(1.0-model_frac)
!          else
!            aot = tau_x470sv_96
!          endif
!        endif
!      endif

      ! Europe
      if (regid == 8 .and. gzflg < 0) then
        aot = tau_x470_new_995
        if (aot < 0.7) then
          model_frac = 1.0-aot/0.7
          aot = tau_x470_new_995*model_frac + tau_x470_new_96*(1.0-model_frac)
        else
          aot = tau_x470_new_96
        endif
      endif

      ! Korea and Japan
      if (regid == 9 .and. gzflg < 0) then 
        aot = tau_x470_new_96
        if (aot < 0.7) then
          model_frac = 1.0-aot/0.7
          aot = tau_x470_new_96*model_frac + tau_x470_new_94*(1.0-model_frac)
        else
          aot = tau_x470_new_94
        endif
      endif
      tau550 = aot

      ! -- force Tinga_Tingana scheme for AOT @ 500nm over barren surfaces even when gzflag is undefined.
      if (lc == 6 .AND. (gzflg /= 16 .AND. gzflg /= 2) .AND. gzflg /= 14 .AND. gzflg /= 21 .AND.  
     1   gzflg /= 10 .AND. gzflg /= 20 .AND. gzflg /= 30 .and. gzflg/=31) then  
        tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 12, lc, stdv,       
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
     1                     tau_x412_dust_91, tau_x412_dust_93, tau_x412_dust_94,
     1                     tau_x412_dust_96, tau_x412_dust_995,tau_x470_dust_91,
     1                     tau_x470_dust_92, tau_x470_dust_93, tau_x470_dust_94,
     1                     tau_x470_dust_95, tau_x470_dust_96, tau_x470_dust_995,
     1                     alpha, status, (lprint > 0),platform)
        if (status /= 0) then     
          tau550 = aot
        endif
      endif
       
      if (gzflg > 0) then
!       -- get AOT @ 500nm based on AERONET regions and case studies 
        tau550 = get_aot500(xlat, xlong, px_elev, scat_ang, season, toa_ndvi, gzflg, lc, stdv, 
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,        
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,     
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,      
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,     
     1                     tau_x412_dust_91, tau_x412_dust_93, tau_x412_dust_94,
     1                     tau_x412_dust_96, tau_x412_dust_995,tau_x470_dust_91,
     1                     tau_x470_dust_92, tau_x470_dust_93, tau_x470_dust_94,
     1                     tau_x470_dust_95, tau_x470_dust_96, tau_x470_dust_995,
     1                     alpha, status, (lprint > 0),platform)
        if (status /= 0) then! revert to manual assignment       
          tau550 = aot
        endif

!       -- force Kanpur (gzflg=15, lc=3) AOT models over gzflg 20 (Thar Desert)
C       if (gzflg == 20) then
C         tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 15, 3, stdv,       
C    1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
C    1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
C    1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
C    1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
C    1                     alpha, status, (lprint > 0),platform)
C         if (status /= 0) then     
C           tau550 = aot
C         endif
C       endif
        
!       -- force Beijing (gzflg=16, low elev) AOT models over gzflg 21 (Asian desert)
        if (gzflg == 21) then
          tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 16, 4, stdv,       
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
     1                     tau_x412_dust_91, tau_x412_dust_93, tau_x412_dust_94,
     1                     tau_x412_dust_96, tau_x412_dust_995,tau_x470_dust_91,
     1                     tau_x470_dust_92, tau_x470_dust_93, tau_x470_dust_94,
     1                     tau_x470_dust_95, tau_x470_dust_96, tau_x470_dust_995,
     1                     alpha, status, (lprint > 0),platform)
          if (status /= 0) then     
            tau550 = aot
          endif
        endif
        
!       -- force SW Asia (Pakistan, Iraq, etc) (gzflg=23) to use Tinga_Tingana
        if (gzflg == 23) then
          tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 12, 6, stdv,       
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
     1                     tau_x412_dust_91, tau_x412_dust_93, tau_x412_dust_94,
     1                     tau_x412_dust_96, tau_x412_dust_995,tau_x470_dust_91,
     1                     tau_x470_dust_92, tau_x470_dust_93, tau_x470_dust_94,
     1                     tau_x470_dust_95, tau_x470_dust_96, tau_x470_dust_995,
     1                     alpha, status, (lprint > 0),platform)
          if (status /= 0) then     
            tau550 = aot
          endif
        endif
        
!       -- force Tinga_Tingana scheme for AOT @ 500nm over barren surfaces
        if (lc == 6 .AND. (gzflg /= 16 .AND. gzflg /= 2) .AND. gzflg /= 14 .AND. gzflg /= 20 .AND. 
     1      gzflg /= 21 .AND. gzflg /= 23 .AND. gzflg /= 24 .AND. gzflg /= 10 .AND. gzflg /= 30 .and. gzflg/=31) then
          tau550 = get_aot500(xlat, xlong, 0.0, scat_ang, season, toa_ndvi, 12, lc, stdv,       
     1                     tau_x412_new_91, tau_x412_new_93, tau_x412_new_94,          
     1                     tau_x412_new_96, tau_x412_new_995, tau_x470_new_91,         
     1                     tau_x470_new_92, tau_x470_new_93, tau_x470_new_94,          
     1                     tau_x470_new_95, tau_x470_new_96, tau_x470_new_995,         
     1                     tau_x412_dust_91, tau_x412_dust_93, tau_x412_dust_94,
     1                     tau_x412_dust_96, tau_x412_dust_995,tau_x470_dust_91,
     1                     tau_x470_dust_92, tau_x470_dust_93, tau_x470_dust_94,
     1                     tau_x470_dust_95, tau_x470_dust_96, tau_x470_dust_995,
     1                     alpha, status, (lprint > 0),platform)
          if (status /= 0) then     
            tau550 = aot
          endif
          
          if ((gzflg <= 6 .OR. gzflg == 27) .and. (gzflg > 0 .and. gzflg /= 2)) then
          tau550 = tau_x412ss    
          if (Dstar1 > 1.08 .and. tau_x470ss > tau550) tau550 = tau_x470ss
          if (xday > 151.0.and.xday < 244.0) then  ! Jun, Jul, Aug
          tau550 = tau_x412ss
          if (Dstar1 > 1.09.and. wv <= 1.5 .and. tau_x470ss > tau550) tau550 = tau_x470ss
            if (Dstar1 > 1.07 .and. tau_x470ss/tau_x412ss > 3.5) tau550 = tau_x470ss
          endif
          ddx = 0.0
          if (xday > 243.0.and.xday < 274.0) ddx = 1.0  ! Sep
          if (xday > 120.0.and.xday < 152.0) ddx = 1.0  ! May
          if (xday > 151.0.and.xday < 244.0) ddx = 1.0  ! Jun, Jul, Aug
          if (ddx > 0.0) then 
            if (xday > 243.0.and.xday < 274.0) then
            if (Dstar1 > 1.01 .and. r412_tbl <12.0 .and. tau_x470ss > tau550) tau550 = tau_x470ss
            if (Dstar1 > 0.98 .and. wv > 1.7 .and. tau_x470ss > tau550) tau550 = tau_x470ss
            endif
            if (tau_x470ss > tau550) then
                if (xlong <= -5.0) tau550 = tau_x470ss
                if (xlong > -5.0 .and. xlong < 0.0) then
                dd = (5.+xlong) / 5.
                tau550 = (1.-dd)* tau_x470ss + dd * tau_x412ss
                endif
                if (xday > 120.0.and.xday < 152.0) then
                if (xlong <= -5.0) tau550 = tau_x470ss
                if (xlong > -5.0 .and. xlong < 10.0) then
                dd = (5.+xlong) / 15.
                tau550 = (1.-dd)* tau_x470ss + dd * tau_x412ss
                endif
                endif
            endif
          endif
            if (tau_x470ss > tau550) then
                if (xlong >= 35.0) tau550 = tau_x470ss
                if (xlong > 27.0 .and. xlong < 35.0) then
                dd = (xlong-27.0) / 8.
                tau550 = dd* tau_x470ss + (1.-dd) * tau_x412ss
                endif
          endif

          
          if (tau550 < 0.1) tau550 = tau550 + 0.05    ! check for geo zone

          if (Dstar1 > 1.1 .and. tau_x470ss > tau_x412ss) tau550 = tau_x470ss  !new

          if (xday > 151.0.and.xday < 258.0.and.
     1        Dstar1 > 1.2.and.tau_x650 < 0.0.and.tau550 < 0.8)
     1        tau550 = tau_x470ss * 2.
 
          if (Dstar1 > 1.2.and.tau_x412ss < 0.5 .and.tau550< 0.7)
     1        tau550 = Dstar1
 
          if (wv < 0.45.and.tau_x470ss >0.9.and.tau_x412ss < 0.7) tau550 = tau_x412ss  !new
          if (Dstar1 > 1.08 .and. tau_x470ss > tau550) tau550 = tau_x470ss

          end if
 
          if (gzflg >= 6 .and. gzflg <= 11 .AND. gzflg /= 10) then
               tau550 = tau_x470ss
               if (tau550 < 0.45) tau550 = (tau_x470ss+tau_x412ss2)/2.
          if (xday > 243.0 .and. xday < 274.0) then
               tau550 = tau_x470ss
               if (tau550 < 0.45.and. tau_x412ss2 > tau550) tau550 = (tau_x470ss+tau_x412ss2)/2.
          endif
               if (xday > 120.0 .and. xday < 151.0) then                   ! May
                tau550 = tau_x470ss
               if (tau550 < 0.45.and. tau_x412ss2 > tau550) tau550 = (tau_x470ss+tau_x412ss2)/2.
               end if
               if (xday > 90.0 .and. xday < 121.0) then                    ! Apr
                   tau550 = tau_x412ss2
                   if (tau550 < 0.5 .and. tau_x470ss > tau550) tau550 = tau_x470ss
               end if 
          if (xlat > 29.5 .and. Dstar1 < 0.97) then
             if (tau550 > 0.6) tau550 = tau_x412ss2_995
          end if 

          if (xday < 182.0) then
          if (xlat > 20.0.and. xlat < 28.0.and.xlong > 42.0.and. xlong < 50.0) then
             if (Dstar1 < 1.06.and.r412_tbl > 14.0.and.tau550 > 0.8) go to 10
          end if
          end if

          if (Dstar1 > 1.1 .and. tau_x470ss <0.45) tau550 = Dstar1

          if (xday > 243.0 .and. xday < 274.0) then    ! Sept
          if (Dstar1 < 0.98.and.r412_tbl > 12.0.and.tau550 > 0.8) go to 10
          endif
          if (xday > 273.0 .and. xday < 305.0) then    ! Oct
          if (Dstar1 < 1.04.and.r412_tbl > 12.0.and.tau550 > 0.75) go to 10
          endif

          if (xday > 150.0 .and. xday < 258.0) then    ! Jun, Jul, Aug
          if (xlat > 20.0.and. xlat < 28.0.and.xlong > 42.0.and. xlong < 50.0) then
             if (Dstar1 < 1.06.and.r412_tbl > 12.8.and.tau550 > 0.78) go to 10
          end if
           if (Dstar1 < 1.04.and.r412_tbl > 12.0.and.tau550 > 0.78.and. xlat>21.0.and.xlong<53.0) go to 10
           if (Dstar1 < 1.04.and.r412_tbl > 12.0.and.tau550 > 0.78.and. xlat>23.0.and.xlong>=53.0) go to 10
          end if

          end if 
           
        end if    ! end of barren
        if (gzflg < 6 .and. gzflg >0 .and. Dstar1 > 1.08 .and. tau_x470ss > tau550) tau550 = tau_x470ss
 
      end if  
      
!     -- for Arabian peninsula, just use AOT412 w/ w0=0.96. Override all of that above.
      if (gzflg >= 6 .and. gzflg <= 11 .AND. gzflg /= 10) then
        if ((r650_135 > 32.0) .AND. (r650_135/r412_135) > 3.7) then
           if (xlat > 20.5 .and. xlong < 50.0) go to 861
           if (xlat > 19.0.and. xlat < 20.5 .and. xlong > 43.0 .and. xlong < 48.0) go to 861
            if ((r412_135 < 7.25) .OR. r650_135 > 42.0 .OR. Dstar1 >1.03 .OR. use_alternate_brdf) then
            continue
           
          else if (xday >= 152.0 .AND. xday <= 243.0) then     ! summer
            if (Dstar1 < 0.99) tau550 = (tau_x412_new_94+tau_x412_new_96)/2.0    ! orig = 0.96

          else if (xday >= 244.0 .AND. xday < 335.0) then
            tau550 = (tau_x412_new_94+tau_x412_new_96)/2.0      ! orig = 0.995

          else
            tau550 = tau_x412_new_94
            if (xday >= 120.0 .AND. xday < 152.0.AND.tau550<0.3.AND.tau_x470ss>tau550) tau550 = tau_x470ss
          end if 
          
        end if
      end if
 861  continue
      
      if (lprint > 0) print *, 'after get_aot500: xlat, xlong, tau550: ', gzflg, xlat, xlong, tau550
!     Weakly absorbing dust
!      if (Dstar1 < 1.01.and.tau550 > 0.8.and.rat1 < 1.6.and.
!     1    gzflg <6 .and. gzflg >0)! .and. lc .eq. 6) 
!     1    tau550 = tau_x412ss_98
!      if (lprint > 0) print *, 'after absorbing dust check: tau550, Dstar1, rat1: ',  
!     1  tau550, Dstar1, rat1, rr470_mod, rr412_mod
      

!     Smoke pixels
      if (((gzflg <6 .and. gzflg >0) .OR. gzflg == 27) .and. xlat <20.0 .and. lc <5) then
        if (xday > 335.0 .or. xday < 32.0) then
          if (tau550 < 0.065 .or. tau_x470 <0.05) then
            dd412 = rr412_mod - r412_tbl
            dd470 = rr470_mod - r470_tbl
            dd650 = ref650    - r650_tbl
            if (debug) print *, 'smoke check, dd412, dd470, dd650: ', dd412, dd470, dd650
            !if (dd412 < 0.0 .AND. dd470 < 0.0 .AND. dd650 < 0.0) then
            if (dd412 < 1.0) then
              call smoke_mod(xthet,scat_ang,dd412,dd470,dd650,  
     1                       tau_x412ss,tau_x412ss_91,tau_smoke)
              if (debug) print *, 'smoke_mod, dd412, dd470, dd650,t412ss, t412_91, tsmoke: ',  
     1                              dd412, dd470, dd650, tau_x412ss, tau_x412ss_91, tau_smoke
              tau550 = tau_smoke
              if (tau_smoke < 0.25) then
                tau550 = tau_smoke *2.0
                if (debug) print *, 'smoke detected, aot550 reset: ', tau550
              end if
            end if        
          end if
        end if
      end if
      if (lprint > 0) print *, 'after smoke check: tau550: ', tau550
      
      if (lprint > 0) print *,                                  
     1   'gzflg,lc, tau_x412ss,tau_x470ss,tau550 =',   
     1    gzflg,lc,tau_x412ss,tau_x470ss,tau550

864   if (lprint > 0) then
        print *, 'lat,lon,gzflg,r650,tau_x650,tau_x470,tau550: ', xlat, xlong, gzflg, r650, 
     1   tau_x650, tau_x470, tau550
        print *, 'abs_aero_flag: ', abs_aero_flag
      end if
      
c     -- Try to detect and remove smoke over clouds. For example, see RGB:
c     --  MYD021KM.A2007081.0640*.hdf
      if (debug) print *, 'smoke detection: ', xlat, xlong, r650, tau_x650, refl6, w0_x470
      if (r650 <8.0 .and. tau_x650 > 3.49 .and. refl6 > 0.1 .and.w0_x470 > 0.999) go to 10
!      if (xday > 0.0 .and. xday < 121.0) then
! 				if (gzflg==16 .and. (w0_x470 > 0.96 .and. alpha < 1.0) .and. Dstar1 <1.0 .and. tau550 >1.0) then
! 					go to 10
! 				end if
!      endif
      
!     -- set tau550 to tau_x470ss over high elevation to minimize the effects of 
!     -- Rayleigh scattering. 
!     -- Morocco region
      if (px_elev > 750.0 .AND. (xlat > 28.0 .AND. xlat < 37.0 .AND. xlong > -12.0 .AND. 
     1    xlong < 10.0)) then
        tau550 = tau_x470
      endif
      
!     -- West Sudan region
      if (px_elev > 750.0 .AND. (xlat > 10.5 .AND. xlat < 19.5 .AND. xlong > 20.5 .AND. 
     1    xlong < 29.0)) then
        tau550 = tau_x470
      endif
      
!     -- Somalia
      if (lc == 6 .AND. gzflg < 6 .AND. (xlat > -5.0 .AND. xlat < 18.0 .AND. xlong > 35.0 .AND. 
     1    xlong < 52.0)) then
        tau550 = tau_x470ss
      endif

!     -- filter out dark Iraq oil fire plumes     
!     1/31/2017
      if ((rr412_mod - r412_135) < -1.0.and.px_elev < 1000.0 .and.
     1     gzflg >= 6 .and. gzflg <= 11) go to 10

!     -- filter out noisy high AOD over US desert   
      if (gzflg == 13 .and. lc == 6 .and. tau550 > 1.2 .and. Dstar1 < 1.0 .and.
     1    realbuf(6)/realbuf(22) > 1.65) go to 10

      if (gzflg_sav == 31 .and. tau550 > 0.6 .and. Dstar1 < 1.1 .and.
     1    realbuf(22)/realbuf(6) < 0.5) go to 10

c     -- check D* for dust plumes globally except for N.Africa (1-5, 26, 27), 
c     -- Arabian Peninsula (6-11), Taklimakan Desert (24), and Beijing (East China, 16) 
c     -- zones. Limit to locations at less than 4500 m elevation. Otherwise,
c			-- abnormally high AOT's are seen over Tibetan Plateau.
c     ----------------------------------------------------------------------------
c     ----------------------------------------------------------------------------
c     -- Increased threshold to 1.15 from 1.05 over globe. Created exception for 
c     -- ConUS (1.05). Needs to be tested and results analysed after next SIPS 
c     -- large-scale test! 
c     ----------------------------------------------------------------------------
c     ----------------------------------------------------------------------------
      if ((gzflg < 0 .OR. gzflg > 11) .AND. (gzflg /= 24 .AND. gzflg /= 26 .AND. 
     & gzflg /= 27 .AND. gzflg /= 16)) then
                if (px_elev < 4750.0) then
          dda1 = 1.15
          if (gzflg_sav == 13 .and. realbuf(22)/realbuf(6) > 0.3) dda1 = 1.05
          if (gzflg_sav == 31 .and. realbuf(22)/realbuf(6) > 0.5) dda1 = 1.05
                if (Dstar1 > dda1) tau550 = tau_x650
                if (Dstar1 > dda1 .and. tau_x470 > tau550) tau550 = tau_x470
          if (gzflg_sav == 31 .and. gzflg == 13 .and.  realbuf(22)/realbuf(6) > 0.3
     1        .and. Dstar1 > 1.05 .and. tau_x470 > tau550) tau550 = tau_x470
        end if
      end if
      
      xtau(1) = tau_x412
      xtau(2) = tau_x470
      xtau(3) = tau_x650
      
      if (w0_x412 > 0.0) w0_x = w0_x412
      if (w0_x470 > 0.0) w0_int_470 = w0_x470
      
      if (alpha < 1.0)  ssa(1) = w0_x
      if (alpha < 1.0)  ssa(2) = w0_int_470
      if (ssa(1) > 0.0) ssa(3) = 0.976
c-- Set Surface Type

      sfc_typ = -999.

      if (terrain_flag_new5.gt.0.0) sfc_typ = 7.
      if (terrain_flag_new5.lt.5.5.and.terrain_flag_new5.gt.0.0)
     1    sfc_typ = 10.

C      dd = xsfc650(ilon,ilat) - xsfc650_bk2(ilon,ilat)
C      dd1 = xsfc470b_bk(ilon,ilat) - xsfc470_bk(ilon,ilat)
 
      if (terrain_flag_new5.lt.5.5.and.terrain_flag_new5.gt.0.0.and.
     1    xlat.le.15.0.and.xlong.lt.26.0) then
          sfc_typ = 1.
C          if (xlong.gt.10.0.and.abs(dd1).ge.0.2) then
C          if (xsfc470_bk(ilon,ilat).gt.10.0) sfc_typ = 2.
C          if (xsfc470_bk(ilon,ilat).le.10.0.and.
C     1        xsfc470_bk(ilon,ilat).gt.9.0.and.
C     1        xsfc650_bk(ilon,ilat).lt.22.0) sfc_typ = 2.
C          endif
C          if (xlong.gt.10.0) then
C          if (xsfc470_bk(ilon,ilat).gt.10.0.and.
C     1        xsfc650_bk(ilon,ilat).lt.23.0) sfc_typ = 2.
C          endif
      endif
C      if (terrain_flag_new5.lt.5.5.and.terrain_flag_new5.gt.0.0.and.
C     1    xlat.gt.15.0.and.xlat.le.20.0.and.dd.ge.1.6) then
C          sfc_typ = 1.
C          if (xlong.gt.10.0.and.abs(dd1).ge.0.2) then
C          if (xsfc470_bk(ilon,ilat).gt.10.0) sfc_typ = 2.
C          if (xsfc470_bk(ilon,ilat).gt.10.0.and.
C     1        xsfc650_bk(ilon,ilat).lt.23.0) sfc_typ = 2.
C          if (xsfc470_bk(ilon,ilat).le.10.0.and.
C     1        xsfc470_bk(ilon,ilat).gt.9.0.and.
C     1        xsfc650_bk(ilon,ilat).lt.22.0) sfc_typ = 2.
C          endif
C      endif
C      if (terrain_flag_new5.lt.5.5.and.terrain_flag_new5.gt.0.0.and.
C     1    xlong.ge.26.0.and.xlat.le.20.0.and.dd.ge.1.6) then
C          sfc_typ = 1.
C          if (xlong.gt.10.0.and.abs(dd1).ge.0.2) then
C          if (xsfc470_bk(ilon,ilat).gt.10.0) sfc_typ = 2.
C          if (xsfc470_bk(ilon,ilat).gt.10.0.and.
C     1        xsfc650_bk(ilon,ilat).lt.23.0) sfc_typ = 2.
C          if (xsfc470_bk(ilon,ilat).le.10.0.and.
C     1        xsfc470_bk(ilon,ilat).gt.9.0.and.
C     1        xsfc650_bk(ilon,ilat).lt.22.0) sfc_typ = 2.
C          endif
C      endif

      if (terrain_flag_new5.gt.2.5.and.terrain_flag_new5.lt.3.5)
     1    sfc_typ = 3.

      if (terrain_flag_new5.gt.1.5.and.terrain_flag_new5.lt.2.5) then
          sfc_typ = 5.
C          if (xsfc470_bk(ilon,ilat).gt.6.7.and.xlong.lt.-7.9)
C     1        sfc_typ = 4.
      endif

C      if (terrain_flag_new5.gt.4.5.and.terrain_flag_new5.lt.5.5.and.
C     1    xlat.gt.31.3.and.xlong.lt.-7.9.and.xsfc470_bk(ilon,ilat)
C     1    .gt.6.7.and.xsfc470_bk(ilon,ilat).lt.9.5)
C     1    sfc_typ = 4.

      if (terrain_flag_new5.gt.3.5.and.terrain_flag_new5.lt.4.5)
     1    sfc_typ = 6.

      if (terrain_flag_new5.gt.8.5.and.terrain_flag_new5.lt.9.5)
     1    sfc_typ = 8.

      if (terrain_flag_new5.gt.9.5.and.terrain_flag_new5.lt.10.5)
     1    sfc_typ = 9.

      if (terrain_flag_new5.gt.5.5.and.terrain_flag_new5.lt.6.5.and.
     1    xlat.gt.24.6.and.xlat.lt.25.2.and.xlong.gt.46.1.and.
     1    xlong.lt.46.8.and.r412_tbl.ge.7.0.and.
     1    r412_tbl.le.9.9)
     1    sfc_typ = 9.

c     -- Set Flags
      if (tau_x650_flag.gt.0.or.tau_x470_flag.gt.0) qa_flag(4)= 3
      if (tau_x412.lt.0.0.and.tau_x470.lt.0.0.and.
     1    tau_x650.lt.0.0) then
      qa_flag(1)= 0
      qa_flag(2)= 0
      qa_flag(3)= 0
      qa_flag(4)= 2
      endif
      if (tau_x412.gt.0.0.or.tau_x470.gt.0.0.or.
     1    tau_x650.gt.0.0) qa_flag(1)= 1
      if (tau_x412.lt.0.05.or.tau_x470.gt.0.05)
     1    qa_flag(2)= 1
      if (tau_x412.ge.0.05.and.tau_x412.lt.0.3)
     1    qa_flag(2)= 2
      if (tau_x412.ge.0.3)
     1    qa_flag(2)= 3
      if (alpha.lt.0.5.and.alpha.gt.0.0)
     1    qa_flag(3)= 1
      if (alpha.ge.0.5.and.alpha.lt.1.4)
     1    qa_flag(3)= 2
    
 865  continue

!     -- in zone 1, in summer, sometimes the BRDF for low NDVI pixels is too high 
!     -- resulting in an abnormally low AOT.  In these cases, perform a second retrieval
!     -- using a constant fit (an alternate fit) for 470nm that will bring up 
!     -- those low areas.
      if (tau550 < 0.1 .AND. (gdt1%month >= 6 .AND. gdt1%month <= 8) .AND. gzflg == 1) then
        if (.NOT. use_alternate_brdf) then
          use_alternate_brdf = .true.
          goto 5
        end if
      end if
      
!     -- over Arabian Peninsula, trouble with some very low AOT's in certain areas. Redo 
!     -- the retrieval if in summer.
      if (.NOT. use_alternate_brdf .AND. (tau550 < 0.25 .AND.  
     &    (gdt1%month >= 6 .AND. gdt1%month< = 8) .AND. gzflg >= 6 .AND. gzflg <= 11 .AND. gzflg /= 10)) then
        if (r650_135 > 32.0 .AND. (r650_135/r412_135) > 3.7) then
          use_alternate_brdf = .true.
          goto 5
        end if
      end if

! from MODIS code for swir-vis surface database
c     -- use tau_x470sv96_dust over high elevation
c
      if (gzflg /= 15 .and. gzflg /= 19 .and. xlong> -20.0 .and. xlong< 70.0 .and.
     &    xlat> 5.0 .and. xlat< 45.0) then

         if (gzflg < 6 .and. gzflg > 0) then
         if (Dstar1 < 1.1) then
          if (r412_tbl > 12.0 .and. Dstar1 < 0.95 .and. tau_x470sv96_dust > 0.6) tau_x470sv96_dust = -999.
          dd = 999.
          if (xday >= 91.0 .and. xday < 244.0) dd = 10.0
          if (xlat >= 20.0 .and. xlong <= dd .and. tau_x470sv96_dust > 0.0) then 
          	tau550 = tau_x470sv96_dust
          	alg_flag = 12
          	old_alg_flag = 0
       	  endif

          if (xday < 91.0 .or. xday >= 244.0) go to 871
          if (xlong >= 10.0 .and. xlong <15.0 .and. xlat > 15.0) then
          dd = (xlong - 10.0) / 5.
            if (tau_x470sv96_dust > 0.0 .and. tau550 > 0.0) then
             tau550 = tau_x470sv96_dust + (tau550 - tau_x470sv96_dust) * dd
             alg_flag = 12
             old_alg_flag = 0
    	    endif
         endif
871      continue

          dd = 999.
          if (xday >= 91.0 .and. xday < 244.0) dd = 10.0
          if (xlat >= 15.0 .and. xlat <20.0 .and. xlong < dd) then
          dd = (xlat - 15.0) / 5.
            if (tau_x470sv96_dust > 0.0 .and. tau550 > 0.0) then
             tau550 = tau550 + (tau_x470sv96_dust - tau550) * dd
             alg_flag = 12
             old_alg_flag = 0
    	    endif
          endif

          tau550_sav = tau550

          if (tau_x470sv96_dust > tau550) then
          	tau550 = tau_x470sv96_dust
          	alg_flag = 12
          	old_alg_flag = 0
          endif            

          if (xlat >= 25.0 .and. xlong >= 20) then
                tau550 = tau550_sav
          endif

          if (xlong >= 10.0 .and. xlong <20.0 .and. xlat > 20.0) then
          dd = (xlong - 10.0) / 10.
            if (tau550_sav > 0.0 .and. tau550 > 0.0) then
             tau550 = tau550 + (tau550_sav - tau550) * dd
            endif
         endif

          if (xlat >= 20.0 .and. xlat <25.0 .and. xlong > 10) then
          dd = (xlat - 20.0) / 5.
            if (tau550_sav > 0.0 .and. tau550 > 0.0) then
             tau550 = tau550 + (tau550_sav - tau550) * dd
            endif
          endif


          if (r412_tbl > 12.0 .and. tau_x470sv96_dust < 0.0) then
             if (realbuf(23) > 25.0 .and. realbuf(23)/realbuf(6) > 0.7
     &           .and. realbuf(23)/realbuf(6) < 0.85) go to 867
             tau550 = -999.0
             alg_flag = -999.0
             old_alg_flag = -999.0             
867          continue
          endif
         if (xday >= 152.0 .and. xday < 244.0 .and. Dstar1 >= 1.04 .and. tau_x470sv96_dust > tau550) then
           tau550 = tau_x470sv96_dust
           alg_flag = 12
           old_alg_flag = 0
         endif
         endif
         if (Dstar1 >= 1.1 .and. tau_x470sv96_dust > tau550) then
         	tau550 = tau_x470sv96_dust
         	alg_flag = 12
         	old_alg_flag = 0
         endif
         if (Dstar1 >= 1.1 .and. tau_x470sv96_dust <0.0 .and. tau550 < 0.3) then
         	tau550 = tau_x470sv96_dust
         	alg_flag = 12
         	old_alg_flag = 0
         endif
         endif

        if (gzflg == 2 .or. (gzflg > 5 .and. gzflg < 12)) then
          tau550 = tau_x470sv96_dust
          alg_flag = 12
          old_alg_flag = 0
				endif
				
        if (gzflg == 30) then
          if (xlong > 65.0 .and. xlong< 70.0) then
            dd = (xlong - 65.0) / 5.
            if (tau_x470sv96_dust > 0.0 .and. tau550 > 0.0) then
              tau550 = tau550 + (tau_x470sv96_dust - tau550) * (1.-dd)
              alg_flag = 12
              old_alg_flag = 0
						endif
            if (tau_x470sv96_dust < 0.0) then
            	tau550 = tau_x470sv96_dust
            	alg_flag = 12
            	old_alg_flag = 0
						endif
           else
           tau550 = tau_x470sv96_dust
           alg_flag = 12
           old_alg_flag = 0
          endif
        endif

        if (gzflg /= 30 .and. tau_x470sv96_dust > tau550 .and.
     &     (px_elev > 750.0 .and. (xlat> 14.0. or. xlong< 5.5 .or. xlong> 19.0)) .and.
     &     (px_elev > 750.0 .and. (xlat> 14.0 .or. xlong > -6.37))) then
          tau550 = tau_x470sv96_dust
          alg_flag = 12
          old_alg_flag = 0
				endif

        if (gzflg == 32) then! Portugal and Spain 
          if (season == 1) then
            if (tau_x470sv96 < 0.3) then
              tau550 = tau_x470sv96
              alg_flag = 12
              old_alg_flag = 0
            elseif (tau_x470sv96 < 0.6) then
              aod_frac = (tau_x470sv96-0.3)/(0.6-0.3)
              tau550 = (1.0-aod_frac)*tau_x470sv96+aod_frac*tau_x470sv94
              alg_flag = 12
              old_alg_flag = 0
            else
              tau550 = tau_x470sv94
              alg_flag = 12
              old_alg_flag = 0
            endif
          else
            if (tau_x470sv995 < 0.3) then 
              tau550 = tau_x470sv995
              alg_flag = 12
              old_alg_flag = 0
            elseif (tau_x470sv995 < 0.6) then
              aod_frac = (tau_x470sv995-0.3)/(0.6-0.3)
              tau550 = (1.0-aod_frac)*tau_x470sv995+aod_frac*tau_x470sv96
              alg_flag = 12
              old_alg_flag = 0
            else
              tau550 = tau_x470sv96
              alg_flag = 12
              old_alg_flag = 0
            endif
          endif
        endif

        if (gzflg == 33) then! Palencia
          if (season == 1) then 
            if (tau_x412sv96 < 0.3) then
              tau550 = tau_x412sv96
              alg_flag = 12
              old_alg_flag = 0
            elseif (tau_x412sv96 < 0.6) then
              aod_frac = (tau_x412sv96-0.3)/(0.6-0.3)
              tau550 = (1.0-aod_frac)*tau_x412sv96+aod_frac*tau_x412sv94
              alg_flag = 12
              old_alg_flag = 0
            else
              tau550 = tau_x412sv94
              alg_flag = 12
              old_alg_flag = 0
            endif
          else
            if (tau_x412sv98 < 0.3) then
              tau550 = tau_x412sv98
              alg_flag = 12
              old_alg_flag = 0
            elseif (tau_x412sv98 < 0.6) then
              aod_frac = (tau_x412sv98-0.3)/(0.6-0.3)
              tau550 = (1.0-aod_frac)*tau_x412sv98+aod_frac*tau_x412sv96
              alg_flag = 12
              old_alg_flag = 0
            else
              tau550 = tau_x412sv96
              alg_flag = 12
              old_alg_flag = 0
            endif
          endif
        endif
      endif

!       24 January 2018 JLee TEST
!      print *, gzflg, px_elev
!      if (gzflg == 15 .and. px_elev > 600.0) then! India
!        tau550 = tau_x412sv94
!      endif

!end from MODIS code for swir-vis surface database

! jlee test 2.2 um AOD everywhere (global) 20170725
!      tau550 = tau_x470sv
! end jlee test

!     -- We only really want to extrapolate aot550 past 3.5. All other AOT values
!     -- should be truncated at 3.5 to maintain current AE and SSA values.
      if (xtau(1) > 3.5) xtau(1) = 3.5
      if (xtau(2) > 3.5) xtau(2) = 3.5
      if (xtau(3) > 3.5) xtau(3) = 3.5

!     -- To detect the smoke plumes for aerosol typing
      dda1 = 999.0
      dda2 = 0.0
      if (realbuf(22) > 0.0 .and. realbuf(22) < 50.0) then
          dda1 = (realbuf(22)*realbuf(6)) / (realbuf(23)*realbuf(23))
      endif

      ddx = 0.90
      if (xlong < -20.0 .or. (xlat <30.0 .and. xlat >-40.0 .and.
     1    xlong >-30.0 .and. xlong <60.0 .and. gzflg < 1) .or.
     1    gzflg == 26 .or. gzflg == 27 .or. (gzflg > 0 .and. gzflg < 6)) go to 868
      if (xphi >= 90.0 .and. xthet > 50.0 .and. realbuf(23) < 33.0) then
        ddx = 0.90 - (0.90-0.83)*(xthet-50.0)/15.0
        if (xthet > 65.0) ddx = 0.83
      endif
868   continue

      if (gzflg == 15) then    ! India
       ddx33 = dda1
       if (xphi < 90.0 .and. xthet > 35.0) then
          if (xthet > 35.0 .and. xthet < 60.0) then
          ddx33 = dda1 - (xthet - 35.0) * (1.00-0.90)/ 25.0
          endif
          if (xthet >= 60.0) ddx33 = dda1 - (1.00-0.90)
       endif
      dda1 = ddx33
      ddx  = 0.90
      endif

      if (tau550 > 0.4 .and. dda1 < ddx .and. bt11 > 286.0) dda2 = 1.0

      dda3 = 999.0
      if (realbuf(22) > 0.0 .and. realbuf(22) < 50.0) then
      dda3 = realbuf(23)-realbuf(22) - (realbuf(6)-realbuf(23))*(488.-412.)/(670.-488.)
      endif
      if (dda1 < 0.90 .and. dda1 > 0.86 .and. dda3 < 1.4) dda2 = 0.0

      !if (gzflg == 15 .and. realbuf(22) > 12.0 .and. dda3 > 1.2 .and. r650 < 11.0) dda2 = 1.0

!     -- check consistency when tau550 is undefined
      if (tau550 < -900.0) then
        xtau(1)       = -999.0
        xtau(2)       = -999.0
        xtau(3)       = -999.0
        ssa(1)        = -999.0
        ssa(2)        = -999.0
        ssa(3)        = -999.0
        tau550        = -999.0
        alpha         = -999.0
        r412          = -999.0
        r470          = -999.0
        r650          = -999.0
        tau_x470_flag = -999
        tau_x650_flag = -999
        alg_flag      = -999
        old_alg_flag  = -999
      end if
c-------------------------------------------------------------
c Set output buffer
c-------------------------------------------------------------
      if (platform .eq. 'AHI') then 
        !-- apply scattering angle correction to India regions (15,19,20,30) 
        !   and Australia (12) and New Zealand (lat/lon limits)      
      if (gzflg == 19 .OR. gzflg == 15 .OR. gzflg == 20 .OR. gzflg == 30 .OR. gzflg == 12 .OR. 
     &  (xlat > -50.0 .AND. xlat < -30.0 .AND. xlong > 160.0 .AND. xlong < 180)) then
        tau550 = tau550 * (2.2 / (1.765 + 0.018*scat_ang))
      end if
      end if 

      if (debug) print *, 'final spectra aot: ', xtau(1), xtau(2), xtau(3)
      if (debug) print *, 'final tau550, ae: ', tau550, alpha
      if (debug) print *, 'final SSA: ', ssa(1), ssa(2), ssa(3)
      do i=1,3
         outbuf(i)    = xtau(i)
         outbuf(i+3)  = ssa(i)
      enddo

      outbuf(7)   = tau550
      outbuf(8)   = alpha
      outbuf(9)   = r412
      outbuf(10)  = -999.0
      if (tau_x470_flag_veg == 1 .or. tau_x470_flag > 0) outbuf(10) = 1.0  ! out-of-bound flag
      !tau_x470_flag_veg : veg code outbound flag, tau_x470_flag : DB code outbound flag
!      outbuf(10)  = 1.0*tau_x650_flag
      outbuf(11)  = r470
      outbuf(12)  = r650
      outbuf(13)  = xthet
      outbuf(14)  = dda2
      outbuf(15)  = sfc_typ
      outbuf(18)  = float(alg_flag)
      outbuf(19)  = outbufvg(21) ! veg bright surface flag
      outbuf(20)  = outbufvg(7)  ! veg AOD to be used for rcndvi_avg
      outbuf(22)  = float(old_alg_flag)
      if (realbuf(22) > 0.0 .and. realbuf(22) < 50.0)
     1    outbuf(21)  = realbuf(22)/realbuf(23)        ! LER412/488nm ratio
      
10    continue
      return
      end

c     --------------
c
      subroutine search(dflag,xbar,x,n,i)
c
c     purpose
c       locate position in table of point at which interpolation is
c       required
c
c     usage
c       call  search (xbar, x, n, i)
c
c     description of parameters
c       xbar   - point at which interpolation is required
c       x      - array containing independent variable
c       n      - number of points in x array
c       i      - index specifying segment containing xbar
c
      logical dflag
      dimension x(n)
      data b/.69314718/
      icnt = 0
      if (n.lt.2) go to 15
      if(x(1).gt.x(2)) go to 17
      m = int((log(float(n)))/b)
      i=2**m
      if (i .ge. n) i = n-1
      k=i
   10 k=k/2
      if (k .eq. 0) icnt = icnt + 1
      if (icnt .ge. 2) goto 27
      if (xbar.ge.x(i).and.xbar.lt.x(i+1))return
      if (xbar.gt.x(i)) go to 12
      i = i-k
      go to 10
   12 i = i+k
      if (i.lt.n) go to 10
      i=n-1
      go to 10
   15 print *, "Search n is less than 2."
      return
   17 print *, "Search table is not in increasing order."
      return

   27 continue
      do 22 i=1,n-1
         if (xbar.ge.x(i).and.xbar.le.x(i+1)) return
   22 continue           
      write(6,*) "setting dflag = true"
      dflag = .true.
      return

      end

      subroutine search2(dflag,xbar,x,n,i,fx)
c
c       call  search (xbar, x, n, i)
c
c     description of parameters
c       xbar   - point at which interpolation is required
c       x      - array containing independent variable
c       n      - number of points in x array
c       i      - index specifying segment containing xbar
c       fx     - fraction from i (between i and i+1)
c
      dimension x(n)
      logical dflag

      call search(dflag,xbar,x,n,i)

      fx   = (xbar-x(i))/ (x(i+1)-x(i))

      end



c------------------------------------------------------------
      subroutine aero_mod (tau_x412,tau_x470,tau_x650,
     1           tau_x412_91,aot_mod)

      real aot_mod(6)

      aot_mod(1) = tau_x650
      aot_mod(2) = tau_x470
      aot_mod(3) = tau_x412_91
      aot_mod(4) = tau_x412*1.9
      aot_mod(5) = tau_x470*2.
      aot_mod(6) = tau_x650*2.8
      return
      end
      
			integer function handle_lut_out_of_bounds(geo_zone, retr_flag, aot) result(status)
				implicit none
				
				integer, intent(in)	:: geo_zone
				integer, intent(in)	:: retr_flag
				real, intent(inout)	:: aot
				
				status = 0
												
				if (retr_flag == -10) then      ! -10 = return value from aero_412, etc indicating
					select case (geo_zone)        ! refl < yy(1) scenario.
						case (5:11, 15, 16, 19, 20, 21, 23, 24, 26, 27)
							aot = 0.06
						case (1:4, 12:14, 17, 18, 22, 25, 28, 29, 30, 31, 32, 33, 34)
							aot = 0.02
						case (-999)
							aot = 0.02
						case default
              aot = 0.02
							print *, "ERROR: Invalid geographical zone in handle_lut_out_of_bounds: ", geo_zone
              print *,  "Use default aot of 0.02 rather than return without output"
							!status = -1
							!return
					end select
				end if
				
				return
				
			end function handle_lut_out_of_bounds
			
			
			subroutine smoke_mod(view,scat_ang,dd412,dd470,dd650,  
     1                     tau_x412ss,tau_x412ss_91,tau_smoke)
        implicit none
        
        real, intent(in)        ::  view
        real, intent(in)        ::  scat_ang
        real, intent(in)        ::  dd412
        real, intent(in)        ::  dd470
        real, intent(in)        ::  dd650
        real, intent(in)        ::  tau_x412ss, tau_x412ss_91
        real                    ::  tau_smoke
     
        if (dd412 <1.0 .and. dd470 <1.0) then
            tau_smoke = tau_x412ss * 2.85  
          if (dd470 > dd412 .and. dd650 > dd412)   
     1      tau_smoke = tau_x412ss_91 * 1.77
          if (abs(dd412-dd470) < 0.25) then
            if (dd412 <0.35) tau_smoke = tau_x412ss * 4.3
            if (scat_ang <158.0 .and. dd412 <0.8) then
            tau_smoke = tau_x412ss * 3.5
            if (view >45.0 .and. dd470 <0.2 .and. abs(dd412-dd470) < 0.24)  
     1        tau_smoke = tau_x412ss * 6.0
            endif
          endif
          if (abs(dd412-dd650) < 0.1 .and. dd412 >0.7)    
     1        tau_smoke = tau_x412ss_91
          if (dd650 > 1.0 .and. abs(dd412-dd470) < 0.25) then
              tau_smoke = tau_x412ss * 4.5
          if (scat_ang <158.0 .and. dd412 >0.8)    
     1        tau_smoke = tau_x412ss * 3.5
          endif
        endif
     
        if (dd412 >1.0 .and. dd470 <1.0) then
             tau_smoke = tau_x412ss * 2.7
          if (abs(dd412-dd470) < 0.25 .and. dd650 <0.)   
     1      tau_smoke = tau_x412ss * 2.4
        endif
     
        if (dd412 <1.0 .and. dd470 >1.0) then
             tau_smoke = tau_x412ss * 2.7
          if (abs(dd412-dd470) < 0.25 .and. dd650 <0.)    
     1        tau_smoke = tau_x412ss * 2.4
          if (dd650 - dd470 > -0.2) then
             tau_smoke = tau_x412ss * 2.1
          if (view < 45.0) tau_smoke = tau_x412ss_91
          endif
        endif
     
        if (dd412 >1.0 .and. dd470 >1.0) then
          if (dd412 > dd470) tau_smoke = tau_x412ss_91 * 1.15
          if (dd470 > dd650 .and. dd650 > dd412 .and. dd650 >1.0) then
            if (view >= 45.0) tau_smoke = tau_x412ss * 2.1
            if (view <  45.0) tau_smoke = tau_x412ss_91
          endif
          if (dd412 <2.0 .and. dd470 >2.0)   
     1      tau_smoke = tau_x412ss * 1.4
          if (dd650 > dd470 .and. dd470 > dd412) then
            if (view <= 25.0) tau_smoke = (tau_x412ss+tau_x412ss_91)/2.
            if (view > 25.0 .and. abs(dd412-dd470) < 0.2)  
     1        tau_smoke = (tau_x412ss+tau_x412ss_91)/2.
          endif
          if (dd650 > 3.5) tau_smoke = tau_x412ss
        endif
     
        if (dd650 <-1.0 .and. abs(dd412-dd650)>1.9 .and. view >45.) then
          if (dd412 >= 0.88) tau_smoke = tau_x412ss * 4.0
          if (dd412 <0.88 .and. dd412 >= 0.6) tau_smoke = tau_x412ss *6.0
          if (dd412 <0.6 .and. dd412 >= 0.35) tau_smoke = tau_x412ss *8.0
          if (dd412 <0.35) tau_smoke = tau_x412ss_91 *8.0
          if (abs(dd412-dd470) < 0.14 .and. dd650 > -1.7 .and. dd412 <0.4)  
     1        tau_smoke = tau_x412ss * 3.5
          if (dd412 <0.) tau_smoke = tau_x412ss_91 *9.0
        endif
     
        return
 
      end subroutine smoke_mod
			
