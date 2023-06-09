      subroutine find_v_veg(month,season,realbuf,tmpvg,
     1           r412sv,r470sv,gzflg,outbufvg,tau_x470_flag,platform) !modified by CH 11/18/2016
c.    subroutine myd_aero_veg(imod,season,itmp,jtmp,realbuf,tmpvg,
c.   1           outbufvg)

c...  Main subroutines for aerosol retrieval over vegetated surfaces.
c...     Written by Myeong-Jae Jeong (MJ)
c...     Last modified Aug 08, 2011
c...
c...  Inputs:
c...     season:   1:DJF; 2:MAM; 3:JJA; 4:SON
c...     imod:     aerosol model index for 470nm --> no longer an input
c...     nvalx470: TOA refl LUT for 470nm
c...     nvalx650: TOA refl LUT for 650nm
c...     xlcvr_2:  MODIS land cover (IGBP)
c...     regid_2:  Region Index (to choose a set of aerosol models
c...               depending on regions and land cover)
c...     realbuf:  see below for definitions
c...     tmpvg:
c...
c...  Output:
c...     outbufvg: see the end of this subroutine for definitions
c...

      use viirs_aerosol_luts, only: aero_470, aero_650

c      include 'aottbl.inc'
      include 'newaottbl.inc'

      parameter(nx2=3600, ny2=1800)  ! landcover data dimension
      character(len=*)    ::   platform
      logical dflag, debug, do_sv, do_nir

      integer bflag ! jlee added
      integer month, season, gzflg

      real realbuf(26), tmpvg(7), outbufvg(21), xnvalm6(6)
      real r412sv, r470sv

c      real nvalx470(10,46,30,10,4,24),nvalx650(10,46,30,10,24)
      real nval(10,46,30), yy(10), yyw(8) !tau(10), 
c.    real xnvalm6(6), realbuf(13), outbuf(20)
c      real*4 xlcvr_2(nx2,ny2)
      real    tau_x470_1, tau_x470_2, tau_x470_3
      real    tau_x470sv_1, tau_x470sv_2, tau_x470sv_3
      integer tau_x470_flag, tau_x650_flag, tau_ini_flag, tau_x470_flag_ini
      integer tau_x470_flag1, tau_x470_flag2, tau_x470_flag3
      integer tau_x470sv_flag1, tau_x470sv_flag2, tau_x470sv_flag3 
     
c      real theta0(10), theta(46), phi(30)
c      real sfc_ref412(20), sfc_ref470(24), sfc_ref650(24)
      real*4 r412db,r470db,r650db,cl_flag
      real xtau(3),ssa(3),qa_flag(4),aot_mod(6) !,w0_470(4)
      character*4 w0_name470(4)
      character*12 aer_tab(10)
      real*4 ctharr(3)     ! 08/05/2011
      integer imodarr(3)   ! 08/05/2011
      real modfrac(3)      ! 16 January 2018 JLee

c      common /angle_node/ theta0, theta, phi
c      common /sfcref_node/ sfc_ref412, sfc_ref470, sfc_ref650
      common /fname_node/ aer_tab, w0_name470
      data pi      /3.14159/
C      data theta0  /0.0,8.0,16.0,24.0,32.0,40.0,48.0,56.0,64.0,72.0/
C      data tau     /0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5/
C      data w0_470  /0.91, 0.94, 0.96, 0.99/
C      data sfc_ref412 /1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,
C     1                 11.,12.,13.,14.,15.,16.,17.,18.,19.,20./
C      data sfc_ref470 /1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,
C     1                 11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,
C     2                 21.,22.,23.,24./
C      data sfc_ref650 /1.,3.,5.,7.,9.,11.,13.,15.,17.,19.,21.,
C     1                 23.,25.,27.,29.,31.,33.,35.,37.,39.,41.,
C     2                 43.,45.,47./

c.    data aer_tab /'table_0.0aot', 'table_0.1aot', 'table_0.3aot',
c.   1              'table_0.5aot', 'table_1.0aot', 'table_1.5aot',
c.   2              'table_2.0aot', 'table_2.5aot', 'table_3.0aot',
c.   3              'table_3.5aot'/
      data w0_name470 /'0.91','0.94','0.96','0.99'/
     
      bflag = 0 ! jlee added 
      debug = .false.

c-------------------------------------------
c... define nodes for angles
c      do i = 1, 46
c         theta(i) = 2.*float(i-1)
c      enddo

c      do i = 1, 30
c         phi(i) = 6. + 6.*float(i-1)
c      enddo
      
      mm = 12     ! solar zenith
      nn = 46     ! satellite zenith
      ll = 30     ! rel. azimuth
      ma = 10     ! tau
      mw =  8     ! ssa
c-------------------------------------------

c-------------------------------------------
c... IGBP Landcover data
      xlatbeg=-90.0
      xlonbeg=-180.0
      xyintv2=0.10
c-------------------------------------------

c... ************************************************
c... ************************************************
c... come back later to set these threshold values, if necessary...(MJ) Jun 15, 2010
      xmnsfc1=1.0  ! allowed min. sfc refl. for B1
      xmxsfc1=47.0 ! allowed max. sfc refl. for B1
      xmnsfc3=1.0  ! allowed min. sfc refl. for B3
      xmxsfc3=24.0 ! allowed max. sfc refl. for B3
c... ************************************************
c... ************************************************

c.?   xlat5  = realbuf(1)
c.?   xlong5 = realbuf(2)
      xlat   = realbuf(1)
      xlong  = realbuf(2)
      sza    = realbuf(3)
      xthet  = realbuf(4)
      xphi   = realbuf(5)
      xnvalm6(1)=tmpvg(1)     ! (B1)
      xnvalm6(2)=tmpvg(2)     ! (B3)
      xnvalm6(3)=tmpvg(3)     ! (B8)
      xnvalm6(4)=tmpvg(4)     ! 2.1um (B7)
      xnvalm6(5)=tmpvg(5)     ! 850nm (B2)
      xnvalm6(6)=tmpvg(6)     ! 1.2um (B5)
c      band26    =realbuf(10)  ! 1.38um(B26)
c      band02    =xnvalm6(5)
c      band05    =xnvalm6(6)
c      band07    =tmpvg(7)     ! jlee commented out
c      band08    =xnvalm6(3)   ! (B8)
c      band03    =xnvalm6(2)   ! (B3)
c      band01    =xnvalm6(1)   ! (B1)
c      cl_flag   =realbuf(13)
      r412db    =realbuf(24)*100.
      r470db    =realbuf(25)*100.
      r650db    =realbuf(26)*100.

c      print *, itmp, jtmp, realbuf, lsf   ! test

c... ============================================================
c     -- sun glint mask
c
      cc     = 3.14159/180.
      psi    = acos(cos(sza*cc)*cos(xthet*cc) +
     1         sin(sza*cc)*sin(xthet*cc)*cos(xphi*cc))
      glint_ang = psi/cc

c      if (abs(psi/cc).lt.35.0) go to 10

c     -- scattering angle (scat_ang)
c
      psi    = acos(cos(sza*cc)*cos(xthet*cc) -
     1         sin(sza*cc)*sin(xthet*cc)*cos(xphi*cc))
      scat_ang = 180. - psi/cc

c      if (scat_ang .gt. 158.) go to 10
c      if (scat_ang .gt. 175.) go to 10

c... ######################################
      lcvr=1  ! initialization
      ioprg=5 ! ititialization
c...  get lcvr (IGBP landcover type; integer) data here
      idx=int((xlong-xlonbeg)/xyintv2) + 1
      idy=int((xlat-xlatbeg)/xyintv2) + 1
            
      if(idx.ge.1.and.idx.le.nx2.and.idy.ge.1.and.idy.le.ny2) then
         sfc_typ = xlcvr_2(idx,idy)
         xreg_id = regid_2(idx,idy)  ! 08/05/2011
      else
c        sfc_typ = 1.0*lcvr  ! @@@@@ just for test
         print *, 'lcvr data out of bound: idx,nx2,idy,ny2: ',idx,nx2,idy,ny2,xlat,xlong
c        noob=noob+1
         stop ! @@@@@
      endif
      lcvr=int(sfc_typ)
      ioprg=int(xreg_id)  ! 08/05/2011

c.    if(itmp.ge.100.and.itmp.lt.110.and.jtmp.ge.200.and.
c.   1   jtmp.lt.210) print *, sfc_typ, lcvr, xreg_id, ioprg
c... ######################################

      !xmu = cos(sza * 3.14159/180.)

c... refl unit conversion (pi*L/F/xmu --> L/F)
      !do i = 1, 6
      !   xnvalm6(i) = xnvalm6(i) * xmu/3.14159
      !enddo

c     if (xphi.gt.179.99) go to 10
c     if (xphi.gt.179.99) return   ! choose this or below (MJ)
      if (xphi.gt.179.99) xphi=179.99  ! choose this or above (MJ)
      if (xphi.lt.6.0) xphi = 6.

      x1 = sza
      x2 = xthet
      x3 = xphi

c... ============================================================

      refl1= xnvalm6(3)    ! 412 nm
      refl6= xnvalm6(1)    ! 650 nm
      refl3= xnvalm6(2)    ! 470 nm
      !refl21=band07        ! 2.1um ! jlee commneted out
      refl21 = xnvalm6(4)  !2.1um ! jlee added
      refl865 = xnvalm6(5) !865 nm !jlee added
      refln21 = xnvalm6(4)*3.14159/cos(sza*cc) ! 2.1 um, normalized ! jlee added
      refln865 = xnvalm6(5)*3.14159/cos(sza*cc)! 865 nm, normalized ! jlee added

      refn672_rc = -999.0
      refn865_rc = -999.0

      call calc_rc_ref672(sza,xthet,xphi,refl6,refn672_rc)
      call calc_rc_ref865(sza,xthet,xphi,refl865,refn865_rc)

      sirndvi=(xnvalm6(6)-xnvalm6(4))/(xnvalm6(6)+xnvalm6(4)) !jlee added, ndvi_swir
      rc_ndvi =(refn865_rc-refn672_rc)/(refn865_rc+refn672_rc) !jlee added, Rayleigh-corrected ndvi_vis
      toa_ndvi = (refl865-refl6)/(refl865+refl6)
      airmass = 1.0/cos(sza*cc)*1.0/cos(vza*cc) !jlee added
      !print *, toa_ndvi, rc_ndvi
c... ============================================================

!===================================================================
c---------------------------------------------------------------
c     Initialization
c---------------------------------------------------------------
c     --intermediate parameters
      w0_x       = -999.
      w0_int     = -999.
      tau_x412   = -999.
      tau_x412_91 = -999.
      tau_x470   = -999.
      tau_x650   = -999.
      xxrat      = -999.
      xxrat2     = -999.
      aot        = -999.
c     -- output parameters
      tau550     = -999.
      alpha      = -999.

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

      do i = 1,21
        outbufvg(i) = -999.
      enddo
c-----------------------------------------------------------------------
c...  estimate surface reflectance at 470nm & 650nm
c-----------------------------------------------------------------------
      xeAs_B1 = -999.0  ! initialization
      xeAs_B3 = -999.0
!==============================================================
!jlee added
      ! do_sv = false and do_nir = false: use swir method
      ! do_sv = false and do_nir = true: use nir method
      ! do_sv = true and do_nir = false: use 2.2 um database, alpha from swir method 
      ! do_sv = true and do_nir = true: use 2.2 um database, alpha from nir method
      do_sv = .false.
      do_nir = .false.
      if (gzflg.eq.15.or.gzflg.eq.19) then !swir method only over india
        do_sv = .false. 
        do_nir = .false.
      elseif (lcvr.eq.12.and.rc_ndvi.lt.0.4.and.(ioprg.eq.1.or.ioprg.eq.3)) then !use nir method over bright croplands in North and South America
        do_sv = .false. 
        do_nir = .true.
      else !default, swir method only
        do_sv = .false.
        do_nir = .false.
      endif

      if (do_nir) then !this condition is also used other parts of the code below 
        call get_sfcrfl_nir_vis(season,lcvr,refn865_rc,rc_ndvi,
     1                          xeAs_B1,xeAs_B3)
      else
        call get_sfcrfl_swir_vis(season,lcvr,refln21,rc_ndvi,
     1                           xeAs_B1,xeAs_B3)
      endif

      if (xeAs_B1.lt.xmnsfc1.and.xeAs_B1.gt.-900.0) xeAs_B1=xmnsfc1
      if (xeAs_B3.lt.xmnsfc3.and.xeAs_B3.gt.-900.0) xeAs_B3=xmnsfc3

      refl = refl6
      r650 = xeAs_B1
      dflag = .false.
      trflg = 1
      call aero_650(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1     r650,tau_ini,tau_ini_flag,0,0.0,0.0,0,trflg)!tau_x470_flag,tau_x412,tau_x470,tau_x412_flag_91

!      call aero_650veg(refl,x1,x2,x3,mm,nn,ll,ma,
!     1                 r650,tau_ini,tau_ini_flag,
!     2                 tau_x470_flag_ini,tau_x470_ini)
      
!      imod=imodarr(1)
!      refl = refl3
!      r470 = xeAs_B3 
!
!      call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,imod,r470,
!     1                 tau_ini,tau_ini_flag)

      if (tau_ini.gt.0.3) then 
        xeAs_B1 = -999.0
        xeAs_B3 = -999.0

        rc_ndvi_c = rc_ndvi+0.257*tau_ini*rc_ndvi
     1                     +0.172*(airmass-1.0)*tau_ini*rc_ndvi
     
        if (scat_ang.lt.120) then
          rc_ndvi_c = rc_ndvi+0.257*tau_ini*rc_ndvi
     1               +0.172*(airmass-1.0)*tau_ini*rc_ndvi
     2               +0.025*(30.0-(scat_ang-90.0))*tau_ini*rc_ndvi
        endif

        if (do_nir) then
          call get_sfcrfl_nir_vis(season,lcvr,refn865_rc,rc_ndvi_c,
     1                            xeAs_B1,xeAs_B3)

!! swir sfc for high aod over cropland
!          refl = refl6
!          r650 = xeAs_B1
!          call aero_650veg(refl,x1,x2,x3,mm,nn,ll,ma,
!     1                     r650,tau_ini,tau_ini_flag,
!     2                     tau_x470_flag_ini,tau_x470_ini)
!
!          if (tau_ini.gt.0.4) then
!            call get_sfcrfl_swir_vis(season,lcvr,refln21,rc_ndvi_c,
!     1                               xeAs_B1,xeAs_B3)
!          endif
        else
          call get_sfcrfl_swir_vis(season,lcvr,refln21,rc_ndvi_c,
     1                             xeAs_B1,xeAs_B3)
        endif

        if (xeAs_B1.lt.xmnsfc1.and.xeAs_B1.gt.-900.0) xeAs_B1=xmnsfc1
        if (xeAs_B3.lt.xmnsfc3.and.xeAs_B3.gt.-900.0) xeAs_B3=xmnsfc3
      endif! if (tau_ini.gt.0.3) then

      if (xeAs_B3.gt.24.0) return
      
!end jlee added
!===================================================================

c     sfc_typ = -999.

c      print *, 'xeAs_B1,xmnsfc1,xeAs_B3,xmnsfc3',xeAs_B1,xmnsfc1,xeAs_B3,xmnsfc3

c      if(xeAs_B1.lt.xmnsfc1.or.xeAs_B1.gt.xmxsfc1.or.
c     1   xeAs_B3.lt.xmnsfc3.or.xeAs_B3.gt.xmxsfc3) return  ! return to main w/o retrieval
c     print *, 'xeAs_B1,xeAs_B3',xeAs_B1,xeAs_B3,refl21,lcvr

c-----------------------------------------------------------------------

c           if(itmp.ge.100.and.itmp.lt.200.and.
c    1         jtmp.ge.200.and.jtmp.lt.300)
c      if(xeAs_B1.gt.0.and.xeAs_B3.gt.0)
c    2 print *, 'sfc refl= ',xeAs_B1,xeAs_B3,r650db,r470db
c... write outputs on fort.23 for tests !MJ@@@@@
c.    if(xeAs_B1.gt.0.and.xeAs_B1.lt.47.and.xeAs_B3.gt.0.and.
c.   1   xeAs_B3.lt.24)
c.   2write(23,2323) itmp,jtmp,xlat,xlong,sza,xthet,xphi,lcvr,
c.   3  (tmpvg(kk),kk=1,6),xeAs_B1,xeAs_B3,r650db,r470db,r412db,
c.   4  sirndvi,toa_ndvi,ioprg,lcvr
c2323  format(2(i4,2x),2(f9.4,2x),3(f8.3,2x),i3,6(2x,e14.6),
c     1       5(2x,f8.3),2(2x,f9.4),2(2x,i3))

c... Begin aerosol retrieval .....
c... getting thresholds for selecting a set of aerosol models
c... for retrieval depending on regions (08/05/2011)
!... below 37 lines --> @@new@@
      ctharr(:) = -999.0
      imodarr(:) = -999 !1:0.91, 2:0.94, 3:0.96, 4:0.995
      modfrac(:) = 0.0 ! aerosol model fraction of imodarr()+1 
      if(ioprg.eq.1) then      ! N. America
        if (season.eq.1) then
          ctharr(1)=0.2
          ctharr(2)=0.5 
          ctharr(3)=1.0
          imodarr(1)=3   
          imodarr(2)=2   
          imodarr(3)=2   
          modfrac(:)=0.5
        else 
          ctharr(1)=0.2   
          ctharr(2)=0.5
          ctharr(3)=0.8 
          imodarr(1)=4   
          imodarr(2)=4   
          imodarr(3)=2   
          modfrac(:)=0.0
          modfrac(3)=0.5
        endif
      elseif(ioprg.eq.2) then  ! China
        if (season.eq.3) then
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=3
          imodarr(2)=3
          imodarr(3)=2
          modfrac(:)=0.0
        else
          ctharr(1)=0.2  
          ctharr(2)=0.5  
          ctharr(3)=1.0  
          imodarr(1)=3   
          imodarr(2)=2   
          imodarr(3)=1   
          modfrac(1)=0.0
          modfrac(2:3)=0.5
        endif
      elseif(ioprg.eq.9) then  ! Korea and Japan
        if (season.eq.3) then
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=3
          imodarr(2)=3
          imodarr(3)=2
          modfrac(:)=0.0
        else
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.5
          imodarr(1)=3
          imodarr(2)=2
          imodarr(3)=1
          modfrac(1)=0.0
          modfrac(2:3)=0.5
        endif
      elseif(ioprg.eq.3) then  ! S. America
        if (season.eq.1.) then 
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=3   
          imodarr(2)=3   
          imodarr(3)=3   
          modfrac(:)=0.0
        else
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=2   
          imodarr(2)=2   
          imodarr(3)=2   
          modfrac(:)=0.5
        endif
      elseif(ioprg.eq.4) then  !India
        if (season.eq.3) then   ! Summer
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=3   
          imodarr(2)=2  
          imodarr(3)=2  
          modfrac(:)=0.5
        else
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=3   
          imodarr(2)=2   
          imodarr(3)=1   
          modfrac(:)=0.0
          modfrac(3)=0.5
        end if
      elseif(ioprg.eq.5) then  ! N. Africa
        ctharr(1)=0.2
        ctharr(2)=0.5
        ctharr(3)=1.0
        imodarr(1)=2   
        imodarr(2)=1   
        imodarr(3)=1   
        modfrac(1)=0.0
        modfrac(2:3)=0.5
      elseif(ioprg.eq.6) then  ! S. Africa
        if (season.eq.3.or.season.eq.4) then 
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=1   
          imodarr(2)=1   
          imodarr(3)=1   
          modfrac(:)=0.5
        elseif (season.eq.1) then
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=2
          imodarr(2)=2
          imodarr(3)=1
          modfrac(:)=0.5
          modfrac(3)=0.2
        else     
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=2   
          imodarr(2)=2   
          imodarr(3)=1   
          modfrac(:)=0.5
        end if
      elseif(ioprg.eq.7) then  ! SE Asia
        if (season.eq.1) then   ! Winter
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.0
          imodarr(1)=1   
          imodarr(2)=1   
          imodarr(3)=1   
          modfrac(:)=0.5
        else
          ctharr(1)=0.2
          ctharr(2)=0.5
          ctharr(3)=1.5
          imodarr(1)=2   
          imodarr(2)=1   
          imodarr(3)=1   
          modfrac(1:2)=0.5
          modfrac(3)=0.0
        end if
      elseif(ioprg.eq.8.or.ioprg.eq.10) then ! Europe, Central Asia         
        if (season.eq.3) then
          ctharr(1)=0.2
          ctharr(2)=0.5 
          ctharr(3)=1.0
          imodarr(1)=3   
          imodarr(2)=3   
          imodarr(3)=3   
          modfrac(:)=0.0
        else
          ctharr(1)=0.2
          ctharr(2)=0.5 
          ctharr(3)=1.0 
          imodarr(1)=3   
          imodarr(2)=3   
          imodarr(3)=2  
          modfrac(:)=0.0 
          modfrac(3)=0.5
        endif
      elseif(ioprg.eq.13) then  ! C. America
        ctharr(1)=0.2
        ctharr(2)=0.5
        ctharr(3)=1.0
        imodarr(1)=2
        imodarr(2)=1
        imodarr(3)=1
        modfrac(1)=0.0
        modfrac(2:3)=0.5
      else !default
        ctharr(1)=0.2
        ctharr(2)=0.5
        ctharr(3)=1.0
        imodarr(1)=2   
        imodarr(2)=2   
        imodarr(3)=1   
        modfrac(:)=0.0
        modfrac(3)=0.5
      endif
      
c---------------------------------------------------------------
c     Screen for pixels outside reasonable ranges of reflectance
c---------------------------------------------------------------
c... will need to re-open the following lines .... (MJ)
      if (platform .eq. 'VIIRS' .and. refl1.gt.0.0.and.refl1.lt.0.09.and.
     1    refl6.gt.0.0.and.refl6.lt.0.14) go to 11
      if (platform .eq. 'AHI' .and. refl6.gt.0.0.and.refl6.lt.0.14) go to 11          
      if (toa_ndvi.gt.0.1.and.
     1    refln21.gt.0.01.and.refln21.le.0.25) go to 11 ! remove water/bright sfc contamination
      return   ! return to main w/o retrieval
    
11    continue

c--------------------------------------------------------
c   Input surface reflectance at 470 nm
c--------------------------------------------------------
c   Input surface reflectance at 470 nm
      r470 = xeAs_B3
c
c   Input toa reflectance (L/F) at 470 nm
      refl = refl3
c     x3 = xphi                 ! ori. position (moved up)

c     Retrieving 470 nm AOT
c     imod = 1                  ! w0 = 0.91
c     imod = 2                  ! w0 = 0.94
c     imod = 3                  ! w0 = 0.96
c     imod = 4                  ! w0 = 0.995
!c... below 40 lines --> AOT-dependent aerosol model selection (new@@@)
      As_21=-999.0
      tau_x470    = -999.
      tau_x470_1   = -999.
      tau_x470_2   = -999.
      tau_x470_3   = -999.
      tau_x470_flag=0
      tau_x650_flag=0
      tau_x470_flag1=0
      tau_x470_flag2=0
      tau_x470_flag3=0
      cth0=ctharr(1)
      cth1=ctharr(2)
      cth2=ctharr(3)
c     print *, ctharr
c     print *, cth0, cth1, cth2
      imod=imodarr(1)

      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1     imod,r470,tau_x470_1,tau_x470_flag1,trflg,0.0,debug)!model_frac

      go to 333 ! skip 2.1 um SFC re-calculation !jlee added

c... **********************************************
c... **********************************************
c.... new addition: calculate 2.1um sfc refl. using
c...  first guess of 470nm AOT
c...  Aug. 12, 2011

      tau_x4701_00=tau_x470_1
      xeAs_B3_00=xeAs_B3
      tref21in=xnvalm6(4)
      alpest=1.5
      aot21est=tau_x470_1*(466.0/2110.)**alpest
      if(aot21est.lt.0.0.or.aot21est.gt.0.5) return  ! no retr. for too-heavy aerosol
      if(aot21est.ge.0.0.and.aot21est.lt.0.06) go to 333

c.
      call calc_sfc21(x1,x2,x3,tref21in,aot21est,As_21)
c.    
c...  refl21r=3.14159*As_21/100.0/xmu
      refl21r=As_21/100.0
c      print *, 'refl21r: ', refl21r, As_21
      if(refl21r.lt.0) return

      call get_sfcrfl_veg(season,lcvr,refl21r,rc_ndvi,xeAs_B1,xeAs_B3)
      if(xeAs_B1.lt.xmnsfc1.or.xeAs_B1.gt.xmxsfc1.or.
     1   xeAs_B3.lt.xmnsfc3.or.xeAs_B3.gt.xmxsfc3) return  ! return to main w/o retrieval
      r470=xeAs_B3
c... **********************************************
c... **********************************************
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1     imod,r470,tau_x470_1,tau_x470_flag1,trflg,0.0,debug)

333   continue
c... ********
c      if(itmp.ge.800.and.itmp.lt.820.and.
c     1         jtmp.ge.900.and.jtmp.lt.920)
c     2  print *,100.*refl21,As_21,aot21est,tau_x4701_00,tau_x4701,
c     3          xeAs_B3_00, xeAs_B3 
c... ********

      if(xeAs_B1.lt.xmnsfc1.or.xeAs_B1.gt.xmxsfc1.or.
     1   xeAs_B3.lt.xmnsfc3.or.xeAs_B3.gt.xmxsfc3) return  ! return to main w/o retrieval

! 16 January 2018 JLee
      dflag = .false.
      trflg = 1 ! default AOD = 0.02 if TOA reflectance < min(lut)
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1              imodarr(1),r470,tau_x470_1,tau_x470_flag1,trflg,modfrac(1),debug)
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1              imodarr(2),r470,tau_x470_2,tau_x470_flag2,trflg,modfrac(2),debug)
      call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1              imodarr(3),r470,tau_x470_3,tau_x470_flag3,trflg,modfrac(3),debug)

      tau_x470 = tau_x470_1
      tau_x470_flag = tau_x470_flag1

      if(tau_x470_1.ge.cth0.and.tau_x470_1.lt.cth1) then
        f_1=(tau_x470_1-cth0)/(cth1-cth0)
        tau_x470=(1.0-f_1)*tau_x470_1+f_1*tau_x470_2
        tau_x470_flag=tau_x470_flag2
      endif

      if(tau_x470_1.ge.cth1.and.tau_x470_1.lt.cth2) then
        f_1=(tau_x470_2-cth1)/(cth2-cth1)
        tau_x470=(1.0-f_1)*tau_x470_2+f_1*tau_x470_3
        tau_x470_flag=tau_x470_flag3
      endif

      if(tau_x470_1.ge.cth2) then
        tau_x470=tau_x470_3
        tau_x470_flag=tau_x470_flag3
      endif

      ! aod using 2.2 um database
      if (do_sv.and.r470sv.gt.0.0.and.r470sv.lt.24.0) then
        if (r470sv.lt.1.0) r470sv = 1.0
        call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1       imodarr(1),r470sv,tau_x470sv_1,tau_x470sv_flag1,trflg,modfrac(1),debug)
        call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1       imodarr(2),r470sv,tau_x470sv_2,tau_x470sv_flag2,trflg,modfrac(2),debug)
        call aero_470(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1       imodarr(3),r470sv,tau_x470sv_3,tau_x470sv_flag3,trflg,modfrac(3),debug)
  
        tau_x470sv = tau_x470sv_1
        tau_x470sv_flag = tau_x470sv_flag1
  
        if(tau_x470sv_1.ge.cth0.and.tau_x470sv_1.lt.cth1) then
          f_1=(tau_x470sv_1-cth0)/(cth1-cth0)
          tau_x470sv=(1.0-f_1)*tau_x470sv_1+f_1*tau_x470sv_2
          tau_x470sv_flag=tau_x470sv_flag2
        endif
  
        if(tau_x470sv_1.ge.cth1.and.tau_x470sv_1.lt.cth2) then
          f_1=(tau_x470sv_2-cth1)/(cth2-cth1)
          tau_x470sv=(1.0-f_1)*tau_x470sv_2+f_1*tau_x470sv_3
          tau_x470sv_flag=tau_x470sv_flag3
        endif
  
        if(tau_x470sv_1.ge.cth2) then
          tau_x470sv=tau_x470sv_3
          tau_x470sv_flag=tau_x470sv_flag3
        endif
      endif
! end JLee added

!      if(tau_x4701.ge.0.0.and.tau_x4701.lt.cth0) then
!         tau_x470=tau_x4701
!         tau_x470_flag=tau_x470_flag1
!      elseif(tau_x4701.ge.cth0.and.tau_x4701.lt.cth1) then
!         imod=imodarr(2)
!         call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
!     1      imod,r470,tau_x4702,tau_x470_flag2)
!         f_1=(tau_x4701-cth0)/(cth1-cth0)
!         tau_x470=f_1*tau_x4702+(1.0-f_1)*tau_x4701
!         tau_x470_flag=tau_x470_flag2
!      elseif(tau_x4701.ge.cth1.and.tau_x4701.lt.cth2) then
!         imod=imodarr(2)
!         call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
!     1      imod,r470,tau_x4702,tau_x470_flag2)
!         imod=imodarr(3)
!         call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
!     1      imod,r470,tau_x4703,tau_x470_flag3)
!         f_2=(tau_x4701-cth1)/(cth2-cth1)
!         tau_x470=f_2*tau_x4703+(1.0-f_2)*tau_x4702
!         tau_x470_flag=tau_x470_flag3
!      elseif(tau_x4701.ge.cth2) then
!         imod=imodarr(3)
!         call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
!     1      imod,r470,tau_x4703,tau_x470_flag3)
!         tau_x470=tau_x4703
!         tau_x470_flag=tau_x470_flag3
!      else
!         tau_x470=tau_x4701
!         tau_x470_flag=tau_x470_flag1
!c        tau_x470=-999.0
!      endif

c;--------------------------------------------<
c.   original call for 470nm AOT retrieval (before 08/05/2011)
c.       imod=3
c.       call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma, ! commented out
c.   1       imod,r470,tau_x470,tau_x470_flag)       ! 08/05/2011

!c. @@@test@@@--> 5lines
      if(tau_x470.gt.5.or.tau_x470.lt.0) then
          tau_x470=-999.
          tau_x650=-999.
          return
      endif
c     print *, 'hereC'
c      print *,'tau_x470 =', tau_x470


c--------------------------------------------------------
c          Retrieving 650 nm AOT --> added Feb 26, 2010
c--------------------------------------------------------
c   Input surface reflectance at 650 nm
      r650 = xeAs_B1

!jlee added
!      if (test_type.eq.1) then 
!        r650 = xeAs_B1
!      elseif (test_type.eq.2) then
!        r650 = xeAs_B1 - 4.67*tau_x470
!      elseif (test_type.eq.3) then
!        r650 = xeAs_B1 
!        if (rc_ndvi.lt.0.4) then
!          r650 = xeAs_B1 - 4.67*tau_x470
!        endif
!      elseif (test_type.eq.4) then
!        r650 = xeAs_B1
!        if (rc_ndvi.lt.0.4.and.lcvr.eq.12) then
!          r650 = xeAs_B1 - 4.67*tau_x470
!        endif
!      else
!      endif
!end jlee added
c
c   Input toa reflectance (L/F) at 650 nm
      refl = refl6

c     Retrieving 650 nm AOT
      call aero_650(dflag,refl,x1,x2,x3,mm,nn,ll,ma,
     1     r650,tau_x650,tau_x650_flag,0,0.0,0.0,0,trflg)!tau_x470_flag,tau_x412,tau_x470,tau_x412_flag_91

!      call aero_650veg(refl,x1,x2,x3,mm,nn,ll,ma,
!     1    r650,tau_x650,tau_x650_flag,
!     2    tau_x470_flag,tau_x470)

c      print *, 'hereD'
c      print *,'tau_x650 =', tau_x650, r650

!c. @@@test@@@--> 5lines
      if(tau_x650.gt.5.or.tau_x650.lt.0) then
          tau_x470=-999.
          tau_x650=-999.
          return
      endif

!     -- only set bright surface flag over N.America (ioprg==1) and 
!     -- croplands (lcvr == 12)
      if (ioprg .eq. 1 .and. lcvr .eq. 12) then
        if (xeAs_B1.gt.17.or.xeAs_B3.gt.17) then
          bflag = 1
        endif
      endif

c--------------------------------------------------------

c   -- tweak AOT over S.Africa, Mongu site.
      if (ioprg .eq. 6) then
        if (season .eq .3) then
          tau_x470 = tau_x470 * 1.25
          tau_x650 = tau_x650 * 1.25
        end if
!        if (month .ge. 9 .AND. month .le. 11) then
!          tau_x470 = tau_x470 * 1.2
!          tau_x650 = tau_x650 * 1.2
!        end if
      end if

c   -- tweak over Mukdahan, SE Asia site.
      if (ioprg .eq. 7) then
!        if (month .ge. 3 .AND. month .le. 5) then
!          if (tau_x470 .gt. 0.7) then
!            tau_x470 = tau_x470 * 1.1
!            tau_x650 = tau_x650 * 1.1
!          end if
!        end if
        
!        if (season .eq. 1) then !DJF
!          tau_x470 = tau_x470 * 1.1
!          tau_x650 = tau_x650 * 1.1
!        end if

!        if (month .ge. 9 .AND. month .le. 11) then
!          tau_x470 = tau_x470 * 1.1
!          tau_x650 = tau_x650 * 1.1
!        end if
      end if

! jlee modified 650->672, 466->488      
c... original set up (temporary)
c...==============================================
c... MJ re-defines output

      alpha=alog10(tau_x470/tau_x650)/alog10(672./488.) 

c     if(tau_x470.lt.0.0) tau_x470=-999.0
c     if(tau_x650.lt.0.0) tau_x650=-999.0
c     if(tau_x470.lt.0.0.and.tau_x650.lt.0.0) alpha=-999.0
c     if(alpha.ge.-0.5.and.alpha.le.3.5.and.tau_x470.gt.0.and.
c    &   tau_x650.gt.0.) then

      if (tau_x470.lt.1.0) then 
        if(alpha.ge.-0.5.and.alpha.le.3.5) then
           tau550=tau_x470*(488.0/500.0)**alpha
c        xtau(1)=tau_x470*(466.0/412.0)**alpha  ! AOT 412nm (extrapol.)
        else
          if((tau_x470_flag.eq.1.or.tau_x470_flag.eq.-10).and.
     1       (tau_x650_flag.eq.1.or.tau_x650_flag.eq.-10)) then
            alpha=1.0
            tau550=(tau_x470+tau_x650)/2.0
            tau_x470=tau550*(500./488.)**alpha
            tau_x650=tau550*(500./672.)**alpha
          elseif((tau_x470_flag.eq.1.or.tau_x470_flag.eq.-10).and.
     1            tau_x650_flag.eq.0) then
            alpha=1.0
            tau550=tau_x650*(672./500.)**alpha
            tau_x470=tau_x650*(672./488.)**alpha
          elseif(tau_x470_flag.eq.0.and.
     1          (tau_x650_flag.eq.1.or.tau_x650_flag.eq.-10)) then
            alpha=1.0
            tau550=tau_x470*(488./500.)**alpha
            tau_x650=tau_x470*(488./672.)**alpha
          else
            alpha=1.0
            tau550=(tau_x470+tau_x650)/2.0
            tau_x470=tau550*(500./488.)**alpha
            tau_x650=tau550*(500./672.)**alpha
          endif
        endif
      else ! if (tau_x470.lt.1.0) !start jlee added 
        if(alpha.ge.-0.5.and.alpha.le.3.5.and.
     1     tau_x470_flag.eq.0.and.tau_x650_flag.eq.0) then 
           tau550=tau_x470*(488.0/500.0)**alpha
c        xtau(1)=tau_x470*(466.0/412.0)**alpha  ! AOT 412nm (extrapol.)
        else
          if((tau_x470_flag.eq.1.or.tau_x470_flag.eq.-10).and.
     1       (tau_x650_flag.eq.1.or.tau_x650_flag.eq.-10)) then
            alpha=1.0
            tau550=(tau_x470+tau_x650)/2.0
            tau_x470=tau550*(500./488.)**alpha
            tau_x650=tau550*(500./672.)**alpha
          elseif(tau_x470_flag.eq.-10.and.tau_x650_flag.eq.0) then
            alpha=1.0
            tau550=tau_x650*(672./500.)**alpha
            tau_x470=tau_x650*(672./488.)**alpha
          elseif((tau_x470_flag.eq.1.or.tau_x470_flag.eq.0).and.
     1            tau_x650_flag.eq.0) then 
            tau550=tau_x470 
            alpha = 1.8 !thick smoke 
            tau_x470=tau550*(500./488.)**alpha
            tau_x650=tau550*(500./672.)**alpha                      
          elseif(tau_x470_flag.eq.0.and.
     1          (tau_x650_flag.eq.1.or.tau_x650_flag.eq.-10)) then
            alpha=1.0
            tau550=tau_x470*(488./500.)**alpha
            tau_x650=tau_x470*(488./672.)**alpha
          else
            alpha=1.0
            tau550=(tau_x470+tau_x650)/2.0
            tau_x470=tau550*(500./488.)**alpha
            tau_x650=tau550*(500./672.)**alpha
          endif
        endif
      endif !if (tau_x470.lt.1.5) ! end jlee added

      if (do_sv) then
        tau550=tau_x470sv*(488./500.)**alpha
        tau_x470 = tau_x470sv
        tau_x470_flag = tau_x470sv_flag
      endif
c...==============================================
!end jlee modified

c...==========================================================
c... MJ re-defines output
c     if(tau_x470.lt.0.0) tau_x470=-999.0
c     if(tau_x650.lt.0.0) tau_x650=-999.0
c     if(tau_x470.lt.0.0.or.tau_x650.lt.0.0) return
c-------------------------------------------------------------
c Set output buffer --> based on find_v.f
c-------------------------------------------------------------
!     22 December 2017 decided not to report Veg SSA
!     To do: report SSA as a diagnostic variable in the extended output
!      read(w0_name470(imod),'(f4.2)') ssatmp
!      ssa(2)=ssatmp
!      ssa(3)=0.976
      xtau(1)=-999.0  ! fillvalue at 412nm for over-veg retr.
      xtau(2)=tau_x470
      xtau(3)=tau_x650
      do i=1,3
        outbufvg(i) = xtau(i)
        outbufvg(i+3) = ssa(i)
      enddo
      outbufvg(7) = tau550
      outbufvg(8) = alpha
      outbufvg(9) =  1.0*tau_x470_flag  ! used to be "r412"
      outbufvg(10) = 1.0*tau_x650_flag
      outbufvg(11) = r470   ! over-veg. sfc. refl.
      outbufvg(12) = r650   ! over-veg. sfc. refl.
      outbufvg(13) = xthet
      outbufvg(14) = scat_ang
      outbufvg(15) = sfc_typ

      outbufvg(16) = xlat
      outbufvg(17) = xlong
      outbufvg(18) = sirndvi
      outbufvg(19) = rc_ndvi
      outbufvg(20) = xreg_id
      outbufvg(21) = bflag

c     print *,'outbufvg(1,2,3) ',outbufvg(1),outbufvg(2),outbufvg(3)
      return
      end
c... --------------------------------------
c... The end of subroutine myd_aero_veg
c... (main subroutine for aerosol retr. over vegetation)
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_sfcrfl_veg(season,lcvr,refln21,rc_ndvi,
     1                          xeAs_B1,xeAs_B3)
c...
c...  subroutine to estimate surface reflectance at 470nm & 650nm over
c...  vegetated surfaces using 2.1um TOA reflectance
c...     program written by Myeong-Jae Jeong (MJ)
c...     last modified Jun 15, 2010
c...
c...     Inputs
c...        season:   season ID (1: MAM; 2: JJA; 3: SON; 4: DJF)
c...        refln21:  2.1um TOA reflectance (normalized, Ipi/fcos(sza))
c...        xnvalm6: TOA reflectances at other bands
c...        lcvr:    MODIS land cover (IGBP)
c...
c...     Output
c...        xeAs_B1: estimated surface refl. for band 1 (LER in %)
c...        xeAs_B3: estimated surface refl. for band 3 (LER in %)
c...
c... 02/11/2014, VIIRS surface parameterization for general vegetation areas, modified by J. Lee

      integer season

      sb7=100.0*refln21 

      if(season.eq.2) then  ! MAM 2012

        if((lcvr.gt.0.and.lcvr.le.12).or.lcvr.eq.14) then
           xeAs_B1 = 0.283 + 0.509*sb7 + 0.0038*sb7**2.0
           sb1 = xeAs_B1
           xeAs_B3 = 0.810 + 0.535*sb1
        elseif(lcvr.eq.13) then
           xeAs_B1 = 0.762*sb7
           sb1 = xeAs_B1
           xeAs_B3 = 0.718*sb1
        else
           return
        endif

      elseif(season.eq.3) then  ! JJA 2012

        if((lcvr.gt.0.and.lcvr.le.12).or.lcvr.eq.14) then
           xeAs_B1 = 0.295 + 0.405*sb7 + 0.0095*sb7**2.0
           sb1 = xeAs_B1
           xeAs_B3 = 0.605 + 0.507*sb1
        elseif(lcvr.eq.13) then
           xeAs_B1 = 0.765*sb7
           sb1 = xeAs_B1
           xeAs_B3 = 0.664*sb1
        else
           return
        endif

      elseif(season.eq.4) then  ! SON 2012

        if((lcvr.gt.0.and.lcvr.le.12).or.lcvr.eq.14) then
           xeAs_B1 =  0.298 + 0.430*sb7 + 0.0084*sb7**2.0
           sb1 = xeAs_B1
           xeAs_B3 = 0.494 + 0.524*sb1
        elseif(lcvr.eq.13) then
           xeAs_B1 = 0.784*sb7
           sb1 = xeAs_B1
           xeAs_B3 = 0.665*sb1
        else
           return
        endif

      elseif(season.eq.1) then  ! DJF 2013

        if((lcvr.gt.0.and.lcvr.le.12).or.lcvr.eq.14) then
           xeAs_B1 = 0.283 + 0.509*sb7 + 0.0038*sb7**2.0
           sb1 = xeAs_B1
           xeAs_B3 = 0.810 + 0.535*sb1
        elseif(lcvr.eq.13) then
           xeAs_B1 = 0.762*sb7
           sb1 = xeAs_B1
           xeAs_B3 = 0.718*sb1
        else
           return
        endif

      else
         return
      endif

      return
      end
c... --------------------------------------
c... The end of subroutine get_sfcrfl_veg
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_sfcrfl_veg_high(season,lcvr,refln21,rc_ndvi,
     1                               xeAs_B1,xeAs_B3)
c...
c...  subroutine to estimate surface reflectance at 470nm & 650nm over
c...  vegetated surfaces using 2.1um TOA reflectance
c...     program written by Myeong-Jae Jeong (MJ)
c...     last modified Jun 15, 2010
c...
c...     Inputs
c...        season:   season ID (1: DJF; 2: MAM; 3: JJA; 4: SON)
c...        refln21:  2.1um TOA reflectance
c...        xnvalm6: TOA reflectances at other bands
c...        lcvr:    MODIS land cover (IGBP)
c...
c...     Output
c...        xeAs_B1: estimated surface refl. for band 1 (LER in %)
c...        xeAs_B3: estimated surface refl. for band 3 (LER in %)
c...
c... 02/11/2014, VIIRS surface parameterization for general vegetation areas, modified by J. Lee

      integer season

      sb7=100.0*refln21

      if(season.eq.2) then  ! MAM 2012

        if((lcvr.gt.0.and.lcvr.le.12).or.lcvr.eq.14) then
           xeAs_B1 = -0.542 +0.575*sb7 +0.0022*sb7**2 
           sb1 = xeAs_B1
           xeAs_B3 = +0.342 +0.598*sb1
        elseif(lcvr.eq.13) then
           xeAs_B1 = -1.201 +0.785*sb7
           sb1 = xeAs_B1
           xeAs_B3 = +0.401 +0.679*sb1
        else
           return
        endif

      elseif(season.eq.3) then  ! JJA 2012

        if((lcvr.gt.0.and.lcvr.le.12).or.lcvr.eq.14) then
           xeAs_B1 = -0.039 +0.401*sb7 +0.0108*sb7**2
           sb1 = xeAs_B1
           xeAs_B3 = +0.568 +0.572*sb1
        elseif(lcvr.eq.13) then
           xeAs_B1 = -2.147 +0.976*sb7
           sb1 = xeAs_B1
           xeAs_B3 = +0.713 +0.620*sb1
        else
           return
        endif

      elseif(season.eq.4) then  ! SON 2012

        if((lcvr.gt.0.and.lcvr.le.12).or.lcvr.eq.14) then
           xeAs_B1 = -0.125 +0.466*sb7 +0.0076*sb7**2
           sb1 = xeAs_B1
           xeAs_B3 = +0.438 +0.518*sb1
        elseif(lcvr.eq.13) then
           xeAs_B1 = -1.428 +0.839*sb7
           sb1 = xeAs_B1
           xeAs_B3 = +0.233 +0.648*sb1
        else
           return
        endif

      elseif(season.eq.1) then  ! DJF 2013

        if((lcvr.gt.0.and.lcvr.le.12).or.lcvr.eq.14) then
           xeAs_B1 = -0.542 +0.575*sb7 +0.0022*sb7**2
           sb1 = xeAs_B1
           xeAs_B3 = +0.342 +0.598*sb1
        elseif(lcvr.eq.13) then
           xeAs_B1 = -1.201 +0.785*sb7
           sb1 = xeAs_B1
           xeAs_B3 = +0.401 +0.679*sb1
        else
           return
        endif

      else
         return
      endif

      return
      end
c... --------------------------------------
c... The end of subroutine get_sfcrfl_veg_high
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_sfcrfl_swir_vis(season,lcvr,refln21,rc_ndvi,
     1                               xeAs_B1,xeAs_B3)
c...
c...  subroutine to estimate surface reflectance at 470nm & 650nm over
c...  vegetated surfaces using 2.25 um TOA reflectance and Rayleigh-corrected NDVI
c...     program written by Jaehwa Lee 
c...     last modified March 24, 2014
c...
c...     Inputs
c...        season:   season ID (1: DJF; 2: MAM; 3: JJA; 4: SON)
c...        lcvr:    MODIS land cover (IGBP)
c...        refln21: 2.1um TOA reflectance
c...        rc_ndvi: Rayleigh-corrected NDVI
c...
c...     Output
c...        xeAs_B1: estimated surface refl. for band 1 (LER in %)
c...        xeAs_B3: estimated surface refl. for band 3 (LER in %)
c...
c...  03/24/2014, Coeffs. w.r.t. Rayleigh-corrected NDVI are applied.
      integer season

      real, dimension(20) :: coeff0_red
      real, dimension(20) :: coeff1_red
      real, dimension(20) :: coeff0_blue
      real, dimension(20) :: coeff1_blue

      rc_ndvi = max(rc_ndvi, 0.1)
      rc_ndvi = min(rc_ndvi, 1.0)
      sb7=100.0*refln21

      if(season.eq.1.or.season.eq.2) then  ! DJF, MAM

        if(lcvr.gt.0.and.lcvr.ne.7.and.lcvr.ne.13) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.688,  0.688,  0.688,  0.679,  0.632,
     1                          0.604,  0.592,  0.578,  0.561,  0.538,
     1                          0.519,  0.508,  0.502,  0.491,  0.468,
     1                          0.426,  0.382,  0.337,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.522,  0.522,  0.522,  0.534,  0.530,
     1                          0.551,  0.567,  0.596,  0.618,  0.633,
     1                          0.651,  0.676,  0.703,  0.734,  0.771,
     1                          0.818,  0.885,  0.906,  0.000,  0.000/)

          ndvi_idx = min(floor(rc_ndvi*20)+1, 18)

          xeAs_B1 = coeff0_red(ndvi_idx) + coeff1_red(ndvi_idx)*sb7
          sb1 = xeAs_B1
          xeAs_B3 = coeff0_blue(ndvi_idx) + coeff1_blue(ndvi_idx)*sb1

!          ndvi_idx = min(floor((rc_ndvi-0.025)*20+1), 18)
!          coeff0_red_fin = coeff0_red(ndvi_idx)
!     1    +(coeff0_red(ndvi_idx+1)-coeff0_red(ndvi_idx))/0.05
!     1    *(rc_ndvi-(0.025+(ndvi_idx-1)*0.05))
!          coeff1_red_fin = coeff1_red(ndvi_idx)
!     1    +(coeff1_red(ndvi_idx+1)-coeff1_red(ndvi_idx))/0.05
!     1    *(rc_ndvi-(0.025+(ndvi_idx-1)*0.05))
!          coeff0_blue_fin = coeff0_blue(ndvi_idx)
!     1    +(coeff0_blue(ndvi_idx+1)-coeff0_blue(ndvi_idx))/0.05
!     1    *(rc_ndvi-(0.025+(ndvi_idx-1)*0.05))
!          coeff1_blue_fin = coeff1_blue(ndvi_idx)
!     1   +(coeff1_blue(ndvi_idx+1)-coeff1_blue(ndvi_idx))/0.05
!     1   *(rc_ndvi-(0.025+(ndvi_idx-1)*0.05))
!
!          xeAs_B1 = coeff0_red_fin + coeff1_red_fin*sb7
!          sb1 = xeAs_B1
!          xeAs_B3 = coeff0_blue_fin + coeff1_blue_fin*sb1

        elseif(lcvr.eq.7) then !open shrublands
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.735,  0.735,  0.735,  0.724,  0.716,
     1                          0.709,  0.692,  0.647,  0.601,  0.565,
     1                          0.555,  0.535,  0.537,  0.513,  0.498,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.527,  0.527,  0.527,  0.508,  0.526,
     1                          0.540,  0.555,  0.555,  0.606,  0.631,
     1                          0.744,  0.800,  0.809,  0.810,  0.813,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)

          ndvi_idx = min(floor(rc_ndvi*20)+1, 15)

          xeAs_B1 = coeff0_red(ndvi_idx) + coeff1_red(ndvi_idx)*sb7
          sb1 = xeAs_B1
          xeAs_B3 = coeff0_blue(ndvi_idx) + coeff1_blue(ndvi_idx)*sb1

        elseif(lcvr.eq.13) then !urban and built-up
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.903,  0.903,  0.903,  0.860,  0.818,
     1                          0.751,  0.721,  0.677,  0.642,  0.627,
     1                          0.604,  0.581,  0.561,  0.530,  0.483,
     1                          0.436,  0.404,  0.000,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.693,  0.693,  0.693,  0.729,  0.735,
     1                          0.725,  0.728,  0.725,  0.729,  0.738,
     1                          0.749,  0.768,  0.793,  0.805,  0.818,
     1                          0.838,  0.862,  0.000,  0.000,  0.000/)

          ndvi_idx = min(floor(rc_ndvi*20)+1, 17)

          xeAs_B1 = coeff0_red(ndvi_idx) + coeff1_red(ndvi_idx)*sb7
          sb1 = xeAs_B1
          xeAs_B3 = coeff0_blue(ndvi_idx) + coeff1_blue(ndvi_idx)*sb1

        else
           return
        endif

      elseif(season.eq.3) then  ! JJA

        if(lcvr.gt.0.and.lcvr.ne.7.and.lcvr.ne.13) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.692,  0.692,  0.692,  0.681,  0.684,
     1                          0.653,  0.626,  0.602,  0.593,  0.580,
     1                          0.564,  0.540,  0.519,  0.496,  0.471,
     1                          0.441,  0.395,  0.332,  0.260,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.522,  0.522,  0.522,  0.516,  0.490,
     1                          0.487,  0.499,  0.516,  0.530,  0.541,
     1                          0.556,  0.583,  0.620,  0.676,  0.767,
     1                          0.827,  0.899,  0.926,  0.665,  0.000/)

          ndvi_idx = min(floor(rc_ndvi*20)+1, 18)

          xeAs_B1 = coeff0_red(ndvi_idx) + coeff1_red(ndvi_idx)*sb7
          sb1 = xeAs_B1
          xeAs_B3 = coeff0_blue(ndvi_idx) + coeff1_blue(ndvi_idx)*sb1

        elseif(lcvr.eq.7) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.735,  0.735,  0.735,  0.724,  0.716,
     1                          0.709,  0.692,  0.647,  0.601,  0.565,
     1                          0.555,  0.535,  0.537,  0.513,  0.498,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.527,  0.527,  0.527,  0.508,  0.526,
     1                          0.540,  0.555,  0.555,  0.606,  0.631,
     1                          0.744,  0.800,  0.809,  0.810,  0.813,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)

          ndvi_idx = min(floor(rc_ndvi*20)+1, 15)

          xeAs_B1 = coeff0_red(ndvi_idx) + coeff1_red(ndvi_idx)*sb7
          sb1 = xeAs_B1
          xeAs_B3 = coeff0_blue(ndvi_idx) + coeff1_blue(ndvi_idx)*sb1

        elseif(lcvr.eq.13) then !urban and built-up
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.901,  0.901,  0.901,  0.874,  0.844,
     1                          0.811,  0.802,  0.778,  0.742,  0.708,
     1                          0.677,  0.646,  0.615,  0.574,  0.526,
     1                          0.481,  0.427,  0.362,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.653,  0.653,  0.653,  0.688,  0.698,
     1                          0.692,  0.699,  0.700,  0.701,  0.695,
     1                          0.710,  0.770,  0.819,  0.847,  0.865,
     1                          0.895,  0.951,  1.045,  0.000,  0.000/)

          ndvi_idx = min(floor(rc_ndvi*20)+1, 18)

          xeAs_B1 = coeff0_red(ndvi_idx) + coeff1_red(ndvi_idx)*sb7
          sb1 = xeAs_B1
          xeAs_B3 = coeff0_blue(ndvi_idx) + coeff1_blue(ndvi_idx)*sb1

        else
           return
        endif

      elseif(season.eq.4) then  ! SON

        if(lcvr.gt.0.and.lcvr.ne.7.and.lcvr.ne.13) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.686,  0.686,  0.686,  0.677,  0.610,
     1                          0.572,  0.564,  0.559,  0.549,  0.534,
     1                          0.526,  0.517,  0.511,  0.500,  0.475,
     1                          0.439,  0.397,  0.337,  0.238,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.506,  0.506,  0.506,  0.489,  0.485,
     1                          0.490,  0.497,  0.519,  0.560,  0.587,
     1                          0.605,  0.634,  0.668,  0.703,  0.736,
     1                          0.782,  0.843,  0.859,  0.622,  0.000/)

          ndvi_idx = min(floor(rc_ndvi*20)+1, 18)

          xeAs_B1 = coeff0_red(ndvi_idx) + coeff1_red(ndvi_idx)*sb7
          sb1 = xeAs_B1
          xeAs_B3 = coeff0_blue(ndvi_idx) + coeff1_blue(ndvi_idx)*sb1

        elseif(lcvr.eq.7) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.735,  0.735,  0.735,  0.724,  0.716,
     1                          0.709,  0.692,  0.647,  0.601,  0.565,
     1                          0.555,  0.535,  0.537,  0.513,  0.498,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.527,  0.527,  0.527,  0.508,  0.526,
     1                          0.540,  0.555,  0.555,  0.606,  0.631,
     1                          0.744,  0.800,  0.809,  0.810,  0.813,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)

          ndvi_idx = min(floor(rc_ndvi*20)+1, 15)

          xeAs_B1 = coeff0_red(ndvi_idx) + coeff1_red(ndvi_idx)*sb7
          sb1 = xeAs_B1
          xeAs_B3 = coeff0_blue(ndvi_idx) + coeff1_blue(ndvi_idx)*sb1

        elseif(lcvr.eq.13) then !urban and built-up
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.854,  0.854,  0.854,  0.846,  0.824,
     1                          0.780,  0.746,  0.699,  0.664,  0.647,
     1                          0.632,  0.606,  0.576,  0.540,  0.510,
     1                          0.481,  0.425,  0.000,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.700,  0.700,  0.700,  0.729,  0.722,
     1                          0.709,  0.709,  0.711,  0.711,  0.732,
     1                          0.751,  0.769,  0.789,  0.808,  0.829,
     1                          0.858,  0.917,  0.000,  0.000,  0.000/)

          ndvi_idx = min(floor(rc_ndvi*20)+1, 17)

          xeAs_B1 = coeff0_red(ndvi_idx) + coeff1_red(ndvi_idx)*sb7
          sb1 = xeAs_B1
          xeAs_B3 = coeff0_blue(ndvi_idx) + coeff1_blue(ndvi_idx)*sb1

        else
           return
        endif

      else
         return
      endif

      return
      end

c... --------------------------------------
c... The end of subroutine get_sfcrfl_swir_vis
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine get_sfcrfl_nir_vis(season,lcvr,refn865_rc,rc_ndvi,
     1                              xeAs_B1,xeAs_B3)
c...
c...  subroutine to estimate surface reflectance at 470nm & 650nm over
c...  vegetated surfaces using 0.86 um TOA reflectance
c...     program written by Jaehwa Lee 
c...     last modified March 24, 2014
c...
c...     Inputs
c...        season:      season ID (1: DJF; 2: MAM; 3: JJA; 4: SON)
c...        lcvr:       MODIS land cover (IGBP)
c...        refn865_rc: Rayleigh-corrected 865 nm reflectance
c...        rc_ndvi:    Rayleigh-corrected NDVI
c...
c...     Output
c...        xeAs_B1: estimated surface refl. for band 1 (LER in %)
c...        xeAs_B3: estimated surface refl. for band 3 (LER in %)
c...
c...  03/24/2014, Coeffs. w.r.t. Rayleigh-corrected NDVI are applied.
      integer season

      real, dimension(20) :: coeff0_red
      real, dimension(20) :: coeff1_red
      real, dimension(20) :: coeff0_blue
      real, dimension(20) :: coeff1_blue

      rc_ndvi = max(rc_ndvi, 0.1)
      rc_ndvi = min(rc_ndvi, 1.0)
      sb2=100.0*refn865_rc

      if(season.eq.1.or.season.eq.2) then  ! MAM, DJF

        if(lcvr.gt.0.and.lcvr.ne.13) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.770,  0.770,  0.770,  0.681,  0.617,
     1                          0.548,  0.487,  0.430,  0.378,  0.330,
     1                          0.286,  0.244,  0.207,  0.171,  0.135,
     1                          0.102,  0.075,  0.056,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.406,  0.406,  0.406,  0.379,  0.344,
     1                          0.308,  0.277,  0.251,  0.226,  0.204,
     1                          0.183,  0.163,  0.144,  0.124,  0.102,
     1                          0.081,  0.064,  0.050,  0.000,  0.000/)

          ndvi_idx = min(floor((rc_ndvi-0.025)*20+1), 18)
          coeff0_red_fin = coeff0_red(ndvi_idx)
     1    +(coeff0_red(ndvi_idx+1)-coeff0_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_red_fin = coeff1_red(ndvi_idx)
     1    +(coeff1_red(ndvi_idx+1)-coeff1_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025)) 
          coeff0_blue_fin = coeff0_blue(ndvi_idx)
     1    +(coeff0_blue(ndvi_idx+1)-coeff0_blue(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025)) 
          coeff1_blue_fin = coeff1_blue(ndvi_idx)
     1   +(coeff1_blue(ndvi_idx+1)-coeff1_blue(ndvi_idx))/0.05
     1   *(rc_ndvi-(ndvi_idx*0.05-0.025))

          xeAs_B1 = coeff0_red_fin + coeff1_red_fin*sb2
          xeAs_B3 = coeff0_blue_fin + coeff1_blue_fin*sb2
        elseif(lcvr.eq.13) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.739,  0.739,  0.739,  0.673,  0.604,
     1                          0.537,  0.478,  0.422,  0.372,  0.326,
     1                          0.282,  0.236,  0.197,  0.162,  0.129,
     1                          0.094,  0.000,  0.000,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.516,  0.516,  0.516,  0.482,  0.438,
     1                          0.386,  0.344,  0.304,  0.269,  0.240,
     1                          0.212,  0.184,  0.157,  0.130,  0.103,
     1                          0.080,  0.000,  0.000,  0.000,  0.000/)

          ndvi_idx = min(floor((rc_ndvi-0.025)*20+1), 16)
          coeff0_red_fin = coeff0_red(ndvi_idx)
     1    +(coeff0_red(ndvi_idx+1)-coeff0_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_red_fin = coeff1_red(ndvi_idx)
     1    +(coeff1_red(ndvi_idx+1)-coeff1_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff0_blue_fin = coeff0_blue(ndvi_idx)
     1    +(coeff0_blue(ndvi_idx+1)-coeff0_blue(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_blue_fin = coeff1_blue(ndvi_idx)
     1   +(coeff1_blue(ndvi_idx+1)-coeff1_blue(ndvi_idx))/0.05
     1   *(rc_ndvi-(ndvi_idx*0.05-0.025))

          xeAs_B1 = coeff0_red_fin + coeff1_red_fin*sb2
          xeAs_B3 = coeff0_blue_fin + coeff1_blue_fin*sb2
        else
           return
        endif

      elseif(season.eq.3) then  ! JJA

        if(lcvr.gt.0.and.lcvr.ne.13) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.764,  0.764,  0.764,  0.703,  0.619,
     1                          0.552,  0.489,  0.431,  0.380,  0.331,
     1                          0.288,  0.247,  0.207,  0.171,  0.137,
     1                          0.107,  0.080,  0.058,  0.041,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.398,  0.398,  0.398,  0.358,  0.300,
     1                          0.269,  0.245,  0.223,  0.202,  0.180,
     1                          0.161,  0.144,  0.128,  0.113,  0.099,
     1                          0.084,  0.069,  0.052,  0.037,  0.000/)

          ndvi_idx = min(floor((rc_ndvi-0.025)*20+1), 19)
          coeff0_red_fin = coeff0_red(ndvi_idx)
     1    +(coeff0_red(ndvi_idx+1)-coeff0_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_red_fin = coeff1_red(ndvi_idx)
     1    +(coeff1_red(ndvi_idx+1)-coeff1_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff0_blue_fin = coeff0_blue(ndvi_idx)
     1    +(coeff0_blue(ndvi_idx+1)-coeff0_blue(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_blue_fin = coeff1_blue(ndvi_idx)
     1   +(coeff1_blue(ndvi_idx+1)-coeff1_blue(ndvi_idx))/0.05
     1   *(rc_ndvi-(ndvi_idx*0.05-0.025))

          xeAs_B1 = coeff0_red_fin + coeff1_red_fin*sb2
          xeAs_B3 = coeff0_blue_fin + coeff1_blue_fin*sb2
        elseif(lcvr.eq.13) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.745,  0.745,  0.745,  0.681,  0.610,
     1                          0.546,  0.485,  0.429,  0.373,  0.325,
     1                          0.281,  0.241,  0.204,  0.170,  0.139,
     1                          0.111,  0.087,  0.000,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.492,  0.492,  0.492,  0.471,  0.422,
     1                          0.376,  0.338,  0.300,  0.260,  0.226,
     1                          0.197,  0.183,  0.169,  0.144,  0.121,
     1                          0.100,  0.083,  0.000,  0.000,  0.000/)

          ndvi_idx = min(floor((rc_ndvi-0.025)*20+1), 17)
          coeff0_red_fin = coeff0_red(ndvi_idx)
     1    +(coeff0_red(ndvi_idx+1)-coeff0_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_red_fin = coeff1_red(ndvi_idx)
     1    +(coeff1_red(ndvi_idx+1)-coeff1_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff0_blue_fin = coeff0_blue(ndvi_idx)
     1    +(coeff0_blue(ndvi_idx+1)-coeff0_blue(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_blue_fin = coeff1_blue(ndvi_idx)
     1   +(coeff1_blue(ndvi_idx+1)-coeff1_blue(ndvi_idx))/0.05
     1   *(rc_ndvi-(ndvi_idx*0.05-0.025))

          xeAs_B1 = coeff0_red_fin + coeff1_red_fin*sb2
          xeAs_B3 = coeff0_blue_fin + coeff1_blue_fin*sb2
        else
           return
        endif

      elseif(season.eq.4) then  ! SON

        if(lcvr.gt.0.and.lcvr.ne.13) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.762,  0.762,  0.762,  0.691,  0.617,
     1                          0.549,  0.488,  0.432,  0.380,  0.332,
     1                          0.288,  0.248,  0.210,  0.175,  0.141,
     1                          0.110,  0.083,  0.063,  0.043,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.395,  0.395,  0.395,  0.342,  0.309,
     1                          0.276,  0.247,  0.223,  0.204,  0.185,
     1                          0.166,  0.150,  0.136,  0.119,  0.101,
     1                          0.085,  0.070,  0.054,  0.028,  0.000/)

          ndvi_idx = min(floor((rc_ndvi-0.025)*20+1), 19)
          coeff0_red_fin = coeff0_red(ndvi_idx)
     1    +(coeff0_red(ndvi_idx+1)-coeff0_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_red_fin = coeff1_red(ndvi_idx)
     1    +(coeff1_red(ndvi_idx+1)-coeff1_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff0_blue_fin = coeff0_blue(ndvi_idx)
     1    +(coeff0_blue(ndvi_idx+1)-coeff0_blue(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_blue_fin = coeff1_blue(ndvi_idx)
     1   +(coeff1_blue(ndvi_idx+1)-coeff1_blue(ndvi_idx))/0.05
     1   *(rc_ndvi-(ndvi_idx*0.05-0.025))

          xeAs_B1 = coeff0_red_fin + coeff1_red_fin*sb2
          xeAs_B3 = coeff0_blue_fin + coeff1_blue_fin*sb2
        elseif(lcvr.eq.13) then
          coeff0_red(1:20) = (/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_red(1:20) = (/ 0.735,  0.735,  0.735,  0.672,  0.601,
     1                          0.532,  0.473,  0.419,  0.369,  0.324,
     1                          0.283,  0.243,  0.204,  0.169,  0.141,
     1                          0.119,  0.093,  0.000,  0.000,  0.000/)
          coeff0_blue(1:20) =(/ 0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000,
     1                          0.000,  0.000,  0.000,  0.000,  0.000/)
          coeff1_blue(1:20) =(/ 0.492,  0.492,  0.492,  0.482,  0.429,
     1                          0.371,  0.330,  0.293,  0.258,  0.232,
     1                          0.209,  0.185,  0.159,  0.136,  0.119,
     1                          0.103,  0.083,  0.000,  0.000,  0.000/)

          ndvi_idx = min(floor((rc_ndvi-0.025)*20+1), 17)
          coeff0_red_fin = coeff0_red(ndvi_idx)
     1    +(coeff0_red(ndvi_idx+1)-coeff0_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_red_fin = coeff1_red(ndvi_idx)
     1    +(coeff1_red(ndvi_idx+1)-coeff1_red(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff0_blue_fin = coeff0_blue(ndvi_idx)
     1    +(coeff0_blue(ndvi_idx+1)-coeff0_blue(ndvi_idx))/0.05
     1    *(rc_ndvi-(ndvi_idx*0.05-0.025))
          coeff1_blue_fin = coeff1_blue(ndvi_idx)
     1   +(coeff1_blue(ndvi_idx+1)-coeff1_blue(ndvi_idx))/0.05
     1   *(rc_ndvi-(ndvi_idx*0.05-0.025))

          xeAs_B1 = coeff0_red_fin + coeff1_red_fin*sb2
          xeAs_B3 = coeff0_blue_fin + coeff1_blue_fin*sb2
        else
           return
        endif

      else
         return
      endif

      return
      end

c... --------------------------------------
c... The end of subroutine get_sfcrfl_nir_vis
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1    imod,r470,tau_x470,tau_x470_flag)

      use viirs_aerosol_luts, only: new_intep  

      include 'aottbl.inc'
      include 'newaottbl.inc'

c      common /angle_node/ theta0, theta, phi
c      common /sfcref_node/ sfc_ref412, sfc_ref470, sfc_ref650
c      real theta0(10), theta(46), phi(30), tau(10)
c      real sfc_ref412(20), sfc_ref470(24), sfc_ref650(24)
      real nnvalx(4,4,2,10), yy(10), yy2(8)
c      real nvalx470(10,46,30,10,4,24)
      integer tau_x470_flag, imod
      data pi   /3.14159/
c      data tau  /0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5/
      
c      print *, 'aero_470veg in: ', refl,x1,x2,x3,mm,nn,ll,ma,imod,r470,tau_x470,tau_x470_flag
      index_ii = r470
      if (index_ii.lt.0) return  ! orig.

      frac = (r470-sfc_ref470(index_ii))/
     1       (sfc_ref470(index_ii+1)-sfc_ref470(index_ii))

      if (index_ii.lt.1.or.index_ii.gt.24)
     1    print *,'index_iir470 = ', index_ii,xlat,xlong
      if (frac.lt.0.0.or.frac.gt.1.0)
     1    print *,'frac on sfc470=', frac
      if (index_ii.lt.1.or.index_ii.gt.24) then  ! MJ added
        tau_x470_flag = 9
        tau_x470=-999.0
        return             
      endif


      call search(dflag2,x3,phi,ll,ii)
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))

      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif

      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif

c      write(6,*) "frac =",frac
      do ia = 1, 10
       do i = 1, 4
        do j = 1, 4
       nnvalx(i,j,1,ia) = nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*
     1  (1.-frac) + nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac
       nnvalx(i,j,2,ia) = nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*
     1  (1.-frac) + nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac
        enddo
       enddo
      enddo

c---     interpolating AOT tables

      do 105 ia = 1, 10

      call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,
     1          x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

      yy(ia) = y/pi
c      print *,'tau, i/f=', tau(ia), y/pi, dy
 105  continue
      if (refl.le.yy(1)) then
       tau_x470 = 0.02
       tau_x470_flag = 2 ! MJ added
       return
      endif

      if (refl.ge.yy(10)) then
       xxrat = 0.8
       tau_x470 = 3.5
       tau_x470_flag = 1
       return
      endif
c
c
c     Check if the reflectance increase with AOT
c

      if (yy(1).lt.yy(2)) go to 650

      if (refl.lt.yy(4)) return

      do i = 1, 7
      yy2(i) = yy(i+3)
      enddo

      if (yy2(2).lt.yy2(1)) return
      call search2(dflag2,refl,yy2,7,index_ii,frac)
c      print *,'after 1st search'
c      if (yy2(1).gt.yy2(2)) print *,'yy=',refl,(yy(i),I=1,10)
      tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
      tau_x470 = tau_x
      return

650   continue
c
c     Pass the monotonic order check
c
      call search2(dflag2,refl,yy,10,index_ii,frac)
c      print *,'after 2nd search'
      tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii)

      tau_x470 = tau_x

c      print *,'tau_x470 =', tau_x470

      return
      end
c... --------------------------------------
c... The end of subroutine aero_470veg
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c--------------------------------------------------------
      subroutine aero_650veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1    r650,tau_x650,tau_x650_flag,
     2    tau_x470_flag,tau_x470)

      use viirs_aerosol_luts, only: new_intep   

      include 'aottbl.inc'
      include 'newaottbl.inc'

c      common /angle_node/ theta0, theta, phi
c      common /sfcref_node/ sfc_ref412, sfc_ref470, sfc_ref650

c      real theta0(10), theta(46), phi(30), tau(10)
c      real sfc_ref412(20), sfc_ref470(24), sfc_ref650(24)
      real nnvalx(4,4,2,10), yy(10), yy2(8), yy3(3), yy5(6)
      real tau_x650, tau_x412, tau_x470
      real refl,x1,x2,x3,r650
c      real nvalx650(10,46,30,10,24)
      integer tau_x470_flag, tau_x650_flag, tau_x412_flag_91
      logical dflag2
      data pi   /3.14159/
c      data tau  /0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5/
      
      index_ii = (r650+1.)/2.

      frac = (r650-sfc_ref650(index_ii))/
     1       (sfc_ref650(index_ii+1)-sfc_ref650(index_ii))

c      if (index_ii.lt.1.or.index_ii.gt.24)
c     1    print *,'index_iir650 = ', index_ii,xlat,xlong
c      if (frac.lt.0.0.or.frac.gt.1.0)
c     1    print *,'frac on sfc=', frac

c... --------------------------------------------------
c...  MJ added 12Feb2010 @@@@@
c... --------------------------------------------------
c.    if(index_ii.lt.2.or.index_ii.gt.12) then ! ok ver.
      if(index_ii.lt.1.or.index_ii.gt.23) then
         tau_x650_flag = 9 
         tau_x650 = -999.0
         return
      endif ! if(index_ii.lt.1.or.index_ii.gt.23) then
c... --------------------------------------------------

      if (index_ii.lt.1) then
        index_ii = 1
        frac = 0.
      endif

      call search(dflag2,x3,phi,ll,ii)
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))

      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif

      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif

c      write(6,*) "frac =",frac
      do ia = 1, 10
       do i = 1, 4
        do j = 1, 4
          nnvalx(i,j,1,ia) = nvalx650(mbeg+i,nbeg+j,ii,ia,index_ii)*
     1       (1.-frac) + nvalx650(mbeg+i,nbeg+j,ii,ia,index_ii+1)*frac
          nnvalx(i,j,2,ia) = nvalx650(mbeg+i,nbeg+j,ii+1,ia,index_ii)*
     1       (1.-frac) + nvalx650(mbeg+i,nbeg+j,ii+1,ia,index_ii+1)*frac
        enddo
       enddo
      enddo

c---     interpolating AOT tables

      do 600 ia = 1, 10

      call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,
     1          x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

      yy(ia) = y/pi
c      print *,'tau, i/f=', tau(ia), y/pi, dy
 600  continue


      if (refl.le.yy(1).and.yy(1).lt.yy(2)) then
       tau_x650_flag = 2  ! MJ added
       tau_x650 = 0.02
       return
      endif

      if (refl.ge.yy(10)) then
      tau_x650 = 3.5
      w0_x = -999.
      tau_x650_flag = 1
      return
      endif

c!!!  if (tau_x470_flag.gt.0) go to 670 ! ori open 20100216 @@@@@
c      if (tau_x412_flag_91.gt.0) go to 680
c
c
c     Check if the reflectance increase with AOT
c
      if (yy(1).lt.yy(2)) go to 650

c... *****************************
c... *****************************
c... *****************************
        tau_x650_flag = 4  ! mj added 08/05/2011
        return             ! mj added 20100216 @@@@@
                           ! no retrieval for non-monotonic
                           ! changes in LUT reflectance
c... needs investigation to comment this out!
c... *****************************
c... *****************************
c... *****************************

      if (refl.lt.yy(4)) return

      do i = 1, 7
      yy2(i) = yy(i+3)
      enddo

      if (yy2(2).lt.yy2(1)) return
      call search2(dflag2,refl,yy2,7,index_ii,frac)
c      print *,'after 1st search'
c      if (yy2(1).gt.yy2(2)) print *,'yy=',refl,(yy(i),I=1,10)
      tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
      tau_x650 = tau_x
      return

670   continue
      if (refl.lt.yy(8)) return

      do i = 1, 3
      yy3(i) = yy(i+7)
      enddo

      if (yy3(2).lt.yy3(1)) return
      call search2(dflag2,refl,yy3,3,index_ii,frac)
c      print *,'after 1st search'
c      if (yy3(1).gt.yy3(2)) print *,'yy=',refl,(yy(i),I=1,10)
      tau_x = frac*tau(index_ii+1+7) + (1.-frac)*tau(index_ii+7)
      tau_x650 = tau_x
      return

680   continue
      if (refl.lt.yy(5)) return

      do i = 1, 6
      yy5(i) = yy(i+4)
      enddo

      if (yy5(2).lt.yy5(1)) return
      call search2(dflag2,refl,yy5,6,index_ii,frac)
c      print *,'after 1st search'
c      if (yy3(1).gt.yy3(2)) print *,'yy=',refl,(yy(i),I=1,10)
      tau_x = frac*tau(index_ii+1+4) + (1.-frac)*tau(index_ii+4)
      tau_x650 = tau_x
      return

650   continue
c
c     Pass the monotonic order check
c
      call search2(dflag2,refl,yy,10,index_ii,frac)
c      print *,'after 2nd search'
      tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii)

      tau_x650 = tau_x

      return
      end
c... --------------------------------------
c... The end of subroutine aero_650veg
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c **********************************************************************************************************


c **********************************************************************************************************
      subroutine cmb_stddb_veg(outbuf, outbufvg, tmpvg, idalg)
c... 
c...  A subroutine to combine standard Deep Blue and over-vegetation
c...  retrievals, depending on ndvi, ndvi_swir, refl at 2.1um, etc.
c...     Written by Myeong-Jae Jeong (MJ)
c...     Ver. 0.5 Aug 8, 2011
c...
c.... Inputs:
c...     outbuf(20)   - standard DeepBlue output
c...     outbufvg(20) - over-vegetation retr. output
c...     tmpvg(6)     - TOA reflectances
c...
c...  Output:
c...     outfub(20)   - combined output
c...     idalg        - index for algorithm in effect
c...
c...     idalg=0 : standard Deep Blue 
c...     idalg=1 : over-vegetation retrieval
c...     idalg=2 : mixture of standard DB and over-veg. retrievals

      real outbuf(20), outbufvg(20), tmpvg(6)
      real outbuftmp(20), xfrdb
      real c_ndsir1,c_ndvi1,c_rf21_1,c_rf21_2,c_rf21_3
      integer idalg

      c_ndsir1=0.1  ! NDVI_swir minimum threshold  
      c_ndvi1=0.1   ! NDVI minimum threshold 
      c_rf21_1=0.01 ! 2.1um refl. threshold 1
      c_rf21_2=0.25 ! 2.1um refl. threshold 2
      c_rf21_3=0.35 ! 2.1um refl. threshold 3

c...  Note. Also try to employ NDVI_swir thresholds later to combine
c...        standard DeepBlue and over-vegeation retrievals... (MJ)

c...  Initialization ... set to standard DeepBlue retrieval results
      do i=1,20
         outbuftmp(i)=outbuf(i)
      enddo
      idalg=0 ! standard DeepBlue 


c      if(outbufvg(18).ge.c_ndsir1.and.outbufvg(19).ge.c_ndvi1) then
c         if(tmpvg(4).ge.c_rf21_1.and.tmpvg(4).lt.c_rf21_2) then  ! veg.
c            do i=1,8
c               outbuftmp(i)=outbufvg(i)
c            enddo
c            do i=11,14
c               outbuftmp(i)=outbufvg(i)
c            enddo
c            outbuftmp(15) = 100.0
c            idalg=1
c         elseif(tmpvg(4).ge.c_rf21_2.and.tmpvg(4).lt.c_rf21_3) then
cc            xfrdb=(tmpvg(4)-c_rf21_2)/(c_rf21_3-c_rf21_2)     
cc            if(outbufvg(7).ge.0.and.outbuf(7).ge.0) then      ! veg+DB
cc               do i=2,3
cc                  outbuftmp(i)=(1.0-xfrdb)*outbufvg(i)+xfrdb*outbuf(i)
cc                  outbuftmp(i+3)=((1.0-xfrdb)*outbufvg(i+3)*outbufvg(i) 
cc     1                           + xfrdb*outbuf(i+3)*outbuf(i)) 
cc     2                           /(outbuftmp(i))
cc               enddo
cc               outbuftmp(7)=(1.0-xfrdb)*outbufvg(7)+xfrdb*outbuf(7)
cc               outbuftmp(1)=-999.0  ! set aot412nm to -999.0
cc               outbuftmp(4)=-999.0  ! set ssa412nm to -999.0
cc               outbuftmp(8)=alog(outbuftmp(2)/outbuftmp(3))/
cc     1                           alog(650./466.)
cc               if(outbufvg(11).gt.0.and.outbuf(11).gt.0) then
cc                  outbuftmp(11)=((1.0-xfrdb)*outbufvg(11)*outbufvg(2)
cc     1                          +xfrdb*outbuf(11)*outbuf(2))
cc     2                          /(outbuftmp(2))       
cc                  outbuftmp(12)=((1.0-xfrdb)*outbufvg(12)*outbufvg(3)
cc     1                          +xfrdb*outbuf(12)*outbuf(3))
cc     2                          /(outbuftmp(3))       
cc               else 
cc                  outbuftmp(11)=-999.0
cc                  outbuftmp(12)=-999.0
cc               endif
cc               idalg=2
c            if(outbufvg(7).ge.0.and.outbuf(7).lt.0) then
c               do i=1,8
c                  outbuftmp(i)=outbufvg(i)
c               enddo
c               do i=11,14
c                  outbuftmp(i)=outbufvg(i)
c               enddo
c               outbuftmp(15) = 100.0
c               idalg=1
c            else
c               do i=1,8
c                  outbuftmp(i)=outbuf(i)
c               enddo
c               do i=11,14
c                  outbuftmp(i)=outbuf(i)
c               enddo
c               idalg=0
c            endif
c         endif
c      endif 

c...  Finalization ... push the results to outbuf
      do i=1,20
         outbuf(i)=outbuftmp(i)
      enddo

      return
      end
c... --------------------------------------
c... The end of subroutine cmb_stddb_veg
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c **********************************************************************************************************


c **********************************************************************************************************
      subroutine calc_sfc21(sza,xthet,xphi,refl21,aot21,As_21)
c...
c...  ----------------------------------------------------------
c...  written by Myeong-Jae Jeong (MJ) at GWNU, Korea
c...  ver. 1.0 Aug 12, 2011
c...
c...  subroutine to calculate surface reflectance at 2.1um.
c...  code inherited from "sfc470_5wav.f" by N.C. Hsu 
c...
c...  Inputs
c...     sza   : solar zenith angle
c...     xthet : viewing zenith angle
c...     xphi  : relative azimuth angle
c...     refl21: 2.1um TOA reflectance
c...     aot21 : 2.1um AOT
c... 
c...  Output
c...     As_21 : 2.1 surface reflectance (%)
c...
c...  ----------------------------------------------------------
c...     refl21  = refl21  * xmu/3.14159

      include 'aottbl.inc'
      include 'sfc21tbl.inc'

c     real tau2(4), theta0(10), theta(46), phi(30)
      real tau2(4)
      real nnvalx(4,4,2), rr0x(4,4,2)
      real ttx(4,4,2), ssx(4,4,2)
      logical dflag3
      integer index_ii

      data tau2    /0.0, 0.1, 0.3, 0.5/
      data pi      /3.14159/
c     data theta0  /0.0,8.0,16.0,24.0,32.0,40.0,48.0,56.0,64.0,72.0/

      x1 = sza
      x2 = xthet
      x3 = xphi
      refl7 = refl21        ! 2.1 micron

c     do i = 1, 46
c        theta(i) = 2.*float(i-1)
c     enddo

c     do i = 1, 30
c        phi(i) = 6. + 6.*float(i-1)
c     enddo

      mm = 10     ! solar zenith
      nn = 46     ! satellite zenith
      ll = 30     ! rel. azimuth
      ma2=  4     ! tau2 (aot at 2.1um)

c     -- Derive SFC

c--   interpolating for 2.1 micron channel

      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif
 
      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif
c... 
      aot2=aot21
      call search(dflag3,aot2,tau2,ma2,index_ii)
      frac  = (aot2-tau2(index_ii))
     1        /(tau2(index_ii+1)-tau2(index_ii))
      if(frac.lt.0.or.frac.gt.1)
     1   print *,'aot2, frac, index_ii=',aot2, frac, index_ii
 
      call search(dflag3,x3,phi,ll,ii)
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))

c     print *, 'nsm, nsn, ii, aot2, frac, index_ii'
c     print *, nsm, nsn, ii, aot2, frac, index_ii
       do 123 i = 1, 4
        do 124 j = 1, 4
          nnvalx(i,j,1) = nvalx21(mbeg+i,nbeg+j,ii,index_ii)*
     1       (1.-frac) + nvalx21(mbeg+i,nbeg+j,ii,index_ii+1)*frac   
          nnvalx(i,j,2) = nvalx21(mbeg+i,nbeg+j,ii+1,index_ii)*
     1       (1.-frac) + nvalx21(mbeg+i,nbeg+j,ii+1,index_ii+1)*frac   
 
          rr0x(i,j,1) = r0x_21(mbeg+i,nbeg+j,ii,index_ii)*
     1       (1.-frac) + r0x_21(mbeg+i,nbeg+j,ii,index_ii+1)*frac
          rr0x(i,j,2) = r0x_21(mbeg+i,nbeg+j,ii+1,index_ii)*
     1       (1.-frac) + r0x_21(mbeg+i,nbeg+j,ii+1,index_ii+1)*frac
 
          ttx(i,j,1) = tx_21(mbeg+i,nbeg+j,ii,index_ii)*
     1       (1.-frac) + tx_21(mbeg+i,nbeg+j,ii,index_ii+1)*frac
          ttx(i,j,2) = tx_21(mbeg+i,nbeg+j,ii+1,index_ii)*
     1       (1.-frac) + tx_21(mbeg+i,nbeg+j,ii+1,index_ii+1)*frac
 
          ssx(i,j,1) = sx_21(mbeg+i,nbeg+j,ii,index_ii)*
     1       (1.-frac) + sx_21(mbeg+i,nbeg+j,ii,index_ii+1)*frac
          ssx(i,j,2) = sx_21(mbeg+i,nbeg+j,ii+1,index_ii)*
     1       (1.-frac) + sx_21(mbeg+i,nbeg+j,ii+1,index_ii+1)*frac

124    continue
123   continue
c       enddo
c      enddo
 
c---     
 
        call new_intepsf(theta0, theta, phi, nnvalx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
        RR = y/pi
 
        call new_intepsf(theta0, theta, phi, rr0x, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
        RR0 = y/pi
        
        call new_intepsf(theta0, theta, phi, ttx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
        
        TT = y/pi
 
        call new_intepsf(theta0, theta, phi, ssx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
c!      SS = y/pi
        SS = y    ! chg20090818
 
        dff = refl7 - RR0
        As_21  = 100.* dff / (TT + SS*dff)

c        print *,nday,sza,xthet,xphi,refl3,refl7,As,As_21
c        print *,'refl7,RR0,dff,TT,SS,As_21= ',
c     1     refl7,RR0,dff,TT,SS,As_21

      return
      end
c... --------------------------------------
c... The end of subroutine calc_sfc21
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c **********************************************************************************************************

c **********************************************************************************************************
      subroutine calc_rc_ref672(sza,xthet,xphi,ref_toa,ref_rc)
c...
c...  ----------------------------------------------------------
c...  written by Jaehwa Lee
c...  ver. 1.0 Mar. 21, 2014
c...
c...  subroutine to calculate Rayleigh-corrected reflectance
c...  code inherited from "calc_sfc21.f" by MJ
c...
c...  Inputs
c...     sza   : solar zenith angle
c...     xthet : viewing zenith angle
c...     xphi  : relative azimuth angle
c...     ref_toa: TOA reflectance (I/F)
c...
c...  Output
c...     ref_rc : Rayleigh-corrected reflectance (normalized) 
c...
c...  ----------------------------------------------------------

      include 'aottbl.inc'
      include 'sfc21tbl.inc'

c     real tau2(4), theta0(10), theta(46), phi(30)
      real tau2(4)
      real nnvalx(4,4,2), rr0x(4,4,2)
      real ttx(4,4,2), ssx(4,4,2)
      logical dflag3

      data pi      /3.14159/

      x1 = sza
      x2 = xthet
      x3 = xphi

      mm = 10     ! solar zenith
      nn = 46     ! satellite zenith
      ll = 30     ! rel. azimuth
      ma2=  4     ! tau2 (aot at 2.1um)

c--   interpolating for 2.1 micron channel

      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif

      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif
      
      call search(dflag3,x3,phi,ll,ii)
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))

       do 123 i = 1, 4
        do 124 j = 1, 4
          nnvalx(i,j,1) = nvalx672(mbeg+i,nbeg+j,ii,1)
          nnvalx(i,j,2) = nvalx672(mbeg+i,nbeg+j,ii+1,1)

          rr0x(i,j,1) = r0x_672(mbeg+i,nbeg+j,ii,1)
          rr0x(i,j,2) = r0x_672(mbeg+i,nbeg+j,ii+1,1)

          ttx(i,j,1) = tx_672(mbeg+i,nbeg+j,ii,1)
          ttx(i,j,2) = tx_672(mbeg+i,nbeg+j,ii+1,1)

          ssx(i,j,1) = sx_672(mbeg+i,nbeg+j,ii,1)
          ssx(i,j,2) = sx_672(mbeg+i,nbeg+j,ii+1,1)

124    continue
123   continue

        call new_intepsf(theta0, theta, phi, nnvalx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

        RR = y/pi

        call new_intepsf(theta0, theta, phi, rr0x, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

        RR0 = y/pi

        call new_intepsf(theta0, theta, phi, ttx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

        TT = y/pi

        call new_intepsf(theta0, theta, phi, ssx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

        SS = y  

        dff = ref_toa - RR0
        ref_rc  = dff / (TT + SS*dff)

      return
      end
c... --------------------------------------
c... The end of subroutine calc_rc_ref672
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c **********************************************************************************************************

c **********************************************************************************************************
      subroutine calc_rc_ref865(sza,xthet,xphi,ref_toa,ref_rc)
c...
c...  ----------------------------------------------------------
c...  written by Jaehwa Lee
c...  ver. 1.0 Mar. 21, 2014
c...
c...  subroutine to calculate Rayleigh-corrected reflectance
c...  code inherited from "calc_sfc21.f" by MJ
c...
c...  Inputs
c...     sza   : solar zenith angle
c...     xthet : viewing zenith angle
c...     xphi  : relative azimuth angle
c...     ref_toa: TOA reflectance (I/F)
c...
c...  Output
c...     ref_rc : Rayleigh-corrected reflectance (normalized)
c...
c...  ----------------------------------------------------------

      include 'aottbl.inc'
      include 'sfc21tbl.inc'

c     real tau2(4), theta0(10), theta(46), phi(30)
      real tau2(4)
      real nnvalx(4,4,2), rr0x(4,4,2)
      real ttx(4,4,2), ssx(4,4,2)
      logical dflag3

      data pi      /3.14159/

      x1 = sza
      x2 = xthet
      x3 = xphi

      mm = 10     ! solar zenith
      nn = 46     ! satellite zenith
      ll = 30     ! rel. azimuth
      ma2=  4     ! tau2 (aot at 2.1um)

c--   interpolating for 2.1 micron channel

      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif

      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif

      call search(dflag3,x3,phi,ll,ii)
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))

       do 123 i = 1, 4
        do 124 j = 1, 4
          nnvalx(i,j,1) = nvalx865(mbeg+i,nbeg+j,ii,1)
          nnvalx(i,j,2) = nvalx865(mbeg+i,nbeg+j,ii+1,1)

          rr0x(i,j,1) = r0x_865(mbeg+i,nbeg+j,ii,1)
          rr0x(i,j,2) = r0x_865(mbeg+i,nbeg+j,ii+1,1)

          ttx(i,j,1) = tx_865(mbeg+i,nbeg+j,ii,1)
          ttx(i,j,2) = tx_865(mbeg+i,nbeg+j,ii+1,1)

          ssx(i,j,1) = sx_865(mbeg+i,nbeg+j,ii,1)
          ssx(i,j,2) = sx_865(mbeg+i,nbeg+j,ii+1,1)

124    continue
123   continue

        call new_intepsf(theta0, theta, phi, nnvalx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

        RR = y/pi

        call new_intepsf(theta0, theta, phi, rr0x, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

        RR0 = y/pi

        call new_intepsf(theta0, theta, phi, ttx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

        TT = y/pi

        call new_intepsf(theta0, theta, phi, ssx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

        SS = y

        dff = ref_toa - RR0
        ref_rc  = dff / (TT + SS*dff)

      return
      end
c... --------------------------------------
c... The end of subroutine calc_rc_ref865
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c **********************************************************************************************************

      subroutine new_intepsf(x1a,x2a,x3a,ya,m,n,l,x1,x2,x3,y,dy,
     1                     mbeg,nbeg,frac)
      
      use viirs_aerosol_luts, only: polint
      
      parameter (nmax=46,mmax=10,lmax=30)
      dimension x1a(m),x2a(n),x3a(l),ya(4,4,2)
      dimension xx2a(4), xx1a(4)
      dimension yntmp(4),ymtmp(4),yltmp(2)
 
      do 12 j=1,4
        do 11 k=1,4
          yltmp(1)=ya(j,k,1)
          yltmp(2)=ya(j,k,2)
          yntmp(k) = yltmp(1)*(1.-frac) + yltmp(2)*frac
          xx2a(k) = x2a(k+nbeg)
11      continue
        call polint(xx2a,yntmp,4,x2,ymtmp(j),dy)
        xx1a(j) = x1a(j+mbeg)
12    continue
      call polint(xx1a,ymtmp,4,x1,y,dy)
      return
      end

