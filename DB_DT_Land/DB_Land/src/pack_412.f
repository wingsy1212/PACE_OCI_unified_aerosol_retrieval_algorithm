c-------------------------------------------------------------
      subroutine rsfc412(dflag2,mod_sfc,xday,xthet,xphi,scat_ang,r412)
c
      integer mod_sfc, nc_wav(5), index_iw
      logical dflag2
      real xday,r412p,xsfc412p,frac_iw
      real xthet,xphi,scat_ang,r412
      real wav1(16),sfc1(16),wav2(16),sfc2(16)
      real wav3(16),sfc3(16),wav4(15),sfc4(15)
      real wav5(15),sfc5(15)
      real wav1_ond(15),sfc1_ond(15)
      real wav1_jf(15), sfc1_jf(15)
      real wav1_mam(15),sfc1_mam(15)
      character *60 xname_wav(5)
      data xname_wav /'Harmim','Mezaira','Al_Khaznah',
     1  'Sai_Salam','SMART'/
      data nc_wav /16,16,16,15,15/
      data nc_wav_ond /15/
      data nc_wav_jf  /15/
      data nc_wav_mam /15/
      data wav1 /-62.2694,-57.1807,-51.0789,-43.4577,-34.2380,
     1 -23.0985, -10.1670,-9.66400, 3.88359, 17.4518, 29.4005,
     1  39.5023,  47.6345, 54.2590, 59.7262, 64.2073/
      data sfc1 /8.99260,9.73873,9.93104,9.84897,9.45576,
     1   9.23501,9.13252,9.12330,8.55226,8.60710,8.95064,
     1   9.00080,9.14428,9.15092,9.09600,9.07312/
      data wav2 /-65.2144,-60.7404,-55.3780,-48.8869,-40.6035,
     1 -30.8647, -19.0751,-5.59166, 8.29930, 21.4201, 32.7246,
     1  42.2137,  49.7960, 56.0005, 61.1576, 65.3304/
      data sfc2 /10.5924,10.4384,10.6512,10.7640,10.2053,
     1   10.0751,10.0022,9.87149,9.03432,9.31455,9.63064,
     1   9.61609,10.0515,10.2903,9.67356,9.55387/
      data wav3 /-65.0665,-60.5814,-55.2427,-48.7424,-40.7624,
     1 -30.8889, -19.1950,-5.84152, 7.86039, 21.0191, 32.2845,
     1  41.6669,  49.4183, 55.6076, 60.8122, 65.1024/
      data sfc3 /11.3743,11.5292,12.1709,12.4900,12.1680,
     1   11.8708,11.5987,10.6814,10.5286,10.6991,10.4998,
     1   10.4047,10.5189,10.1996,10.1292,10.0239/
      data wav4 /-61.5411,-56.4792,-50.1862,-42.6219,-33.1955,
     1 -21.9909, -9.07847, 4.54637, 17.8186, 29.6754, 39.4847,
     1  47.5092,  54.1248, 59.5025, 63.9649/
      data sfc4 /13.8277,14.2728,13.9865,14.0596,13.5513,
     1   13.1604,12.3758,12.3556,12.1580,12.2147,12.3659,
     1   12.3154,12.1898,12.1819,12.4827/
      data wav5 /-62.1827,-57.1542,-51.0592,-43.6501,-34.5908,
     1 -23.3993, -10.6227, 3.20328, 16.5169, 28.6361, 38.6162,
     1  46.9379,  53.6876, 59.1098, 63.7301/
      data sfc5 /8.55188,9.39847,9.08119,8.26033,7.64707,
     1   7.60514,7.36482,7.11801,6.97991,7.16096,7.57389,
     1   8.00454,8.02307,7.92663,7.96206/

      data wav1_ond /-62.320,-57.277,-51.116,-43.613,-34.427,
     1       -23.160,-10.255, 3.8275, 17.257, 29.243, 39.410,
     1        47.510, 54.190, 59.630, 64.150/
      data sfc1_ond /10.9800,10.0325, 10.230, 9.6033, 9.160,
     1        8.0350, 8.7350, 8.3525, 8.3100, 8.8467, 8.900,
     1        9.1900, 9.7500, 9.3000, 9.2700/
      data wav1_jf  /-62.367,-57.350,-51.200,-43.760,-34.425,
     1       -23.160,-10.283, 3.4000, 17.147, 29.137, 39.263,
     1        47.500, 54.140, 59.590, 64.160/
      data sfc1_jf  / 9.3475, 8.9250, 9.3300, 8.5600, 8.8600,
     1        9.0400, 8.2333, 8.1650, 8.0100, 8.3450, 8.4967,
     1        8.6400, 8.9500, 8.8000, 7.7700/
      data wav1_mam /-62.410,-57.403,-51.347,-43.770,-34.595,
     1       -23.380,-10.533,  3.510, 16.877, 28.947, 39.235,
     1        47.353, 54.147, 59.533, 64.097/
      data sfc1_mam / 9.5100, 9.9200, 9.3300, 9.5600, 8.9150,
     1        9.5375, 8.6567, 8.5633, 8.5100, 8.7433, 8.6950,
     1        8.8400, 8.9933, 8.9500, 8.6700/
       
      dflag2 = .false.

       xxfac = 8.1650 / 8.55226
       if (xday.gt.59.0.and.xday.le.181.0)
     1     xxfac = 8.5633 / 8.55226
       if (xday.gt.181.0.and.xday.le.273.0)
     1     xxfac = 1.0
       if (xday.gt.273.0.and.xday.le.360.0)
     1     xxfac = 8.3525 / 8.55226

      if (mod_sfc.eq.9) then
c-- JF
       dd = xthet
       if (xphi.gt.90.0) dd = -1. *xthet
       if (dd.le.wav1_jf(1)) then
          r412 = sfc1_jf(1)
          return
       endif
       if (dd.ge.wav1_jf(nc_wav_jf)) then
          r412 = sfc1_jf(nc_wav_jf)
          return
       endif
       call search2(dflag2,dd,wav1_jf,nc_wav_jf,index_iw,frac_iw)
      if (dflag2) return
       r412 =frac_iw*sfc1_jf(index_iw+1)+(1.-frac_iw)*sfc1_jf(index_iw)

c-- MAM
       if (xday.gt.59.0.and.xday.le.181.0) then
       dd = xthet
       if (xphi.gt.90.0) dd = -1. *xthet
       if (dd.le.wav1_mam(1)) then
          r412 = sfc1_mam(1)
          return
       endif
       if (dd.ge.wav1_mam(nc_wav_mam)) then
          r412 = sfc1_mam(nc_wav_mam)
          return
       endif
       call search2(dflag2,dd,wav1_mam,nc_wav_mam,index_iw,frac_iw)
      if (dflag2) return
       r412 =frac_iw*sfc1_mam(index_iw+1)
     1      +(1.-frac_iw)*sfc1_mam(index_iw)
       endif

c-- JAS
       if (xday.gt.181.0.and.xday.le.273.0) then
       dd = xthet
       if (xphi.gt.90.0) dd = -1. *xthet
       if (dd.le.wav1(1)) then
          r412 = sfc1(1)
          return
       endif
       if (dd.ge.wav1(nc_wav(1))) then
          r412 = sfc1(nc_wav(1))
          return
       endif
       call search2(dflag2,dd,wav1,nc_wav(1),index_iw,frac_iw)
      if (dflag2) return
       r412 = frac_iw*sfc1(index_iw+1) + (1.-frac_iw)*sfc1(index_iw)
       endif

c-- OND
       if (xday.gt.273.0.and.xday.le.360.0) then
       dd = xthet
       if (xphi.gt.90.0) dd = -1. *xthet
       if (dd.le.wav1_ond(1)) then
          r412 = sfc1_ond(1)
          return
       endif
       if (dd.ge.wav1_ond(nc_wav_ond)) then
          r412 = sfc1_ond(nc_wav_ond)
          return
       endif
       call search2(dflag2,dd,wav1_ond,nc_wav_ond,index_iw,frac_iw)
      if (dflag2) return
       r412 =frac_iw*sfc1_ond(index_iw+1)
     1      +(1.-frac_iw)*sfc1_ond(index_iw)
       endif

      endif
 
      if (mod_sfc.eq.10) then
       dd = xthet
       if (xphi.gt.90.0) dd = -1. *xthet
       if (dd.le.wav2(1)) then
          r412 = sfc2(1) *xxfac
          return
       endif
       if (dd.ge.wav2(nc_wav(2))) then
          r412 = sfc2(nc_wav(2)) *xxfac
          return
       endif
       call search2(dflag2,dd,wav2,nc_wav(2),index_iw,frac_iw)
      if (dflag2) return
       r412 =(frac_iw*sfc2(index_iw+1) +(1.-frac_iw)*sfc2(index_iw))
     1        *xxfac
      endif
 
      if (mod_sfc.eq.11) then
       dd = xthet
       if (xphi.gt.90.0) dd = -1. *xthet
       if (dd.le.wav3(1)) then
          r412 = sfc3(1) *xxfac
          return
       endif
       if (dd.ge.wav3(nc_wav(3))) then
          r412 = sfc3(nc_wav(3)) *xxfac
          return
       endif
       call search2(dflag2,dd,wav3,nc_wav(3),index_iw,frac_iw)
      if (dflag2) return
       r412 =(frac_iw*sfc3(index_iw+1) +(1.-frac_iw)*sfc3(index_iw))
     1        *xxfac
      endif
 
      if (mod_sfc.eq.12) then
       dd = xthet
       if (xphi.gt.90.0) dd = -1. *xthet
       if (dd.le.wav4(1)) then
          r412 = sfc4(1) *xxfac
          return
       endif
       if (dd.ge.wav4(nc_wav(4))) then
          r412 = sfc4(nc_wav(4)) *xxfac
          return
       endif
       call search2(dflag2,dd,wav4,nc_wav(4),index_iw,frac_iw)
      if (dflag2) return
       r412 =(frac_iw*sfc4(index_iw+1) +(1.-frac_iw)*sfc4(index_iw))
     1        *xxfac
      endif

      if (mod_sfc.eq.15) then
       dd = xthet
       if (xphi.gt.90.0) dd = -1. *xthet
       if (dd.le.wav5(1)) then
          r412 = sfc5(1) *xxfac
          return
       endif
       if (dd.ge.wav5(nc_wav(5))) then
          r412 = sfc5(nc_wav(5)) *xxfac
          return
       endif
       call search2(dflag2,dd,wav5,nc_wav(5),index_iw,frac_iw)
      if (dflag2) return
       r412 =(frac_iw*sfc5(index_iw+1) +(1.-frac_iw)*sfc5(index_iw))
     1        *xxfac
      endif

      return
      end

c-------------------------------------------------------------
      subroutine newsfc412_arab(dflag2,mod_sfc,xday,xlatp,xlonp,xthet,
     1         xphi,scat_ang,terrain_flag_new5,r412_135, r412new) 

c
      include 'aottbl.inc'
      include 'newaottbl.inc'

      integer mod_sfc, jtime
      logical dflag2
      real xday,r412p,xsfc412p
      integer xlatp, xlonp 
      real terrain_flag_new5
      real xthet,xphi,scat_ang,r412new1,r412new2
      real dd1, xnorm_fac1, xnorm_fac2
      real r412new, xx, xnorm_fac
      real xfac9(4), xfac10(4), xfac11(4), xfac12(4)
      real xfac13(4), xfac14(4), xfac15(4), xfac16(4)
      real xcc(8,4), xfacp(8)
      character* 3 name(4)

      data name /'win', 'spr', 'sum', 'fal'/
      data xfac9  / 7.8650,  8.1939,  9.0053, 9.0/    !for viirs
      data xfac10 / 10.4847, 10.6111, 9.78469, 0.0/
      data xfac11 / 11.0870, 10.6779, 10.3870, 0.0/
      data xfac12 / 12.8442, 11.9320, 12.1442, 0.0/
      data xfac13 / 12.8442, 11.9320, 12.1442, 0.0/
      data xfac14 / 12.8442, 11.9320, 12.1442, 0.0/
      data xfac15 / 8.52884, 8.17331, 7.82884, 0.0/
      data xfac16 / 11.2968, 12.1234, 10.5265, 0.0/

c -- Interpolate seasonal tables
      do i = 1,8
        do j = 1,3
          xcc(i,j) = -999.0
        enddo
      enddo

      do i = 1, 8
       if (i.eq.1) then
         do j = 1, 4
         xcc(i,j) = xfac9(j)
         enddo
       endif
       if (i.eq.2) then
         do j = 1, 4
         xcc(i,j) = xfac10(j)
         enddo
       endif
       if (i.eq.3) then
         do j = 1, 4
         xcc(i,j) = xfac11(j)
         enddo
       endif
       if (i.eq.4) then
         do j = 1, 4
         xcc(i,j) = xfac12(j)
         enddo
       endif
       if (i.eq.5) then
         do j = 1, 4
         xcc(i,j) = xfac13(j)
         enddo
       endif
       if (i.eq.6) then
         do j = 1, 4
         xcc(i,j) = xfac14(j)
         enddo
       endif
       if (i.eq.7) then
         do j = 1, 4
         xcc(i,j) = xfac15(j)
         enddo
       endif
       if (i.eq.8) then
         do j = 1, 4
         xcc(i,j) = xfac16(j)
         enddo
       endif
      enddo

      jtime = 1
      if (xday.ge.60.0.and.xday.lt.152.0) jtime = 2
      if (xday.ge.152.0.and.xday.lt.244.0) jtime = 3
      if (xday.ge.244.0.and.xday.lt.335.0) jtime = 4

      do i = 1, 8
      xfacp(i) = xcc(i,jtime)
      enddo

c -- interpolate bi-directional factors

      dflag2 = .false.

       mod_sfc = 9
       call rsfc412(dflag2,mod_sfc,xday,xthet,xphi,scat_ang,r412new2)
       xnorm_fac1 = r412new2 / xfacp(1)
       r412new = r412_135 * xnorm_fac1
      return
      end
