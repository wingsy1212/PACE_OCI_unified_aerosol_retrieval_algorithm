      subroutine reflec(iofset,xlat,isnow,grref,clref,pwtlo,pwthi,
     1 partial,pcloud,pteran,xnvalm,lprint,clfrac,ref,ref10,ref04,
     2 algflg,r412,r470,sfref412,sfref470,sfref650,qdif412,qdif470, 
     3 qdif650)
c
c***********************************************************************
c
c     purpose
c       compute cloud fraction and reflectivity for 650 nm wavelength 
c       and for both table pressures
c
c     method
c       select interp. indices for table lookup.
c       perform table interpolations, then calculate reflectivity.
c
c     variables
c       name      type  i/o   description
c       ----      ----  ---   -----------
c      arguments
c       iofset    i*4    i    offset pointer.
c       xlat      i*4    i    latitude
c       isnow     i*4    i    snow/ice indicator
c       grref     r*4    i    ground reflectivity
c       clref     r*4    i    cloud reflectivity
c       pwtlo     r*4    i    weight calculated using terrain pressure
c       pwthi     r*4    i    weight calculated using cloud pressure
c       partial   i*4    i    partial controls how snow/ice is handled
c       pcloud    r*4    i    cloud top pressure
c       pteran    r*4    i    terrain pressure
c       lprint    l*4    i    if = .true., print
c       clfrac    r*4    o    cloud fraction
c       ref       r*4    o    reflectivity from 650 nm channel
c       ref10     r*4    o    low  reflectivity (needed for clfrac < 0, > 1)
c       ref04     r*4    o    high reflectivity (needed for clfrac < 0, > 1)
c       algflg    i*4    o    algorithm flag (set to 10 for snow/ice, 0 otherwise)
c
c      internal
c       indx1     i*4         table index
c       indx2     i*4         coefficient index
c       iptab1(2) i*4         pressure offset for table index
c       iptab2(2) i*4         pressure offset for coefficient index
c
c***********************************************************************
c
      include 'refl.inc'
c     -- input parameters
      real xlat, pwtlo, pwthi, pcloud, pteran
      real grref, clref, xnvalm(6)
      real sfref412,sfref470,sfref650
      integer iofset, isnow, partial
      logical lprint(20)
c     -- internal parameters
      real(kind=8) qgcalc, qccalc, qsavlo, qsavhi
      real(kind=8) qgcalcx, qccalcx,cl412,cl470
      real reflo, refhi
c
      real(kind=8) ezero(2)
      real t(2), sbx(2)
      real(kind=8) ez412(2)
      real t412(2), sb412(2)
      real(kind=8) ez470(2)
      real t470(2), sb470(2)
      real(kind=8) qgc412x, qcc412x, qgc470x, qcc470x, alb, alb412, alb470
      real(kind=8) den, d412,d470
c
      integer iptab1(2),iptab2(2)
      integer ipt1(2),ipt2(2)
      integer indx1, indx2, ip
c
c     -- output parameters
      integer algflg
      real clfrac, ref
      real lamtb1(5), lamtb2(5)
      real qdif412, qdif470, qdif650
c
c     -- data statements
      data iptab1/1,81/, iptab2/1,2/
c     -- these data statements are for the 412 and 470 reflectivities
      data ipt1/1,10401/, ipt2/1,131/
      data lamtb1/0,2080,4160,6240,8320/
      data lamtb2/0,26,52,78,104/
c
c     -- convert 650, 412, and 470 nm channel n values to albedos
      alb    = 10.**(-xnvalm(6))
      alb412 = 10.**(-xnvalm(5))
      alb470 = 10.**(-xnvalm(3))
      
c      print *,alb412,alb470, alb
c      if (lprint(11)) then
c         write (6,1000) xlat, isnow, iofset, grref, clref, 
c     1    pwtlo, pwthi, pcloud, pteran, alb
c      endif
      if (lprint(1)) then
      endif
c
      qsavhi = -999.0
      qsavlo = -999.0
      qsh412 = -999.0
      qsh470 = -999.0
      qsl412 = -999.0
      qsl470 = -999.0
         
c     -- calculate Rayleigh correction
      do 350 ip = 1,2
c        -- compute table indices for 650 channel
         indx1  = ipt1(ip) + lamtb1(4) + iofset
         indx2  = ipt2(ip) + lamtb2(4)
c        -- compute table indices for 412 channel
         id1  = ipt1(ip) + lamtb1(1) + iofset
         id2  = ipt2(ip) + lamtb2(1)
c        -- compute table indices for 470 channel
         ind1 = ipt1(ip) + lamtb1(3) + iofset 
         ind2 = ipt2(ip) + lamtb2(3)
c        -- perform interpolations
c        -- for 650
         call intero(indx1,indx2,lprint,ezero(ip),t(ip),sbx(ip))
c        -- for 412
         call intero(id1,id2,lprint,ez412(ip),t412(ip),sb412(ip))
c        -- for 470
         call intero(ind1,ind2,lprint,ez470(ip),t470(ip),sb470(ip))
c
c        -- determine calculated q values
c
c        -- 650 --
c
         qgcalc = ezero(ip) + grref*t(ip)/(1-grref*sbx(ip))
         qccalc = ezero(ip) + clref*t(ip)/(1-clref*sbx(ip))
c
c        -- 412 --
c
         qgc412 = ez412(ip) + sfref412*t412(ip)/(1-sfref412*sb412(ip))
         qcc412 = ez412(ip) + clref*t412(ip)/(1-clref*sb412(ip))
c
c        -- 470 --
c
         qgc470 = ez470(ip) + sfref470*t470(ip)/(1-sfref470*sb470(ip))
         qcc470 = ez470(ip) + clref*t470(ip)/(1-clref*sb470(ip))
c
c         if (lprint(11)) write (6,1100) ip, iptab1(ip),
c     1    indx1, iptab2(ip), indx2, sbx(ip), qgcalc, qccalc

         if (ip.eq.1) then
c           -- for ip = 1 (p = 1 atm) save q values
c           -- and denominator used in ref calculation
            qsavlo = qgcalc
            qsavhi = qccalc
            qsl412 = qgc412
            qsh412 = qcc412
            qsl470 = qgc470
            qsh470 = qcc470
         else
c           -- for ip = 2 (p = 0.4 atm) interpolate
c           -- 0.4 and 1 atm q values to obtain
c           -- final q values
            qgcalc = pwtlo*qsavlo + qgcalc*(1.-pwtlo)
            qccalc = pwthi*qsavhi + qccalc*(1.-pwthi)
            qdif650= qsavlo - qgcalc
c           -- 412 --
            qgc412 = pwtlo*qsl412 + qgc412*(1.-pwtlo)
            qcc412 = pwthi*qsh412 + qcc412*(1.-pwthi)
            qdif412= qsl412 - qgc412
c           -- 470 --
            qgc470 = pwtlo*qsl470 + qgc470*(1.-pwtlo)
            qcc470 = pwthi*qsh470 + qcc470*(1.-pwthi)
            qdif470= qsl470 - qgc470
c            if (lprint(11)) write (6,1200) qgcalc,qccalc,
c     1       reflo, refhi, den
         endif
350   continue

c     -- Calculate reflectivity
      rh412 = -999.0
      rh470 = -999.0
      rl412 = -999.0
      rl470 = -999.0
      cl412 = -999.0
      cl470 = -999.0
      
      grref_base = 0.02
      do 351 ip = 1,2
c        -- compute table indices for 650 channel
         indx1  = ipt1(ip) + lamtb1(4) + iofset
         indx2  = ipt2(ip) + lamtb2(4)
c        -- compute table indices for 412 channel
         id1  = ipt1(ip) + lamtb1(1) + iofset 
         id2  = ipt2(ip) + lamtb2(1)
c        -- compute table indices for 470 channel
         ind1 = ipt1(ip) + lamtb1(3) + iofset 
         ind2 = ipt2(ip) + lamtb2(3)
c        -- perform interpolations
c        -- for 650
         call intero(indx1,indx2,lprint,ezero(ip),t(ip),sbx(ip))
c        -- for 412
         call intero(id1,id2,lprint,ez412(ip),t412(ip),sb412(ip))
c        -- for 470
         call intero(ind1,ind2,lprint,ez470(ip),t470(ip),sb470(ip))
c
c        -- determine calculated q values
c
c        -- 650 --
c
         qgcalcx = ezero(ip) + grref_base*t(ip)/(1-grref_base*sbx(ip))
         qccalcx = ezero(ip) + clref*t(ip)/(1-clref*sbx(ip))
c
c        -- 412 --
c
         qgc412x = ez412(ip) + grref_base*t412(ip)
     1                       /(1-grref_base*sb412(ip))
         qcc412x = ez412(ip) + clref*t412(ip)/(1-clref*sb412(ip))
c
c        -- 470 --
c
         qgc470x = ez470(ip) + grref_base*t470(ip)
     1                      /(1-grref_base*sb470(ip))
         qcc470x = ez470(ip) + clref*t470(ip)/(1-clref*sb470(ip))
c
         if (lprint(11)) write (6,1100) ip, iptab1(ip),
     1    indx1, iptab2(ip), indx2, sbx(ip), qgcalcx, qccalcx

         if (ip.eq.1) then
c           -- for ip = 1 (p = 1 atm) save q values
c           -- and denominator used in ref calculation
            qsavlo = qgcalcx
            qsavhi = qccalcx
            qsl412 = qgc412x
            qsh412 = qcc412x
            qsl470 = qgc470x
            qsh470 = qcc470x
            den  = alb     - ezero(ip)
            d412 = alb412  - ez412(ip)
            d470 = alb470  - ez470(ip)
c           -- calculate reflectivity using version
c           -- 6 method and save
            if (den.ne.0. .and. d412.ne.0. .and .d470.ne.0.) then
               ref10 = 1. / (t(ip)/den + sbx(ip))
               r412_10 = 1. / (t412(ip)/d412 + sb412(ip))
               r470_10 = 1. / (t470(ip)/d470 + sb470(ip))
c               print *,'ref10 =', ref10
            else
               ref10 = 0.
               r412_10 = 0.
               r470_10 = 0.
            endif
         else
c           -- for ip = 2 (p = 0.4 atm) interpolate
c           -- 0.4 and 1 atm q values to obtain
c           -- final q values
            pwtlo_base = 1.0
            qgcalcx = pwtlo_base*qsavlo + qgcalcx*(1.-pwtlo_base)
            qccalcx = pwthi*qsavhi + qccalcx*(1.-pwthi)
c           -- 412 --
            qgc412x = pwtlo_base*qsl412 + qgc412x*(1.-pwtlo_base)
            qcc412x = pwthi*qsh412 + qcc412x*(1.-pwthi)
c           -- 470 --
            qgc470x = pwtlo_base*qsl470 + qgc470x*(1.-pwtlo_base)
            qcc470x = pwthi*qsh470 + qcc470x*(1.-pwthi)
c           -- determine denominator used in Version 6
c           -- ref calculation and calculate Version 6
c           -- reflectivity
            den  = alb     - ezero(ip)
            d412 = alb412  - ez412(ip)
            d470 = alb470  - ez470(ip)
            if (den.ne.0. .and. d412.ne.0. .and. d470.ne.0.) then
               ref04 = 1. / (t(ip)/den + sbx(ip))
               reflo = pwtlo_base*ref10 + (1.-pwtlo_base)*ref04
               refhi = pwthi*ref10 + (1.-pwthi)*ref04
c               print *,'pwtlo,ref10,ref04 =', pwtlo_base,ref10,ref04
c              -- 412 --
               r412_04 = 1. / (t412(ip)/d412 + sb412(ip))
               rl412 = pwtlo_base*r412_10 + (1.-pwtlo_base)*r412_04
               rh412 = pwthi*r412_10 + (1.-pwthi)*r412_04
c              -- 470 --
               r470_04 = 1. / (t470(ip)/d470 + sb470(ip))
               rl470 = pwtlo_base*r470_10 + (1.-pwtlo_base)*r470_04
               rh470 = pwthi*r470_10 + (1.-pwthi)*r470_04
            else
               reflo = ref10
               refhi = ref10
c              -- 412 --
               rl412 = r412_10
               rh412 = r412_10
c              -- 470 --
               rl470 = r470_10
               rh470 = r470_10
            endif
            if (lprint(11)) write (6,1200) qgcalcx,qccalcx,
     1       reflo, refhi, den
         endif
351   continue
c     -- Calculate cloud fraction
      if (qccalcx-qgcalcx.ne.0.) then
         clfrac = (alb-qgcalcx) / (qccalcx-qgcalcx)
         cl412  = (alb412-qgc412x) / (qcc412x-qgc412x)
         cl470  = (alb470-qgc470x) / (qcc470x-qgc470x)
      else
         clfrac = 0.
      endif
      if (partial.eq.0) then
        clfrac = -1.
      else
c       -- if snow/ice probability > 50% recalculate clfrac, grref 
c       -- and, if necessary, clref.
        if (isnow.ge.5) then
          algflg = 10
          if (clfrac.gt.0.0 .and. clfrac.lt.1.0) then
            clfrac = clfrac/2.
            alb = (alb-clfrac*qccalcx)/(1.-clfrac)
            do ip = 1,2
              den = alb - ezero(ip)
              if (ip.eq.1) then
                if (den.ne.0.) then
                  ref10 = 1. / (t(ip)/den + sbx(ip))
                else
                  ref10 = 0.
                endif
              else
                if (den.ne.0.) then
                  ref04 = 1. / (t(ip)/den + sbx(ip))
                  grref = pwtlo*ref10 + (1.-pwtlo)*ref04
                else
                  grref = ref10
                endif
              endif
            enddo
          else if (clfrac.ge.1.0) then
            clfrac = 0.5
            grref  = ref10*pwtlo + ref04*(1.-pwtlo)
            clref  = ref10*pwthi + ref04*(1.-pwthi)
          else
            clfrac = 0.
            grref = ref10*pwtlo + ref04*(1.-pwtlo)
          endif
        else
          algflg = 0
        endif
      endif
c      print *,'partial, clfrac =', partial, clfrac
c     -- calculate reflectivity
c     -- use partial cloud algorithm if 0 < clfrac < 1
c     -- otherwise use Version 6 calculation
      if (clfrac.lt.0.) then
         ref  = reflo
         r412 = rl412
         r470 = rl470
      else if (clfrac.gt.1.) then
         ref  = refhi
         r412 = rh412
         r470 = rh470
      else
         ref = grref_base + clfrac*(clref-grref_base)
         r412 = grref_base + cl412*(clref-grref_base)
         r470 = grref_base + cl470*(clref-grref_base)
      endif
      if (lprint(11)) write (6,1300) clfrac, ref
c      write (6,*) clfrac,cl412,cl470,grref_base, ref,r412,
c     1               r470, qdif412,qdif470, qdif650
c      print *,alb412,alb470,alb 
c      print *,r412,ref,r470,grref
c      print *,'cl412,clfrac,cl470=',cl412,clfrac,cl470

      return
1000  format (/'Subroutine reflec'/'Input:   ','lat    = ',f8.3,
     1 ' snow   = ',i8,' iofset = ',i8/
     2 '         grref  = ',f8.4,' clref  = ',f8.4,
     3 ' pwtlo  = ',f8.4,' pwthi  = ',f8.4/
     3 '         pcloud = ',f8.4,' pteran = ',f8.4,' rsf    = ',f8.4)
1100  format ('Intern:  ip     = ',i8,' iptab1 = ',i8,
     1 ' indx1  = ',i8/'         iptab2 = ',i8,
     2 ' indx2  = ',i8,' sbx    = ',f8.5/
     3 '         qgcalc = ',f8.6,' qccalc = ',f8.6)
1200  format ('Final:   qgcalc = ',f8.6,' qccalc = ',f8.6,
     1 ' reflo  = ',f8.4,' refhi  = ',f8.4/
     2 '         den    = ',f8.4)
1300  format ('Output:  clfrac = ',f8.4,' ref    = ',f8.5,
     1 ' qdif412    = ',f9.6,' qdif470    = ',f9.6,
     2 ' qdif650    = ',f9.6)
      end
