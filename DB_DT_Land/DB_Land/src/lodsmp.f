      subroutine lodsmp
c
c***********************************************************************
c
c    lodsmp 
c
c    calling routine           total
c
c***********************************************************************
c
      use viirs_ler_luts
      
      include 'contrl.inc'
c      include 'sample.inc'
      common/longit/xlonlo,xlonhi
c      common/lpoly/ xzlog(10),xlog(8),densol(4,7),denscn(4,5),
c     1              cthet0(4),ctheta(4),cofs(16),indsol,indscn,iofset,
c     2              p1,pr,dum10(10)
      double precision convrt
      data convrt /0.005729577951d0/
c     -- unpack angles and lat, & long
c     -- set latitude flag
c      ilat=1
c      if(abs(xlat).gt.15.) ilat=2
c      if(abs(xlat).gt.60.) ilat=3
100   continue
c     -- unpack thir cloud top pressure, terrain height, surface
c     -- category, percent cloudiness and snow thickness
c     -- compute phase factor terms,lagrange coeffs. and index offset
c     -- for solar and satellite zenith angle interpolations
      xtemp1=sza/(convrt*10000.)
      xtemp2=sin(xtemp1)
      xtemp1=cos(xtemp1)
      xzlog1=0.
      if(xtemp1.gt.0.)xzlog1=alog(1.0/xtemp1)
c     -- phi variables
      ztemp1=xphi/(convrt*10000.)
      cphi=cos(ztemp1)
      c2phi=cos(ztemp1*2.0)
c     -- theta variables
      ytemp1=xthet/(convrt*10000.)
      ytemp2=sin(ytemp1)
      ytemp1=cos(ytemp1)
      xlog1=0.
      if(ytemp1.gt.0.)xlog1=alog(1.0/ytemp1)
      pathl=1./xtemp1+1./ytemp1
c     -- phase factors
      p1=-.375*xtemp1*xtemp2*ytemp2
      pr=2.0*p1/(3.0*ytemp1*xtemp1*xtemp1)
c     -- set up index offset value for theta,theta0
      do 200 i=1,10
         indsol=i
c        print *,i,xzlog(i),xzlog1
         if(xzlog(i).ge.xzlog1) go to 250
200   continue
250   continue
      indsol=indsol-2
c     print *,'indsol = ',indsol
c     -- check for range
      if(indsol.lt.1) indsol=1
      if(indsol.gt.7) indsol=7
c     -- theta
      do 300 i=1,8
         indscn=i
c        print *,i,xlog(i),xlog1
         if(xlog(i).ge.xlog1) go to 350
300   continue
350   continue
      indscn=indscn-2
c     print *,'indscn = ',indscn
c     -- check range
      if(indscn.lt.1) indscn=1
      if(indscn.gt.5) indscn=5
c     -- compute lagrange coeffs for solar z.a.
      indmax=indsol+3
      j=1
      do 450 k=indsol,indmax
         xnom=1.0
         do 400 i=indsol,indmax
            if(i.eq.k) go to 400
            xnom=(xzlog1-xzlog(i))*xnom
400      continue
         cthet0(j)=xnom/densol(j,indsol)
         j=j+1
450   continue
      if (lprint(20)) write (6,*) cthet0
      if (lprint(20)) write (6,*) 
      if (lprint(20)) write (6,*) densol
c     -- compute lagrange coeffs. for theta interp
      j=1
      indmax=indscn+3
      do 550 k=indscn,indmax
         xnom=1.0
         do 500 i=indscn,indmax
            if(i.eq.k) go to 500
            xnom=(xlog1-xlog(i))*xnom
500      continue
         ctheta(j)=xnom/denscn(j,indscn)
         j=j+1
550   continue
c     -- compute offset into tables for theta0,theta block
      iofset=(indsol-1)*8+(indscn-1)
c     -- precompute quantities to be used in interpolations
      f1=xzlog1
      f2=xzlog1*xzlog1
      f3=xlog1
      f4=xzlog1*xlog1
      f5=f2*xlog1
c     -- store products of coeffs
      l=1
      do 600 i=1,4
      do 601 k=1,4
         cofs(l)=ctheta(k)*cthet0(i)
         l=l+1
601   continue         
600   continue
c     -- debug for coeffs and indices computed by lodsmp
      if(lprint(9))write(6,1900) cthet0,ctheta,cofs,indsol,
     1 indscn,iofset,p1,pr
610   continue
      return
1900  format(' debug form subroutine lodsmp',/,
     1       ' cthet0= ',4f8.5,1x,'ctheta= ',4f8.5,/,
     2       ' cofs= ',4(/,1x,4f8.5,1x),/,
     3       ' indsol= ',i4,' indscn= ',i4,' iofset= ',i6,/,
     4       ' p1= ',f8.5,1x,' pr= ',f9.5,//)
      end
