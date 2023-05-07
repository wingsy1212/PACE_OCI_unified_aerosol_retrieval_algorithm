      subroutine intero(index1,index2,lprint,ezero,tr,sbx)
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C    Calculates intensity terms used in computation of
C    reflecitivity for non 650 nm channels.
c    Performs table interpolation for theta and theta0.
c    Uses precomputed lagrangian coefs.
c    Calculates izero (intensity term arising from atmospheric
c    scattering; also calculates transmittance for reflected
c    radiation).
C
C !INPUT PARAMETERS:
C
C    Type       Name             Description
C    ====       ====             ===========
C    INTEGER*4  index1           table pointer for interpolation tables
C    INTEGER*4  index2           pointer for r-coeffs.
C    INTEGER*4  lprint           printout control flags
C
C !OUTPUT PARAMETERS:
C
C    Type       Name             Description
C    ====       ====             ===========
C    REAL*4     t                transmittance for ground reflected
C                                 radiation
C    REAL*4     sb               table values of sb
C
C !REVISION HISTORY:
C
C    Initial Version by Jeremy Warner   12/01/2006
C     Legacy code from C. Seftor.
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the Deep Blue Science Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-02041.
C
C !REFERENCES AND CREDITS
C
C !DESIGN NOTES:
C
C   Externals:
C
c     cphi                       (sample.inc)
c     c2phi                      (sample.inc)
C     li0r                       (table.inc)
C     z1i0r                      (table.inc)
C     z2i0r                      (table.inc)
C     ti0r                       (table.inc)
C
C   Functions:
C
C !END
C-----------------------------------------------------------------------
c     -- input parameters
      logical lprint(20)
      integer index1, index2
c     -- internal parameters
      real(kind=8) inot, ione, itwo, tran
      real(kind=8) zone, ztwo
c     -- output parameters
      real(kind=8) ezero
      real tr, sbx
      common/lpoly/ xzlog(10),xlog(8),densol(4,7),denscn(4,5),
     1              cthet0(4),ctheta(4),cofs(16),indsol,indscn,iofset,
     2              p1,pr,dum10(10)
      include 'table.inc'
      include 'sample.inc'
c
c     -- interp. for i0, i1, i2, t  
      inot = 0.0
      zone = 0.0
      ztwo = 0.0
      tran = 0.0
      m=0
      do 100 i=1,4
         l=index1+8*i-8
         do 50 k=1,4
            m=m+1
            y=cofs(m)
            inot = inot + logi0(l)*y
            zone = zone + z1i0(l)*y
            ztwo = ztwo + z2i0(l)*y
            tran = tran + ti0(l)*y
            l=l+1
50       continue
100   continue
c     -- convert table values into i1, i2, t, and i0
      inot = 10.**inot
      ione = zone * p1 * inot
      itwo = ztwo * pr * p1 * inot
      tran = tran * inot
c     -- compute izero, t, and sbar
      ezero = inot + ione*cphi + itwo*c2phi
      tr    = tran
      sbx   = sb(index2)
c     -- debug print out
      if (lprint(10)) then
         write (6,1000) index1, index2
         write (6,1100) inot,ione,itwo,tran,cphi,c2phi
         write (6,1200) tr,ezero,sbx
      endif
      return
 1000 format('Subroutine intero',/'Input:   index1 = ',i8,' index2 = ',
     1 i8)
 1100 format('Intern:  inot   = ',f8.5,' ione   = ',f8.5,
     1 ' itwo   = ',f8.5/'         tran   = ',f8.5,' cphi   = ',f8.5,
     2 ' c2phi  = ',f8.5)
 1200 format('Output:  t      = ',f8.5,' ezero  = ',f8.5,
     1 ' sbx    = ',f8.5)
      end
