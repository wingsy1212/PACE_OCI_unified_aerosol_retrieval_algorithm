      subroutine read_data(error_file,xnvalm5)
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C !INPUT PARAMETERS:
C
C    Type       Name             Description
C    ====       ====             ===========
C    REAL*4     xnvalm5          MODIS reflectivities to be converted
C                                 to n-values
C
C !OUTPUT PARAMETERS:
C
C    Type       Name             Description
C    ====       ====             ===========
C    LOGICAL*4  error_file       returned status flag
C
C !REVISION HISTORY:
C
C    Initial Version by Jeremy Warner   12/01/2006
C     Modified from legacy code by Christina Hsu.
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
C     xnvalm                     (sample.inc)
C     xlat                       (sample.inc)
C     xlong                      (sample.inc)
C     sfref412                   (sample.inc)
C     sfref470                   (sample.inc)
C     sfref650                   (sample.inc)
C
C   Functions:  none
C
C !END
C-----------------------------------------------------------------------
c
      implicit none
      logical error_file
      integer*4 i,j,k, ilat5, ilon5, itmp
      real*4    xnvalm5(3), psi, cc, s1, s2, ss
      real*4    rr, xmu, xtmp
      
      real xlat
      real xlong
      real sza
      real xthet
      real xphi
      real cphi
      real c2phi
      real pteran
      real xnvalm
      integer ilat
      integer isnow
      integer partial
      real so2ind
      real resn
      real sens
      real rsens
      real ozbst
      real ref
      real estozn
      real ozcld
      real pcloud
      real prfrac
      real clfrac
      real rayval
      real r412
      real r470
      real sfref412
      real sfref470
      real sfref650
      real qdif412
      real qdif470
      real qdif650
      real stdv
      
      common/sample/ xlat,xlong,sza,xthet,xphi,cphi,c2phi,pteran,
     2 xnvalm(6),ilat,isnow,so2ind,resn(5),sens(5),
     3 rsens(6),ozbst,ref,estozn,ozcld,pcloud,prfrac,clfrac,partial,
     4 rayval(5),r412,r470,sfref412,sfref470,sfref650,qdif412,
     5 qdif470,qdif650, stdv

c      include 'sample.inc'
c      include 'aottbl.inc'
 
c     -- first undo sza normalization in the L1B code
 
      xmu = cos(sza * 3.14159/180.)
 
      do i = 1,3
      rr   = xnvalm5(i) !* xmu/3.14159 ! use this for data from Wei
      if (rr .gt. 0.0) then
      xnvalm5(i) =  -100.*alog10(rr)
      	else
c      	  print *, 'err - xvalm5(i) le 0, i= ',i, rr
      	  go to 300
      endif
      enddo
 
      do i = 1,2
      xnvalm(i) = xnvalm5(1)
      enddo
      do i = 3,4
      xnvalm(i) = xnvalm5(2)
      enddo
 
      xnvalm(5) = xnvalm5(3)
      xnvalm(6) = xnvalm5(1)

c     -- find terrain pressure
c      ilat5 = (-1.* xlat + 90.) *2.
c      ilat5 = ilat5 + 1
c      if (ilat5.gt.360) ilat5 = 360
c      if (ilat5.lt.1)   ilat5 = 1

c      if (xlong.ge.0.0) then
c      ilon5 = xlong *2.
c      ilon5 = ilon5 + 1
c        else
c         ilon5 = (xlong + 360.) *2.
c         ilon5 = ilon5 + 1
c      endif
c      if (ilon5.gt.720) ilon5 = 720
c      if (ilon5.lt.1)   ilon5 = 1

c      pteran = sfcprs(ilat5,ilon5)/1013.

      isnow = 0
      do i = 1,6
         xnvalm(i) = xnvalm(i)/100.
      enddo
      error_file = .false.
      return

300   continue
      error_file = .true.
      return
      end
      
