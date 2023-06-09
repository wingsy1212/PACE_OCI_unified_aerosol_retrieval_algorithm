      subroutine total
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C    Controls reflectivity processing.  First gets surface reflectivity,
C    and then computes cloud fraction and reflectivity using 650 nm 
C    channel.
C
C !INPUT PARAMETERS:  none
C
C !OUTPUT PARAMETERS:  none
C
C !REVISION HISTORY:
C
C    Initial Version by Jeremy Warner   12/01/2006
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
C     sfref412       (sample.inc)
C     sfref470       (sample.inc)
C     sfref650       (sample.inc)
C     xlat           (sample.inc)
C     isnow          (sample.inc)
C     pteran         (sample.inc)
C     prfrac         (sample.inc)
C     pcloud         (sample.inc)
C     partial        (sample.inc)
C     xnvalm         (sample.inc)
C     clfrac         (sample.inc)
C     r412           (sample.inc)
C     r470           (sample.inc)
C     qdif412        (sample.inc)
C     qdif470        (sample.inc)
C     qdif650        (sample.inc)
C     lprint         (contrl.inc)
C
C   Functions:
C
C     lodsmp
C     prflec
C     reflec
C     lodbuf
C
C !END
C-----------------------------------------------------------------------
      use viirs_ler_luts

      integer algflg
      include 'contrl.inc'
c      include 'sample.inc'
c      common/lpoly/ xzlog(10),xlog(8),densol(4,7),denscn(4,5),
c     1              cthet0(4),ctheta(4),cofs(16),indsol,indscn,iofset,
c     2              p1,pr,dum10(10)
c     -- load in measurement     
      call lodsmp
c     -- set ground reflecitivity
      grref = sfref650

      call prflec(pteran,pcloud,lprint,grref,clref,pwtlo,pwthi)

c     -- calculate the terrain factor
      fteran = (1.0-pteran)/0.60

c     -- compute cloud fraction and reflectivity using 650 nm channel
      partial = 1
      grref = sfref650

      call reflec(iofset,xlat,isnow,grref,clref,pwtlo,pwthi,partial,
     1 pcloud,pteran,xnvalm,lprint,clfrac,ref,ref10,ref04,algflg,
     2 r412,r470,sfref412,sfref470,sfref650,qdif412,qdif470,qdif650)

c     -- load data into realbuf
      call lodbuf
      return
      end
