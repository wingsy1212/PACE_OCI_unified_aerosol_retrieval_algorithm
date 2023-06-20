      subroutine lodbuf
c***********************************************************************
c
c     purpose
c        loads buffer /bufout/ with data for current measurement
c
c     variables
c       name        type   i/o  description
c       ----        ----   ---  -----------
c
c      common /sample/
c       xlat         r*4    i   ifov latitude
c       xlong        r*4    i   ifov longitude in degs
c       sza          r*4    i   ifov solar zenith angle (degs)
c       xthet        r*4    i   ifov scan angle (degs)
c       xphi         r*4    i   ifov scan plane angle (degs)
c       pteran       r*4    i   terrain height pressure
c       xnvalm(6)    r*4    i   measured n-values
c       resn         r*4    i   n-value residuals
c       sens         r*4    i   n-value sensitivities
c
c     calling routine  total
c
c******************************************************************
c
      use viirs_ler_luts
      
      integer errflg

c     -- store lat, long
      realbuf(1)  = xlat
      realbuf(2)  = xlong
c     -- store angles
      realbuf(3)  = sza
      realbuf(4)  = xthet
      realbuf(5)  = xphi
c     -- store ozone and reflectivity
c     -- store n values
      realbuf(6)  = ref*100.
      realbuf(7)  = xnvalm(1)*100.
      realbuf(8)  = xnvalm(2)*100.
      realbuf(9)  = xnvalm(3)*100.
      realbuf(10) = xnvalm(4)*100.
      realbuf(11) = xnvalm(5)*100.
      realbuf(12) = xnvalm(6)*100.
c     -- random data
      realbuf(13) = stdv  
      realbuf(14) = -999.
      realbuf(15) = -999.
      realbuf(16) = -999.
      realbuf(17) = -999.
c     -- store rayleigh
      realbuf(18) = qdif412
      realbuf(19) = qdif470
      realbuf(20) = qdif650
      realbuf(21) = pteran
c     -- store reflec
      realbuf(22) = r412*100.
      realbuf(23) = r470*100.
      realbuf(24) = sfref412
      realbuf(25) = sfref470
      realbuf(26) = sfref650
      
      return
      end
