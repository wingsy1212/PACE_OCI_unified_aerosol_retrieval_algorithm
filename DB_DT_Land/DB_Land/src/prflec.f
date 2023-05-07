      subroutine prflec(pteran,pcloud,lprint,grref,clref,
     1 pwtlo,pwthi)
c***********************************************************************
c
c     purpose
c        calculate cloud pressure, pwtlo, pwthi, and initialize ground
c        reflectivity and cloud reflectivity for 650 nm reflectivity 
c        channel
c
c     calling sequence
c        call prflec
c
c     variables
c       name      type  i/o  description
c       ----      ----  ---  -----------
c      pteran      i*4   i   terrain pressure
c      pcloud      r*4   i   cloud pressure
c      lprint      l*4   i   if = .true., print
c      grref       r*4   o   ground reflectivity
c      clref       r*4   o   cloud reflectivity
c
c     calling routine       total
c
c
c***********************************************************************
c
c     -- input parameters
      real pteran, pcloud, grref
      logical lprint(20)
c     -- output parameters
      real pwtlo, pwthi, clref
c     -- intialize ground cloud reflectivity
      clref  = 0.80
c     -- set terrain and cloud weighting fractions
      pwtlo = (pteran - .4)/.6
      pwthi = (pcloud - .4)/.6
      if (lprint(15)) then
         write (6,1000) pteran, pcloud
         write (6,1100) pwtlo,pwthi,grref,clref
      endif
      return
1000  format (/'Subroutine prflec'/'Input:   pteran = ',f8.4,
     1 ' pcloud = ',f8.4)
1100  format ('Output:  pwtlo  = ',f8.4,' pwthi  = ',f8.4/
     1 '         grref  = ',f8.4,' clref  = ',f8.4)
      end
