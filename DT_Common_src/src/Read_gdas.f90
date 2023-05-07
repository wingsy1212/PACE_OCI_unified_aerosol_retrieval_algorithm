      SUBROUTINE Get_An_gdas(lat,lon,met_ugrd,met_vgrd,&
           met_pwat,ozone,set_counter_for_anc,RTN_NCEP)
     
         

      implicit none 
      save
      INCLUDE 'read_Sat_MODIS.inc'
 
! --- constants
      real		MISSING
      parameter 	(MISSING = -999.0)

      integer		GRIDSIZE
      parameter		(GRIDSIZE = 721)      

! --- parameters 
      real		lat
      real		lon
      real		met_pres(0:25)
      real		met_temp(0:25)
      real		met_mixr(0:25)
      real		met_land  
      real		met_sfctmp
      real		met_prmsl
      real		met_pwat
      real		met_ugrd
      real		met_vgrd
      real		ozone
      real		icec
      real		sst_bl
      real		sst
      integer		lun_met
      integer		lun_ozone
      integer		lun_icec
      integer		lun_sst
      integer		lun_nisen
      integer		lun_nises
      integer	 	nise,reclen
      
! --- external functions 
      integer 		ezlh_convert
      integer 		read_int
      real 		ppv
      real              read_real
      external 		ezlh_convert
      external          bl_int
      external          blint_met
      external          read_int
      external          read_real

      
! --- internal variables
      character*1       nise_north( GRIDSIZE, GRIDSIZE )
      character*1       nise_south( GRIDSIZE, GRIDSIZE )
      character*160     errmsg
      integer           xsize, ysize
      integer		i
      integer		j
      integer		k
      integer		ret
      integer		iret
      integer		level
      integer		ios
      integer           status,ii,jj,met_flag,ozone_flag
! ... Temporary sst holding variables
      integer  iword,iyrst,imst,idst,iyrnd,imnd,idnd,ndays,index,temp1
      real		x
      real		x0
      real		dx
      real		y
      real		y0
      real		dy
      real		p(0:25)
      real		satmix
      real		xlon
      real		met_grid( 0:359, 0:180, 0:3 )
      real		ozn_grid( 0:359, 0:180 )
      real		ice_grid( 0:719, 0:359 )
      real		sst_grid( 0:359, 0:179 )
      logical 		init
      logical		met_success
      logical		ozn_success
      logical		ice_success
      logical		sst_success
      logical		nisen_success      
      logical		nises_success  
      character * 132 ncepgdas_name,ozone_name
      integer  set_counter_for_anc,RTN_NCEP
      integer iscan,idata
      CHARACTER*255 cwd
      
! ... Temperature and moisture profile pressure levels (hPa)

      data p / 1000.0, 975.0, 950.0, 925.0, 900.0, 850.0, 800.0,&
         750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 450.0, 400.0,&
        350.0, 300.0, 250.0, 200.0, 150.0, 100.0,  70.0,  50.0,&
          30.0,  20.0,  10.0 /
                  
! --- Initialization flag
      data  		init/ .TRUE. /

!-----------------------------------------------------------------------
!    INITIALIZATION
!-----------------------------------------------------------------------

!-----Open and read input data files if this is the first call and
!     set data ingest success/fail flags

   
      
! open and read first time onley

       if( set_counter_for_anc .eq.1) then
         met_flag=0
         ozone_flag=0

!-----------------------------------------------------------------------
!     Get NCEP Meteorological Data
!-----------------------------------------------------------------------
         
          if ( init ) then
        met_success = .false.
        ozn_success = .false.
        ice_success = .false.
       endif
         
        lun_met= 222
      

         
       
        
        
		reclen = 360*181*4*4

        
 
!                 ncepgdas_name =  GDAS_file 
   		       open (unit = lun_met, file = "GDAS_Out", form = 'unformatted',&  
       	       access = 'direct', status = 'old', &  
               recl = reclen, iostat = ios ) 
         
      	        
            
		!  read the unpacked met file 

				read(lun_met, rec = 1, iostat = ios ) met_grid 
               if(ios .eq. 0) met_flag=1 
				close(lun_met)

             
             
             if( ios.ne.0 ) then
               level = 1 
             else
             met_success = .true.
             endif

         endif 

        RTN_NCEP =0 
!-----------------------------------------------------------------------

! --- SET MISSING VALUES
      ozone      = missing
      icec       = missing
      sst        = missing
      sst_bl     = missing
      nise       = int( missing )
      met_ugrd = missing
      met_vgrd = missing
      met_pwat = missing
      
!-----------------------------------------------------------------------
 
! --- MET DATA
        
       
       if( met_flag .ge.1) then   
!       print*,'file name',GDAS_file ,met_flag 
!      if( met_success ) then
        RTN_NCEP=RTN_NCEP+1
! ------ Compute cell coordinates in met and ozn grids
         x = min( max( lon,  -179.99 ), 179.99 )
         if( x .lt. 0.0 ) x = x + 360.0
         x0 = 0.0
         dx = 1.0
         i = int( ( x - x0 + 0.5*dx ) / dx )
         if( i .eq. 360 ) i = 0
 
         y = min( max( lat, -89.99 ), 89.99 )
         y0 = 90.0
         dy = -1.0
         j = int( ( y - y0 + 0.5*dy ) / dy )  

          met_pwat   = met_grid( i, j, 0)
          met_ugrd   = met_grid( i, j, 1 )
          met_vgrd   = met_grid( i, j, 2 )
         ozone    =    met_grid( i, j, 3 )  
       Else
        RTN_NCEP =0
        endif 
      return     
      end