! write_vaerdt.f90

! Draft netCDF writer for VIIRS Dark Target
! Simplified to fewer variables for testing

! Virginia Sawyer
! Created 2017-03-10
! Updated 2017-04-17

 module write_UV

! Declare use of Fortran 90 netCDF package
! If this part doesn't compile: module load sips/gcc

use netcdf

implicit none

contains

! ---------------------------------------------------------------------
!
! Begin subroutine make_nc(ALL VARIABLES GO HERE,XL,YL)
!
! ---------------------------------------------------------------------

      subroutine make_UV_nc(Save_Lat,Save_Lon,save_SolZen,save_View_angle,&
             save_View_phi,save_solar_phi,Save_average_Tau_Ocean_UV,&
             Save_average_Omega_Ocean_UV,IX10KM,IY10KM,&
             ipart,yy_part1,yy_part2,Sat_Flag,&
             Save_mode_F_FUV,Save_mode_C_FUV,Save_Small_weighting_FUV,&
             Save_Index_Albedo,Save_Index_Height,Save_Fitting_Error,&
             Save_ref_allwav_uv)
                   
                   
!     Define input variable dimensions, etc.
!     Define parameters from Main_Driver.f90

      include 'read_Sat_MODIS.inc'
      include 'mod04.inc'

!     Declare input variables

      include 'Save_data_declare.inc'
      include 'Save_data_UV_declare.inc'

!     Declare output variables
       character(len=10) :: Sat_Flag  
       integer yy1,yy2,ipart,xx1,yy_part1,yy_part2,num_joint
       integer yy3,yy4,number_of_parts         
      Real nc_lat(IX10KM,IY10KM),nc_lon(IX10KM,IY10KM)
      real nc_solzen(IX10KM,IY10KM),nc_senzen(IX10KM,IY10KM)
      real nc_solazi(IX10KM,IY10KM),nc_senazi(IX10KM,IY10KM)
      real nc_aod_ocean(IX10KM,IY10KM,NWAV_uv+7) 
      real nc_solindx_lg(IX10KM,IY10KM)
      real nc_solindx_sm(IX10KM,IY10KM)
      real nc_Albedo_indx(IX10KM,IY10KM)
      real nc_Height_indx(IX10KM,IY10KM)
      real nc_W0_Ocean_UV(IX10KM,IY10KM,NWAV_uv)
      real nc_sm_weight(IX10KM,IY10KM)
      real nc_fitting_error(IX10KM,IY10KM)
      real nc_ref_allwav_uv(IX10KM,IY10KM,NWAV_uv+7)
       integer  ixp, IYp,A_Num_part,iscan,idata,ik,ij
      real  new_array(IX10KM,IY10KM)
!     Declare variables for netCDF processing

      integer XL,YL,ZL,li,lj, YY1_new, YY2_new,ZL2,ZL3,ZL4,ZLS
      integer ncid,grpid
      character (len=*), parameter :: nc_name = 'UV_output.nc'
      real fv3,fv4

!     Open the netCDF file already created via ncgen
       
                    
      call check( nf90_open(nc_name,NF90_WRITE, ncid) )
       ipart =1
      XL = IX10KM
      YL = IY10KM
      ZL = NWAV_uv 
      ZL3 = NWAV_uv+7
      ZL4 = Land_Sol4 
      fv3 = -990.
      fv4 = -9990. 
    
!     Fill output variables with data
!     Apply scale factors and offsets, excluding fill values

          
            XX1   =    XL
           A_Num_part = YL  
           
           if( ipart.eq.1) then 
                 YY1=1
                 YY2=yy_part2/ILINE   
           elseif( ipart.gt.1) then
                 YY1=(yy_part1/ILINE) +1
                 YY2= yy_part2/ILINE  
           endif 
       
           print*,'XL,YL', XL,YL,ipart 
           
        
      nc_lat(1:XL,1:YL) = Save_Lat(1:XL,1:YL) 
      nc_lat = scale_offset(nc_lat, XX1,SCALE1,OFFSET1,fv3,XL,YL) 
      nc_lon(1:XL,1:YL) = Save_Lon(1:XL,1:YL)
      nc_lon = scale_offset(nc_lon, XX1,SCALE1,OFFSET1,fv3,XL,YL)  
      nc_solzen(1:XL,1:YL) = save_SolZen(1:XL,1:YL)
      nc_solzen = scale_offset(nc_solzen, XX1,SCALE2,OFFSET2,fv4,XL,YL)
      nc_senzen(1:XL,1:YL) = save_View_angle(1:XL,1:YL)
      nc_senzen = scale_offset(nc_senzen, XX1,SCALE2,OFFSET2,fv4,XL,YL) 
      nc_senazi(1:XL,1:YL) = save_View_phi(1:XL,1:YL)
      nc_senazi = scale_offset(nc_senazi, XX1,SCALE2,OFFSET2,fv4,XL,YL)
      nc_solazi(1:XL,1:YL) = save_solar_phi(1:XL,1:YL)
      nc_solazi = scale_offset(nc_solazi, XX1,SCALE2,OFFSET2,fv4,XL,YL)
      nc_solindx_lg(1:XL,1:YL) = Save_mode_C_FUV(1:XL,1:YL)
      nc_solindx_lg = scale_offset(nc_solindx_lg,XX1,SCALE1,OFFSET1,fv4,XL,YL)
      nc_solindx_sm(1:XL,1:YL) = Save_mode_F_FUV(1:XL,1:YL)
      nc_solindx_sm = scale_offset(nc_solindx_sm,XX1,SCALE1,OFFSET1,fv4,XL,YL)
     
      nc_Albedo_indx(1:XL,1:YL) = Save_Index_Albedo(1:XL,1:YL)
      nc_Albedo_indx = scale_offset(nc_Albedo_indx,XX1,SCALE1,OFFSET1,fv4,XL,YL) 
      
      nc_Height_indx(1:XL,1:YL) = Save_Index_Height(1:XL,1:YL)
      nc_Height_indx = scale_offset(nc_Height_indx,XX1,SCALE1,OFFSET1,fv4,XL,YL)
       
      nc_sm_weight(1:XL,1:YL) = Save_Small_weighting_FUV(1:XL,1:YL) 
      nc_sm_weight=scale_offset(nc_sm_weight, XX1,SCALE3,OFFSET3,fv4,XL,YL) 
      nc_fitting_error(1:XL,1:YL) = Save_fitting_error(1:XL,1:YL)
      nc_fitting_error=scale_offset(nc_fitting_error,XX1,SCALE3,OFFSET3,fv4,XL,YL)        
      nc_W0_Ocean_UV(1:XL,1:YL,1:NWAV_uv)=Save_average_Omega_Ocean_UV(1:XL,1:YL,1:NWAV_uv)
      nc_ref_allwav_uv(1:XL,1:YL,1:NWAV_uv+7)= Save_ref_allwav_uv(1:XL,1:YL,1:NWAV_uv+7) 
         do li=1,ZL 
        nc_W0_Ocean_UV(1:XL,1:YL,li) = scale_offset(nc_W0_Ocean_UV(1:XL,1:YL,li),&
         XX1,SCALE3,OFFSET3,fv4,XL,YL)
         enddo 
         
         do li=1,ZL3 
         nc_ref_allwav_uv(1:XL,1:YL,li) = scale_offset(nc_ref_allwav_uv(1:XL,1:YL,li),&
         XX1,SCALE4,OFFSET4,fv4,XL,YL)
         enddo
         
        
      
      nc_aod_ocean(1:XL,1:YL,1:NWAV_uv+7) = &
      Save_average_Tau_Ocean_UV(1:XL,1:YL,1:NWAV_uv+7)
    
      do li=1,Zl3
      nc_aod_ocean(1:XL,1:YL,li) = scale_offset(nc_aod_ocean(1:XL,1:YL,li),&
                                XX1,SCALE3,OFFSET3,fv4,XL,YL)  
      enddo

      
!     Open netCDF file group geolocation_data
1000    continue
      call check( nf90_inq_ncid(ncid, 'geolocation_data', grpid) )

!     Write variables to netCDF group geolocation_data
         
      call write_nc_2d(nc_lat,'latitude', XX1,YY1,YY2,grpid,0,A_Num_part) 
      call write_nc_2d(nc_lon,'longitude',XX1,YY1,YY2,grpid,0,A_Num_part) 
      call write_nc_2d(nc_solzen,'solar_zenith_angle',XX1,YY1,YY2,grpid,1,A_Num_part)  
      call write_nc_2d(nc_senzen,'sensor_zenith_angle',XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_senazi,'sensor_azimuth_angle',XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_solazi,'solar_azimuth_angle',XX1,YY1,YY2,grpid,1,A_Num_part)
       
        
!     Open netCDF file group geophysical_data

      call check( nf90_inq_ncid(ncid, 'geophysical_data', grpid) )
         
!     Write 2D variables to netCDF group geophysical_data 
!     Write 2D variables to netCDF group geophysical_data
       call write_nc_2d(nc_solindx_sm,'Fine_MODE_used',&
            XX1,YY1,YY2,grpid,1,A_Num_part)
       call write_nc_2d(nc_solindx_lg,'Coarse_MODE_used',&
            XX1,YY1,YY2,grpid,1,A_Num_part) 
        call write_nc_2d(nc_Albedo_indx,'Index_Albedo',&
            XX1,YY1,YY2,grpid,1,A_Num_part) 
        call write_nc_2d(nc_Height_indx,'Height',&
            XX1,YY1,YY2,grpid,1,A_Num_part)   
             
           
       call write_nc_2d(nc_sm_weight,'Optical_Depth_Ratio_Small_Ocean_used',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
        call write_nc_2d(nc_fitting_error,'Fitting_Error_Height',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
          call check( nf90_inq_ncid(ncid, 'geophysical_data', grpid) ) 
      call write_nc_3d(nc_aod_ocean,'Effective_Optical_Depth_Average_Ocean_UV',&
           XX1,YY1,YY2,ZL3,grpid,A_Num_part)
      call write_nc_3d(nc_W0_Ocean_UV,'Single_Scattering_Albedo_Average_Ocean_UV',&
           XX1,YY1,YY2,ZL,grpid,A_Num_part)
      call write_nc_3d(nc_ref_allwav_uv,'Mean_Reflectance_UV_ocean',&
           XX1,YY1,YY2,ZL3,grpid,A_Num_part)
      
!     Close the netCDF file
3000   continue
       if( ipart .eq.1)call check( nf90_close(ncid) )
      print *, 'Wrote variables to ',nc_name

      return
       end subroutine make_UV_nc

! ---------------------------------------------------------------------
!
! Begin function scale_offset(var,sox,soy,scl,ofs,fv)
!
!     var = data array to be scaled
!     sox,soy = dimensions of var
!     scl,ofs = scale factor and offset to be applied
!     fv = fill value (assumed to be a large negative real)
!
! ---------------------------------------------------------------------

      function scale_offset(var,sox,scl,ofs,fv,xx,yy)
        include 'mod04.inc'
      integer sox,i,j,ex,xx,yy
      real, dimension(sox,yy) :: var,scale_offset
      real*8 scl,ofs
      real fv
        
       do i=1,sox
         do j=1,yy
            if (var(i,j) .gt. fv) then
               ex = exponent(var(i,j))
               if (ex .lt. -10) then
                  var(i,j) = 0.
               endif 
               var(i,j) = var(i,j)*scl + ofs
              
            end if
         enddo
      enddo

      scale_offset = var

      end function scale_offset
 
! ---------------------------------------------------------------------
!
! Begin subroutine check(status)
!
!     Checks for errors in netCDF handling
!
! ---------------------------------------------------------------------

      subroutine check(status)

      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if

      end subroutine check
 
! ---------------------------------------------------------------------
!
! Begin subroutine write_nc_2d(var_real,var_name,nx,ny,grpid,dty,A_Num_part)
!
!     var_real = data to be written, type real
!     var_name = name of variable to write in netCDF
!     nx,ny = dimensions of var_real
!     ncid,grpid = IDs for currently open netCDF file and group
!     dty = 0 to write variable as real, 1 to write as integer
!
! ---------------------------------------------------------------------

      subroutine write_nc_2d(var_real,var_name,nx,ny0,ny,grpid,dty,A_Num_part)
       
      include 'mod04.inc'
      character(len=*),intent(in) :: var_name
      integer start,count,A_Num_part
      integer grpid,varid,nx,ny,dty,ny0,xx,yy
      real, dimension(nx,A_Num_part) :: var_real
      integer*2, dimension(nx,A_Num_part):: var_int
      integer ii,ij,ikk
       
!     Get the variable ID within the group ID passed from make_nc
      call check( nf90_inq_varid(grpid, var_name, varid) )
       
      if (dty .gt. 0) then
        
!        If indicated, convert to integer before writing to file
      
           do ii=1,nx 
            do ij=1,A_Num_part
             var_int(ii,ij) = nint(var_real(ii,ij))
            enddo
         enddo 
!        if( A_Num_part .gt. 406 )then
!         call check( nf90_put_var(grpid,varid,var_int))
!           print*,var_name,count,A_Num_part,nx,ny0
!        else
        call check(nf90_put_var(grpid,varid,var_int,start=(/1,ny0/),&
                            count = (/nx,A_Num_part /)))
!        Endif
                            
      else
         
!        Otherwise, write variable to file as real
         call check(nf90_put_var(grpid,varid,var_real,start=(/1,ny0/),&
                            count = (/nx,A_Num_part /)))

      end if
       

      return
      end subroutine write_nc_2d
      ! ---------------------------------------------------------------------
!
! Begin subroutine write_nc_3d
!                  (var_real,var_name,nx,ny,nz,grpid)

!       nx,ny0,ny,grpid
!
!     var_real = data to be written, type real
!     var_name = name of variable to write in netCDF
!     nx,ny,nz = dimensions of var_real
!     ncid,grpid = IDs for currently open netCDF file and group
!
! ---------------------------------------------------------------------

      subroutine write_nc_3d(var_real,var_name,nx,ny0,ny,nz,grpid,A_Num_part)
       include 'mod04.inc' 
      character(len=*),intent(in) :: var_name
      integer ncid,grpid,varid,nx,ny,nz,ny0,A_Num_part
      real, dimension(nx,A_Num_part,nz) :: var_real 
      integer*2, dimension(nz,nx,A_Num_part) :: var_int 
      integer ii,ij,ik,iz,ikk,xx,yy
      integer st(3),ct(3)

 
         
!     Get the variable ID within the group ID passed from make_nc
      call check( nf90_inq_varid(grpid, var_name, varid) )
     
      do ii=1,nx 
         do ij=1,A_Num_part
             do ik=1,nz  
               var_int(ik,ii,ij) = nint(var_real(ii,ij,ik))
            enddo
         enddo
      enddo
       
      
!     print*,'3D',ny0,ny,ikk,nx,A_Num_part
!      call check( nf90_put_var(grpid,varid,var_int))
       call check( nf90_put_var(grpid,varid,var_int,start=(/1,1,ny0/),&
                             count = (/nz,nx,A_Num_part/)))
!
      return
      end subroutine write_nc_3d
      
end module write_UV


