! write_vaerdt.f90

! Draft netCDF writer for VIIRS Dark Target
! Simplified to fewer variables for testing

! Virginia Sawyer
! Created 2017-03-10
! Updated 2017-04-17

 module write_intrim_forUVtau
 
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

       
        subroutine write_Output_forUVtau(Ret_Xtrack,Ret_Lines,Ret_Lat,Ret_Lon,&
                   Ret_SolZen,Ret_View_angle,Ret_View_phi,Ret_solar_phi,uvdbdtaod,&
                   Ret_tau_ocean,Month,nc_name)
                   
!     Define input variable dimensions, etc.
!     Define parameters from Main_Driver.f90

      include 'output_Variables.inc' 
!     Declare input variables
 
       CHARACTER(255) :: nc_name

!     Declare output variables
       character(len=10) :: Sat_Flag  
       integer yy1,yy2,ipart,xx1,yy_part1,yy_part2,num_joint
       integer yy3,yy4,number_of_parts,month
       Real nc_lat(Ret_Xtrack,Ret_Lines),nc_lon(Ret_Xtrack,Ret_Lines)
       real nc_solzen(Ret_Xtrack,Ret_Lines),nc_senzen(Ret_Xtrack,Ret_Lines)
       real nc_solazi(Ret_Xtrack,Ret_Lines),nc_senazi(Ret_Xtrack,Ret_Lines) 
       real nc_aod_ocean(Ret_Xtrack,Ret_Lines,9)
       real uvdbdtaod(Ret_Xtrack,Ret_Lines,5)  
       real nc_aod_land (Ret_Xtrack,Ret_Lines,5) 
       integer  ixp, IYp,A_Num_part,iscan,idata,ik,ij 
!     Declare variables for netCDF processing

      integer XL,YL,ZL,li,lj, YY1_new, YY2_new,ZL2,ZL3,ZL4,ZLS
      integer IL,IX,IY
      integer ncid,grpid
!      character (len=*), parameter :: nc_name = 'Interm_file.nc'
      real fv3,fv4

!     Open the netCDF file already created via ncgen
       
                    
      call check( nf90_open(nc_name,NF90_WRITE, ncid) )
       ipart =1
    
      XL = Ret_Xtrack
      YL = Ret_Lines
      ZL3 = 5
      ZL4 = 9
      fv3 = -990.
      fv4 = -9990. 
    
!     Fill output variables with data
!     Apply scale factors and offsets, excluding fill values
! Adapted from a previous written subroutine that used to write in chunks
! A little confusing  but works          
                 XX1   =    XL
                 A_Num_part = YL   
                 YY1=1
                 YY2=YL  
!           
      nc_lat(1:XL,1:YL) = Ret_Lat (1:XL,1:YL) 
      nc_lat = scale_offset(nc_lat, XX1,SCALE1,OFFSET1,fv3,XL,YL) 
      nc_lon(1:XL,1:YL) = Ret_Lon (1:XL,1:YL)
      nc_lon = scale_offset(nc_lon, XX1,SCALE1,OFFSET1,fv3,XL,YL) 
      
      nc_solzen(1:XL,1:YL) = Ret_SolZen (1:XL,1:YL)
      nc_solzen = scale_offset(nc_solzen, XX1,SCALE2,OFFSET2,fv4,XL,YL)
      nc_senzen(1:XL,1:YL) = Ret_View_angle(1:XL,1:YL)
      nc_senzen = scale_offset(nc_senzen, XX1,SCALE2,OFFSET2,fv4,XL,YL) 
      nc_senazi(1:XL,1:YL) = Ret_View_phi(1:XL,1:YL)
      nc_senazi = scale_offset(nc_senazi, XX1,SCALE2,OFFSET2,fv4,XL,YL)
      nc_solazi(1:XL,1:YL) = Ret_solar_phi(1:XL,1:YL)
      nc_solazi = scale_offset(nc_solazi, XX1,SCALE2,OFFSET2,fv4,XL,YL) 
          
      nc_aod_ocean(1:XL,1:YL,1:zl4)  = Ret_tau_ocean(1:XL,1:YL,1:zl4) 
      nc_aod_land(1:XL,1:YL,1:zl3)   =  uvdbdtaod(1:XL,1:YL,1:zl3)
      
    
      do li=1,Zl4
      nc_aod_ocean(1:XL,1:YL,li) = scale_offset(nc_aod_ocean(1:XL,1:YL,li),&
                                   XX1,SCALE3,OFFSET3,fv4,XL,YL) 
     enddo 
      do li=1,Zl3                              
      nc_aod_land(1:XL,1:YL,li) = scale_offset(nc_aod_land(1:XL,1:YL,li),&
                                 XX1,SCALE3,OFFSET3,fv4,XL,YL)                             
      enddo

      
!     Open netCDF file group geolocation_data
1000    continue
      call check( nf90_inq_ncid(ncid, 'geolocation_data', grpid) )

!     Write variables to netCDF group geolocation_data
      call write_nc_one(month,'month',grpid,1)
      call write_nc_2d(nc_lat,'latitude', XX1,YY1,YY2,grpid,0,A_Num_part) 
      call write_nc_2d(nc_lon,'longitude',XX1,YY1,YY2,grpid,0,A_Num_part) 
      call write_nc_2d(nc_solzen,'solar_zenith_angle',XX1,YY1,YY2,grpid,1,A_Num_part)  
      call write_nc_2d(nc_senzen,'sensor_zenith_angle',XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_senazi,'sensor_azimuth_angle',XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_solazi,'solar_azimuth_angle',XX1,YY1,YY2,grpid,1,A_Num_part) 
!     Open netCDF file group geophysical_data

       call check( nf90_inq_ncid(ncid, 'geophysical_data', grpid) )   
          call check( nf90_inq_ncid(ncid, 'geophysical_data', grpid) ) 
          call write_nc_3d(nc_aod_land,'Optical_Depth_Land',&
           XX1,YY1,YY2,ZL3,grpid,1,A_Num_part)  
       call write_nc_3d(nc_aod_ocean,'Optical_Depth_Ocean',&
           XX1,YY1,YY2,ZL4,grpid,1,A_Num_part) 
        
!     Close the netCDF file
3000   continue
       if( ipart .eq.1)call check( nf90_close(ncid) )
      print *, 'Wrote variables to ',nc_name 
      return
       end subroutine write_Output_forUVtau

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
       include 'output_Variables.inc' 
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
  !             print*,var(i,j),scl,ofs
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
! Begin subroutine write_nc_one(var_real,var_name,grpid,dty)
!
!     var_real = data to be written, type real
!     var_name = name of variable to write in netCDF
!     ncid,grpid = IDs for currently open netCDF file and group
!     dty = 0 to write variable as real, 1 to write as integer
!
! ---------------------------------------------------------------------

      subroutine write_nc_one(var_int,var_name,grpid,dty)

      include 'output_Variables.inc'
      character(len=*),intent(in) :: var_name
      integer start,count,A_Num_part
      integer grpid,varid,nx,dty,var_int
      real var_real

!     Get the variable ID within the group ID passed from make_nc
      call check( nf90_inq_varid(grpid, var_name, varid) )

      if (dty .gt. 0) then
!        If indicated, convert to integer before writing to file
         call check(nf90_put_var(grpid,varid,var_int))
!        Endif

      else
        var_real=float(var_int)
!        Otherwise, write variable to file as real
         call check(nf90_put_var(grpid,varid,var_real))

      end if


      return
      end subroutine write_nc_one
 
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
       
      include 'output_Variables.inc' 
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

      subroutine write_nc_3d(var_real,var_name,nx,ny0,ny,nz,grpid,dty,A_Num_part)
       include 'output_Variables.inc' 
      character(len=*),intent(in) :: var_name
      integer ncid,grpid,varid,nx,ny,nz,ny0,A_Num_part
      real, dimension(nx,A_Num_part,nz) :: var_real 
      integer*2, dimension(nz,nx,A_Num_part) :: var_int 
      real, dimension(nz,nx,A_Num_part) :: dummy_real
      integer ii,ij,ik,iz,ikk,xx,yy,dty
      integer st(3),ct(3)

 
         
!     Get the variable ID within the group ID passed from make_nc
      call check( nf90_inq_varid(grpid, var_name, varid) )
        
 
     if (dty .gt. 0) then 
      do ii=1,nx 
         do ij=1,A_Num_part
             do ik=1,nz  
               var_int(ik,ii,ij) = nint(var_real(ii,ij,ik))
            enddo
         enddo
      enddo 
       call check( nf90_put_var(grpid,varid,var_int,start=(/1,1,ny0/),&
                             count = (/nz,nx,A_Num_part/)))
    Else  
           do ii=1,nx 
           do ij=1,A_Num_part
             do ik=1,nz  
             dummy_real(ik,ii,ij) = (var_real(ii,ij,ik))
            enddo
         enddo
      enddo                          
        call check( nf90_put_var(grpid,varid,dummy_real,start=(/1,1,ny0/),&
                             count = (/nz,nx,A_Num_part/)))                           
!
    end if
      return
      end subroutine write_nc_3d
      
end module write_intrim_forUVtau






