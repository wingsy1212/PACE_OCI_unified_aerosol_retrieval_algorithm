! write_vaerdt.f90

! Draft netCDF writer for VIIRS Dark Target
! Simplified to fewer variables for testing

! Virginia Sawyer
! Created 2017-03-10
! Updated 2017-04-17

 module write_Pace_merged_1KM
 
! Declare use of Fortran 90 netCDF package
! If this part doesn't compile: module load sips/gcc

use netcdf
USE OCIUAAER_Config_Module

implicit none

contains

! ---------------------------------------------------------------------
!
! Begin subroutine make_nc(ALL VARIABLES GO HERE,XL,YL)
!
! ---------------------------------------------------------------------

                  
Subroutine write_Output_merged_1KM(l1b_nXTrack, l1b_nLines, Latitude, Longitude, &
	       UVAI, Residue_1km, Reflectivity_1km, out_file_1KM) 
	      

	   
	       
!     Define input variable dimensions, etc.
!     Define parameters from Main_Driver.f90
      Include 'common_l1b_var.inc'
      include 'output_Variables.inc' 

!     Declare input variables
      CHARACTER(255) :: out_file_1KM

!     Declare output variables 
      character(len=10) :: Sat_Flag  
      integer yy1,yy2,ipart,xx1,yy_part1,yy_part2,num_joint,A_Num_part
      integer yy3,yy4,number_of_parts   
!     Declare variables for netCDF processing
      real nc_lat_1km(l1b_nXTrack, l1b_nLines)
      real nc_lon_1km(l1b_nXTrack, l1b_nLines)
      real nc_uvai_1km(l1b_nXTrack, l1b_nLines)
      real nc_reflectivity_1km(l1b_nXTrack, l1b_nLines,2)
      real nc_residue_1km(l1b_nXTrack, l1b_nLines) 
      
      integer XL,YL,ZL,li,lj, YY1_new, YY2_new,ZL1,ZL2,ZL3,ZL4,ZLS
      integer IL,IX,IY
      integer ncid,grpid, retval
!      character (len=*), parameter :: nc_name = 'PACE_output.nc'
      real fv3,fv4

!     Open the netCDF file already created via ncgen
       
                    
!       call check( nf90_open(nc_name,NF90_WRITE, ncid) )
       retval= nf90_open(trim(out_file_1KM),NF90_WRITE, ncid)
       IF (retval /= nf90_noerr) THEN
            PRINT *, 'Error: Unable to open 1KM file for writing'
            STOP
      END IF

      ! Write the attributes :time_coverage_start and :time_coverage_end
       retval = nf90_put_att(ncid, NF90_GLOBAL, 'time_coverage_start', cfg%coverage_start)
       retval = nf90_put_att(ncid, NF90_GLOBAL, 'time_coverage_end', cfg%coverage_end)
       retval = nf90_put_att(ncid, NF90_GLOBAL, 'history', cfg%history)

       ipart =1
    
      ZL1 = 3
      ZL2 = 2
      ZL3 = 5
      ZL4 = 9
      fv3 = -9999
      fv4 = -9999 
 
!     Fill output variables with data
!     Apply scale factors and offsets, excluding fill values
! Adapted from a previous written subroutine that used to write in chunks
! A little confusing  but works          
                 XX1   =    XL
                 A_Num_part = YL   
                 YY1=1
                 YY2=YL 
           
           print*,'XX1,YY2', XX1,YY2
              

! Native 1km resolution variables   
        
      nc_lat_1km(1:l1b_nXTrack, 1:l1b_nLines) = Latitude(1:l1b_nXTrack, 1:l1b_nLines) 
      nc_lat_1km = scale_offset(nc_lat_1km, l1b_nXTrack,SCALE1,OFFSET1,fv3,l1b_nXTrack, l1b_nLines) 

      nc_lon_1km(1:l1b_nXTrack, 1:l1b_nLines) = Longitude(1:l1b_nXTrack, 1:l1b_nLines)
      nc_lon_1km = scale_offset(nc_lon_1km, l1b_nXTrack,SCALE1,OFFSET1,fv3,l1b_nXTrack, l1b_nLines)  
      
      nc_uvai_1km(1:l1b_nXTrack, 1:l1b_nLines) = UVAI(1:l1b_nXTrack, 1:l1b_nLines)
      nc_uvai_1km = scale_offset(nc_uvai_1km, l1b_nXTrack,SCALE3,OFFSET3,fv4,l1b_nXTrack, l1b_nLines)  

      nc_residue_1km(1:l1b_nXTrack, 1:l1b_nLines) = Residue_1km(1:l1b_nXTrack, 1:l1b_nLines)
      nc_residue_1km = scale_offset(nc_residue_1km, l1b_nXTrack,SCALE3,OFFSET3,fv4,l1b_nXTrack, l1b_nLines)  

      nc_reflectivity_1km(1:l1b_nXTrack, 1:l1b_nLines, 1) = Reflectivity_1km(1, 1:l1b_nXTrack, 1:l1b_nLines)
      nc_reflectivity_1km(1:l1b_nXTrack, 1:l1b_nLines, 2) = Reflectivity_1km(2, 1:l1b_nXTrack, 1:l1b_nLines)
      do li=1,2
      nc_reflectivity_1km(1:l1b_nXTrack, 1:l1b_nLines, li) = &
          scale_offset(nc_reflectivity_1km(1:l1b_nXTrack, 1:l1b_nLines, li), &
                    l1b_nXTrack,SCALE3,OFFSET3,fv4,l1b_nXTrack, l1b_nLines)  
      enddo
      

 
!     Open netCDF file group geolocation_data

      call check( nf90_inq_ncid(ncid, 'geolocation_data', grpid) )

!     Write variables to netCDF group geolocation_data 
      call write_nc_2d(nc_lat_1km,'Latitude_1km',l1b_nXTrack,1,l1b_nLines,grpid,0,l1b_nLines)  
      call write_nc_2d(nc_lon_1km,'Longitude_1km',l1b_nXTrack,1,l1b_nLines,grpid,0,l1b_nLines) 
      
        
!     Open netCDF file group geophysical_data

      call check( nf90_inq_ncid(ncid, 'geophysical_data', grpid) )  
      
       
!  NUV variables    

      call write_nc_2d(nc_uvai_1km,'NUV_AerosolIndex_1km',&
           l1b_nXTrack,1,l1b_nLines,grpid,1,l1b_nLines)  

      call write_nc_2d(nc_residue_1km,'NUV_Residue_1km',&
           l1b_nXTrack,1,l1b_nLines,grpid,1,l1b_nLines)  

      call write_nc_3d(nc_reflectivity_1km,'NUV_Reflectivity_1km',&
           l1b_nXTrack,1,l1b_nLines,ZL2,grpid,1,l1b_nLines)  
	                
!  Close the netCDF file
if( ipart .eq.1)call check( nf90_close(ncid) )
print *, 'Wrote variables to ', out_file_1KM

return

end subroutine write_Output_merged_1KM

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
   	     if(var_int(ii,ij) .lt. -9999) var_int(ii,ij)= -9999
            enddo
         enddo 
        call check(nf90_put_var(grpid,varid,var_int,start=(/1,ny0/),&
                            count = (/nx,A_Num_part /)))
      else
!        Otherwise, write variable to file as real
         call check(nf90_put_var(grpid,varid,var_real,start=(/1,ny0/),&
                            count = (/nx,A_Num_part /)))

      end if
       
      Print *, 'Written NETCDF-4 variable = ', var_name
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
              if(var_int(ik,ii,ij) .lt. -9999) var_int(ik,ii,ij)= -9999
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
    
    Print *, 'Written NETCDF-4 variable = ', var_name
    
    return
    end subroutine write_nc_3d


      subroutine write_nc_4d(var_real,var_name,nx,ny0,ny,nz,nl,grpid,dty,A_Num_part)
       include 'output_Variables.inc' 
      character(len=*),intent(in) :: var_name
      integer ncid,grpid,varid,nx,ny,nz, nl, ny0,A_Num_part
      real, dimension(nx,A_Num_part,nz,nl) :: var_real 
      integer*2, dimension(nl,nz,nx,A_Num_part) :: var_int 
      real, dimension(nl,nz,nx,A_Num_part) :: dummy_real
      integer ii,ij,ik,il,iz,ikk,xx,yy,dty
      integer st(3),ct(3)

 
         
!     Get the variable ID within the group ID passed from make_nc
      call check( nf90_inq_varid(grpid, var_name, varid) )
        
 
     if (dty .gt. 0) then 
      do ii=1,nx 
         do ij=1,A_Num_part
             do ik=1,nz  
	       do il=1,nl
               var_int(il,ik,ii,ij) = nint(var_real(ii,ij,ik,il))
              if(var_int(il,ik,ii,ij) .lt. -9999) var_int(il,ik,ii,ij)= -9999
	      enddo
            enddo
         enddo
      enddo 
       call check( nf90_put_var(grpid,varid,var_int,start=(/1,1,1,ny0/),&
                             count = (/nl,nz,nx,A_Num_part/)))
    Else  
           do ii=1,nx 
           do ij=1,A_Num_part
             do ik=1,nz  
	     do il=1,nl
             dummy_real(il,ik,ii,ij) = (var_real(ii,ij,ik,il)) 
	     enddo
            enddo
         enddo
      enddo                          
        call check( nf90_put_var(grpid,varid,dummy_real,start=(/1,1,1,ny0/),&
                             count = (/nl,nz,nx,A_Num_part/)))                           
!
    end if
    
    Print *, 'Written NETCDF-4 variable = ', var_name
    
    return
    end subroutine write_nc_4d

      
end module write_Pace_merged_1KM

 



