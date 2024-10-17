! write_vaerdt.f90

! Draft netCDF writer for VIIRS Dark Target
! Simplified to fewer variables for testing

! Virginia Sawyer
! Created 2017-03-10
! Updated 2017-04-17

 module write_Pace_merged
 
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

                  
Subroutine write_Output_merged(Ret_Lat,Ret_Lon,Ret_SolZen,Ret_View_angle,&
              Ret_View_phi,Ret_solar_phi,Ret_Xtrack,Ret_Lines,Ret_Small_weighting,&
	      Land_sea_flag,Ret_ref_LandOceanwOutUV,&
              Ret_ref_LandOcean_UV,Ret_Tau_LandOcean,&
              Ret_average_Omega_Ocean_UV,Ret_Index_Height,Cloud_Frac_LandOcean,&
	      Ret_Quality_LandOcean,Ret_Quality_LandOcean_W0,& 
              NUV_AI, NUV_COD, NUV_CldFrac, UVReflectivity, UVResidue, &
	      NUV_SSA, NUV_ALH, NUV_ACAOD, NUV_AerCorrCOD,& 
	      NUV_ACAODVsHeight, NUV_AerCorrCODVsHeight, &
	      NUV_FinalAlgorithmFlagsACA, NUV_UncertaintyACAODToSSA, NUV_UncertaintyCODToSSA,&
	      out_file)

	   
	       
!     Define input variable dimensions, etc.
!     Define parameters from Main_Driver.f90
      Include 'common_l1b_var.inc'
      include 'output_Variables.inc' 

!     Declare input variables
      CHARACTER(255) :: out_file

!     Declare output variables 
      character(len=10) :: Sat_Flag  
      integer yy1,yy2,ipart,xx1,yy_part1,yy_part2,num_joint
      integer yy3,yy4,number_of_parts         
      Real nc_lat(Ret_Xtrack,Ret_Lines),nc_lon(Ret_Xtrack,Ret_Lines)
      real nc_solzen(Ret_Xtrack,Ret_Lines),nc_senzen(Ret_Xtrack,Ret_Lines)
      real nc_solazi(Ret_Xtrack,Ret_Lines),nc_senazi(Ret_Xtrack,Ret_Lines)
      real nc_aod_ocean_land(Ret_Xtrack,Ret_Lines,9)  
      real nc_sm_weight(Ret_Xtrack,Ret_Lines)
      real nc_ref_allwav(Ret_Xtrack,Ret_Lines,7)
      real nc_ref_allwav_UV(Ret_Xtrack,Ret_Lines,2)
      real nc_ls_flag(Ret_Xtrack,Ret_Lines) 
      real nc_W0_Ocean_UV(Ret_Xtrack,Ret_Lines,3)
      real nc_Height_indx(Ret_Xtrack,Ret_Lines) 
      real nc_Cloud_Frac_LandOcean(Ret_Xtrack,Ret_Lines)
      real nc_Quality_LandOcean(Ret_Xtrack,Ret_Lines)
      real nc_Quality_LandOcean_W0(Ret_Xtrack,Ret_Lines)
      integer  ixp, IYp,A_Num_part,iscan,idata,ik,ij 

!     Declare variables for netCDF processing
      
      real nc_uvai(Ret_Xtrack,Ret_Lines)
      real nc_reflectivity(Ret_Xtrack,Ret_Lines,2)
      real nc_residue(Ret_Xtrack,Ret_Lines)
      real nc_cod(Ret_Xtrack,Ret_Lines)
      real nc_cldfrac(Ret_Xtrack,Ret_Lines)
      real nc_nuvalh(Ret_Xtrack,Ret_Lines)
      real nc_nuvssa(Ret_Xtrack,Ret_Lines,5)
      real nc_nuvacaod(Ret_Xtrack,Ret_Lines,3)
      real nc_nuvaercorrcod(Ret_Xtrack,Ret_Lines)
      real nc_nuvacaodvshgt(Ret_Xtrack,Ret_Lines,3,5)
      real nc_nuvaercorrcodvshgt(Ret_Xtrack,Ret_Lines,5)
      real nc_nuvacqf(Ret_Xtrack,Ret_Lines)
      real nc_nuvacuncer1(Ret_Xtrack,Ret_Lines,2)
      real nc_nuvacuncer2(Ret_Xtrack,Ret_Lines,2)
      
      
      integer XL,YL,ZL,li,lj, YY1_new, YY2_new,ZL1,ZL2,ZL3,ZL4,ZLS
      integer IL,IX,IY
      integer ncid,grpid
!      character (len=*), parameter :: nc_name = 'PACE_output.nc'
      real fv3,fv4

!     Open the netCDF file already created via ncgen
       
                    
!       call check( nf90_open(nc_name,NF90_WRITE, ncid) )
       call check( nf90_open(trim(out_file),NF90_WRITE, ncid) )
       ipart =1
    
      XL = Ret_Xtrack
      YL = Ret_Lines
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
              


! Aggregate resolution variables

      nc_lat(1:XL,1:YL) = Ret_Lat(1:XL,1:YL) 
      nc_lat = scale_offset(nc_lat, XX1,SCALE1,OFFSET1,fv3,XL,YL) 

      nc_lon(1:XL,1:YL) = Ret_Lon(1:XL,1:YL)
      nc_lon = scale_offset(nc_lon, XX1,SCALE1,OFFSET1,fv3,XL,YL)  

      nc_solzen(1:XL,1:YL) = Ret_SolZen(1:XL,1:YL)
      nc_solzen = scale_offset(nc_solzen, XX1,SCALE2,OFFSET2,fv4,XL,YL)

      nc_senzen(1:XL,1:YL) = Ret_View_angle(1:XL,1:YL)
      nc_senzen = scale_offset(nc_senzen, XX1,SCALE2,OFFSET2,fv4,XL,YL) 

      nc_senazi(1:XL,1:YL) = Ret_View_phi(1:XL,1:YL)
      nc_senazi = scale_offset(nc_senazi, XX1,SCALE2,OFFSET2,fv4,XL,YL)

      nc_solazi(1:XL,1:YL) = Ret_solar_phi(1:XL,1:YL)
      nc_solazi = scale_offset(nc_solazi, XX1,SCALE2,OFFSET2,fv4,XL,YL)

      nc_sm_weight(1:XL,1:YL) = Ret_Small_weighting(1:XL,1:YL) 
      nc_sm_weight=scale_offset(nc_sm_weight, XX1,SCALE3,OFFSET3,fv4,XL,YL)  

!     
      nc_Height_indx(1:XL,1:YL) =  Ret_Index_Height(1:XL,1:YL)
      nc_Height_indx = scale_offset(nc_Height_indx,XX1,SCALE1,OFFSET1,fv4,XL,YL)
       
      nc_Cloud_Frac_LandOcean(1:XL,1:YL) = Cloud_Frac_LandOcean(1:XL,1:YL) 
      nc_Cloud_Frac_LandOcean= &
        scale_offset(nc_Cloud_Frac_LandOcean,XX1,SCALE3,OFFSET3,fv4,XL,YL)
        
      nc_Quality_LandOcean(1:XL,1:YL) = Ret_Quality_LandOcean(1:XL,1:YL) 
      nc_Quality_LandOcean = &
        scale_offset(nc_Quality_LandOcean,XX1,SCALE3,OFFSET3,fv4,XL,YL)
        
      nc_Quality_LandOcean_W0(1:XL,1:YL) = Ret_Quality_LandOcean_W0(1:XL,1:YL) 
      nc_Quality_LandOcean_W0= &
        scale_offset(nc_Quality_LandOcean_W0,XX1,SCALE3,OFFSET3,fv4,XL,YL)
          
      nc_W0_Ocean_UV(1:XL,1:YL,1:ZL1)= Ret_average_Omega_Ocean_UV(1:XL,1:YL,1:ZL1)
      do li=1,ZL1
      nc_W0_Ocean_UV(1:XL,1:YL,li) = scale_offset(nc_W0_Ocean_UV(1:XL,1:YL,li),&
         XX1,SCALE3,OFFSET3,fv4,XL,YL)
      enddo
       
!  Saperated Reflectance into gas corrected /not corrected
     
      nc_ref_allwav(1:XL,1:YL,1:ZL4-2)= Ret_ref_LandOceanwOutUV(1:XL,1:YL,1:ZL4-2)  
    
      do li=1,ZL4-2
      nc_ref_allwav(1:XL,1:YL,li) = scale_offset(nc_ref_allwav(1:XL,1:YL,li),&
         XX1,SCALE4,OFFSET4,fv4,XL,YL)
      enddo  
 !  UV channels are not gas corrected corrected 
       nc_ref_allwav_UV(1:XL,1:YL,1:ZL2)= Ret_ref_LandOcean_UV(1:XL,1:YL,1:ZL2) 
    
      do li=1,ZL2 
      nc_ref_allwav_uv(1:XL,1:YL,li) = scale_offset(nc_ref_allwav_uv(1:XL,1:YL,li),&
         XX1,SCALE4,OFFSET4,fv4,XL,YL)
      enddo  
      
       
      nc_ls_flag(1:XL,1:YL) =  Land_sea_flag(1:XL,1:YL)  
      nc_aod_ocean_land(1:XL,1:YL,1:ZL4) = Ret_Tau_LandOcean(1:XL,1:YL,1:ZL4) 
      do li=1,ZL4
      nc_aod_ocean_land(1:XL,1:YL,li) = scale_offset(nc_aod_ocean_land(1:XL,1:YL,li),&
                                XX1,SCALE3,OFFSET3,fv4,XL,YL)   
      enddo
         

!  NUV Variables  

        nc_uvai(1:XL,1:YL)= NUV_AI(1:XL,1:YL)
        nc_uvai = scale_offset(nc_uvai,XX1,SCALE3,OFFSET3,fv4,XL,YL) 

        nc_cod(1:XL,1:YL) = NUV_COD(1:XL,1:YL)
        nc_cod = scale_offset(nc_cod, XX1,SCALE3,OFFSET3,fv4,XL,YL)  

        nc_residue(1:XL,1:YL) = UVResidue(1:XL,1:YL)
        nc_residue = scale_offset(nc_residue, XX1,SCALE3,OFFSET3,fv4,XL,YL)  

        nc_reflectivity(1:XL,1:YL,1:2) = UVReflectivity(1:XL,1:YL,1:2)
	do li=1,2
        nc_reflectivity(1:XL,1:YL,li) = scale_offset(nc_reflectivity(1:XL,1:YL,li), &
	                              XX1,SCALE3,OFFSET3,fv4,XL,YL)  
        enddo

        nc_cldfrac(1:XL,1:YL) = NUV_CldFrac(1:XL,1:YL)
        nc_cldfrac = scale_offset(nc_cldfrac, XX1,SCALE3,OFFSET3,fv4,XL,YL)  

        nc_nuvalh(1:XL,1:YL) = NUV_ALH(1:XL,1:YL)
        nc_nuvalh = scale_offset(nc_nuvalh, XX1,SCALE3,OFFSET3,fv4,XL,YL)  
         
        nc_nuvssa(1:XL,1:YL,1:5)= NUV_SSA(1:XL,1:YL,1:5) 
        do li=1,5
        nc_nuvssa(1:XL,1:YL,li) = scale_offset(nc_nuvssa(1:XL,1:YL,li),&
                                XX1,SCALE3,OFFSET3,fv4,XL,YL)  				
        enddo  

          
! NUV  Above-cloud variables   
  
        nc_nuvaercorrcod(1:XL,1:YL) = NUV_AerCorrCOD(1:XL,1:YL)
        nc_nuvaercorrcod = scale_offset(nc_nuvaercorrcod, XX1,SCALE_10,OFFSET_10,fv4,XL,YL)  

        nc_nuvacaod(1:XL,1:YL,1:3)= NUV_ACAOD(1:XL,1:YL,1:3) 
        do li=1,3
        nc_nuvacaod(1:XL,1:YL,li) = scale_offset(nc_nuvacaod(1:XL,1:YL,li),&
                                XX1,SCALE3,OFFSET3,fv4,XL,YL)  				
        enddo  

        nc_nuvaercorrcodvshgt(1:XL,1:YL,1:5) = NUV_AerCorrCODVsHeight(1:XL,1:YL,1:5)
	do li=1,5
        nc_nuvaercorrcodvshgt(1:XL,1:YL,li) = scale_offset(nc_nuvaercorrcodvshgt(1:XL,1:YL,li), &
						XX1,SCALE_10,OFFSET_10,fv4,XL,YL)  
	enddo
	
        nc_nuvacaodvshgt(1:XL,1:YL,1:3,1:5)= NUV_ACAODVsHeight(1:XL,1:YL,1:3,1:5) 
        do li=1,3
	do lj=1,5
        nc_nuvacaodvshgt(1:XL,1:YL,li,lj) = scale_offset(nc_nuvacaodvshgt(1:XL,1:YL,li,lj),&
                                XX1,SCALE3,OFFSET3,fv4,XL,YL)  				
        enddo  
        enddo
	
        nc_nuvacqf(1:XL,1:YL) =  NUV_FinalAlgorithmFlagsACA(1:XL,1:YL) 
        nc_nuvacqf = scale_offset(nc_nuvacqf, XX1,SCALE3,OFFSET3,fv4,XL,YL)  

        nc_nuvacuncer1(1:XL,1:YL,1:2)= NUV_UncertaintyACAODToSSA(1:XL,1:YL,1:2) 
        do li=1,2
        nc_nuvacuncer1(1:XL,1:YL,li) = scale_offset(nc_nuvacuncer1(1:XL,1:YL,li),&
                                XX1,SCALE_10,OFFSET_10,fv4,XL,YL)  				
        enddo  

        nc_nuvacuncer2(1:XL,1:YL,1:2)= NUV_UncertaintyCODToSSA(1:XL,1:YL,1:2) 
        do li=1,2
        nc_nuvacuncer2(1:XL,1:YL,li) = scale_offset(nc_nuvacuncer2(1:XL,1:YL,li),&
                                XX1,SCALE_10,OFFSET_10,fv4,XL,YL)  				
        enddo  

         
!     Open netCDF file group geolocation_data

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
      
      call write_nc_2d(nc_ls_flag,'Land_Sea_Flag',&
           XX1,YY1,YY2,grpid,1,A_Num_part) 
      
      call write_nc_2d(nc_sm_weight,'Optical_Depth_Ratio_Small_Ocean_used',&
           XX1,YY1,YY2,grpid,1,A_Num_part)  
             
      call write_nc_3d(nc_aod_ocean_land,'Aerosol_Optical_Depth',&
           XX1,YY1,YY2,ZL4,grpid,1,A_Num_part) 
            
      call write_nc_3d(nc_ref_allwav,'Mean_Gas_Corrected_Reflectance',&
           XX1,YY1,YY2,ZL4-2,grpid,1,A_Num_part) 
      call write_nc_3d(nc_ref_allwav_uv,'Mean_Reflectance',& 
           XX1,YY1,YY2,ZL2,grpid,1,A_Num_part) 
                  
            
      call write_nc_3d(nc_W0_Ocean_UV,'DT_AerosolSingleScattAlbedo',&
            XX1,YY1,YY2,ZL1,grpid,1,A_Num_part)  

            
! Diagonostic for DT>>>>>>>>>          

      call write_nc_2d(nc_Height_indx,'DT_AerosolLayerHeight',&
            XX1,YY1,YY2,grpid,1,A_Num_part) 
              
      call write_nc_2d(nc_Cloud_Frac_LandOcean,'Aerosol_Cld_Fraction_Land_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part) 
           
      call write_nc_2d(nc_Quality_LandOcean ,'Quality_flag_Aerosol_Optical_Depth',&
           XX1,YY1,YY2,grpid,1,A_Num_part)    
         
      call write_nc_2d(nc_Quality_LandOcean_W0,'Quality_flag_SingleScattAlbedo',&
           XX1,YY1,YY2,grpid,1,A_Num_part)                   


!  NUV variables    

     
      call write_nc_2d(nc_uvai,'NUV_AerosolIndex',&
           XX1,YY1,YY2,grpid,1,A_Num_part)  
            
      call write_nc_2d(nc_residue,'NUV_Residue',&
           XX1,YY1,YY2,grpid,1,A_Num_part)  

      call write_nc_3d(nc_reflectivity,'NUV_Reflectivity',&
            XX1,YY1,YY2,ZL2,grpid,1,A_Num_part)  

      call write_nc_2d(nc_cod,'NUV_CloudOpticalDepth',&
           XX1,YY1,YY2,grpid,1,A_Num_part) 
            
      call write_nc_2d(nc_cldfrac,'NUV_RadiativeCloudFraction',&
          XX1,YY1,YY2,grpid,1,A_Num_part) 
             
      call write_nc_2d(nc_nuvalh,'NUV_AerosolLayerHeight',&
          XX1,YY1,YY2,grpid,1,A_Num_part)                   
     
      call write_nc_3d(nc_nuvssa,'NUV_AerosolSingleScattAlbedo',&
            XX1,YY1,YY2,ZL3,grpid,1,A_Num_part)  


! NUV  Above-cloud variables

      call write_nc_2d(nc_nuvaercorrcod,'NUV_AerosolCorrCloudOpticalDepth',&
           XX1,YY1,YY2,grpid,1,A_Num_part)  
	   	   
      call write_nc_3d(nc_nuvacaod,'NUV_AerosolOpticalDepthOverCloud',&
            XX1,YY1,YY2,ZL1,grpid,1,A_Num_part)  
	    
      call write_nc_3d(nc_nuvaercorrcodvshgt,'NUV_AerosolCorrCloudOpticalDepthVsHeight',&
           XX1,YY1,YY2,ZL3,grpid,1,A_Num_part)  
	   	   
      call write_nc_4d(nc_nuvacaodvshgt,'NUV_AerosolOpticalDepthOverCloudVsHeight',&
            XX1,YY1,YY2,ZL1,ZL3,grpid,1,A_Num_part)  

      call write_nc_3d(nc_nuvacuncer1,'NUV_UncertaintyACAODToSSA',&
            XX1,YY1,YY2,ZL2,grpid,0,A_Num_part)  
	    
      call write_nc_3d(nc_nuvacuncer2,'NUV_UncertaintyCODToSSA',&
            XX1,YY1,YY2,ZL2,grpid,0,A_Num_part)  

      call write_nc_2d(nc_nuvacqf,'NUV_FinalAlgorithmFlagsACA',&
           XX1,YY1,YY2,grpid,1,A_Num_part) 
	                
!  Close the netCDF file
if( ipart .eq.1)call check( nf90_close(ncid) )
print *, 'Wrote variables to ', out_file

return

end subroutine write_Output_merged

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

      
end module write_Pace_merged

 



