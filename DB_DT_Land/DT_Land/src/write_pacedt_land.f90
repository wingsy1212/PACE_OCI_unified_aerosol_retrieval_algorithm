! write_pace_dt

! Draft netCDF writer for VIIRS Dark Target
! Simplified to fewer variables for testing

! Virginia Sawyer
! Created 2017-03-10
! Updated  for PACE Shana Mattoo

module write_pace_dt

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

      subroutine make_nc_Land(Save_Lat,Save_Lon,save_SolZen,save_View_angle,&
                 save_View_phi,save_solar_phi,Save_Scatt_Angle,&
                 Save_Glint_angle,save_Sea_Land_Flag,Save_Cld_mask_500,&
                 Save_cloud_distance_500,Save_land_ocean_Quality,&
                 Save_land_ocean_tau,Save_Image_land_ocean_tau,& 
                  Save_Surface_Reflectance_Land,Save_corrected,&
                  Save_dust_weighting,Save_NUMPIXELS_land,Save_ref_Land,&
                  Save_SDS_Land,Save_SDS_mass_conc_land,&
                  Save_CLDFRC_land,Save_average_Tau_Ocean,&
                  Save_SDSTAUS_average,Save_SDSTAUL_average,&
                  Save_SDS_mass_conc_ocean,Save_CLDFRC_Ocean,&
                  Save_effrad,Save_PSML003_Ocean,Save_ASSY,Save_BACK,&
                  Save_SDS_angs_coeff1,Save_SDS_angs_coeff2,&
                  Save_Least_error,Save_small_weighting,&
                  Save_SDS_AOT_model,Save_NUMPIXELS_Ocean,&
                  Save_ref_Ocean,Save_STD_Ocean,Save_Wspeed_Ocean,&
                  Save_Topo_Height,Save_land_ocean_Q_Ret,&
                  Save_Aerosol_Type,Save_Fitting_Error_Land,Save_ave_clddis_land_ocean,&
                  Ret_Xtrack,Ret_Lines,ipart,num_joint)
                   
                   
!     Define input variable dimensions, etc.
!     Define parameters from Main_Driver.f90

      include 'read_Sat_MODIS.inc'
      include 'mod04.inc'

!     Declare input variables

      include 'Save_data_declare.inc'

!     Declare output variables
       character(len=10) :: Sat_Flag  
       integer yy1,yy2,ipart,xx1,yy_part1,yy_part2,num_joint
       integer yy3,yy4,number_of_parts         
      Real nc_lat(Ret_Xtrack,Ret_Lines),nc_lon(Ret_Xtrack,Ret_Lines)
      real nc_solzen(Ret_Xtrack,Ret_Lines),nc_senzen(Ret_Xtrack,Ret_Lines)
      real nc_solazi(Ret_Xtrack,Ret_Lines),nc_senazi(Ret_Xtrack,Ret_Lines)
      real nc_scatt(Ret_Xtrack,Ret_Lines),nc_glint(Ret_Xtrack,Ret_Lines)
      real nc_ls_flag(Ret_Xtrack,Ret_Lines)
      real nc_cldmask(Ret_Xtrack*ILINE,Ret_Lines*ILINE)
      real nc_clddist(Ret_Xtrack*ILINE,Ret_Lines*ILINE)
      real nc_ls_qual(Ret_Xtrack,Ret_Lines),nc_ls_aod(Ret_Xtrack,Ret_Lines)
      real nc_img_ls_aod(Ret_Xtrack,Ret_Lines),nc_aero_type(Ret_Xtrack,Ret_Lines)
      real nc_error_land(Ret_Xtrack,Ret_Lines),nc_sfc_land(Ret_Xtrack,Ret_Lines,Land_Sol3)
      real nc_aod_land(Ret_Xtrack,Ret_Lines,Land_Sol4),nc_aodrat_sm_land(Ret_Xtrack,Ret_Lines)
      real nc_numpix_land(Ret_Xtrack,Ret_Lines,Land_Sol3)
      real nc_ref_land(Ret_Xtrack,Ret_Lines,Land_Sol3)
      real nc_refsd_land(Ret_Xtrack,Ret_Lines,Land_Sol3)
      real nc_massconc_land(Ret_Xtrack,Ret_Lines),nc_cldfrc_land(Ret_Xtrack,Ret_Lines)
      real nc_aod_ocean(Ret_Xtrack,Ret_Lines,NWAVN)
      real nc_aod_sm_ocean(Ret_Xtrack,Ret_Lines,NWAVN)
      real nc_aod_lg_ocean(Ret_Xtrack,Ret_Lines,NWAVN)
      real nc_massconc_ocean(Ret_Xtrack,Ret_Lines)
      real nc_cldfrc_ocean(Ret_Xtrack,Ret_Lines)
      real nc_ave_clddis_land_ocean(Ret_Xtrack,Ret_Lines)
      real nc_effrad(Ret_Xtrack,Ret_Lines),nc_psml003_ocean(Ret_Xtrack,Ret_Lines)
      real nc_assy(Ret_Xtrack,Ret_Lines,NWAVN),nc_back(Ret_Xtrack,Ret_Lines,NWAVN)
      real nc_angstrom1(Ret_Xtrack,Ret_Lines),nc_angstrom2(Ret_Xtrack,Ret_Lines)
      real nc_least_err(Ret_Xtrack,Ret_Lines),nc_sm_weight(Ret_Xtrack,Ret_Lines)
      real nc_aod_model(Ret_Xtrack,Ret_Lines,num_model_index)
      real nc_numpix_ocean(Ret_Xtrack,Ret_Lines,NWAVN)
      real nc_ref_ocean(Ret_Xtrack,Ret_Lines,NWAVN)
      real nc_refsd_ocean(Ret_Xtrack,Ret_Lines,NWAVN)
      real nc_ws_ocean(Ret_Xtrack,Ret_Lines)
      real nc_topo_alt(Ret_Xtrack,Ret_Lines)
      real nc_ls_qret(Ret_Xtrack,Ret_Lines,Land_Sol1)
      real nc_solindx_lg(Ret_Xtrack,Ret_Lines)
      real nc_solindx_sm(Ret_Xtrack,Ret_Lines)
       integer  ixp, IYp,A_Num_part,iscan,idata,ik,ij
      real  new_array(Ret_Xtrack,Ret_Lines)
!     Declare variables for netCDF processing

      integer XL,YL,ZL,li,lj, YY1_new, YY2_new,ZL2,ZL3,ZL4,ZLS
      integer ncid,grpid
      character (len=*), parameter :: nc_name = 'vnpaerdt_output_Land.nc'
      real fv3,fv4

!     Open the netCDF file already created via ncgen
       
         ipart =1            
      call check( nf90_open(nc_name,NF90_WRITE, ncid) )

      XL = Ret_Xtrack
      YL = Ret_Lines
      ZL = NWAVN
      ZL2 = Land_Sol1
      ZL3 = Land_Sol3
      ZL4 = Land_Sol4
      ZLS = num_model_index
      fv3 = -990.
      fv4 = -9990. 
    
!     Fill output variables with data
!     Apply scale factors and offsets, excluding fill values

          
                  XX1   =    XL
                  A_Num_part = YL    
                  YY1=1
                  YY2=Ret_Lines
            
       
             print*,'starting write',trim(nc_name)
           
        
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
      nc_scatt(1:XL,1:YL) = Save_Scatt_angle(1:XL,1:YL)
      nc_scatt = scale_offset(nc_scatt, XX1,SCALE2,OFFSET2,fv4,XL,YL)
      nc_glint(1:XL,1:YL) = Save_Glint_angle(1:XL,1:YL)
      nc_glint = scale_offset(nc_glint, XX1,SCALE2,OFFSET2,fv4,XL,YL)
      nc_ls_flag(1:XL,1:YL) = save_Sea_Land_Flag(1:XL,1:YL) 
      nc_ls_flag = scale_offset(nc_ls_flag, XX1,SCALE1,OFFSET1,fv4,XL,YL)
      nc_ls_qual(1:XL,1:YL) = Save_land_ocean_Quality(1:XL,1:YL)
      nc_ls_qual = scale_offset(nc_ls_qual, XX1,SCALE1,OFFSET1,fv4,XL,YL)
      nc_ls_aod(1:XL,1:YL) = Save_land_ocean_tau(1:XL,1:YL)
      nc_ls_aod = scale_offset(nc_ls_aod, XX1,SCALE3,OFFSET3,fv4,XL,YL)
      nc_img_ls_aod(1:XL,1:YL) = Save_Image_land_ocean_tau(1:XL,1:YL)
      nc_img_ls_aod = scale_offset(nc_img_ls_aod, XX1,SCALE3,OFFSET3,fv4,XL,YL)
      nc_aero_type(1:XL,1:YL) = Save_Aerosol_Type(1:XL,1:YL) 
      nc_aero_type = scale_offset(nc_aero_type, XX1,SCALE1,OFFSET1,fv4,XL,YL)
      nc_error_land(1:XL,1:YL) = Save_Fitting_Error_Land(1:XL,1:YL)
      nc_error_land = scale_offset(nc_error_land,XX1,SCALE3,OFFSET3,fv4,XL,YL)  
      nc_sfc_land(1:XL,1:YL,:) = Save_Surface_Reflectance_Land(1:XL,1:YL,:)
       
 !       nc_ave_clddis_land_ocean(1:XL,1:YL) = Save_ave_clddis_land_ocean(1:XL,1:YL)
 !       nc_ave_clddis_land_ocean = scale_offset(nc_ave_clddis_land_ocean,XX1,SCALE1,OFFSET1,fv4,XL,YL) 
 !       nc_cldmask(1:XX1*Iline,1:YL*ILINE) = Save_Cld_mask_500(1:XL*ILINE,1:YL*ILINE)
 !       nc_cldmask = scale_offset(nc_cldmask,(XX1*ILINE),SCALE1,OFFSET1,fv4,(XL*ILINE),(YL*ILINE)) 
 !       nc_clddist(1:XX1*Iline,1:YL*ILINE) = Save_cloud_distance_500(1:XL*ILINE,1:YL*ILINE)
 !       nc_clddist = scale_offset(nc_clddist,(XX1*ILINE),SCALE1,OFFSET1,fv4,(XL*ILINE),(YL*ILINE))
        
! 
         
      YL = Ret_Lines
      ZL = NWAVN
      ZL2 = Land_Sol1
      ZL3 = Land_Sol3
      ZL4 = Land_Sol4
      
      do li=1,ZL3
         nc_sfc_land(1:XL,1:YL,li) = scale_offset(nc_sfc_land(1:XL,1:YL,li),&
         XX1,SCALE3,OFFSET3,fv4,XL,YL)
      enddo 
      nc_aod_land(1:XL,1:YL,:) = save_corrected(1:XL,1:YL,:)
      do li=1,ZL4
         nc_aod_land(1:XL,1:YL,li) =scale_offset(nc_aod_land(1:XL,1:YL,li),&
         XX1,SCALE3,OFFSET3,fv4,XL,YL)
      enddo
       
      nc_aodrat_sm_land(1:XL,1:YL) = Save_dust_weighting(1:XL,1:YL)
      nc_aodrat_sm_land = scale_offset(nc_aodrat_sm_land, XX1,SCALE3,&
                          OFFSET3,fv4,XL,YL)
      nc_numpix_land(1:XL,1:YL,:) = save_NUMPIXELS_land(1:XL,1:YL,:)
      nc_ref_land(1:XL,1:YL,:) = save_ref_Land(1:XL,1:YL,:)
      nc_refsd_land(1:XL,1:YL,:) = save_SDS_Land(1:XL,1:YL,:)
      do li=1,ZL3
         nc_numpix_land(1:XL,1:YL,li) = scale_offset(nc_numpix_land(1:XL,1:YL,li),&
                                XX1,SCALE1,OFFSET1,fv4,XL,YL)
         nc_ref_land(1:XL,1:YL,li) = scale_offset(nc_ref_land(1:XL,1:YL,li),&
                              XX1,SCALE4,OFFSET4,fv4,XL,YL)
         nc_refsd_land(1:XL,1:YL,li) = scale_offset(nc_refsd_land(1:XL,1:YL,li),&
                               XX1,SCALE4,OFFSET4,fv4,XL,YL)
      enddo
       
      nc_massconc_land(1:XL,1:YL) = Save_SDS_mass_conc_land(1:XL,1:YL)
      nc_massconc_land = scale_offset(nc_massconc_land, XX1,SCALE1,OFFSET1,fv3,XL,YL)
      nc_cldfrc_land(1:XL,1:YL) = Save_CLDFRC_land(1:XL,1:YL)
      nc_cldfrc_land = scale_offset(nc_cldfrc_land, XX1,SCALE3,OFFSET3,fv4,XL,YL)
      nc_aod_ocean(1:XL,1:YL,:)    = Save_average_Tau_Ocean(1:XL,1:YL,:) 
      nc_aod_sm_ocean(1:XL,1:YL,:) = Save_SDSTAUS_average(1:XL,1:YL,:)
      nc_aod_lg_ocean(1:XL,1:YL,:) = Save_SDSTAUL_average(1:XL,1:YL,:)
       
                                 
      do li=1,ZL
         nc_aod_ocean(1:XL,1:YL,li) = scale_offset(nc_aod_ocean(1:XL,1:YL,li),&
                                XX1,SCALE3,OFFSET3,fv4,XL,YL) 
         nc_aod_sm_ocean(1:XL,1:YL,li) = scale_offset(nc_aod_sm_ocean(1:XL,1:YL,li),&
                                XX1,SCALE3,OFFSET3,fv4,XL,YL)
         nc_aod_lg_ocean(1:XL,1:YL,li) = scale_offset(nc_aod_lg_ocean(1:XL,1:YL,li),&
                                XX1,SCALE3,OFFSET3,fv4,XL,YL)
      enddo

      nc_massconc_ocean(1:XL,1:YL) = Save_SDS_mass_conc_ocean(1:XL,1:YL)
      nc_massconc_ocean = scale_offset(nc_massconc_ocean, XX1,SCALE1,OFFSET1,fv3,XL,YL)
      nc_cldfrc_ocean(1:XL,1:YL) = Save_CLDFRC_Ocean(1:XL,1:YL)
      nc_cldfrc_ocean = scale_offset(nc_cldfrc_ocean, XX1,SCALE3,OFFSET3,fv4,XL,YL)
      nc_effrad(1:XL,1:YL) = Save_effrad(1:XL,1:YL)
      nc_effrad = scale_offset(nc_effrad, XX1,SCALE3,OFFSET3,fv4,XL,YL)
      nc_psml003_ocean(1:XL,1:YL) = Save_PSML003_Ocean(1:XL,1:YL)
      nc_psml003_ocean = scale_offset(nc_psml003_ocean, XX1,SCALE1,OFFSET1,fv3,XL,YL)

      nc_assy(1:XL,1:YL,:) = Save_ASSY(1:XL,1:YL,:)
      nc_back(1:XL,1:YL,:) = Save_BACK(1:XL,1:YL,:)
      do li=1,ZL
         nc_assy(1:XL,1:YL,li) = scale_offset(nc_assy(1:XL,1:YL,li),&
         XX1,SCALE3,OFFSET3,fv4,XL,YL)
         nc_back(1:XL,1:YL,li) = scale_offset(nc_back(1:XL,1:YL,li),&
         XX1,SCALE3,OFFSET3,fv4,XL,YL)
      enddo
             
      nc_angstrom1(1:XL,1:YL) = Save_SDS_angs_coeff1(1:XL,1:YL)
      nc_angstrom1 =scale_offset(nc_angstrom1, XX1,SCALE3,OFFSET3,fv4,XL,YL)
      nc_angstrom2(1:XL,1:YL) = Save_SDS_angs_coeff2(1:XL,1:YL)
      nc_angstrom2 = scale_offset(nc_angstrom2, XX1,SCALE3,OFFSET3,fv4,XL,YL)
      nc_least_err(1:XL,1:YL) = Save_Least_error(1:XL,1:YL)
      nc_least_err =scale_offset(nc_least_err, XX1,SCALE3,OFFSET3,fv4,XL,YL)
      nc_sm_weight(1:XL,1:YL) = Save_small_weighting(1:XL,1:YL)
      nc_sm_weight = scale_offset(nc_sm_weight, XX1,SCALE3,OFFSET3,fv4,XL,YL)
      nc_aod_model(1:XL,1:YL,:) = Save_SDS_AOT_model(1:XL,1:YL,:)
      do li=1,9
         nc_aod_model(1:XL,1:YL,li) = scale_offset(nc_aod_model(1:XL,1:YL,li),&
                      XX1,SCALE3,OFFSET3,fv4,XL,YL)
      enddo
      nc_numpix_ocean(1:XL,1:YL,:) = Save_NUMPIXELS_Ocean(1:XL,1:YL,:)
      nc_ref_ocean(1:XL,1:YL,:) = Save_ref_Ocean(1:XL,1:YL,:)
      nc_refsd_ocean(1:XL,1:YL,:) = Save_STD_Ocean(1:XL,1:YL,:)
      do li=1,ZL
         nc_numpix_ocean(1:XL,1:YL,li) = scale_offset(nc_numpix_ocean(1:XL,1:YL,li),&
                  XX1,SCALE1,OFFSET1,fv4,XL,YL)
         nc_ref_ocean(1:XL,1:YL,li) = scale_offset(nc_ref_ocean(1:XL,1:YL,li),&
                  XX1,SCALE4,OFFSET4,fv4,XL,YL)
         nc_refsd_ocean(1:XL,1:YL,li) =scale_offset(nc_refsd_ocean(1:XL,1:YL,li),&
                  XX1,SCALE4,OFFSET4,fv4,XL,YL)
      enddo

      nc_ws_ocean(1:XL,1:YL) = Save_Wspeed_Ocean(1:XL,1:YL)
      nc_ws_ocean = scale_offset(nc_ws_ocean, XX1,SCALE2,OFFSET2,fv4,XL,YL)
      nc_topo_alt(1:XL,1:YL) = Save_Topo_Height(1:XL,1:YL)
      nc_topo_alt = scale_offset(nc_topo_alt, XX1,SCALE2,OFFSET2,fv4,XL,YL)
      nc_ls_qret(1:XL,1:YL,:) = Save_land_ocean_Q_Ret(1:XL,1:YL,:)
      nc_ls_qret(1:XL,1:YL,1) =scale_offset(nc_ls_qret(1:XL,1:YL,1), XX1,SCALE1,OFFSET1,fv4,XL,YL)
      nc_ls_qret(1:XL,1:YL,2) =scale_offset(nc_ls_qret(1:XL,1:YL,2), XX1,SCALE1,OFFSET1,fv4,XL,YL)
      
            print*,'done1' 
              
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
      call write_nc_2d(nc_scatt,'Scattering_Angle',XX1,YY1,YY2,grpid,1,A_Num_part)
      
      call write_nc_2d(nc_glint,'Glint_Angle',XX1,YY1,YY2,grpid,1,A_Num_part) 
        
!     Open netCDF file group geophysical_data

       print*,'done2' 

      call check( nf90_inq_ncid(ncid, 'geophysical_data', grpid) )
         
!     Write 2D variables to netCDF group geophysical_data
           
      call write_nc_2d(nc_ls_flag,'Land_Sea_Flag',&
           XX1,YY1,YY2,grpid,1,A_Num_part) 
!      call write_nc_2d(nc_ave_clddis_land_ocean,'Average_Cloud_Pixel_Distance_Land_Ocean',&
!           XX1,YY1,YY2,grpid,1,A_Num_part)
! !       call write_nc_2d(nc_cldmask,'Aerosol_Cldmask_Land_Ocean',&
!             (XX1*ILINE),YY1,(YY2*ILINE),grpid,1,(A_Num_part*ILINE))  
!        call write_nc_2d(nc_clddist,'Cloud_Pixel_Distance_Land_Ocean',&
!             (XX1*ILINE),YY1,(YY2*ILINE),grpid,1,(A_Num_part*ILINE)) 

 
      call write_nc_2d(nc_ls_qual,'Land_Ocean_Quality_Flag',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_ls_aod,'Optical_Depth_Land_And_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_img_ls_aod,'Image_Optical_Depth_Land_And_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_aero_type,'Aerosol_Type_Land',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_error_land,'Fitting_Error_Land',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_aodrat_sm_land,'Optical_Depth_Ratio_Small_Land',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_massconc_land,'Mass_Concentration_Land',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_cldfrc_land,'Aerosol_Cloud_Fraction_Land',&
           XX1,YY1,YY2,grpid,1,A_Num_part) 
      call write_nc_2d(nc_massconc_ocean,'Mass_Concentration_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_cldfrc_ocean,'Aerosol_Cloud_Fraction_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part) 
      call write_nc_2d(nc_effrad,'Effective_Radius_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_psml003_ocean,'PSML003_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_angstrom1,'Angstrom_Exponent_1_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_angstrom2,'Angstrom_Exponent_2_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_least_err,'Least_Squares_Error_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_sm_weight,'Optical_Depth_Ratio_Small_Ocean_0p55micron',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_ws_ocean,'Wind_Speed_Ncep_Ocean',&
           XX1,YY1,YY2,grpid,1,A_Num_part)
      call write_nc_2d(nc_topo_alt,'Topographic_Altitude_Land',&
           XX1,YY1,YY2,grpid,1,A_Num_part)  
          call check( nf90_inq_ncid(ncid, 'geophysical_data', grpid) )
      call write_nc_3d(nc_sfc_land,'Surface_Reflectance_Land',&
           XX1,YY1,YY2,ZL3,grpid,A_Num_part)  
      call write_nc_3d(nc_aod_land,'Corrected_Optical_Depth_Land',&
           XX1,YY1,YY2,ZL4,grpid,A_Num_part)  
       call write_nc_3d(nc_numpix_land,'Number_Pixels_Used_Land',&
            XX1,YY1,YY2,ZL3,grpid,A_Num_part) 
       call write_nc_3d(nc_ref_land,'Mean_Reflectance_Land',&
            XX1,YY1,YY2,ZL3,grpid,A_Num_part)  
       call write_nc_3d(nc_refsd_land,'STD_Reflectance_Land',&
            XX1,YY1,YY2,ZL3,grpid,A_Num_part)  
      call write_nc_3d(nc_aod_ocean,'Effective_Optical_Depth_Average_Ocean',&
           XX1,YY1,YY2,ZL,grpid,A_Num_part)
            
      call write_nc_3d(nc_aod_sm_ocean,'Optical_Depth_Small_Average_Ocean',&
           XX1,YY1,YY2,ZL,grpid,A_Num_part)
            
      call write_nc_3d(nc_aod_lg_ocean,'Optical_Depth_Large_Average_Ocean',&
           XX1,YY1,YY2,ZL,grpid,A_Num_part)
           
      call write_nc_3d(nc_assy,'Asymmetry_Factor_Average_Ocean',&
           XX1,YY1,YY2,ZL,grpid,A_Num_part)
        
      
!     Close the netCDF file
3000   continue
       if( ipart .eq.1)call check( nf90_close(ncid) )
      print *, 'Wrote variables to ',nc_name

      return
      end subroutine make_nc_Land

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
      
end module write_pace_dt

