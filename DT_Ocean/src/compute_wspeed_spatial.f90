 module compute_wspeed_spatial
 
      implicit none

contains   
! ---------------------------------------------------------------------
!     
!!SUBROUTINE interpol_wspeed_spatial
!
! This routine reads in the wind speeds of each grid cell surrounding the 
! current grid cell, and bilinearly interpolates wind speed between grid cell.
!
! This routine is intended to eliminate the problem of having discontinuities
!in AOD at the edge of a 1x1 lat/lon grid cell due to sharp changes in the 
! wind speed and therefore wind speed correction.
!
! Bilinear interpolation copied from Numerical Recipes in FORTRAN 77 
! Formula 3.6.5 
!
! Introduced by L. Munchak on March 27, 2013. leigh.a.munchak@nasa.gov
!
!================================================================================

       SUBROUTINE interpol_wspeed_spatial(set_counter_for_anc,&
        Lat_center,Lon_center,ugrd,vgrd)  
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

!     Arguments passed from main
      REAL ugrd,vgrd,lat_center, lon_center       

!     Dummy variables for met I don't want to read
      REAL temp_null, pwat_null, oz_null
!     FLOOR of lat_center & lon_center
      REAL fl_lat, fl_lon                       
      integer R_temp,set_counter_for_anc        

!     CEILING of lat_center & lon_center
      REAL ci_lon, ci_lat                        

!     u-wind component for upper lefthand corner, lower-right, etc.
      REAL ul_u, ur_u, ll_u, lr_u                
 
!     v-wind component for upper lefthand corner, lower-right, etc.
      REAL ul_v, ur_v, ll_v, lr_v

!     Weighting for bilnear interpolation
      REAL bi_t, bi_u         
     
!    Get integer part of lat/lon
      fl_lat = INT(lat_center)
      fl_lon = INT(lon_center) 

      IF(fl_lat .LT. 0 ) THEN
         fl_lat = fl_lat -1
      ENDIF

      IF(fl_lon .LT. 0) THEN
         fl_lon = fl_lon - 1
      ENDIF

!     Get ceiling of lat/lon
      ci_lat = fl_lat + 1
      ci_lon = fl_lon + 1      

!     Read in wind speeds at each corner, put temp/oz/pwat into dummy variables

        set_counter_for_anc =set_counter_for_anc+1
     
      CALL Get_An_gdas(ci_lat,ci_lon,ur_u,ur_v,pwat_null,oz_null,set_counter_for_anc,R_temp)
      CALL Get_An_gdas(ci_lat,fl_lon,ul_u,ul_v,pwat_null,oz_null,set_counter_for_anc,R_temp)
      CALL Get_An_gdas(fl_lat,ci_lon,lr_u,lr_v,pwat_null,oz_null,set_counter_for_anc,R_temp)
      CALL Get_An_gdas(fl_lat,fl_lon,ll_u,ll_v,pwat_null,oz_null,set_counter_for_anc,R_temp)
!     Compute weighting factors for interpolation
      bi_t  = (lon_center-fl_lon)/(ci_lon-fl_lon)    
      bi_u  = (lat_center-fl_lat)/(ci_lat-fl_lat)   

!     Perform interpolation on u & v
      ugrd  = ((1-bi_t)*(1-bi_u)*ll_u)+(bi_t*(1-bi_u)*lr_u)&
       +(bi_t*bi_u*ur_u)+((1-bi_t)*bi_u*ul_u)
      vgrd  = ((1-bi_t)*(1-bi_u)*ll_v)+(bi_t*(1-bi_u)*lr_v)&
       +(bi_t*bi_u*ur_v)+((1-bi_t)*bi_u*ul_v) 

      RETURN
      end  subroutine interpol_wspeed_spatial
           end module  compute_wspeed_spatial
           
           
           
