!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:  none
!
! !OUTPUT PARAMETERS:  none
!
! !REVISION HISTORY:
!
!
! !TEAM-UNIQUE HEADER:
!
! !REFERENCES AND CREDITS
!
! !DESIGN NOTES:
!
!
!   Externals:
!
!
! !END
!-----------------------------------------------------------------------
program l2_merge
  use netcdf
  
  use viirs_config, only:                              &
                            viirs_merge_type,          &
                            load_merge_config
  use viirs_db_utils, only:   create_viirs_l2  

  type (viirs_merge_type)        ::  config
  character(len=255)              ::  config_file 
  integer                         ::  args, sd
  character(len=255)              ::  usage
  character(len=5)                ::  platform  
  character(len=255)              ::  l2_file, output_file,segf
  character(len=255),dimension(10)::  seg_file
  integer                         ::  status
  logical ::  hires
  integer                               ::  nc_id,nseg,addf,st_inx,fi_inx,st_inx2,fi_inx2
  integer                               ::  dim_id, scan, xscan
  integer                               ::  dset_id    
  integer                               ::  grp_id  
  character(len=255)                    ::  dset_name ,group_name
  
  real(kind=8), dimension(:), allocatable   ::  time_avg          
   
   
!   real, dimension(:,:), allocatable   ::  lat, lon, sza_sav, vza_sav, raa_sav, &
!   & sca_sav, aot550_avg, aot550_best, aot412_avg, aot488_avg, aot670_avg, &
!   & sd550, ae_avg, ae_best, ssa_avg, oaot, oaot550, oaot550_best, oae, oae_best, ofmf,&
!   & ofmf_best, oaot550_stdv, oerror, oss, caot550_avg, caot550_best, cae_avg, cae_best, &
!   & ws_avg, oz_avg, wv_avg, wd_avg, ndvi_avg, rcndvi_avg, ndvi, sr_avg, dstar_avg, &
!   & btd11_avg, model_results, turbid_res, elev_avg, chl_avg, rcond, elev_avg_land, &
!   & elev_avg_ocean, smoke_count,oreflc_mean,lreflc_mean    
!   
!   integer, dimension(:,:), allocatable   :: conf_flag, alg_flag, oalg_flag, onpixels, &
!     & omodel_flag, oconf_flag, n_total_pixels, aerosol_type, combined_type, alg_flag_old, &
!     & alg_flag_old2, naot550_avg
  real, dimension(:,:), allocatable       ::  lat,lon,sza,vza,raa,sca,aot550,aot550_best
  real, dimension(:,:), allocatable       ::  aot412,aot488,aot670,aot550_sd,ae,ae_best
  integer, dimension(:,:), allocatable    ::  naot550,qa_flag,alg, alg_old,alg_old2,oalg
  integer, dimension(:,:), allocatable    ::  onaot550,oqa_flag,omodel_flag, ltype, ctype
  integer, dimension(:,:), allocatable    ::  n_valid_pixels,smoke_count  
  real, dimension(:,:), allocatable       ::  oaot550_best,oaot550,oae, oae_best,ofmf, ofmf_best
  real, dimension(:,:), allocatable       ::  osd550, oss,caot550,caot550_best,cae,cae_best
  real, dimension(:,:), allocatable       ::  ws,oz,wv ,wd,ndvi_avg,rcndvi_avg
  real, dimension(:,:), allocatable       ::  elev_avg,elev_avg_land,elev_avg_ocean,chl_avg        
  real, dimension(:,:), allocatable       ::  dstar_avg,btd11_avg,turbid_res ,ndvi ,rcond, dum
  real, dimension(:,:,:), allocatable     ::  oaot,oreflc_mean, lreflc_mean,oerr,sr_avg
  real, dimension(:,:,:), allocatable     :: ssa,dum3
                     
! -- config        
  usage = "./deep_blue <CONFIG_FILE>"
  args = command_argument_count()
  if ( args < 1 ) then
   write(*,*) usage
   stop
  end if

  call get_command_argument(1, value=config_file, status=status)
  if (status /= 0) then
    print *, "ERROR: Failed to get config file from command arguments: ", status
    stop
  end if
  
  config = load_merge_config(config_file, status)
  if (status /= 0) then
    print *, "ERROR: Failed to read in VIIRS configuration file: ", status
    stop
  end if
  platform  =  config%platform    !'VIIRS' or 'AHI' or 'GOES' 
  output_file = config%output_l2
!   l2_file = config%l2_file
  seg_file(1) = config%seg01_file
  seg_file(2) = config%seg02_file
  seg_file(3) = config%seg03_file
  seg_file(4) = config%seg04_file
  seg_file(5) = config%seg05_file
  seg_file(6) = config%seg06_file
  seg_file(7) = config%seg07_file
  seg_file(8) = config%seg08_file
  seg_file(9) = config%seg09_file
  seg_file(10) = config%seg10_file
  
  if (platform .eq. 'AHI') sd  = 1375
  if (platform .eq. 'GOES') sd  = 1356
  
  
! -- setting & allocate
  hires  = .false.        
  
 allocate(time_avg(sd), stat=status)
 if (status /= 0) then 
  	print *, "ERROR: Failed to allocate output dataset array: ", status
  	stop
 end if  

 allocate(lat(sd,sd), lon(sd,sd), sza(sd,sd), vza(sd,sd), raa(sd,sd), &
&         sca(sd,sd), aot550(sd,sd), aot550_best(sd,sd), naot550(sd,sd), aot412(sd,sd),&
&         aot488(sd,sd), aot670(sd,sd), aot550_sd(sd,sd), &
&         ae(sd,sd), ae_best(sd,sd), qa_flag(sd,sd), alg(sd,sd), oalg(sd,sd),&
&         oaot550(sd,sd), oaot550_best(sd,sd), &
&         oae(sd,sd), oae_best(sd,sd), ofmf(sd,sd), ofmf_best(sd,sd), osd550(sd,sd),&
&         onaot550(sd,sd), &
&         oss(sd,sd), oqa_flag(sd,sd), omodel_flag(sd,sd), caot550(sd,sd), &
&         caot550_best(sd,sd), cae(sd,sd), cae_best(sd,sd), n_valid_pixels(sd,sd), &
&         ws(sd,sd), oz(sd,sd), wv(sd,sd), wd(sd,sd), ndvi_avg(sd,sd), rcndvi_avg(sd,sd),&
&         ndvi(sd,sd), dstar_avg(sd,sd), btd11_avg(sd,sd), &
&         turbid_res(sd,sd), elev_avg(sd,sd), chl_avg(sd,sd), &
&         ltype(sd,sd), ctype(sd,sd),elev_avg_land(sd,sd),&
&         elev_avg_ocean(sd,sd), smoke_count(sd,sd),&
&         alg_old(sd,sd),alg_old2(sd,sd),rcond(sd,sd), stat=status)
 if (status /= 0) then 
  	print *, "ERROR: Failed to allocate output dataset array: ", status
  	stop
 end if  

 allocate(oaot(sd,sd,7),oreflc_mean(8,sd,sd), &
     &lreflc_mean(8,sd,sd), sr_avg(sd,sd,3), ssa(sd,sd,3), stat=status)
  if (status /= 0) then 
  	print *, "ERROR: Failed to allocate output dataset array: ", status
  	stop
 end if  
!   allocate(ssa(sd,sd,3), stat=status)
!   if (status /= 0) then 
!   	print *, "ERROR: Failed to allocate output dataset array: ", status
!   	stop
!  end if   
 
lat(:,:)=-999.
lon(:,:)=-999. 
sza(:,:)=-999. 
vza(:,:)=-999. 
raa(:,:)=-999.
sca(:,:)=-999. 
aot550(:,:)=-999.
aot550_best(:,:)=-999. 
naot550(:,:)=-999. 
aot412(:,:)=-999.
aot488(:,:)=-999.
aot670(:,:)=-999.
aot550_sd(:,:)=-999.
ae(:,:)=-999. 
ae_best(:,:)=-999. 
qa_flag(:,:)=-999. 
alg(:,:)=-999.
oalg(:,:)=-999.
oaot550(:,:)=-999. 
oaot550_best(:,:)=-999.
oae(:,:)=-999. 
oae_best(:,:)=-999. 
ofmf(:,:)=-999.
ofmf_best(:,:)=-999. 
osd550(:,:)=-999.
onaot550(:,:)=-999.
oss(:,:)=-999. 
oqa_flag(:,:)=-999.
omodel_flag(:,:)=-999.
caot550(:,:)=-999. 
caot550_best(:,:)=-999.
cae(:,:)=-999. 
cae_best(:,:)=-999. 
n_valid_pixels(:,:)=-999.
ws(:,:)=-999. 
oz(:,:)=-999.
wv(:,:)=-999. 
wd(:,:)=-999. 
ndvi_avg(:,:)=-999. 
rcndvi_avg(:,:)=-999.        
ndvi(:,:)=-999.
dstar_avg(:,:)=-999.
btd11_avg(:,:)=-999.
turbid_res(:,:)=-999.
elev_avg(:,:)=-999.
chl_avg(:,:)=-999. 
ltype(:,:)=-999.
ctype(:,:)=-999.
elev_avg_land(:,:)=-999.
elev_avg_ocean(:,:)=-999.
smoke_count(:,:)=-999.
alg_old(:,:)=-999.
alg_old2(:,:)=-999.
rcond(:,:)=-999.

oaot(:,:,:)=-999.
oreflc_mean(:,:,:)=-999.
lreflc_mean(:,:,:)=-999.
sr_avg(:,:,:)=-999.
ssa(:,:,:)=-999.0


! read segment files
  do nseg = 1, 10
    segf=seg_file(nseg)
    if (segf .eq. '-1') cycle 
    status = nf90_open(segf, nf90_nowrite, nc_id)
    call check(status,"ERROR: Failed to open geolocation file: ") 

    !define indx
     if (platform .eq. 'AHI') then
       addf    = nseg/2
       st_inx2  = 138*(nseg-1)-addf+1
       fi_inx2  = 138*(nseg)-addf   
       st_inx   = 1    
       fi_inx   = fi_inx2 -st_inx2+1   
     endif !'AHI'  
     if (platform .eq. 'GOES') then
      st_inx  = 2
      fi_inx  = 137 
      if (nseg .eq. 1) then
        st_inx  = 1
        fi_inx  = 136
      endif
      if (nseg .eq. 10) then
        fi_inx  = 133    
      endif
      
      st_inx2 = 136*((nseg)-1)+1
      fi_inx2 = (136*(nseg))
      if (nseg .eq. 10) then 
        st_inx2 = 1225     
        fi_inx2 = 1356  
      endif 
     endif !'GOES'     
      
     dset_name = 'Idx_Atrack'
     status = nf90_inq_dimid(nc_id, dset_name, dim_id)
     call check(status,"ERROR: Failed to get ID of dimension "//trim(dset_name)//": ")
     status = nf90_inquire_dimension(nc_id, dim_id, len=scan)
     call check(status,"ERROR: Failed to get size of dimension "//trim(dset_name)//": ")

     dset_name = 'Idx_Xtrack'
     status = nf90_inq_dimid(nc_id, dset_name, dim_id)
     call check(status,"ERROR: Failed to get ID of dimension "//trim(dset_name)//": ")
     status = nf90_inquire_dimension(nc_id, dim_id, len=xscan)
     call check(status,"ERROR: Failed to get size of dimension "//trim(dset_name)//": ")  
     
     allocate(dum(xscan, scan), stat=status)
     dum(:,:)=-999.


     dset_name = 'Longitude'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     lon(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Latitude'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     lat(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Aerosol_Optical_Thickness_550_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     aot550(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Aerosol_Optical_Thickness_550_Land_Best_Estimate'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     aot550_best(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Aerosol_Optical_Thickness_550_Land_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     caot550(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     caot550_best(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Aerosol_Optical_Thickness_550_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oaot550(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)     

     dset_name = 'Aerosol_Optical_Thickness_550_Ocean_Best_Estimate'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oaot550_best(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)      

     dset_name = 'Aerosol_Optical_Thickness_550_STDV_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     aot550_sd(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Aerosol_Optical_Thickness_550_STDV_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     osd550(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)     

     dset_name = 'Aerosol_Optical_Thickness_QA_Flag_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     qa_flag(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Aerosol_Optical_Thickness_QA_Flag_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oqa_flag(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Aerosol_Type_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     ltype(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)
     
     dset_name = 'Aerosol_Type_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     omodel_flag(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Aerosol_Type_Land_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     ctype(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Algorithm_Flag_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     alg_old(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Algorithm_Flag_Land_for_test'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     alg(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Algorithm_Flag_Land_for_test2'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     alg_old2(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Algorithm_Flag_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oalg(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)    

     dset_name = 'Angstrom_Exponent_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     ae(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)    

     dset_name = 'Angstrom_Exponent_Land_Best_Estimate'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     ae_best(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)    

     dset_name = 'Angstrom_Exponent_Land_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     cae(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)    

     dset_name = 'Angstrom_Exponent_Land_Ocean_Best_Estimate'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     cae_best(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)    

     dset_name = 'Angstrom_Exponent_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oae(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)        

     dset_name = 'Angstrom_Exponent_Ocean_Best_Estimate'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oae_best(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)        

     dset_name = 'Cell_Average_Elevation_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     elev_avg_land(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)     

     dset_name = 'Cell_Average_Elevation_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     elev_avg_ocean(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)       

     dset_name = 'Cell_Average_Elevation_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     elev_avg(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)       

     dset_name = 'Fine_Mode_Fraction_550_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     ofmf(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)    

     dset_name = 'Fine_Mode_Fraction_550_Ocean_Best_Estimate'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     ofmf_best(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)    

     dset_name = 'Number_Of_Pixels_Used_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     naot550(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Number_Of_Pixels_Used_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     onaot550(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)       

     dset_name = 'Number_Valid_Pixels'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     n_valid_pixels(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)     

     dset_name = 'Ocean_Sum_Squares'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oss(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)      

     dset_name = 'Precipitable_Water'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     wv(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)               

     dset_name = 'Relative_Azimuth_Angle'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     raa(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Scan_Start_Time'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     time_avg(st_inx2:fi_inx2)=dum(1,st_inx:fi_inx)

     dset_name = 'Scattering_Angle'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     sca(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Solar_Zenith_Angle'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     sza(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'TOA_NDVI'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     ndvi_avg(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)  

     dset_name = 'Total_Column_Ozone'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oz(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)       

     dset_name = 'Viewing_Zenith_Angle'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     vza(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)

     dset_name = 'Wind_Direction'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     wd(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)  

     dset_name = 'Wind_Speed'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     ws(:,st_inx2:fi_inx2)=dum(:,st_inx:fi_inx)  
  
     allocate(dum3(xscan, scan,3), stat=status)
     dum3(:,:,:)=-999.               
     dset_name = 'Spectral_Single_Scattering_Albedo_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum3)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     ssa(:,st_inx2:fi_inx2,:)=dum3(:,st_inx:fi_inx,:)  

     dset_name = 'Spectral_Surface_Reflectance'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum3)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     sr_avg(:,st_inx2:fi_inx2,:)=dum3(:,st_inx:fi_inx,:)              
     deallocate(dum3, stat=status)

     allocate(dum3(xscan, scan,7), stat=status)
     dum3(:,:,:)=-999.               
     dset_name = 'Spectral_Aerosol_Optical_Thickness_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum3)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oaot(:,st_inx2:fi_inx2,:)=dum3(:,st_inx:fi_inx,:)          
     deallocate(dum3, stat=status)

     allocate(dum3(8,xscan, scan), stat=status)
     dum3(:,:,:)=-999.               
     dset_name = 'Spectral_TOA_Reflectance_Land'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum3)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     lreflc_mean(:,:,st_inx2:fi_inx2)=dum3(:,:,st_inx:fi_inx)  
     
     dset_name = 'Spectral_TOA_Reflectance_Ocean'
     status = nf90_inq_varid(nc_id, dset_name, dset_id)
     call check(status,"ERROR: Failed to get ID of dataset "//trim(dset_name)//": ")   
     status = nf90_get_var(nc_id, dset_id, dum3)
     call check(status,"ERROR: Failed to read dataset "//trim(dset_name)//": ")  
     oreflc_mean(:,:,st_inx2:fi_inx2)=dum3(:,:,st_inx:fi_inx) 
     deallocate(dum3, stat=status)           
     
     !combine land_ocean product     
!      where (oaot550 >= 0)
!        caot550 = oaot550
!      end where          
!      where (aot550 >= 0)
!        caot550 = aot550
!      end where  
               
    deallocate(dum, stat=status)  
    status = nf90_close(nc_id)
    call check(status,"ERROR: Failed to close M-band L1B file: ")       
      
  enddo

  
  status = create_viirs_l2(output_file, time_avg, lat, lon, sza, vza, raa, &
&         sca, aot550, aot550_best, naot550, aot412, aot488, aot670, aot550_sd, &
&         ae, ae_best, ssa, qa_flag, alg, oalg, oaot,oaot550, oaot550_best, &
&         oae, oae_best, ofmf, ofmf_best, osd550, onaot550, &
&         oss, oqa_flag, omodel_flag, caot550, caot550_best, cae, cae_best, n_valid_pixels, &
&         ws, oz, wv, wd, ndvi_avg, rcndvi_avg, ndvi, sr_avg, dstar_avg, btd11_avg, &
&         turbid_res, elev_avg, chl_avg, ltype, ctype,elev_avg_land,&
&         elev_avg_ocean, smoke_count,oreflc_mean, lreflc_mean,hires,alg_old,alg_old2,platform)
  
                              
end program l2_merge

subroutine check(status, str)
  integer, intent ( in) :: status
  character(len=*)    :: str
  
  if(status /= NF90_NOERR) then 
    print *, str, status
  end if
end subroutine check
