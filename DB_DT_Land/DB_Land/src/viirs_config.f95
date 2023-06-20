module viirs_config
!
!	viirs_config.f95
!
!	Description: 	  Configuration module for Deep Blue aerosol retrievals.
!
!	Modules:			  None
!
!	Author:			    Corey Bettenhausen
!						      Science Systems and Applications, Inc.
!						      NASA Goddard Space Flight Center
!						      corey_bettenhausen@ssaihq.com
!
!	Last Changed:	  
!
!	History: 		    2012-08	Original version
!
!------------------------------------------------------------------------------
	implicit none

! Restrict access to module components. 
  private 

  public  ::  load_viirs_config
  public  ::  viirs_config_type
  public  ::  load_merge_config
  public  ::  viirs_merge_type
  
  type  ::  viirs_config_type
    integer                 ::  year
    integer                 ::  month
    integer                 ::  day
    integer                 ::  hour
    integer                 ::  minute
    integer                 ::  second
    character(len=255)      ::  aerosol_land_file
    character(len=255)      ::  aerosol_dust_file
    character(len=255)      ::  aerosol_ocean_dust_file
    character(len=255)      ::  aerosol_ocean_fine_file
    character(len=255)      ::  aerosol_ocean_mari_file
    character(len=255)      ::  aerosol_ocean_mix_file
    character(len=255)      ::  bathymetry_lut_file
    character(len=255)      ::  chl_lut_file
    character(len=255)      ::  ler_lut_file
    character(len=255)      ::  landcover_file
    character(len=255)      ::  surfpressure_file
    character(len=255)      ::  geozone_file
    character(len=255)      ::  seasonal_deserts_file
    character(len=255)      ::  brdfbase_file
    character(len=255)      ::  modis_surfdb_file
    character(len=255)      ::  viirs_surfdb_file
    character(len=255)      ::  surfcoeffs_file
    character(len=255)      ::  swir_vis_surfcoeffs_file
    character(len=255)      ::  veg_landcover_file
    character(len=255)      ::  veg_sfc21_file
    
    character(len=255)      ::  rayl412_file
    character(len=255)      ::  rayl488_file
    character(len=255)      ::  rayl670_file
    character(len=255)      ::  xcal412_file
    character(len=255)      ::  xcal488_file
    character(len=255)      ::  xcal670_file
    
    character(len=255)      ::  gmtco_file
    character(len=255)      ::  iicmo_file
    character(len=255)      ::  svm01_file
    character(len=255)      ::  svm02_file
    character(len=255)      ::  svm03_file
    character(len=255)      ::  svm04_file
    character(len=255)      ::  svm05_file
    character(len=255)      ::  svm07_file
    character(len=255)      ::  svm08_file
    character(len=255)      ::  svm09_file
    character(len=255)      ::  svm10_file
    character(len=255)      ::  svm11_file
    character(len=255)      ::  svm14_file
    character(len=255)      ::  svm15_file
    character(len=255)      ::  svm16_file
	  character(len=255)      ::  gdas_file1
	  character(len=255)      ::  gdas_file2
    
    character(len=255)      ::  l1b_m
    character(len=255)      ::  geo_m
    character(len=255)      ::  platform
    character(len=255)      ::  segment
    character(len=255)      ::  goes_1
    character(len=255)      ::  goes_2
    character(len=255)      ::  goes_3
    character(len=255)      ::  goes_4
    character(len=255)      ::  goes_5
    character(len=255)      ::  goes_6
    character(len=255)      ::  goes_7
    character(len=255)      ::  goes_8
    character(len=255)      ::  goes_9
    character(len=255)      ::  goes_10
    character(len=255)      ::  goes_11
    character(len=255)      ::  goes_12
    character(len=255)      ::  goes_13
    character(len=255)      ::  goes_14
    character(len=255)      ::  goes_15
       
    character(len=255)      ::  output_l2
    character(len=255)      ::  dbdt_file
    
  end type viirs_config_type

  type  ::  viirs_merge_type
    
    character(len=255)      ::  l2_file
    character(len=255)      ::  platform
    character(len=255)      ::  output_l2
    character(len=255)      ::  seg01_file
    character(len=255)      ::  seg02_file
    character(len=255)      ::  seg03_file
    character(len=255)      ::  seg04_file
    character(len=255)      ::  seg05_file
    character(len=255)      ::  seg06_file
    character(len=255)      ::  seg07_file
    character(len=255)      ::  seg08_file
    character(len=255)      ::  seg09_file
    character(len=255)      ::  seg10_file
    
  end type viirs_merge_type
  
	contains


	type(viirs_merge_type) function load_merge_config(config_file, status) result(vcfg)
		implicit none
		
		character(len=255), intent(in)	:: config_file
		integer, intent(inout)          :: status
		
		character(255)					        :: c_line, var, val
		integer							            :: eq_index
		
!   -- initialize our config variables.
    vcfg%l2_file                  = ''
  	vcfg%output_l2                = ''
  	vcfg%platform                 = ''        
		vcfg%seg01_file               = ''
		vcfg%seg02_file               = ''
		vcfg%seg03_file               = ''
		vcfg%seg04_file               = ''
		vcfg%seg05_file               = ''
		vcfg%seg07_file               = ''
		vcfg%seg08_file               = ''
		vcfg%seg09_file               = ''
		vcfg%seg10_file               = ''

		
		open(100, file=config_file, status='OLD', FORM='FORMATTED',  ACTION='READ', iostat=status)
    if (status /= 0) then
      print *, 'ERROR: failed to open configuration file: ', status
      return
    end if 
	  	  
		do 
		! Read a line from the config file.  Die on error.
			read(100, fmt='(A)', iostat=status) c_line
			if (status /= 0) then
				if (status < 0) then	
				  status = 0	
					return
				else if (status > 0) then
					print *, "ERROR: Failed to read configuration file: ", status
					status = -1
					return
				end if
			end if
						
		! Clean up whitespace.
			c_line =rm_whitespace(c_line)
			
		! Detect comment lines or empty lines and skip.
			if (index(c_line, '!') /= 0) then
				cycle
			end if

			if (len_trim(c_line) == 0) then
				cycle
			end if
			
		! Begin parsing for variable and value.  
		! If no '=' is detected, we have a funky line, cycle.
			eq_index = index(c_line,'=')
			if (eq_index /= 0) then
				var = c_line(1:eq_index-1)
				val = c_line(eq_index+1:len(c_line))
				
				select case (trim(var))

				case ('l2')
				  vcfg%l2_file = trim(val)
				case ('output_file')
				  vcfg%output_l2 = trim(val)		
				case ('sensor')
				  vcfg%platform = trim(val)	
				case ('seg01')
				  vcfg%seg01_file = trim(val)
				case ('seg02')
				  vcfg%seg02_file = trim(val)
				case ('seg03')
				  vcfg%seg03_file = trim(val)
				case ('seg04')
			    vcfg%seg04_file = trim(val)
			  case ('seg05')
			    vcfg%seg05_file = trim(val)
			  case ('seg06')
			    vcfg%seg06_file = trim(val)			    
				case ('seg07')
				  vcfg%seg07_file = trim(val)
				case ('seg08')
				  vcfg%seg08_file = trim(val)
				case ('seg09')
				  vcfg%seg09_file = trim(val)
				case ('seg10')
			    vcfg%seg10_file = trim(val)

									  		  		  
				case default
					cycle
				end select
			else 
				cycle
			end if
		  
	
		end do
		
		return
		
	end function load_merge_config

!
!	load_viirs_config()
!	
!		Loads configuration file, config_file, into viirs_config object.
!	
!-------------------------------------------------------------------
	type(viirs_config_type) function load_viirs_config(config_file, status) result(vcfg)
		implicit none
		
		character(len=255), intent(in)	:: config_file
		integer, intent(inout)          :: status
		
		character(255)					        :: c_line, var, val
		integer							            :: eq_index
		
!   -- initialize our config variables.
    vcfg%year     = 1900
    vcfg%month    = 1
    vcfg%day      = 1
    vcfg%hour     = 0
    vcfg%minute   = 0
    vcfg%second   = 0
    
    vcfg%aerosol_land_file        = ''
    vcfg%aerosol_dust_file        = ''
    vcfg%aerosol_ocean_dust_file  = ''
    vcfg%aerosol_ocean_fine_file  = ''
    vcfg%aerosol_ocean_mari_file  = ''
    vcfg%aerosol_ocean_mix_file   = ''
    vcfg%bathymetry_lut_file      = ''
    vcfg%chl_lut_file             = ''
    vcfg%ler_lut_file             = ''
    vcfg%landcover_file           = ''
    vcfg%surfpressure_file        = ''
    vcfg%geozone_file             = ''
    vcfg%seasonal_deserts_file    = ''
    vcfg%brdfbase_file            = '' 
    vcfg%modis_surfdb_file        = ''
    vcfg%viirs_surfdb_file        = ''
    vcfg%surfcoeffs_file          = ''  
    vcfg%swir_vis_surfcoeffs_file = ''
    vcfg%veg_landcover_file       = ''
    vcfg%veg_sfc21_file           = ''
        
		vcfg%gmtco_file               = ''
		vcfg%svm01_file               = ''
		vcfg%svm02_file               = ''
		vcfg%svm03_file               = ''
		vcfg%svm04_file               = ''
		vcfg%svm05_file               = ''
		vcfg%svm07_file               = ''
		vcfg%svm08_file               = ''
		vcfg%svm09_file               = ''
		vcfg%svm10_file               = ''
		vcfg%svm11_file               = ''
		vcfg%svm14_file               = ''
		vcfg%svm15_file               = ''
		vcfg%svm16_file               = ''
		vcfg%iicmo_file               = ''
		
		vcfg%l1b_m                    = ''
		vcfg%geo_m                    = ''
		
		vcfg%gdas_file1               = ''
  	vcfg%gdas_file2               = ''
    
    vcfg%rayl412_file             = ''
    vcfg%rayl488_file             = ''
	  vcfg%rayl670_file             = ''
	  vcfg%xcal412_file             = ''
	  vcfg%xcal488_file             = ''
	  vcfg%xcal670_file             = ''
  	vcfg%output_l2                = ''
  	vcfg%platform                 = ''
  	vcfg%segment                  = ''
  	vcfg%goes_1                   = ''  	
  	vcfg%goes_2                   = '' 
  	vcfg%goes_3                   = '' 
  	vcfg%goes_4                   = '' 
  	vcfg%goes_5                   = '' 
  	vcfg%goes_6                   = '' 
  	vcfg%goes_7                   = '' 
  	vcfg%goes_8                   = '' 
  	vcfg%goes_9                   = '' 
  	vcfg%goes_10                  = '' 
  	vcfg%goes_11                  = '' 
  	vcfg%goes_12                  = '' 
  	vcfg%goes_13                  = '' 
  	vcfg%goes_14                  = '' 
  	vcfg%goes_15                  = '' 
  	vcfg%dbdt_file                = '' 
		
		open(100, file=config_file, status='OLD', FORM='FORMATTED',  ACTION='READ', iostat=status)
    if (status /= 0) then
      print *, 'ERROR: failed to open configuration file: ', status
      return
    end if 
	  	  
		do 
		! Read a line from the config file.  Die on error.
			read(100, fmt='(A)', iostat=status) c_line
			if (status /= 0) then
				if (status < 0) then	
				  status = 0	
					return
				else if (status > 0) then
					print *, "ERROR: Failed to read configuration file: ", status
					status = -1
					return
				end if
			end if
						
		! Clean up whitespace.
			c_line =rm_whitespace(c_line)
			
		! Detect comment lines or empty lines and skip.
			if (index(c_line, '!') /= 0) then
				cycle
			end if

			if (len_trim(c_line) == 0) then
				cycle
			end if
			
		! Begin parsing for variable and value.  
		! If no '=' is detected, we have a funky line, cycle.
			eq_index = index(c_line,'=')
			if (eq_index /= 0) then
				var = c_line(1:eq_index-1)
				val = c_line(eq_index+1:len(c_line))
				
				select case (trim(var))
				case ('year')
				  read(val, fmt='(I4)') vcfg%year
				case ('month')
				  read(val, fmt='(I2)') vcfg%month
				case ('day')
				  read(val, fmt='(I2)') vcfg%day
				case ('hour')
				  read(val, fmt='(I2)') vcfg%hour
				case ('minute')
				  read(val, fmt='(I2)') vcfg%minute
				case ('second')
				  read(val, fmt='(I2)') vcfg%second
				case ('output_file')
				  vcfg%output_l2 = trim(val)
				case ('aero_lut_land')
				  vcfg%aerosol_land_file = trim(val)  
				case ('aero_lut_dust')
				  vcfg%aerosol_dust_file = trim(val)  
				case ('aero_lut_ocean_dust')
				  vcfg%aerosol_ocean_dust_file = trim(val)  
				case ('aero_lut_ocean_fine')
				  vcfg%aerosol_ocean_fine_file = trim(val)  
				case ('aero_lut_ocean_mari')
				  vcfg%aerosol_ocean_mari_file = trim(val)  
				case ('aero_lut_ocean_mix')
				  vcfg%aerosol_ocean_mix_file = trim(val)
				case ('bathymetry_lut')
				  vcfg%bathymetry_lut_file = trim(val)
				case ('chl_lut')
				  vcfg%chl_lut_file = trim(val)
                                case ('ler_lut')
				  vcfg%ler_lut_file = trim(val)
				case ('land_cover')
				  vcfg%landcover_file = trim(val)
        case ('surface_pressure')
          vcfg%surfpressure_file = trim(val)				  
				case ('geozone')
				  vcfg%geozone_file = trim(val)
				case ('seasonal_deserts')
				  vcfg%seasonal_deserts_file = trim(val)
				case ('brdfbase')
				  vcfg%brdfbase_file = trim(val)
				case ('modis_surf_db')
				  vcfg%modis_surfdb_file = trim(val)
				case ('viirs_surf_db')
				  vcfg%viirs_surfdb_file = trim(val)
				case ('surf_coeffs')
				  vcfg%surfcoeffs_file = trim(val)
        case ('swir_vis_surf_coeffs')
          vcfg%swir_vis_surfcoeffs_file = trim(val)
				case ('veg_land_cover')
				  vcfg%veg_landcover_file = trim(val)
				case ('veg_21sfc')
				  vcfg%veg_sfc21_file = trim(val)  
				
				case ('rayl_412')
				  vcfg%rayl412_file = trim(val)
				case ('rayl_488')
				  vcfg%rayl488_file = trim(val)
				case ('rayl_670')
				  vcfg%rayl670_file = trim(val)
				case ('xcal_412')
				  vcfg%xcal412_file = trim(val)
				case ('xcal_488')
				  vcfg%xcal488_file = trim(val)
				case ('xcal_670')
				  vcfg%xcal670_file = trim(val)
				
				case ('gmtco')
				  vcfg%gmtco_file = trim(val)
				case ('iicmo')
				  vcfg%iicmo_file = trim(val)
				case ('svm01')
				  vcfg%svm01_file = trim(val)
				case ('svm02')
				  vcfg%svm02_file = trim(val)
				case ('svm03')
				  vcfg%svm03_file = trim(val)
				case ('svm04')
			    vcfg%svm04_file = trim(val)
			  case ('svm05')
			    vcfg%svm05_file = trim(val)
				case ('svm07')
				  vcfg%svm07_file = trim(val)
				case ('svm08')
				  vcfg%svm08_file = trim(val)
				case ('svm09')
				  vcfg%svm09_file = trim(val)
				case ('svm10')
			    vcfg%svm10_file = trim(val)
				case ('svm11')
				  vcfg%svm11_file = trim(val)
				case ('svm14')
				  vcfg%svm14_file = trim(val)
				case ('svm15')
				  vcfg%svm15_file = trim(val)
				case ('svm16')
				  vcfg%svm16_file = trim(val)
				case ('gdas1')
				  vcfg%gdas_file1 = trim(val)
				case ('gdas2')
				  vcfg%gdas_file2 = trim(val)
				case ('l1b_m')
				  vcfg%l1b_m = trim(val)
				case ('geo_m')
				  vcfg%geo_m = trim(val)
				case ('platform')
				  vcfg%platform = trim(val)		
				case ('segment')
				  vcfg%segment = trim(val)			
				case ('goes_01')
				  vcfg%goes_1 = trim(val)
				case ('goes_02')
				  vcfg%goes_2 = trim(val)
				case ('goes_03')
				  vcfg%goes_3 = trim(val)
				case ('goes_04')
				  vcfg%goes_4 = trim(val)
				case ('goes_05')
				  vcfg%goes_5 = trim(val)
				case ('goes_06')
				  vcfg%goes_6 = trim(val)
				case ('goes_07')
				  vcfg%goes_7 = trim(val)
				case ('goes_08')
				  vcfg%goes_8 = trim(val)
				case ('goes_09')
				  vcfg%goes_9 = trim(val)
				case ('goes_10')
				  vcfg%goes_10 = trim(val)
				case ('goes_11')
				  vcfg%goes_11 = trim(val)
				case ('goes_12')
				  vcfg%goes_12 = trim(val)
				case ('goes_13')
				  vcfg%goes_13 = trim(val)
				case ('goes_14')
				  vcfg%goes_14 = trim(val)
				case ('goes_15')
				  vcfg%goes_15 = trim(val)    
				case ('dbdt')  
				  vcfg%dbdt_file = trim(val)   
									  		  		  
				case default
					cycle
				end select
			else 
				cycle
			end if
		  
	
		end do
		
		return
		
	end function load_viirs_config

!
! rm_whitespace()
!
!	  Removes whitespace from a specified string.
!
!	  Inputs:		      input_string		string to be cleared of whitespace
!	  Outputs:		    Returns character array cleared of all whitespace.
!   Side effects:   None
!
!----------------------------------------------------------------------
	function rm_whitespace(input_string)
		implicit none
		
		character(len=*), intent(in)		:: input_string
		
		character(len=len(input_string))	:: rm_whitespace
		integer									:: space_index, end_index
		
		rm_whitespace = input_string

		! Empty string?
		if (len_trim(rm_whitespace) == 0) then
			return
		end if

		! Remove frontal spaces.
		do 
			space_index = index(rm_whitespace,' ')
			if ( space_index == 1) then
				rm_whitespace = rm_whitespace(2:len(rm_whitespace))
			else
				exit
			end if
		end do
		
		! Remove any internal spaces
		do
			space_index = index(rm_whitespace,' ')
			end_index = len_trim(rm_whitespace)
			if (space_index > 0 .AND. space_index < end_index) then
				rm_whitespace = rm_whitespace(1:space_index-1) //           &
		  &	rm_whitespace(space_index+1:end_index)
			else
				exit
			end if
		end do
		
		return 
	end function rm_whitespace

end module viirs_config
