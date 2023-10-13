MODULE OCIUAAER_Config_Module


IMPLICIT NONE

!---------------------------------------------------------------
!---------------------------------------------------------------
TYPE, PUBLIC  :: ociuaaer_config_type
  CHARACTER(LEN=255)      ::  read_nc_landonly = 'NULL'
  CHARACTER(LEN=255)      ::  input_l1file = 'NULL'
  CHARACTER(LEN=255)      ::  proxy_l1file = 'NULL'
  CHARACTER(LEN=255)      ::  output_dir = 'NULL'
!
  CHARACTER(LEN=255)      ::  landwater_mask = 'NULL'
!
  CHARACTER(LEN=255)	  ::  uv_ai_mielut = 'NULL'
  CHARACTER(LEN=255)      ::  uv_ai_lerlut = 'NULL'
  CHARACTER(LEN=255)      ::  uv_surfalbfile = 'NULL'
  CHARACTER(LEN=255)      ::  uv_surfprsfile = 'NULL'
  CHARACTER(LEN=255)      ::  uv_snowicefile = 'NULL'
  CHARACTER(LEN=255)      ::  uv_ocncorr_lut = 'NULL'
  CHARACTER(LEN=255)      ::  uv_aer_lut = 'NULL'
  CHARACTER(LEN=255)      ::  uv_sfc_file = 'NULL'
  CHARACTER(LEN=255)      ::  uv_zaer_file = 'NULL'
  CHARACTER(LEN=255)      ::  uv_zaer2_file = 'NULL'
  CHARACTER(LEN=255)      ::  uv_AIRSCO_clm_file = 'NULL'
  CHARACTER(LEN=255)      ::  uv_omacalut_file = 'NULL'
  CHARACTER(LEN=255)      ::  uv_ssa388_clim_file = 'NULL'
 
 ! DT OCEAN
  CHARACTER(LEN=255)      ::  dt_ext_coeff = 'NULL'
  CHARACTER(LEN=255)      ::  dt_gas_coeff = 'NULL'
  CHARACTER(LEN=255)      ::  VIIRS_ocean  = 'NULL'  
  CHARACTER(LEN=255)      ::  Pace_ocean   = 'NULL' 
   
 ! DT LAND 
  CHARACTER(LEN=255)      ::  VIIRS_land  = 'NULL'   
  CHARACTER(LEN=255)      ::  db_aero_lut_land = 'NULL'
  CHARACTER(LEN=255)      ::  db_aero_lut_dust = 'NULL'
  CHARACTER(LEN=255)      ::  db_config = 'NULL'
END TYPE ociuaaer_config_type
!
TYPE(ociuaaer_config_type), PUBLIC :: cfg
!---------------------------------------------------------------
!---------------------------------------------------------------

PUBLIC :: read_ociuaaer_config

CONTAINS

!---------------------------------------------------------------
!---------------------------------------------------------------
!
FUNCTION read_ociuaaer_config(par_file, cfg) RESULT(STATUS)

IMPLICIT NONE

CHARACTER(255), INTENT(IN) :: par_file

TYPE(ociuaaer_config_type), INTENT(INOUT) :: cfg
CHARACTER(255), PARAMETER :: config_file = 'OCIUAAER_Config'
CHARACTER*255 , Filename
CHARACTER(255)  :: c_line, var, val
INTEGER         :: STATUS, eq_index, nline

STATUS = 0
nline = 0

PRINT *, 'Now reading OCIUAAER Configuration file : ',  config_file	 
 OPEN(1101, FILE=par_file, STATUS='OLD', FORM='FORMATTED',  ACTION='READ', iostat=STATUS)
    IF (STATUS /= 0) THEN
      PRINT *, 'ERROR: failed to OPEN configuration file: ', STATUS
      RETURN
    ENDIF 
	  	  
    DO 
    ! Read a line from the config file.  Die on error.
    READ(1101, fmt='(A)', iostat=STATUS) c_line
    IF (STATUS /= 0) THEN
	IF (STATUS < 0) THEN	
        STATUS = 0	
        RETURN
	ELSE IF (STATUS > 0) THEN
	    PRINT *, "ERROR: Failed to read configuration file: ", STATUS
	    STATUS = -1
	    RETURN
	ENDIF
    ENDIF
						
    ! Clean up whitespace.
    c_line =    rm_whitespace(c_line)
			
    ! Detect comment lines or empty lines and skip.
    IF (index(c_line, '!') /= 0) THEN
	CYCLE
    ENDIF
    IF (LEN_trim(c_line) == 0) THEN
	CYCLE
    ENDIF
			
    ! Begin parsing for variable and value.  
    ! If no '=' is detected, we have a funky line, CYCLE.
    eq_index = index(c_line,'=')
    IF (eq_index /= 0) THEN
	var = c_line(1:eq_index-1)
	val = c_line(eq_index+1:LEN(c_line))
	!
	!print *, trim(var), trim(val)
	select case (trim(var))
        case ('read_nc_landonly')
      cfg%read_nc_landonly = trim(val)
        case ('input_l1file')
		cfg%input_l1file = trim(val)
        case ('proxy_l1file')
		cfg%proxy_l1file = trim(val)
        case ('output_dir')
		cfg%output_dir = trim(val)
        case ('landwater_mask')
      cfg%landwater_mask = trim(val)

	! UV-LUTs
        case ('uv_ai_mielut')
		cfg%uv_ai_mielut = trim(val)
        case ('uv_ai_lerlut')
		cfg%uv_ai_lerlut = trim(val)
        case ('uv_surfalbfile')
                cfg%uv_surfalbfile = trim(val)
        case ('uv_surfprsfile')
                cfg%uv_surfprsfile = trim(val)
        case ('uv_snowicefile')
                cfg%uv_snowicefile = trim(val)
        case ('uv_ocncorr_lut')
                cfg%uv_ocncorr_lut = trim(val)
        case ('uv_aer_lut')
                cfg%uv_aer_lut = trim(val)
        case ('uv_sfc_file')
                cfg%uv_sfc_file = trim(val)
        case ('uv_zaer_file')
                cfg%uv_zaer_file = trim(val)
        case ('uv_zaer2_file')
                cfg%uv_zaer2_file = trim(val)
        case ('uv_AIRSCO_clm_file')
                cfg%uv_AIRSCO_clm_file = trim(val)
        case ('uv_omacalut_file')
                cfg%uv_omacalut_file = trim(val)
        case ('uv_ssa388_clim_file')
                cfg%uv_ssa388_clim_file = trim(val)

	! DT-LUTs
	    case ('dt_ext_coeff')
		cfg%dt_ext_coeff = trim(val)  
	    case ('dt_gas_coeff')
		cfg%dt_gas_coeff = trim(val) 
	    case ('VIIRS_ocean')
		cfg%VIIRS_ocean = trim(val)   
	    case ('Pace_ocean')
		cfg%Pace_ocean = trim(val)  
	! DT-Land_LUTs  
	    case ('VIIRS_land')
		cfg%VIIRS_land = trim(val)    		
					
	! DB-LUTs
	case ('db_aero_lut_land')
		cfg%db_aero_lut_land = trim(val)  
	case ('db_aero_lut_dust')
		cfg%db_aero_lut_dust = trim(val) 
	case ('db_config')
		cfg%db_config = trim(val) 	
		print*,'db_config',cfg%db_config
	case default
		cycle
	end select
	!
    ELSE 
	CYCLE
    ENDIF
    ENDDO


END FUNCTION read_ociuaaer_config
!
!---------------------------------------------------------------
!---------------------------------------------------------------
!
FUNCTION rm_whitespace(input_string)

    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)	:: input_string
    CHARACTER(LEN=LEN(input_string))	:: rm_whitespace
    INTEGER				:: space_index, end_index
            
    rm_whitespace = input_string
    
    ! Empty string?
    IF (LEN_trim(rm_whitespace) == 0) THEN
        RETURN
    ENDIF
    
    ! Remove frontal spaces.
     DO 
      space_index = index(rm_whitespace,' ')
      IF ( space_index == 1) THEN
        rm_whitespace = rm_whitespace(2:LEN(rm_whitespace))
      ELSE
        EXIT
      ENDIF
     ENDDO
            
    ! Remove any internal spaces
     DO
      space_index = index(rm_whitespace,' ')
      end_index = LEN_trim(rm_whitespace)
      IF (space_index > 0 .AND. space_index < end_index) THEN
                rm_whitespace = rm_whitespace(1:space_index-1) // &
              & rm_whitespace(space_index+1:end_index)
      ELSE
       EXIT
      ENDIF
     ENDDO
            
    RETURN 
    
END FUNCTION rm_whitespace
!---------------------------------------------------------------
!---------------------------------------------------------------

End Module OCIUAAER_Config_Module    
