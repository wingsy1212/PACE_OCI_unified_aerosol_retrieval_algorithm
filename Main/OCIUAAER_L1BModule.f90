Module OCIUAAER_L1BModule

USE H5Util_class
USE HDF5
USE DataTypeDef
USE MyConstants

IMPLICIT NONE


!===============================================================
! OCI-L1B SDS PATHS
!===============================================================
CHARACTER(LEN=*), PARAMETER :: &
    path_lat = '/geolocation_data/latitude', &
    path_lon = '/geolocation_data/longitude', &
    path_saz = '/geolocation_data/solar_azimuth', &
    path_sza = '/geolocation_data/solar_zenith', &
    path_vza = '/geolocation_data/sensor_zenith', &
    path_vaz = '/geolocation_data/sensor_azimuth', &
    path_dem = '/geolocation_data/height', &

    path_time = '/scan_line_attributes/time', &

    path_swir_wav = '/sensor_band_parameters/SWIR_wavelength', &
    path_blue_wav = '/sensor_band_parameters/blue_wavelength', &
    path_red_wav  = '/sensor_band_parameters/red_wavelength', &

    path_swir_irrad = '/sensor_band_parameters/SWIR_solar_irradiance', &
    path_blue_irrad = '/sensor_band_parameters/blue_solar_irradiance', &
    path_red_irrad  = '/sensor_band_parameters/red_solar_irradiance', &

!    path_swir_rad = '/observation_data/Lt_SWIR', &
!    path_blue_rad = '/observation_data/Lt_blue', &
!    path_red_rad  = '/observation_data/Lt_red'

     path_swir_rad = '/observation_data/rhot_SWIR', &
     path_blue_rad = '/observation_data/rhot_blue', &
     path_red_rad  = '/observation_data/rhot_red'


TYPE, PUBLIC :: l1b_rad_type
  INTEGER(HID_T) 	:: file_id
  CHARACTER(LEN=22) 	:: l1b_date_time_str
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE   :: l1b_lat, l1b_lon, &
                                                 l1b_sza, l1b_vza, &
                                                 l1b_saz, l1b_vaz, l1b_dem
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: l1b_swir_rad, l1b_red_rad, l1b_blue_rad
  REAL(KIND=4), DIMENSION(:),	  ALLOCATABLE :: l1b_time, l1b_swir_wav, l1b_red_wav, l1b_blue_wav
  INTEGER(KIND=4) :: l1b_nLines, l1b_nXTrack, l1b_swir_nWav, l1b_red_nWav, l1b_blue_nWav
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: l1b_wv_idx
  INTEGER(KIND=4) :: dims_geo(2), dims_swir_rad(3), dims_red_rad(3), dims_blue_rad(3), &
                                  dims_swir_wav(1), dims_red_wav(1), dims_blue_wav(1)
END TYPE l1b_rad_type
TYPE(l1b_rad_type), PUBLIC :: L1B
!!
!!
TYPE, PUBLIC :: lwmask_type
  INTEGER(KIND=4), DIMENSION(:,:),ALLOCATABLE :: val
  INTEGER(KIND=4) :: dims(2), nlon, nlat
  REAL(KIND=4), DIMENSION(:),  ALLOCATABLE :: lat, lon
END TYPE lwmask_type
TYPE(lwmask_type), PUBLIC :: lwmask
!!
!! 
INTEGER(KIND=4) ,DIMENSION(:,:),  ALLOCATABLE ::  arr_2di4
!
!======================
!   PUBLIC METHODS
!======================
PUBLIC :: OCI_l1b_get_reflect
PUBLIC :: OCI_l1b_DEALLOCATE
PUBLIC :: OCI_l1b_read_date


CONTAINS
!
!===============================================================
!===============================================================
! BELOW ARE THE FUNCTIONS, SUBROUTINES DEFINED IN THIS MODULE.
!===============================================================
!===============================================================
FUNCTION OCI_l1b_get_alldims(l1b_rad_file, L1B) RESULT(STATUS)

USE H5read_module
  
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN)    :: l1b_rad_file
INTEGER(HID_T) :: rad_file_id, dataspace_id, dataset_id, datatype_id
INTEGER(HSIZE_T), DIMENSION(7) :: maxdims
INTEGER :: STATUS, tmp_rank, errorStatus
INTEGER(KIND=8),DIMENSION(7) :: tmp_dims
TYPE(l1b_rad_type), INTENT(INOUT) :: L1B
!

! Open HDF-5 FORTRAN interface
CALL H5fopen_f(l1b_rad_file, H5F_ACC_RDONLY_F, rad_file_id, errorStatus)
!
! Get DIMENSIONs from Latitude and (SWIR/Red/Blue) wavelengths
! PATH_LAT
 CALL get_rank_dims(rad_file_id, path_lat, tmp_rank, tmp_dims)
L1B%dims_geo(:) = tmp_dims(1:tmp_rank)
L1B%l1b_nXTrack= L1B%dims_geo(1)
L1B%l1b_nLines = L1B%dims_geo(2)

! PATH_SWIR_WAV
 CALL get_rank_dims(rad_file_id, path_swir_wav, tmp_rank, tmp_dims)
L1B%dims_swir_wav(:) = tmp_dims(1:tmp_rank)
L1B%l1b_swir_nWav = tmp_dims(1)

! PATH_RED_WAV
 CALL get_rank_dims(rad_file_id, path_red_wav, tmp_rank, tmp_dims)
L1B%dims_red_wav(:) = tmp_dims(1:tmp_rank)
L1B%l1b_red_nWav = tmp_dims(1)

! PATH_BLUE_WAV
 CALL get_rank_dims(rad_file_id, path_blue_wav, tmp_rank, tmp_dims)
L1B%dims_blue_wav(:) = tmp_dims(1:tmp_rank)
L1B%l1b_blue_nWav = tmp_dims(1)

L1B%dims_swir_rad = (/L1B%l1b_nXTrack, L1B%l1b_nLines, L1B%l1b_swir_nWav/)
L1B%dims_red_rad  = (/L1B%l1b_nXTrack, L1B%l1b_nLines, L1B%l1b_red_nWav /)
L1B%dims_blue_rad = (/L1B%l1b_nXTrack, L1B%l1b_nLines, L1B%l1b_blue_nWav/)

!PRINT*,'Number of Bands from 310-600 nm at 2.5 intervals = ',L1B%l1b_blue_nWav
!PRINT*,'Number of Bands from 600-890 nm at 2.5 intervals = ',L1B%l1b_red_nWav
!PRINT*,'Number of SWIR Discrete Bands from 900-2260 nm   = ',L1B%l1b_SWIR_nWav
!PRINT*,'Number of Scanlines in the granule (AlongTrack) = ' ,L1B%l1b_nLines
!PRINT*,'Number of Pixels per scan (AcrossTrack)  =        ' ,L1B%l1b_nXTrack

CALL H5fclose_f(rad_file_id, errorStatus)

STATUS=1

END FUNCTION OCI_l1b_get_alldims
!
!===============================================================
!===============================================================
!
SUBROUTINE OCI_l1b_init(L1B) 

IMPLICIT NONE

INTEGER :: STATUS
TYPE(l1b_rad_type), INTENT(INOUT) :: L1B

 STATUS = 1
 L1B%l1b_lat = FILLVALUE_SP
 L1B%l1b_lon = FILLVALUE_SP
 L1B%l1b_sza = FILLVALUE_SP
 L1B%l1b_vza = FILLVALUE_SP
 L1B%l1b_saz = FILLVALUE_SP
 L1B%l1b_vaz = FILLVALUE_SP
 L1B%l1b_dem = FILLVALUE_SP

 L1B%l1b_Blue_rad = FILLVALUE_SP
 L1B%l1b_Red_rad  = FILLVALUE_SP
 L1B%l1b_SWIR_rad = FILLVALUE_SP

END SUBROUTINE OCI_l1b_init
!
!===============================================================
!===============================================================
!
SUBROUTINE OCI_l1b_ALLOCATE(L1B) 

IMPLICIT NONE

INTEGER :: STATUS

TYPE(l1b_rad_type), INTENT(INOUT) :: L1B
STATUS =1

ALLOCATE(L1B%l1b_lat(L1B%l1b_nXTrack, L1B%l1b_nLines))
ALLOCATE(L1B%l1b_lon(L1B%l1b_nXTrack, L1B%l1b_nLines))
ALLOCATE(L1B%l1b_sza(L1B%l1b_nXTrack, L1B%l1b_nLines))
ALLOCATE(L1B%l1b_vza(L1B%l1b_nXTrack, L1B%l1b_nLines))
ALLOCATE(L1B%l1b_saz(L1B%l1b_nXTrack, L1B%l1b_nLines))
ALLOCATE(L1B%l1b_vaz(L1B%l1b_nXTrack, L1B%l1b_nLines))
ALLOCATE(L1B%l1b_dem(L1B%l1b_nXTrack, L1B%l1b_nLines))
!
ALLOCATE(L1B%l1b_time(L1B%l1b_nLines) )
!
ALLOCATE(L1B%l1b_Blue_rad(L1B%l1b_nXTrack, L1B%l1b_nLines, L1B%l1b_Blue_nWav))
ALLOCATE(L1B%l1b_Red_rad( L1B%l1b_nXTrack, L1B%l1b_nLines, L1B%l1b_Red_nWav ))
ALLOCATE(L1B%l1b_SWIR_rad(L1B%l1b_nXTrack, L1B%l1b_nLines, L1B%l1b_SWIR_nWav))
ALLOCATE(L1B%l1b_Blue_wav(L1B%l1b_Blue_nWav))
ALLOCATE(L1B%l1b_Red_wav( L1B%l1b_Red_nWav ))
ALLOCATE(L1B%l1b_SWIR_wav(L1B%l1b_SWIR_nWav))

CALL OCI_l1b_init(L1B)

END SUBROUTINE OCI_l1b_ALLOCATE
!
!===============================================================
!===============================================================
!
SUBROUTINE OCI_l1b_DEALLOCATE(L1B) ! RESULT(STATUS)

IMPLICIT NONE

INTEGER :: STATUS

TYPE(l1b_rad_type), INTENT(INOUT) :: L1B
STATUS =1

IF (ALLOCATED(L1B%l1b_lat)) DEALLOCATE(L1B%l1b_lat)
IF (ALLOCATED(L1B%l1b_lon)) DEALLOCATE(L1B%l1b_lon)
IF (ALLOCATED(L1B%l1b_sza)) DEALLOCATE(L1B%l1b_sza)
IF (ALLOCATED(L1B%l1b_vza)) DEALLOCATE(L1B%l1b_vza)
IF (ALLOCATED(L1B%l1b_saz)) DEALLOCATE(L1B%l1b_saz)
IF (ALLOCATED(L1B%l1b_vaz)) DEALLOCATE(L1B%l1b_vaz)
IF (ALLOCATED(L1B%l1b_dem)) DEALLOCATE(L1B%l1b_dem)
!
IF (ALLOCATED(L1B%l1b_time)) DEALLOCATE(L1B%l1b_time)
!
IF (ALLOCATED(L1B%l1b_Blue_rad)) DEALLOCATE(L1B%l1b_Blue_rad)
IF (ALLOCATED(L1B%l1b_Red_rad )) DEALLOCATE(L1B%l1b_Red_rad )
IF (ALLOCATED(L1B%l1b_SWIR_rad)) DEALLOCATE(L1B%l1b_SWIR_rad)
IF (ALLOCATED(L1B%l1b_Blue_wav)) DEALLOCATE(L1B%l1b_Blue_wav)
IF (ALLOCATED(L1B%l1b_Red_wav )) DEALLOCATE(L1B%l1b_Red_wav )
IF (ALLOCATED(L1B%l1b_SWIR_wav)) DEALLOCATE(L1B%l1b_SWIR_wav)

END SUBROUTINE OCI_l1b_DEALLOCATE
!
!===============================================================
!===============================================================
!
FUNCTION OCI_l1b_get_reflect(l1b_rad_file, l1b_nXTrack, l1b_nLines, l1b_date_time_str, &
                             Latitude, Longitude, SolarZenithAngle, ViewingZenithAngle, &
                             SolarAzimuthAngle, ViewingAzimuthAngle, TerrainHeight, &
                             UVtoSWIR_nWavel, UVtoSWIR_wavelengths, & 
                             UVtoSWIR_Reflectances) RESULT(STATUS)
USE H5read_module

IMPLICIT NONE

Include 'common_l1b_var.inc'

CHARACTER(LEN=255), INTENT(IN)  :: l1b_rad_file
INTEGER(KIND=4),    INTENT(IN)  :: UVtoSWIR_nWavel
REAL(KIND=4), DIMENSION(UVtoSWIR_nWavel), INTENT(IN)     :: UVtoSWIR_wavelengths 
REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: UVtoSWIR_Reflectances
!
TYPE(l1b_rad_type) :: L1B
INTEGER :: STATUS, index1
INTEGER(KIND=4) :: errorStatus
CHARACTER(LEN=128) :: str1
REAL(KIND=4), DIMENSION(:), ALLOCATABLE    :: UVtoSWIR_IRRadiances
!

STATUS =-1

! Reads the L1B file to extract DIMENSIONs of the granule
STATUS = OCI_l1b_get_alldims(l1b_rad_file, L1B)

! Allocate memory for all OCI-L1B variables
CALL OCI_l1b_ALLOCATE(L1B)

! Read L1B variables
! First integers
ALLOCATE( arr_2di4(L1B%l1b_nXTrack, L1B%l1b_nLines) )
errorStatus = read_2di4(l1b_rad_file, path_sza, L1B%dims_geo, arr_2di4)
L1B%l1b_sza(:,:) = arr_2di4(:,:)
errorStatus = read_2di4(l1b_rad_file, path_vza, (/L1B%l1b_nXTrack, L1B%l1b_nLines/), arr_2di4, str1)
L1B%l1b_vza(:,:) = arr_2di4(:,:)
errorStatus = read_2di4(l1b_rad_file, path_saz, (/L1B%l1b_nXTrack, L1B%l1b_nLines/), arr_2di4, str1)
L1B%l1b_saz(:,:) = arr_2di4(:,:)
errorStatus = read_2di4(l1b_rad_file, path_vaz, (/L1B%l1b_nXTrack, L1B%l1b_nLines/), arr_2di4, str1)
L1B%l1b_vaz(:,:) = arr_2di4(:,:)
errorStatus = read_2di4(l1b_rad_file, path_dem, (/L1B%l1b_nXTrack, L1B%l1b_nLines/), arr_2di4, str1)
L1B%l1b_dem(:,:) = arr_2di4(:,:)

! Geometry has a scale_factor 
! For now it is hard-coded. Will update the HDF-5 utilities later.
L1B%l1b_sza = L1B%l1b_sza/100.
L1B%l1b_saz = L1B%l1b_saz/100.
L1B%l1b_vza = L1B%l1b_vza/100.
L1B%l1b_vaz = L1B%l1b_vaz/100.

! Now real numbers
errorStatus = read_h5_2dr4(l1b_rad_file, path_lat, L1B%dims_geo, L1B%l1b_lat)
errorStatus = read_h5_2dr4(l1b_rad_file, path_lon, L1B%dims_geo, L1B%l1b_lon)

errorStatus = read_h5_1dr4(l1b_rad_file, path_Blue_wav, L1B%dims_Blue_wav, L1B%l1b_Blue_wav)
errorStatus = read_h5_1dr4(l1b_rad_file, path_Red_wav,  L1B%dims_Red_wav,  L1B%l1b_Red_wav )
errorStatus = read_h5_1dr4(l1b_rad_file, path_SWIR_wav, L1B%dims_SWIR_wav, L1B%l1b_SWIR_wav)

errorStatus = read_h5_3dr4(l1b_rad_file, path_Blue_rad, L1B%dims_Blue_rad, L1B%l1b_Blue_rad)
errorStatus = read_h5_3dr4(l1b_rad_file, path_Red_rad,  L1B%dims_Red_rad,  L1B%l1b_Red_rad )
errorStatus = read_h5_3dr4(l1b_rad_file, path_SWIR_rad, L1B%dims_SWIR_rad, L1B%l1b_SWIR_rad)

! Now get discrete band refelctances.
CALL OCI_l1b_select_wv(L1B, UVtoSWIR_nWavel, &
			    UVtoSWIR_wavelengths, & 
                            UVtoSWIR_Reflectances)

! Assign all required L1B variables to common variables. 
index1 = index(l1b_rad_file, 'PACE_OCI')
l1b_date_time_str = l1b_rad_file(index1+9:index1+9+14)
!PRINT *,'L1B start date-time : ',l1b_date_time_str
l1b_nXTrack = L1B%l1b_nXTrack
l1b_nLines = L1B%l1b_nLines
Latitude = L1B%l1b_lat
Longitude = L1B%l1b_lon
SolarZenithAngle = L1B%l1b_sza
ViewingZenithAngle = L1B%l1b_vza
SolarAzimuthAngle = L1B%l1b_saz
ViewingAzimuthAngle = L1B%l1b_vaz
TerrainHeight = L1B%l1b_dem

CALL OCI_l1b_DEALLOCATE(L1B)

STATUS = 1

RETURN

END FUNCTION OCI_l1b_get_reflect
!===============================================================

!===============================================================
FUNCTION OCI_l1b_read_date(l1b_rad_file,  l1b_date_time_str,l1b_nXTrack, l1b_nLines) RESULT(STATUS)

USE H5read_module

IMPLICIT NONE

CHARACTER(LEN=255), INTENT(IN) :: l1b_rad_file
CHARACTER(LEN=256), INTENT(inout) :: l1b_date_time_str
INTEGER(KIND=4), INTENT(out)  :: l1b_nXTrack, l1b_nLines
INTEGER :: STATUS, index1
INTEGER(KIND=4) :: errorStatus
CHARACTER(LEN=128) :: str1
TYPE(l1b_rad_type) :: L1B

STATUS =1
!print *, 'start OCI_l1b_read_date'

index1 = index(l1b_rad_file, 'PACE_OCI_SIM')
l1b_date_time_str = l1b_rad_file(index1+13:index1+13+14)
!PRINT *,'L1B start date-time : ',l1b_date_time_str

STATUS = OCI_l1b_get_alldims(l1b_rad_file, L1B)
l1b_nXTrack = L1B%l1b_nXTrack
l1b_nLines = L1B%l1b_nLines
RETURN

END FUNCTION OCI_l1b_read_date
!===============================================================

!===============================================================

SUBROUTINE OCI_l1b_select_wv(L1B, UVtoSWIR_nWavel, &
                                  UVtoSWIR_wavelengths, &
                                  UVtoSWIR_Reflectances)

USE InterpolationModule

IMPLICIT NONE

TYPE(l1b_rad_type), INTENT(INOUT) :: L1B
INTEGER(KIND=4), INTENT(IN) :: UVtoSWIR_nWavel
REAL(KIND=4), DIMENSION(UVtoSWIR_nWavel), INTENT(IN)     :: UVtoSWIR_wavelengths 
REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: UVtoSWIR_Reflectances
!
REAL(KIND=4), DIMENSION(:), ALLOCATABLE    :: UVtoSWIR_IRRadiances
REAL(KIND=4) :: FractionalVal, BoundVal1, BoundVal2, rad_pixel, this_wavelen
INTEGER, PARAMETER :: fp_diag = 10061
INTEGER :: ii,jj, idx_wv, index1, index2, STATUS, spectra_nWav
REAL(KIND=4), DIMENSION(:),     ALLOCATABLE :: spectra_wav
REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: spectra_data
REAL(KIND=4), PARAMETER :: DTOR = (PI/180.)
REAL(KIND=4) :: cossza

ALLOCATE(UVtoSWIR_Reflectances(UVtoSWIR_nWavel, L1B%l1b_nXTrack, L1B%l1b_nLines))
ALLOCATE(UVtoSWIR_IRRadiances(UVtoSWIR_nWavel))

! For now Hard-code these values 
! (V8 synthetic files will include solar irradiances directly!!)
! Irradiances (W m-2 um-1 sr-1)
UVtoSWIR_IRRadiances(1) = 10.0099 ! 340 
UVtoSWIR_IRRadiances(2) = 10.8003 ! 354
UVtoSWIR_IRRadiances(3) = 10.8937 ! 388
UVtoSWIR_IRRadiances(4) = 16.9532 ! 412
UVtoSWIR_IRRadiances(5) = 19.1250 ! 488
UVtoSWIR_IRRadiances(6) = 18.7093 ! 550
UVtoSWIR_IRRadiances(7) = 15.1462 ! 670
UVtoSWIR_IRRadiances(8) = 14.7619 ! 680
UVtoSWIR_IRRadiances(9) = 14.6245 ! 688
UVtoSWIR_IRRadiances(10)=  9.7864 ! 860
UVtoSWIR_IRRadiances(11)=  4.4829 ! 1250.15
UVtoSWIR_IRRadiances(12)=  3.6468 ! 1378.04
UVtoSWIR_IRRadiances(13)=  2.4439 ! 1615.98
UVtoSWIR_IRRadiances(14)=  0.7685 ! 2260.64

! Units of Irradiances are wrong in the V7 L1B-data
UVtoSWIR_IRRadiances = UVtoSWIR_IRRadiances * 100.

DO 30 idx_wv = 1 , UVtoSWIR_nWavel  !Wavelength Loop
  this_wavelen = UVtoSWIR_wavelengths(idx_wv)

      IF (this_wavelen.LE.599) THEN
           ALLOCATE(spectra_wav(L1B%l1b_blue_nWav))
           ALLOCATE(spectra_data(L1B%l1b_nXTrack, L1B%l1b_nLines, L1B%l1b_blue_nWav))
           spectra_nWav = L1B%l1b_blue_nWav          
           spectra_wav = L1B%l1b_blue_wav
           spectra_data = L1B%l1b_Blue_rad
      ELSEIF (this_wavelen.GE.600 .AND. this_wavelen.LE.899) THEN 
           ALLOCATE(spectra_wav(L1B%l1b_red_nWav))
           ALLOCATE(spectra_data(L1B%l1b_nXTrack, L1B%l1b_nLines, L1B%l1b_red_nWav))
           spectra_nWav = L1B%l1b_red_nWav          
           spectra_wav = L1B%l1b_red_wav
           spectra_data = L1B%l1b_Red_rad
      ELSEIF (this_wavelen.GT.900) THEN
           ALLOCATE(spectra_wav(L1B%l1b_swir_nWav))
           ALLOCATE(spectra_data(L1B%l1b_nXTrack, L1B%l1b_nLines, L1B%l1b_SWIR_nWav))
           spectra_nWav = L1B%l1b_swir_nWav          
           spectra_wav = L1B%l1b_swir_wav
           spectra_data = L1B%l1b_SWIR_rad
      ENDIF

    !Find indices in the wavelength DIMENSION
    STATUS = FindTableEntry(this_wavelen, spectra_wav, spectra_nWav,&
                            BoundVal1, BoundVal2, Index1,Index2,FractionalVal)
    !WRITE(*, '(F8.3,1X,2(I4,1X),F8.4,1X,2(F11.3))') & 
    !  this_wavelen,Index1,Index2,FractionalVal,BoundVal1,BoundVal2
          
  DO 20 ii =1, L1B%l1b_nXTrack   ! X-track Loop
    DO 10 jj =1, L1B%l1b_nLines  !Long_track looop

 
     IF (spectra_data(ii,jj,Index1) .GT. 1.0e-6 .AND. &
         spectra_data(ii,jj,Index2) .GT. 1.0e-6) THEN	 
     
        !Interpolate radiances to the current wavelength for this pixel
        STATUS = Interp1D( spectra_data(ii,jj,Index1), &
                           spectra_data(ii,jj,Index2), &
                           FractionalVal, rad_pixel ) 
        cossza = COS(DTOR * L1B%l1b_sza(ii,jj))
 !        UVtoSWIR_Reflectances(idx_wv,ii,jj) =(rad_pixel*PI)/(UVtoSWIR_IRRadiances(idx_wv)*cossza)
          UVtoSWIR_Reflectances(idx_wv,ii,jj) = rad_pixel 
         
      ELSE			 
	rad_pixel = -9999.
	UVtoSWIR_Reflectances(idx_wv,ii,jj) = -9999.9		 
      END IF
      
      
10 END DO
20 END DO

IF (ALLOCATED(spectra_wav)) DEALLOCATE(spectra_wav)
IF (ALLOCATED(spectra_data)) DEALLOCATE(spectra_data)

30 END DO
!PRINT *,'Wavelengths of interest from UVtoSWIR are selected.'

! DEALLOCATE THE HUGE ARRAYS TO SAVE MEMORY
!IF (ALLOCATED(L1B%l1b_Blue_rad)) DEALLOCATE(L1B%l1b_Blue_rad)
!IF (ALLOCATED(L1B%l1b_Red_rad )) DEALLOCATE(L1B%l1b_Red_rad )
!IF (ALLOCATED(L1B%l1b_SWIR_rad)) DEALLOCATE(L1B%l1b_SWIR_rad)

END SUBROUTINE


!===============================================================
!===============================================================
FUNCTION get_gmt15arc_lwmask(lwmask_file, l1b_nXTrack, l1b_nLines, &
                                          l1b_lat, l1b_lon, &
					  grn_mask) RESULT(STATUS)

USE H5read_module

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: lwmask_file
INTEGER(KIND=4),  INTENT(IN) :: l1b_nXTrack, l1b_nLines
REAL(KIND=4),DIMENSION(:,:), INTENT(IN) :: l1b_lat, l1b_lon

TYPE(lwmask_type)  :: lwmask
INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT) :: grn_mask

INTEGER(HID_T) :: lwmask_fileid
INTEGER :: STATUS, tmp_rank, errorStatus
INTEGER :: ix, iy
INTEGER(KIND=8),DIMENSION(7) :: tmp_dims
CHARACTER(LEN=128) :: str1
!
STATUS = -1

CALL H5fopen_f(lwmask_file, H5F_ACC_RDONLY_F, lwmask_fileid, errorStatus)
CALL get_rank_dims(lwmask_fileid, 'watermask', tmp_rank, tmp_dims)
CALL H5fclose_f(lwmask_fileid, errorStatus)

lwmask%dims(1:tmp_rank) = tmp_dims(1:tmp_rank)
lwmask%nlon = lwmask%dims(1)
lwmask%nlat = lwmask%dims(2)

  !! Allocate memory to read global LandWater-mask
  ALLOCATE( lwmask%val( lwmask%nlon, lwmask%nlat ), &
            lwmask%lon( lwmask%nlon ), &
            lwmask%lat( lwmask%nlat ), STAT = STATUS)
  IF (STATUS < 0) THEN 
      PRINT *,'Error : Allocation of variables for Land-water mask failed.'
      CALL EXIT(1)
  ENDIF     

errorStatus = read_2di4(lwmask_file, 'watermask', lwmask%dims, lwmask%val, str1)
!errorStatus = read_h5_1dr4(lwmask_file, 'lat', lwmask%dims(2), lwmask%lat)
!errorStatus = read_h5_1dr4(lwmask_file, 'lon', lwmask%dims(1), lwmask%lon)

do ix = 1, lwmask%nlon 
   lwmask%lon(ix) = ix * (360./(lwmask%nlon-1)) - 180.
end do 
do iy = 1, lwmask%nlat
   lwmask%lat(iy) = iy * (180./(lwmask%nlat-1)) - 90.
end do 


ALLOCATE( grn_mask(l1b_nXTrack, l1b_nLines), STAT=STATUS )
  IF (STATUS < 0) THEN
     PRINT *,"Error : Error allocating granule mask data"
     CALL EXIT(1)
  ENDIF 

DO iy = 1, l1b_nLines
    DO ix = 1, l1b_nXTrack

      STATUS = Get_granule_lwmask(lwmask, l1b_lat(ix,iy), &
                                          l1b_lon(ix,iy), &
                                         grn_mask(ix,iy))
  	IF (STATUS < 0) THEN
     	  PRINT *,"Error : Error getting grn_mask"
     	  CALL EXIT(1)
	ENDIF 

    END DO
END DO 

! DEALLOCATE THE HUGE ARRAYS TO SAVE MEMORY
IF (ALLOCATED(lwmask%val)) DEALLOCATE(lwmask%val)
IF (ALLOCATED(lwmask%lat)) DEALLOCATE(lwmask%lat)
IF (ALLOCATED(lwmask%lon)) DEALLOCATE(lwmask%lon)

STATUS = 1

RETURN

END FUNCTION get_gmt15arc_lwmask
!===============================================================
!===============================================================
!
FUNCTION Get_granule_lwmask(lwmask, lat, lon, thismask) RESULT(STATUS)

  IMPLICIT NONE

  REAL(KIND=4),    INTENT(IN)  :: lat, lon
  INTEGER(KIND=4), INTENT(OUT) :: thismask
  INTEGER(KIND=4) :: ilon, ilat, STATUS
  TYPE(lwmask_type), INTENT(IN)  :: lwmask

! -- Extract the mask ---
   IF (lat .lt. 0.0) then 
       ilat = ((lwmask%nlat-1)/2) + int(((lwmask%nlat-1)/180.)*lat + 0.000001)
   ENDIF
   IF (lat .ge. 0.0) then 
       ilat = ((lwmask%nlat-1)/2)+ 1 + int(((lwmask%nlat-1)/180.)*lat - 0.000001)
   ENDIF
   IF (lon .lt. 0.0) then 
       ilon = ((lwmask%nlon-1)/2) + int(((lwmask%nlon-1)/360.)*lon + 0.000001)
   ENDIF
   IF (lon .ge. 0.0) then 
       ilon = ((lwmask%nlon-1)/2)+ 1 + int(((lwmask%nlon-1)/360.)*lon - 0.000001)
   ENDIF
!
   if (ilat .gt. lwmask%nlat)  ilat = lwmask%nlat
   if (ilat .lt. 1) ilat = 1
   if (ilon .gt. lwmask%nlon) ilon = lwmask%nlon
   if (ilon .lt. 1) ilon = 1
!
   thismask = lwmask%val(ilon,ilat)
!
  STATUS = 1

  RETURN

 END FUNCTION Get_granule_lwmask

!===============================================================
!===============================================================

End Module OCIUAAER_L1BModule
!
