Module proxydata_L1BModule

USE H5Util_class
USE HDF5
USE DataTypeDef
USE MyConstants

IMPLICIT NONE


!===============================================================
! Trop-in-Viirs SDS PATHS
!===============================================================
CHARACTER(LEN=*), PARAMETER :: &
    path_lat1 = '/VIIRS_Data/Viirs_Lat', &
    path_lon1 = '/VIIRS_Data/Viirs_Lon', &
    path_saz1 = '/VIIRS_Data/Viirs_SolAzi', &
    path_sza1 = '/VIIRS_Data/Viirs_SolZen', &
    path_vza1 = '/VIIRS_Data/Viirs_SenZen', &
    path_vaz1 = '/VIIRS_Data/Viirs_SenAzi', &
    path_dem1 = '/VIIRS_Data/Viirs_TerrainHeight', &

    path_Ref340 = '/TropOMI_Data/Trop_Ref340', &
    path_Ref354 = '/TropOMI_Data/Trop_Ref354', &
    path_Ref388 = '/TropOMI_Data/Trop_Ref388', &

    path_Ref412 = '/VIIRS_Data/Viirs_Ref415', &
    path_Ref488 = '/VIIRS_Data/Viirs_Ref490', &
    path_Ref550 = '/VIIRS_Data/Viirs_Ref550', &
    path_Ref670 = '/VIIRS_Data/Viirs_Ref670', &
    path_Ref680 = '/VIIRS_Data/Viirs_Ref670', &
    path_Ref688 = '/VIIRS_Data/Viirs_Ref670', &
    path_Ref860 = '/VIIRS_Data/Viirs_Ref865', &
    path_Ref1250= '/VIIRS_Data/Viirs_Ref1240', &
    path_Ref1378= '/VIIRS_Data/Viirs_Ref1380', &
    path_Ref1640= '/VIIRS_Data/Viirs_Ref1640', &
    path_Ref2260= '/VIIRS_Data/Viirs_Ref2250'

TYPE, PUBLIC :: l1b_proxy_type
  INTEGER(HID_T) 	:: file_id
  CHARACTER(LEN=22) 	:: l1b_date_time_str
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE   :: l1b_lat, l1b_lon, &
                                                 l1b_sza, l1b_vza, &
                                                 l1b_saz, l1b_vaz, l1b_dem
  INTEGER(KIND=4) :: l1b_nLines, l1b_nXTrack
  REAL(KIND=4), DIMENSION(:,:),  ALLOCATABLE :: Ref340, Ref354, Ref388, Ref412, Ref488, Ref550, Ref670, &
  						Ref680, Ref688, Ref860, Ref1250, Ref1378, Ref1640, Ref2260
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: l1b_wv_idx
  INTEGER(KIND=4) :: dims_geo(2)
END TYPE l1b_proxy_type
TYPE(l1b_proxy_type), PUBLIC :: proxy
!!
!! 
INTEGER(KIND=4) ,DIMENSION(:,:),  ALLOCATABLE ::  arr_2di4
!
!======================
!   PUBLIC METHODS
!======================
PUBLIC :: read_l1bproxydata
PUBLIC :: read_proxydate


CONTAINS
!
!
FUNCTION read_l1bproxydata(proxy_filename, l1b_nXTrack, l1b_nLines, l1b_date_time_str, &
                             Latitude, Longitude, SolarZenithAngle, ViewingZenithAngle, &
                             SolarAzimuthAngle, ViewingAzimuthAngle, TerrainHeight, &
                             UVtoSWIR_nWavel, UVtoSWIR_wavelengths, & 
                             UVtoSWIR_Reflectances) RESULT(status)
USE H5read_module

IMPLICIT NONE

Include 'common_l1b_var.inc'

CHARACTER(LEN=255), INTENT(IN)  :: proxy_filename
INTEGER(KIND=4),    INTENT(IN)  :: UVtoSWIR_nWavel
REAL(KIND=4), DIMENSION(UVtoSWIR_nWavel), INTENT(IN)     :: UVtoSWIR_wavelengths 
REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: UVtoSWIR_Reflectances
!
TYPE(l1b_proxy_type) :: Proxy
REAL(KIND=4), DIMENSION(:), ALLOCATABLE    :: UVtoSWIR_IRRadiances
INTEGER(KIND=4) :: errorStatus, STATUS, tmp_rank, index1
INTEGER(KIND=8),dimension(7) :: tmp_dims
INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id
CHARACTER(LEN=128) :: str1
REAL(KIND=4), PARAMETER :: DTOR = (PI/180.)
REAL(KIND=4)      :: cossza
INTEGER(KIND=4)   :: iPix, jPix
 
status =1

! Open HDF-5 files
 CALL h5fopen_f(proxy_filename, H5F_ACC_RDONLY_F, file_id, errorStatus)

 CALL get_rank_dims(file_id, path_lat1, tmp_rank, tmp_dims)
 PRINT*,'Proxydata file Lat rank:',tmp_rank
 PRINT*,'Proxydata file Lat shape:',tmp_dims(1:tmp_rank)

 Proxy%dims_geo(:) = tmp_dims(1:tmp_rank)
 Proxy%l1b_nXTrack = tmp_dims(1)
 Proxy%l1b_nLines  = tmp_dims(2)

! Allocate memory for input variables.
ALLOCATE(Proxy%l1b_lat( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%l1b_lon( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%l1b_sza( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%l1b_vza( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%l1b_saz( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%l1b_vaz( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%l1b_dem( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
!
ALLOCATE(Proxy%Ref340( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref354( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref388( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref412( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref488( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref550( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref670( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref680( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref688( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref860( Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref1250(Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref1378(Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref1640(Proxy%l1b_nXTrack, Proxy%l1b_nLines))
ALLOCATE(Proxy%Ref2260(Proxy%l1b_nXTrack, Proxy%l1b_nLines))
!ALLOCATE(cossza(Proxy%l1b_nXTrack, Proxy%l1b_nLines))

!
!Read input variables now.
errorStatus = read_h5_2dr4(proxy_filename, path_lat1, Proxy%dims_geo, Proxy%l1b_lat)
errorStatus = read_h5_2dr4(proxy_filename, path_lon1, Proxy%dims_geo, Proxy%l1b_lon)
errorStatus = read_h5_2dr4(proxy_filename, path_sza1, Proxy%dims_geo, Proxy%l1b_sza)
errorStatus = read_h5_2dr4(proxy_filename, path_vza1, Proxy%dims_geo, Proxy%l1b_vza)
errorStatus = read_h5_2dr4(proxy_filename, path_saz1, Proxy%dims_geo, Proxy%l1b_saz)
errorStatus = read_h5_2dr4(proxy_filename, path_vaz1, Proxy%dims_geo, Proxy%l1b_vaz)
errorStatus = read_h5_2dr4(proxy_filename, path_dem1, Proxy%dims_geo, Proxy%l1b_dem)

errorStatus = read_h5_2dr4(proxy_filename, path_Ref340, Proxy%dims_geo, Proxy%Ref340 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref354, Proxy%dims_geo, Proxy%Ref354 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref388, Proxy%dims_geo, Proxy%Ref388 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref412, Proxy%dims_geo, Proxy%Ref412 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref488, Proxy%dims_geo, Proxy%Ref488 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref550, Proxy%dims_geo, Proxy%Ref550 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref670, Proxy%dims_geo, Proxy%Ref670 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref680, Proxy%dims_geo, Proxy%Ref680 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref688, Proxy%dims_geo, Proxy%Ref688 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref860, Proxy%dims_geo, Proxy%Ref860 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref1250,Proxy%dims_geo, Proxy%Ref1250 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref1378,Proxy%dims_geo, Proxy%Ref1378 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref1640,Proxy%dims_geo, Proxy%Ref1640 )
errorStatus = read_h5_2dr4(proxy_filename, path_Ref2260,Proxy%dims_geo, Proxy%Ref2260 )

CALL h5fclose_f(file_id, errorStatus)

!
! Assign all required Proxy variables to common variables. 
l1b_nXTrack = Proxy%l1b_nXTrack
l1b_nLines = Proxy%l1b_nLines
Latitude = Proxy%l1b_lat
Longitude = Proxy%l1b_lon
SolarZenithAngle = Proxy%l1b_sza
ViewingZenithAngle = Proxy%l1b_vza
SolarAzimuthAngle = Proxy%l1b_saz
ViewingAzimuthAngle = Proxy%l1b_vaz
TerrainHeight = Proxy%l1b_dem
index1 = index(proxy_filename, 'TROP-in-Viirs_V4.1_A')
l1b_date_time_str = proxy_filename(index1+20:index1+20+11)
PRINT *,'Proxy start date-time : ',l1b_date_time_str


! To make the entire array UVtoSWIR_Reflectances consistent, 
! the VIIRS L1-B data needs to be divided by cossza.
ALLOCATE(UVtoSWIR_Reflectances(UVtoSWIR_nWavel, Proxy%l1b_nXTrack, Proxy%l1b_nLines))
UVtoSWIR_Reflectances(:,:,:) = -9999.9

DO jPix = 1, l1b_nLines
    DO iPix = 1, l1b_nXTrack
     IF (Proxy%Ref340(iPix,jPix) .GT. 0) UVtoSWIR_Reflectances(1,iPix,jPix) = PI * Proxy%Ref340(iPix,jPix)
     IF (Proxy%Ref354(iPix,jPix) .GT. 0) UVtoSWIR_Reflectances(2,iPix,jPix) = PI * Proxy%Ref354(iPix,jPix)
     IF (Proxy%Ref388(iPix,jPix) .GT. 0) UVtoSWIR_Reflectances(3,iPix,jPix) = PI * Proxy%Ref388(iPix,jPix)
     IF (SolarZenithAngle(iPix,jPix) .GT. -999.0) THEN
     cossza = COS(DTOR * SolarZenithAngle(iPix,jPix))
     UVtoSWIR_Reflectances( 4,iPix,jPix) = Proxy%Ref412(iPix,jPix) /cossza
     UVtoSWIR_Reflectances( 5,iPix,jPix) = Proxy%Ref488(iPix,jPix) /cossza
     UVtoSWIR_Reflectances( 6,iPix,jPix) = Proxy%Ref550(iPix,jPix) /cossza
     UVtoSWIR_Reflectances( 7,iPix,jPix) = Proxy%Ref670(iPix,jPix) /cossza
     UVtoSWIR_Reflectances( 8,iPix,jPix) = Proxy%Ref680(iPix,jPix) /cossza
     UVtoSWIR_Reflectances( 9,iPix,jPix) = Proxy%Ref688(iPix,jPix) /cossza
     UVtoSWIR_Reflectances(10,iPix,jPix) = Proxy%Ref860(iPix,jPix) /cossza
     UVtoSWIR_Reflectances(11,iPix,jPix) = Proxy%Ref1250(iPix,jPix)/cossza
     UVtoSWIR_Reflectances(12,iPix,jPix) = Proxy%Ref1378(iPix,jPix)/cossza
     UVtoSWIR_Reflectances(13,iPix,jPix) = Proxy%Ref1640(iPix,jPix)/cossza
     UVtoSWIR_Reflectances(14,iPix,jPix) = Proxy%Ref2260(iPix,jPix)/cossza
     ENDIF
    END DO 
END DO
    
!
IF (ALLOCATED(Proxy%l1b_lat)) DEALLOCATE(Proxy%l1b_lat)
IF (ALLOCATED(Proxy%l1b_lon)) DEALLOCATE(Proxy%l1b_lon)
IF (ALLOCATED(Proxy%l1b_sza)) DEALLOCATE(Proxy%l1b_sza)
IF (ALLOCATED(Proxy%l1b_vza)) DEALLOCATE(Proxy%l1b_vza)
IF (ALLOCATED(Proxy%l1b_saz)) DEALLOCATE(Proxy%l1b_saz)
IF (ALLOCATED(Proxy%l1b_vaz)) DEALLOCATE(Proxy%l1b_vaz)
IF (ALLOCATED(Proxy%l1b_dem)) DEALLOCATE(Proxy%l1b_dem)

IF (ALLOCATED(Proxy%Ref340)) DEALLOCATE(Proxy%Ref340)
IF (ALLOCATED(Proxy%Ref354)) DEALLOCATE(Proxy%Ref354)
IF (ALLOCATED(Proxy%Ref388)) DEALLOCATE(Proxy%Ref388)
IF (ALLOCATED(Proxy%Ref412)) DEALLOCATE(Proxy%Ref412)
IF (ALLOCATED(Proxy%Ref488)) DEALLOCATE(Proxy%Ref488)
IF (ALLOCATED(Proxy%Ref550)) DEALLOCATE(Proxy%Ref550)
IF (ALLOCATED(Proxy%Ref670)) DEALLOCATE(Proxy%Ref670)
IF (ALLOCATED(Proxy%Ref680)) DEALLOCATE(Proxy%Ref680)
IF (ALLOCATED(Proxy%Ref688)) DEALLOCATE(Proxy%Ref688)
IF (ALLOCATED(Proxy%Ref860)) DEALLOCATE(Proxy%Ref860)
IF (ALLOCATED(Proxy%Ref1250)) DEALLOCATE(Proxy%Ref1250)
IF (ALLOCATED(Proxy%Ref1378)) DEALLOCATE(Proxy%Ref1378)
IF (ALLOCATED(Proxy%Ref1640)) DEALLOCATE(Proxy%Ref1640)
IF (ALLOCATED(Proxy%Ref2260)) DEALLOCATE(Proxy%Ref2260)

RETURN

END FUNCTION read_l1bproxydata
!
!===============================================================

FUNCTION read_proxydate(proxy_filename,  l1b_date_time_str,l1b_nXTrack, l1b_nLines) RESULT(status)

USE H5read_module

IMPLICIT NONE

Include 'common_l1b_var.inc'

CHARACTER(LEN=255), INTENT(IN)  :: proxy_filename

TYPE(l1b_proxy_type) :: Proxy
INTEGER(KIND=4) :: errorStatus, STATUS, tmp_rank, index1
INTEGER(KIND=8),dimension(7) :: tmp_dims
INTEGER(HID_T) :: file_id, group_id, dataset_id, datatype_id, attribute_id
CHARACTER(LEN=128) :: str1

status =1

! Open HDF-5 files
 CALL h5fopen_f(proxy_filename, H5F_ACC_RDONLY_F, file_id, errorStatus)

 CALL get_rank_dims(file_id, path_lat1, tmp_rank, tmp_dims)
 PRINT*,'Proxydata file Lat rank:',tmp_rank
 PRINT*,'Proxydata file Lat shape:',tmp_dims(1:tmp_rank)

 l1b_nXTrack = tmp_dims(1)
 l1b_nLines = tmp_dims(2)

 index1 = index(proxy_filename, 'TROP-in-Viirs_V4.1_A')
 l1b_date_time_str = proxy_filename(index1+20:index1+20+11)
 PRINT *,'Proxy start date-time : ',l1b_date_time_str

RETURN

END FUNCTION read_proxydate
!
!===============================================================
!===============================================================

End Module proxydata_L1BModule
!
