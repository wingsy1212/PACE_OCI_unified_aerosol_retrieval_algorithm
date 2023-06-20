MODULE H5read_module

  USE HDF5

  PUBLIC :: get_rank_dims
  PUBLIC :: read_h5_1dr4
  PUBLIC :: read_h5_2dr4
  PUBLIC :: read_h5_3dr4
  PUBLIC :: read_h5_4dr4
  PUBLIC :: read_2di4
  PUBLIC :: read_4di4
  PUBLIC :: read_3di4

  CONTAINS
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
SUBROUTINE get_rank_dims(file_id, input_path, rank_out, dims_out)

  USE H5Util_class

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: input_path
  INTEGER(HID_T)  , INTENT(IN) :: file_id
  INTEGER, INTENT(OUT) :: rank_out
  INTEGER(KIND=8), DIMENSION(7), INTENT(OUT) :: dims_out
  TYPE(H5SDS_T) :: sds1
  INTEGER :: errorStatus

  dims_out(:)=0
  sds1 = H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),-999.0)
  sds1%name = input_path
  errorStatus = H5Util_selectDataset(file_id, sds1)

  rank_out = sds1%rank
  dims_out(1:rank_out) = sds1%dims(1:rank_out)
  errorStatus = H5Util_disposeDataset(sds1)

END SUBROUTINE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
FUNCTION read_h5_1dr4(l1b_file, dataset_name, dims1_arr, data_values) RESULT(STATUS)

  CHARACTER(LEN=*), INTENT(IN):: dataset_name , l1b_file
  INTEGER(KIND=4), DIMENSION(1), INTENT(IN) :: dims1_arr
  INTEGER(KIND=4) :: STATUS, hdf_err
  INTEGER(HSIZE_T), DIMENSION(1) :: dims1
  INTEGER(HID_T) :: file_id, group_id, dataset_id, dataTYPE_id, attribute_id
  REAL(KIND=4),DIMENSION(dims1_arr(1)) :: data_values

  STATUS=1
!
  CALL H5Fopen_f(l1b_file, H5F_ACC_RDONLY_F, file_id, hdf_err)
  CALL h5dopen_f(file_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) THEN
    print *,'Unable to open file :',l1b_file
    GO TO 100
  ENDIF
!
  CALL h5dget_TYPE_f(dataset_id, dataTYPE_id, hdf_err)
  IF (hdf_err .NE. 0) THEN 
    print *,'Unable to open dataset : ',dataset_name
    CALL h5fclose_f(file_id, hdf_err)
    GO TO 100
  ENDIF
!
  dims1 = SHAPE(data_values)
  CALL h5dread_f(dataset_id, dataTYPE_id, data_values, dims1, hdf_err)
  IF (hdf_err .NE. 0) THEN 
    print *,'Error reading dataset : ',dataset_name
    CALL h5fclose_f(file_id, hdf_err)
    GO TO 100
  ENDIF
!  PRINT *,'Successful in reading.. ', dataset_name
!
  CALL h5dclose_f(dataset_id, hdf_err)
  CALL h5fclose_f(file_id, hdf_err)
  RETURN
100 STATUS = -1

END FUNCTION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
FUNCTION read_h5_2dr4(l1b_file, dataset_name, dims2_arr, data_values) RESULT(STATUS)

  CHARACTER(LEN=*), INTENT(IN):: dataset_name , l1b_file
  INTEGER(KIND=4), DIMENSION(2), INTENT(IN) :: dims2_arr
  INTEGER(KIND=4) :: STATUS, hdf_err
  INTEGER(HSIZE_T), DIMENSION(2) :: dims2
  INTEGER(HID_T) :: file_id, group_id, dataset_id, dataTYPE_id, attribute_id
  REAL(KIND=4),DIMENSION(dims2_arr(1),dims2_arr(2)) :: data_values

  STATUS=1
!
  CALL H5Fopen_f(l1b_file, H5F_ACC_RDONLY_F, file_id, hdf_err)
  CALL h5dopen_f(file_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) THEN
    print *,'Unable to open file :',l1b_file
    GO TO 200
  ENDIF
!
  CALL h5dget_TYPE_f(dataset_id, dataTYPE_id, hdf_err)
  IF (hdf_err .NE. 0) THEN 
    print *,'Unable to open dataset : ',dataset_name
    CALL h5fclose_f(file_id, hdf_err)
    GO TO 200
  ENDIF
!
  dims2 = SHAPE(data_values)
  CALL h5dread_f(dataset_id, dataTYPE_id, data_values, dims2, hdf_err)
  IF (hdf_err .NE. 0) THEN 
    print *,'Error reading dataset : ',dataset_name
    CALL h5fclose_f(file_id, hdf_err)
    GO TO 200
  ENDIF
!  PRINT *,'Successful in reading.. ', dataset_name
!
  CALL h5dclose_f(dataset_id, hdf_err)
  CALL h5fclose_f(file_id, hdf_err)
  RETURN
200 STATUS = -1

END FUNCTION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
FUNCTION read_h5_3dr4(l1b_file, dataset_name, dims3_arr, data_values) RESULT(STATUS)

  CHARACTER(LEN=*), INTENT(IN):: dataset_name , l1b_file
  INTEGER(KIND=4), DIMENSION(3), INTENT(IN) :: dims3_arr
  INTEGER(HSIZE_T), DIMENSION(3) :: dims3
  INTEGER(KIND=4) :: STATUS, hdf_err
  INTEGER(HID_T) :: file_id, group_id, dataset_id, dataTYPE_id, attribute_id
  REAL(KIND=4),DIMENSION(dims3_arr(1),dims3_arr(2),dims3_arr(3)) :: data_values
  STATUS=1
!
  CALL H5Fopen_f(l1b_file, H5F_ACC_RDONLY_F, file_id, hdf_err)
  CALL h5dopen_f(file_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) THEN
    print *,'Unable to open file :',l1b_file
    GO TO 300
  ENDIF
!
  CALL h5dget_TYPE_f(dataset_id, dataTYPE_id, hdf_err)
  IF (hdf_err .NE. 0) THEN 
    print *,'Unable to open dataset : ',dataset_name
    CALL h5fclose_f(file_id, hdf_err)
    GO TO 300
  ENDIF
!
  dims3 = SHAPE(data_values)
  CALL h5dread_f(dataset_id, dataTYPE_id, data_values, dims3, hdf_err)
  IF (hdf_err .NE. 0) THEN 
    print *,'Error reading dataset : ',dataset_name
    CALL h5fclose_f(file_id, hdf_err)
    GO TO 300
  ENDIF
!  PRINT *,'Successful in reading.. ', dataset_name
!
  CALL h5dclose_f(dataset_id, hdf_err)
  CALL h5fclose_f(file_id, hdf_err)
  RETURN
300 STATUS = -1

END FUNCTION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
FUNCTION read_h5_4dr4(l1b_file, dataset_name, dims4_arr, data_values) RESULT(STATUS)

  CHARACTER(LEN=*), INTENT(IN):: dataset_name , l1b_file
  INTEGER(KIND=4), DIMENSION(4), INTENT(IN) :: dims4_arr
  INTEGER(HSIZE_T), DIMENSION(4) :: dims4
  INTEGER(KIND=4) :: STATUS, hdf_err
  INTEGER(HID_T) :: file_id, group_id, dataset_id, dataTYPE_id, attribute_id
  REAL(KIND=4),DIMENSION(dims4_arr(1),dims4_arr(2),dims4_arr(3), dims4_arr(4)) :: data_values
  STATUS=1
!
  CALL H5Fopen_f(l1b_file, H5F_ACC_RDONLY_F, file_id, hdf_err)
  CALL h5dopen_f(file_id, dataset_name, dataset_id, hdf_err)
  IF (hdf_err .NE. 0) THEN
    print *,'Unable to open file :',l1b_file
    GO TO 400
  ENDIF
!
  CALL h5dget_TYPE_f(dataset_id, dataTYPE_id, hdf_err)
  IF (hdf_err .NE. 0) THEN 
    print *,'Unable to open dataset : ',dataset_name
    CALL h5fclose_f(file_id, hdf_err)
    GO TO 400
  ENDIF
!
  dims4 = SHAPE(data_values)
  CALL h5dread_f(dataset_id, dataTYPE_id, data_values, dims4, hdf_err)
  IF (hdf_err .NE. 0) THEN 
    print *,'Error reading dataset : ',dataset_name
    CALL h5fclose_f(file_id, hdf_err)
    GO TO 400
  ENDIF
  PRINT *,data_values(1,1,1,1)
!
  CALL h5dclose_f(dataset_id, hdf_err)
  CALL h5fclose_f(file_id, hdf_err)
  RETURN
400  STATUS=-1

END FUNCTION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
FUNCTION read_2di4(l1b_file, path_field, dims, data_arr2, units_str) RESULT(STATUS)

USE DataTypeDef
USE MyConstants
USE HDF5
USE H5Util_class

  CHARACTER(LEN=*), INTENT(IN)::  l1b_file, path_field
  INTEGER(KIND=4), DIMENSION(2), INTENT(IN) :: dims
  INTEGER(KIND=4) :: STATUS, hdferr
  TYPE(H5SDS_T)  :: sds1 = & 
       H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)

  INTEGER(KIND=4) :: hdf_err
  INTEGER(HID_T) :: file_id, group_id
  INTEGER(I4B), DIMENSION(dims(1), dims(2)), INTENT(OUT) :: data_arr2
  CHARACTER(LEN=128), INTENT(OUT), optional ::  units_str
  STATUS=1

  CALL h5fopen_f(l1b_file, H5F_ACC_RDONLY_F, file_id, hdf_err)

  sds1%name = path_field
  hdferr = H5Util_selectDataset(file_id, sds1)
  if( sds1%dims(1) /= dims(1) .or. sds1%dims(2) /= dims(2) ) THEN
    PRINT *,'ERROR : DIMENSIONS DOES NOT MATCH FOR ',path_field
    PRINT *,'DIMENSION OF FIELD :',sds1%dims(1:2)
    PRINT *,'INPUT DIMENSION :',dims
    STATUS = -1
  ENDif
  hdferr = H5Util_readDataset(sds1, data_arr2)
  if (present(units_str)) then
    units_str = sds1%units
  ENDif
  hdferr = h5Util_disposeDataset(sds1)

  CALL h5fclose_f(file_id, hdferr)

END FUNCTION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
FUNCTION read_3di4(l1b_file, path_field, dims, data_arr3) RESULT(STATUS)

USE DataTypeDef
USE MyConstants
USE HDF5
USE H5Util_class

CHARACTER(LEN=*), INTENT(IN)::  l1b_file, path_field
INTEGER(KIND=4), DIMENSION(3), INTENT(IN) :: dims
INTEGER(KIND=4) :: STATUS, hdferr
TYPE(H5SDS_T)  :: sds1 = & 
     H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)

  INTEGER(KIND=4) :: hdf_err
  INTEGER(HID_T) :: file_id, group_id
  INTEGER(I4B), DIMENSION(dims(1), dims(2),dims(3)), & 
                INTENT(OUT) :: data_arr3

  STATUS=1
  CALL h5fopen_f(l1b_file, H5F_ACC_RDONLY_F, file_id, hdf_err)

  sds1%name = path_field
  hdferr = H5Util_selectDataset(file_id, sds1)
  if( sds1%dims(1) /= dims(1) .or. sds1%dims(2) /= dims(2)  .or. &
      sds1%dims(3) /= dims(3)  ) THEN
    PRINT *,'ERROR : DIMENSIONS DOES NOT MATCH FOR ',path_field
    PRINT *,'DIMENSION OF FIELD :',sds1%dims(1:3)
    PRINT *,'INPUT DIMENSION :',dims
    STATUS = -1
  ENDif
  hdferr = H5Util_readDataset(sds1, data_arr3)
  hdferr = h5Util_disposeDataset(sds1)

  CALL h5fclose_f(file_id, hdferr)

END FUNCTION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
FUNCTION read_4di4(l1b_file, path_field, dims, data_arr4) RESULT(STATUS)

USE DataTypeDef
USE MyConstants
USE HDF5
USE H5Util_class

CHARACTER(LEN=*), INTENT(IN)::  l1b_file, path_field
INTEGER(KIND=4), DIMENSION(4), INTENT(IN) :: dims
INTEGER(KIND=4) :: STATUS, hdferr
TYPE(H5SDS_T)  :: sds1 = & 
     H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)

  INTEGER(KIND=4) :: hdf_err
  INTEGER(HID_T) :: file_id, group_id
  INTEGER(I4B), DIMENSION(dims(1), dims(2),dims(3),dims(4)), & 
                INTENT(OUT) :: data_arr4

  STATUS=1
  CALL h5fopen_f(l1b_file, H5F_ACC_RDONLY_F, file_id, hdf_err)

  sds1%name = path_field
  hdferr = H5Util_selectDataset(file_id, sds1)
  if( sds1%dims(1) /= dims(1) .or. sds1%dims(2) /= dims(2)  .or. &
      sds1%dims(3) /= dims(3) .or. sds1%dims(4) /= dims(4) ) THEN
    PRINT *,'ERROR : DIMENSIONS DOES NOT MATCH FOR ',path_field
    PRINT *,'DIMENSION OF FIELD :',sds1%dims(1:4)
    PRINT *,'INPUT DIMENSION :',dims
    STATUS = -1
  ENDif
  hdferr = H5Util_readDataset(sds1, data_arr4)
  hdferr = h5Util_disposeDataset(sds1)

  CALL h5fclose_f(file_id, hdferr)

END FUNCTION
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
END MODULE
