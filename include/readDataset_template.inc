status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name 
   call Display_Message(ROUTINE_NAME, msg, status) 
   RETURN  
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
IF(typeClass /=  H5T_INTEGER_F) then
   status = FAILURE_STATE
   write( msg,* ) "class type should be H5T_INTEGER_F: ", ds%name 
   call Display_Message(ROUTINE_NAME, msg, status) 
   RETURN  
ENDIF

rankm = RANK( data_out ) 
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name 
   call Display_Message(ROUTINE_NAME, msg, status) 
   RETURN  
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

