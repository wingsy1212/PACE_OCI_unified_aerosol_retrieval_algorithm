status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ",ds%name
   call Display_Message(ROUTINE_NAME, msg, status) 
   RETURN  
ENDIF

rankm = RANK( data_in )
dimsm = SHAPE( data_in )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call h5screate_simple_f(rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, memtypeId, data_in, dimsm, error, &
                 memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name 
   call Display_Message(ROUTINE_NAME, msg, status) 
   RETURN  
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)
