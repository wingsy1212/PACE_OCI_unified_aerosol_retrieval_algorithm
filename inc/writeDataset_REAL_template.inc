! IF( ds%datatype_id /= H5T_NATIVE_REAL) then
! or 
! IF( ds%datatype_id /= H5T_NATIVE_DOUBLE) then

call h5tequal_f(ds%datatype_id, memtypeId, sameType, error)
if(.not. sameType) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL / H5T_NATIVE_DOUBLE: ", ds%name 
   call Display_Message(ROUTINE_NAME, msg, status) 
   RETURN  
ENDIF

include 'writeDataset_template.inc'
