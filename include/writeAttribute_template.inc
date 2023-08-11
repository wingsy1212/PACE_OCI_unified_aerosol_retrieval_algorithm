errorFlag = 0 

buffer = attrValue

!
! datatype in disk file:
!
call H5Screate_f(H5S_SCALAR_F, dataSpaceId, errorFlag)
if(present(h5ObjectDataTypeId)) then
   call H5Tcopy_f(h5ObjectDataTypeId, dataTypeIdOnDisk, errorFlag)
else 
   call H5Tcopy_f(defaultDatatypeId, dataTypeIdOnDisk, errorFlag)
endif
call H5Acreate_f(h5ObjectId, attrName, &
                 dataTypeIdOnDisk, dataSpaceId, attrId,errorFlag)


!
! datatype in memory:
!
call H5Tcopy_f(defaultDatatypeId, dataTypeIdInMemory, errorFlag)
call H5Awrite_f(attrId, datatypeIdInMemory, buffer, dims, errorFlag)
if(errorFlag /= 0) then 
   errorFlag = FAILURE_STATE
   call Display_Message(ROUTINE_NAME, &
               'error writing attribute ' // trim(attrName), errorFlag)
   goto 9999
endif

call H5Aclose_f(attrId, errorFlag)
call H5Tclose_f(datatypeIdOnDisk, errorFlag)
call H5Tclose_f(datatypeIdInMemory, errorFlag)
call H5Sclose_f(dataSpaceId, errorFlag)

9999 return

