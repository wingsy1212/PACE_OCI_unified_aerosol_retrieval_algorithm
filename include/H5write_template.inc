
	
  PRINT *,"Writing field :",field_path

  if ( .NOT. PRESENT(field_sds)) THEN
    arr_shape = shape(field_values)
    sds1 = sds_const
    sds1%name = field_path
    sds1%datatype_id = memtypeId
    sds1%rank = rank
    sds1%dims(1:sds1%rank) = arr_shape(1:sds1%rank)
    shall_write_attr = .FALSE.
  else
    sds1 = field_sds
  endif

  if(PRESENT(dcpl)) then 
    errorStatus = H5Util_createDataset(group_id, sds1, shall_write_attr, dcpl)
  else
    errorStatus = H5Util_createDataset(group_id, sds1, shall_write_attr )
  endif

  errorStatus = H5Util_writeDataset(sds1, field_values)
  errorstatus = h5util_disposeDataset(sds1)
