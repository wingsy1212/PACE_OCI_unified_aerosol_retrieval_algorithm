
  if (PRESENT(this_kind)) then
    if (this_kind == 1) then 
      memtypeId = H5T_STD_I8LE
    elseif (this_kind ==2) then 
      memtypeId = H5T_STD_I16LE
    else
      memtypeId = H5T_NATIVE_INTEGER
    endif
  else 
    memtypeId = H5T_NATIVE_INTEGER
  endif

