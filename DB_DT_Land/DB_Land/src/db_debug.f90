module db_debug

implicit none

public   :: db_debug_output
public   :: pace_test_output

type :: dataset
    character(len=255)                    ::  dset_name
    real, dimension(:,:), allocatable     ::  values
end type dataset
  
type :: output_file
    character(len=255)                        ::  name
    type(dataset), dimension(:), allocatable  ::  datasets
end type output_file

contains

integer function db_debug_output(of)  result(status)
      
  implicit none
 
 	include 'hdf.f90' 

  type(output_file), intent(in)           ::  of
  
  integer, parameter                      ::  SHORT = selected_int_kind(4)
  
  character(len=255)                      ::  dset_name
  character(len=255)                      ::  attr_name
  integer                                 ::  sd_id
  integer                                 ::  sds_id
  integer, dimension(2)                   ::  dims2, chunk_dims2
  integer, dimension(2)                   ::  start2, stride2, edges2
  
  integer, dimension(4)										::  comp_params
  real                                    ::  min_val, max_val
  
  integer   ::  sfstart, sfcreate, sfwdata, sfendacc, sfend, sfschnk
  
  integer   ::  i

  type(dataset) ::  ds
  integer(kind=SHORT), dimension(:,:), allocatable  ::  scaled
  real                                              ::  scale_factor
  
  print *, "output_file: ", trim(of%name)

	chunk_dims2 = (/128,128/)
	comp_params(:) = 0
	comp_params(1) = 4
	
  sd_id = sfstart(trim(of%name), DFACC_CREATE)
  if (sd_id == -1) then 
    print *, "ERROR: Unable to create VIIRS DB output file: ", sd_id
    print *, "File: ", trim(of%name)
    return
  end if
  
  do i = 1, size(of%datasets)
  	ds = of%datasets(i)
    
    dset_name = trim(ds%dset_name)
    
    dims2 = (/size(ds%values,1), size(ds%values,2)/)
    sds_id = sfcreate(sd_id, trim(dset_name), 5, 2, dims2)
    if (sds_id < 0) then
      print *, "ERROR: Failed to make dataset "//trim(dset_name)//": ", status
      return
    end if
  
  	status = sfschnk(sds_id, chunk_dims2, COMP_CODE_DEFLATE, comp_params)
		if (status /= 0) then
			print *, "ERROR: Failed to set up compression on SDS "//trim(dset_name)//": ", status
			!stop
		end if
  
    start2 = (/0, 0/)
    stride2 = (/1, 1/)
    edges2 = dims2
    status = sfwdata(sds_id, start2, stride2, edges2, ds%values)
    if (status == -1) then
      print *, "ERROR: Failed to write dataset "//trim(dset_name)//": ", status
      return
    end if
  
    status = sfendacc(sds_id)
    if (status /= 0) then
      print *, "ERROR: Failed to close dataset "//trim(dset_name)//": ", status
      return
    end if
  
  end do
  
  status = sfend(sd_id)
  if (status /= 0) then
     print *,"WARNING: Failed to close output file: ", status
  end if

end function db_debug_output


integer function pace_test_output(of)  result(status)
      
  implicit none
 
 	include 'hdf.f90' 

  type(output_file), intent(in)           ::  of
  
  integer, parameter                      ::  SHORT = selected_int_kind(4)
  
  character(len=255)                      ::  dset_name
  character(len=255)                      ::  attr_name
  integer                                 ::  sd_id
  integer                                 ::  sds_id
  integer, dimension(2)                   ::  dims2, chunk_dims2
  integer, dimension(2)                   ::  start2, stride2, edges2
  
  integer, dimension(4)										::  comp_params
  real                                    ::  min_val, max_val
  
  integer   ::  sfstart, sfcreate, sfwdata, sfendacc, sfend, sfschnk
  
  integer   ::  i

  type(dataset) ::  ds
  integer(kind=SHORT), dimension(:,:), allocatable  ::  scaled
  real                                              ::  scale_factor
  
  print *, "output_file: ", trim(of%name)

	chunk_dims2 = (/128,128/)
	comp_params(:) = 0
	comp_params(1) = 4
	
  sd_id = sfstart(trim(of%name), DFACC_CREATE)
  if (sd_id == -1) then 
    print *, "ERROR: Unable to create VIIRS DB output file: ", sd_id
    print *, "File: ", trim(of%name)
    return
  end if
  
  do i = 1, size(of%datasets)
  	ds = of%datasets(i)
    
    dset_name = trim(ds%dset_name)
    
    dims2 = (/size(ds%values,1), size(ds%values,2)/)
    sds_id = sfcreate(sd_id, trim(dset_name), 5, 2, dims2)
    if (sds_id < 0) then
      print *, "ERROR: Failed to make dataset "//trim(dset_name)//": ", status
      return
    end if
  
  	status = sfschnk(sds_id, chunk_dims2, COMP_CODE_DEFLATE, comp_params)
		if (status /= 0) then
			print *, "ERROR: Failed to set up compression on SDS "//trim(dset_name)//": ", status
			!stop
		end if
  
    start2 = (/0, 0/)
    stride2 = (/1, 1/)
    edges2 = dims2
    status = sfwdata(sds_id, start2, stride2, edges2, ds%values)
    if (status == -1) then
      print *, "ERROR: Failed to write dataset "//trim(dset_name)//": ", status
      return
    end if
  
    status = sfendacc(sds_id)
    if (status /= 0) then
      print *, "ERROR: Failed to close dataset "//trim(dset_name)//": ", status
      return
    end if
  
  end do

  sds_id = sfcreate(sd_id, trim('nXTrack'), 5, 2, 1)
  if (sds_id < 0) then
    print *, "ERROR: Failed to make dataset "//trim(dset_name)//": ", status
    return
  end if

  status = sfschnk(sds_id, chunk_dims2, COMP_CODE_DEFLATE, comp_params)
  if (status /= 0) then
    print *, "ERROR: Failed to set up compression on SDS "//trim(dset_name)//": ", status
    !stop
  end if

  status = sfwdata(sds_id, 0, 1, 1, 1271)
  if (status == -1) then
    print *, "ERROR: Failed to write dataset "//trim(dset_name)//": ", status
    return
  end if

  status = sfendacc(sds_id)
  if (status /= 0) then
    print *, "ERROR: Failed to close dataset "//trim(dset_name)//": ", status
    return
  end if 

  sds_id = sfcreate(sd_id, trim('nLines'), 5, 2, 1)
  if (sds_id < 0) then
    print *, "ERROR: Failed to make dataset "//trim(dset_name)//": ", status
    return
  end if

  status = sfschnk(sds_id, chunk_dims2, COMP_CODE_DEFLATE, comp_params)
  if (status /= 0) then
    print *, "ERROR: Failed to set up compression on SDS "//trim(dset_name)//": ", status
    !stop
  end if

  status = sfwdata(sds_id, 0, 1, 1, 1722)
  if (status == -1) then
    print *, "ERROR: Failed to write dataset "//trim(dset_name)//": ", status
    return
  end if

  status = sfendacc(sds_id)
  if (status /= 0) then
    print *, "ERROR: Failed to close dataset "//trim(dset_name)//": ", status
    return
  end if 
  
 
  
  status = sfend(sd_id)
  if (status /= 0) then
     print *,"WARNING: Failed to close output file: ", status
  end if

end function pace_test_output


end module db_debug