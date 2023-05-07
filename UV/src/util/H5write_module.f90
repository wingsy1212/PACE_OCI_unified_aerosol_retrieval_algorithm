module H5write_module
  USE HDF5
  USE H5Util_class
  USE DataTypeDef
  USE MyConstants
  implicit none

  INTEGER(I4B) :: errorStatus, hdferr
  TYPE(H5SDS_T) :: sds1
  TYPE(H5SDS_T), PARAMETER :: sds_const = & 
      H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),-999.0)
  CHARACTER(len=256) :: msg


  public :: H5write_2DRK4
  public :: H5write_3DRK4
  public :: H5write_4DRK4
  public :: H5write_2DIK1
  public :: H5write_2DIK2
  public :: H5write_3DIK2
  public :: H5write_1DIK4
  public :: H5write_2DIK4
  public :: H5write_3DIK4
  public :: H5write_4DIK4
  public :: H5write_2DUI8
  public :: H5write_3DUI8
  public :: H5write_4DUI8
  public :: H5write_2DUI16
  public :: H5write_3DUI16
  public :: H5write_4DUI16
  public :: H5write_utc_time
  public :: H5write_data

  
  private ::  errorStatus, hdferr, msg
  private ::  sds1, sds_const

!function overloading
interface H5write_data
  module procedure H5write_1DRK4
  module procedure H5write_2DRK4
  module procedure H5write_3DRK4
  module procedure H5write_4DRK4
  module procedure H5write_5DRK4
  module procedure H5write_6DRK4
  module procedure H5write_7DRK4

  module procedure H5write_1DRK8
  module procedure H5write_2DRK8
  module procedure H5write_3DRK8
  module procedure H5write_4DRK8
  module procedure H5write_5DRK8
  module procedure H5write_6DRK8
  module procedure H5write_7DRK8

  module procedure H5write_1DIK4
  module procedure H5write_2DIK4
  module procedure H5write_3DIK4
  module procedure H5write_4DIK4
  module procedure H5write_5DIK4
  module procedure H5write_6DIK4
  module procedure H5write_7DIK4


end interface

interface H5write_UNSIGNED_data
  module procedure H5write_1DUI16
  module procedure H5write_2DUI16
  module procedure H5write_3DUI16
  module procedure H5write_4DUI16

  module procedure H5write_1DUI8
  module procedure H5write_2DUI8
  module procedure H5write_3DUI8
  module procedure H5write_4DUI8
end interface


  contains


subroutine H5write_1DRK4(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=4), dimension(:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(1) :: arr_shape
  integer , PARAMETER :: rank=1
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_1DRK4"

  memtypeId = H5T_NATIVE_REAL

  include 'H5write_template.inc'

end subroutine H5write_1DRK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_2DRK4(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=4), dimension(:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(2) :: arr_shape
  integer , PARAMETER :: rank=2
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_2DRK4"

  memtypeId = H5T_NATIVE_REAL

  include 'H5write_template.inc'

end subroutine H5write_2DRK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_3DRK4(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=4), dimension(:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(3) :: arr_shape
  integer , PARAMETER :: rank=3
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK4"

  memtypeId = H5T_NATIVE_REAL

  include 'H5write_template.inc'

end subroutine H5write_3DRK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_4DRK4(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=4), dimension(:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(4) :: arr_shape
  integer , PARAMETER :: rank=4
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_4DRK4"

  memtypeId = H5T_NATIVE_REAL

  include 'H5write_template.inc'

end subroutine H5write_4DRK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_5DRK4(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=4), dimension(:,:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(5) :: arr_shape
  integer , PARAMETER :: rank=5
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_5DRK4"

  memtypeId = H5T_NATIVE_REAL

  include 'H5write_template.inc'

end subroutine H5write_5DRK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_6DRK4(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=4), dimension(:,:,:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(6) :: arr_shape
  integer , PARAMETER :: rank=6
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_6DRK4"

  memtypeId = H5T_NATIVE_REAL

  include 'H5write_template.inc'

end subroutine H5write_6DRK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_7DRK4(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=4), dimension(:,:,:,:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(7) :: arr_shape
  integer , PARAMETER :: rank=7
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_7DRK4"

  memtypeId = H5T_NATIVE_REAL

  include 'H5write_template.inc'

end subroutine H5write_7DRK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!               Real Kind 8
!-----------------------------------------------------------
subroutine H5write_1DRK8(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=8), dimension(:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(1) :: arr_shape
  integer , PARAMETER :: rank=1
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_1DRK8"

  memtypeId = H5T_NATIVE_DOUBLE

  include 'H5write_template.inc'

end subroutine H5write_1DRK8

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_2DRK8(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=8), dimension(:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(2) :: arr_shape
  integer , PARAMETER :: rank=2
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_2DRK8"

  memtypeId = H5T_NATIVE_DOUBLE

  include 'H5write_template.inc'

end subroutine H5write_2DRK8

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_3DRK8(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=8), dimension(:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(3) :: arr_shape
  integer , PARAMETER :: rank=3
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK8"


  memtypeId = H5T_NATIVE_DOUBLE

  !RAM TODO include 'H5write_template.inc'
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


end subroutine H5write_3DRK8

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_4DRK8(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=8), dimension(:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(4) :: arr_shape
  integer , PARAMETER :: rank=4
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_4DRK8"

  memtypeId = H5T_NATIVE_DOUBLE

  include 'H5write_template.inc'

end subroutine H5write_4DRK8

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_5DRK8(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=8), dimension(:,:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(5) :: arr_shape
  integer , PARAMETER :: rank=5
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_5DRK8"

  memtypeId = H5T_NATIVE_DOUBLE

  include 'H5write_template.inc'

end subroutine H5write_5DRK8

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_6DRK8(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=8), dimension(:,:,:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(6) :: arr_shape
  integer , PARAMETER :: rank=6
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_6DRK8"

  memtypeId = H5T_NATIVE_DOUBLE

  include 'H5write_template.inc'

end subroutine H5write_6DRK8

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_7DRK8(group_id, field_path, field_values, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  real(kind=8), dimension(:,:,:,:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(7) :: arr_shape
  integer , PARAMETER :: rank=7
  INTEGER(HID_T) :: memtypeId
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_7DRK8"

  memtypeId = H5T_NATIVE_DOUBLE

  include 'H5write_template.inc'

end subroutine H5write_7DRK8

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_2DIK1(group_id, field_path, field_values )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=1), dimension(:,:), intent(in) :: field_values
  integer(HSIZE_T),dimension(2) :: arr_dims
  integer , PARAMETER :: rank=2
  INTEGER(HID_T) :: memtypeId, dspace_id, dset_id

  memtypeId = H5T_STD_U8LE
  arr_dims = shape(field_values)
  CALL h5screate_simple_f(rank, arr_dims, dspace_id, hdferr)
  CALL h5dcreate_f(group_id, field_path, memtypeId, &
    dspace_id, dset_id, hdferr)
  CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,INT(field_values,4), &
    arr_dims, hdferr)
  !Close Datasets
  CALL h5dclose_f(dset_id , hdferr)
  !Close Dataspace
  CALL h5sclose_f(dspace_id, hdferr)


end subroutine H5write_2DIK1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_2DIK2(group_id, field_path, field_values )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=2), dimension(:,:), intent(in) :: field_values
  integer(HSIZE_T),dimension(2) :: arr_dims
  integer , PARAMETER :: rank=2
  INTEGER(HID_T) :: memtypeId, dspace_id, dset_id

  memtypeId = H5T_STD_U16LE
  arr_dims = shape(field_values)
  CALL h5screate_simple_f(rank, arr_dims, dspace_id, hdferr)
  CALL h5dcreate_f(group_id, field_path, memtypeId, &
    dspace_id, dset_id, hdferr)
  CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,INT(field_values,4), &
    arr_dims, hdferr)
  !Close Datasets
  CALL h5dclose_f(dset_id , hdferr)
  !Close Dataspace
  CALL h5sclose_f(dspace_id, hdferr)


end subroutine H5write_2DIK2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_3DIK2(group_id, field_path, field_values )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=2), dimension(:,:,:), intent(in) :: field_values
  integer(HSIZE_T),dimension(3) :: arr_dims
  integer , PARAMETER :: rank=3
  INTEGER(HID_T) :: memtypeId, dspace_id, dset_id

  memtypeId = H5T_STD_U16LE
  arr_dims = shape(field_values)
  CALL h5screate_simple_f(rank, arr_dims, dspace_id, hdferr)
  CALL h5dcreate_f(group_id, field_path, memtypeId, &
    dspace_id, dset_id, hdferr)
  CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,INT(field_values,4), &
    arr_dims, hdferr)
  !Close Datasets
  CALL h5dclose_f(dset_id , hdferr)
  !Close Dataspace
  CALL h5sclose_f(dspace_id, hdferr)


end subroutine H5write_3DIK2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_1DIK4(group_id, field_path, field_values, this_kind, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=4), dimension(:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(1) :: arr_shape
  integer , PARAMETER :: rank=1
  INTEGER(HID_T) :: memtypeId
  INTEGER ,intent(in), optional :: this_kind
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_1DIK4"

  memtypeId = H5T_NATIVE_INTEGER

  include 'H5write_INT_template.inc'
  include 'H5write_template.inc'

end subroutine H5write_1DIK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_2DIK4(group_id, field_path, field_values, this_kind, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=4), dimension(:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(2) :: arr_shape
  integer , PARAMETER :: rank=2
  INTEGER(HID_T) :: memtypeId
  INTEGER ,intent(in), optional :: this_kind
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_2DIK4"

  include 'H5write_INT_template.inc'
  include 'H5write_template.inc'

end subroutine H5write_2DIK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_3DIK4(group_id, field_path, field_values, this_kind, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=4), dimension(:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(3) :: arr_shape
  integer , PARAMETER :: rank=3
  INTEGER(HID_T) :: memtypeId
  INTEGER ,intent(in), optional :: this_kind
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DIK4"

  memtypeId = H5T_NATIVE_INTEGER

  include 'H5write_INT_template.inc'
  include 'H5write_template.inc'

end subroutine H5write_3DIK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_4DIK4(group_id, field_path, field_values, this_kind, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=4), dimension(:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(4) :: arr_shape
  integer , PARAMETER :: rank=4
  INTEGER(HID_T) :: memtypeId
  INTEGER ,intent(in), optional :: this_kind
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_4DIK4"

  memtypeId = H5T_NATIVE_INTEGER

  include 'H5write_INT_template.inc'
  include 'H5write_template.inc'

end subroutine H5write_4DIK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_5DIK4(group_id, field_path, field_values, this_kind, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=4), dimension(:,:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(5) :: arr_shape
  integer , PARAMETER :: rank=5
  INTEGER(HID_T) :: memtypeId
  INTEGER ,intent(in), optional :: this_kind
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_5DIK4"

  memtypeId = H5T_NATIVE_INTEGER

  include 'H5write_INT_template.inc'
  include 'H5write_template.inc'

end subroutine H5write_5DIK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_6DIK4(group_id, field_path, field_values, this_kind, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=4), dimension(:,:,:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(6) :: arr_shape
  integer , PARAMETER :: rank=6
  INTEGER(HID_T) :: memtypeId
  INTEGER ,intent(in), optional :: this_kind
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_6DIK4"

  memtypeId = H5T_NATIVE_INTEGER

  include 'H5write_INT_template.inc'
  include 'H5write_template.inc'

end subroutine H5write_6DIK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_7DIK4(group_id, field_path, field_values, this_kind, dcpl, field_sds )
  implicit none
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=4), dimension(:,:,:,:,:,:,:), intent(in) :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(7) :: arr_shape
  integer , PARAMETER :: rank=7
  INTEGER(HID_T) :: memtypeId
  INTEGER ,intent(in), optional :: this_kind
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_7DIK4"

  memtypeId = H5T_NATIVE_INTEGER

  include 'H5write_INT_template.inc'
  include 'H5write_template.inc'

end subroutine H5write_7DIK4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!           UNSIGNED INTEGERS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine H5write_1DUI8(group_id, field_path, input_values, dcpl, field_sds )
  implicit none
  type(H5SDS_T) :: this_sds = &
    H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=1), dimension(:), intent(in) :: input_values
  integer(kind=4), dimension(:),allocatable :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(1) :: arr_shape
  integer , PARAMETER :: rank=1
  INTEGER(HID_T) :: memtypeId
  integer, PARAMETER :: fill_val=255
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK8"

  memtypeId = H5T_STD_U8LE
  arr_shape = shape(input_values)
  allocate( field_values(arr_shape(1) ))
  field_values = input_values
  where (field_values == -1) field_values = fill_val

  !write temporary array
  include 'H5write_template.inc'
  
  !de-allocate temporary array
  if(allocated(field_values)) deallocate(field_values)

end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_2DUI8(group_id, field_path, input_values, dcpl, field_sds )
  implicit none
  type(H5SDS_T) :: this_sds = &
    H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=1), dimension(:,:), intent(in) :: input_values
  integer(kind=4), dimension(:,:),allocatable :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(2) :: arr_shape
  integer , PARAMETER :: rank=2
  INTEGER(HID_T) :: memtypeId
  integer, PARAMETER :: fill_val=255
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK8"

  memtypeId = H5T_STD_U8LE
  arr_shape = shape(input_values)
  allocate( field_values(arr_shape(1), arr_shape(2) ))
  field_values = input_values
  where (field_values == -1) field_values = fill_val

  !write temporary array
  include 'H5write_template.inc'
  
  !de-allocate temporary array
  if(allocated(field_values)) deallocate(field_values)

end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_3DUI8(group_id, field_path, input_values, dcpl, field_sds )
  implicit none
  type(H5SDS_T) :: this_sds = &
    H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=1), dimension(:,:,:), intent(in) :: input_values
  integer(kind=4), dimension(:,:,:),allocatable :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(3) :: arr_shape
  integer , PARAMETER :: rank=3
  INTEGER(HID_T) :: memtypeId
  integer, PARAMETER :: fill_val=255
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK8"

  memtypeId = H5T_STD_U8LE
  arr_shape = shape(input_values)
  allocate( field_values(arr_shape(1), arr_shape(2), arr_shape(3) ))
  field_values = input_values
  where (field_values == -1) field_values = fill_val

  !write temporary array
  include 'H5write_template.inc'
  
  !de-allocate temporary array
  if(allocated(field_values)) deallocate(field_values)

end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_4DUI8(group_id, field_path, input_values, dcpl, field_sds )
  implicit none
  type(H5SDS_T) :: this_sds = &
    H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=1), dimension(:,:,:,:), intent(in) :: input_values
  integer(kind=4), dimension(:,:,:,:),allocatable :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(4) :: arr_shape
  integer , PARAMETER :: rank=4
  INTEGER(HID_T) :: memtypeId
  integer, PARAMETER :: fill_val=255
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK8"

  memtypeId = H5T_STD_U8LE
  arr_shape = shape(input_values)
  allocate( field_values(arr_shape(1), arr_shape(2), arr_shape(3), &
    arr_shape(4) ) )
  field_values = input_values
  where (field_values == -1) field_values = fill_val

  !write temporary array
  include 'H5write_template.inc'
  
  !de-allocate temporary array
  if(allocated(field_values)) deallocate(field_values)

end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ------16-bit unsigned integers ------
!
subroutine H5write_1DUI16(group_id, field_path, input_values, dcpl, field_sds )
  implicit none
  type(H5SDS_T) :: this_sds = &
    H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=2), dimension(:), intent(in) :: input_values
  integer(kind=4), dimension(:),allocatable :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(1) :: arr_shape
  integer , PARAMETER :: rank=1
  INTEGER(HID_T) :: memtypeId
  integer, PARAMETER :: fill_val=65535
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK8"

  memtypeId = H5T_STD_U16LE
  arr_shape = shape(input_values)
  allocate( field_values(arr_shape(1) ))
  field_values = input_values
  where (field_values == -1) field_values = fill_val

  !write temporary array
  include 'H5write_template.inc'
  
  !de-allocate temporary array
  if(allocated(field_values)) deallocate(field_values)

end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_2DUI16(group_id, field_path, input_values, dcpl, field_sds )
  implicit none
  type(H5SDS_T) :: this_sds = &
    H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=2), dimension(:,:), intent(in) :: input_values
  integer(kind=4), dimension(:,:),allocatable :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(2) :: arr_shape
  integer , PARAMETER :: rank=2
  INTEGER(HID_T) :: memtypeId
  integer, PARAMETER :: fill_val=65535
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK8"

  memtypeId = H5T_STD_U16LE
  arr_shape = shape(input_values)
  allocate( field_values(arr_shape(1), arr_shape(2) ))
  field_values = input_values
  where (field_values == -1) field_values = fill_val

  !write temporary array
  include 'H5write_template.inc'
  
  !de-allocate temporary array
  if(allocated(field_values)) deallocate(field_values)

end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_3DUI16(group_id, field_path, input_values, dcpl, field_sds )
  implicit none
  type(H5SDS_T) :: this_sds = &
    H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=2), dimension(:,:,:), intent(in) :: input_values
  integer(kind=4), dimension(:,:,:),allocatable :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(3) :: arr_shape
  integer , PARAMETER :: rank=3
  INTEGER(HID_T) :: memtypeId
  integer, PARAMETER :: fill_val=65535
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK8"

  memtypeId = H5T_STD_U16LE
  arr_shape = shape(input_values)
  allocate( field_values(arr_shape(1), arr_shape(2), arr_shape(3) ))
  field_values = input_values
  where (field_values == -1) field_values = fill_val

  !write temporary array
  include 'H5write_template.inc'
  
  !de-allocate temporary array
  if(allocated(field_values)) deallocate(field_values)

end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine H5write_4DUI16(group_id, field_path, input_values, dcpl, field_sds )
  implicit none
  type(H5SDS_T) :: this_sds = &
    H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)
  integer(HID_T), intent(in) :: group_id
  character(len=*), intent(in) :: field_path
  integer(kind=2), dimension(:,:,:,:), intent(in) :: input_values
  integer(kind=4), dimension(:,:,:,:),allocatable :: field_values
  type(H5SDS_T), intent(in), optional :: field_sds
  INTEGER(HID_T), INTENT(IN), optional :: dcpl
  integer,dimension(4) :: arr_shape
  integer , PARAMETER :: rank=4
  INTEGER(HID_T) :: memtypeId
  integer, PARAMETER :: fill_val=65535
  LOGICAL :: shall_write_attr = .FALSE.
  character(len=*), parameter :: routineName = "H5write_3DRK8"

  memtypeId = H5T_STD_U16LE
  arr_shape = shape(input_values)
  allocate( field_values(arr_shape(1), arr_shape(2), arr_shape(3), &
    arr_shape(4) ) )
  field_values = input_values
  where (field_values == -1) field_values = fill_val

  !write temporary array
  include 'H5write_template.inc'
  
  !de-allocate temporary array
  if(allocated(field_values)) deallocate(field_values)

end subroutine

subroutine H5write_utc_time(group_id, StringLength, StringArrayDim, utc_time)
  implicit none
  integer, INTENT(IN) :: StringLength
  integer, INTENT(IN) :: StringArrayDim
  character(LEN=*), dimension(:), intent(in):: utc_time
  type(H5SDS_CHAR_T)  :: sds_utc = & 
    H5SDS_CHAR_T("","","","",-1,-1,-1,0,(/0,0,0,0,0,0,0/) )
  integer(HID_T), intent(in) :: group_id
  integer :: errorStatus

  sds_utc%name = "UTC_Time"
  sds_utc%datatype_id = H5T_NATIVE_CHARACTER
  sds_utc%title = 'Time in UTC'
  sds_utc%rank = 1
  sds_utc%character_length = StringLength
  sds_utc%dims(1) = StringArrayDim
  sds_utc%content_type = "referenceInformation"
  errorStatus = H5Util_createDataset(group_id, sds_utc )
  errorStatus = H5Util_writeDataset(sds_utc, utc_time)
  errorStatus = H5Util_disposeDataset(sds_utc)

end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module H5write_module
