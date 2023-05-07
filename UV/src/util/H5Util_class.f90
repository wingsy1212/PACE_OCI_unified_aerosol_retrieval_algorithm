MODULE H5Util_class
!///////////////////////////////////////////////////////////////////////////////
!
! *** DESCRIPTION *** :
!
!    a set of general purpose high-level HDF5 utilities to handle dataset 
!    creation, writing and reading. It borrows many lines source code and some 
!    ideas from Kai Yang's UTIL_lh5_class.f90 in TO3_CORE. It takes care of all
!    my HDF5 I/O needs for working with OMPS L1B and L2 type of HDF5 files.
!
!
!...............................................................................
! REVISION HISTORY:
!    initial version:  Wed Mar 30 16:11:26 EDT 2011
!    complete history: svn log file:///omps/cm/app/JLi_Fortran_Toolkit
!
! AUTHOR: 
!    Jason Li (SSAI)
!///////////////////////////////////////////////////////////////////////////////

use HDF5
use DataTypeDef, only: SP, DP, I1B, I2B, I4B 
use MyConstants, only: SUCCESS_STATE, FAILURE_STATE, WARNING_STATE
use ErrorHandler_class, only: Display_Message

IMPLICIT NONE
PRIVATE


!...............................................................................
!
! constants:
!
!...............................................................................

INTEGER (I4B), PARAMETER :: MAXRANK = 7 ! max rank supported by Fortran90 

INTEGER (I4B), PARAMETER :: zero = 0

!...............................................................................
!
! public methods and data types:
!
!...............................................................................

! valid_range and fill_value in memory are double precision float (REAL(DP)).
! In the file,  they are stored, rightly so, as the same type as the dataset, 
! in other words it is sds%datatype_id.

! initialization:
!    type(H5SDS_T) :: sds = &
!         H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","","","",(/0,0/),0)

type, public :: H5SDS_T 
      character(len=128) :: name = ""
      integer(HID_T) :: dataset_id  = -1
      integer(HID_T) :: datatype_id = -1
      integer :: rank = 0
      integer(HSIZE_T), dimension(MAXRANK) :: dims = (/0,0,0,0,0,0,0/)
      character(len=128) :: title = ""
      character(len=128) :: units = ""
      character(len=128) :: content_type = ""
      character(len=1024) :: comment = ""
      real(DP), dimension(2) :: valid_range = (/ 0 , 0 /)
      real(DP) :: fill_value = -999.0
end type H5SDS_T


! initialization:
!    type(H5SDS_CHAR_T) :: sds = &
!         H5SDS_CHAR_T("","","","",-1,-1,-1,0,(/0,0,0,0,0,0,0/) )
type, public :: H5SDS_CHAR_T 
      character(len=128) :: name = ""
      character(len=128) :: content_type = ""
      character(len=128) :: title = ""
      character(len=1024) :: comment = ""
      integer(HID_T) :: dataset_id  = -1
      integer(HID_T) :: datatype_id = -1
      integer :: rank = 0 
      integer :: character_length = 0
      integer(HSIZE_T), dimension(MAXRANK) :: dims = (/0,0,0,0,0,0,0/)
end type H5SDS_CHAR_T


public :: H5Util_createDataset
public :: H5Util_selectDataset
public :: H5Util_readDataset
public :: H5Util_writeDataset
public :: H5Util_writeScalarAttribute
public :: H5Util_disposeDataset

!...............................................................................
! 
! function overloading:
!
!...............................................................................

interface H5Util_createDataset
   module procedure createDataset_char
   module procedure createDataset
end interface

interface H5Util_selectDataset
   module procedure selectDataset_char
   module procedure selectDataset
end interface

interface H5Util_disposeDataset
   module procedure disposeDataset_char
   module procedure disposeDataset
end interface



interface H5Util_writeDataset
   MODULE PROCEDURE writeDataset_0DIK4
   MODULE PROCEDURE writeDataset_1DIK4
   MODULE PROCEDURE writeDataset_2DIK4
   MODULE PROCEDURE writeDataset_3DIK4
   MODULE PROCEDURE writeDataset_4DIK4
   MODULE PROCEDURE writeDataset_5DIK4
   MODULE PROCEDURE writeDataset_6DIK4
   MODULE PROCEDURE writeDataset_7DIK4

   MODULE PROCEDURE writeDataset_0DRK4
   MODULE PROCEDURE writeDataset_1DRK4
   MODULE PROCEDURE writeDataset_2DRK4
   MODULE PROCEDURE writeDataset_3DRK4
   MODULE PROCEDURE writeDataset_4DRK4
   MODULE PROCEDURE writeDataset_5DRK4
   MODULE PROCEDURE writeDataset_6DRK4
   MODULE PROCEDURE writeDataset_7DRK4

   MODULE PROCEDURE writeDataset_0DRK8
   MODULE PROCEDURE writeDataset_1DRK8
   MODULE PROCEDURE writeDataset_2DRK8
   MODULE PROCEDURE writeDataset_3DRK8
   MODULE PROCEDURE writeDataset_4DRK8
   MODULE PROCEDURE writeDataset_5DRK8
   MODULE PROCEDURE writeDataset_6DRK8
   MODULE PROCEDURE writeDataset_7DRK8

  MODULE PROCEDURE writeDataset_0DCHAR
  MODULE PROCEDURE writeDataset_1DCHAR
  MODULE PROCEDURE writeDataset_2DCHAR
  MODULE PROCEDURE writeDataset_3DCHAR
  MODULE PROCEDURE writeDataset_4DCHAR
  MODULE PROCEDURE writeDataset_5DCHAR
  MODULE PROCEDURE writeDataset_6DCHAR
  MODULE PROCEDURE writeDataset_7DCHAR
end interface


interface H5Util_writeScalarAttribute
   module procedure writeScalarAttr_char
   module procedure writeScalarAttr_short
   module procedure writeScalarAttr_long
   module procedure writeScalarAttr_float
   module procedure writeScalarAttr_double
end interface


interface H5Util_readDataset
   MODULE PROCEDURE readDataset_0DIK4
   MODULE PROCEDURE readDataset_1DIK4
   MODULE PROCEDURE readDataset_2DIK4
   MODULE PROCEDURE readDataset_3DIK4
   MODULE PROCEDURE readDataset_4DIK4
   MODULE PROCEDURE readDataset_5DIK4
   MODULE PROCEDURE readDataset_6DIK4
   MODULE PROCEDURE readDataset_7DIK4

   MODULE PROCEDURE readDataset_0DRK4
   MODULE PROCEDURE readDataset_1DRK4
   MODULE PROCEDURE readDataset_2DRK4
   MODULE PROCEDURE readDataset_3DRK4
   MODULE PROCEDURE readDataset_4DRK4
   MODULE PROCEDURE readDataset_5DRK4
   MODULE PROCEDURE readDataset_6DRK4
   MODULE PROCEDURE readDataset_7DRK4

   MODULE PROCEDURE readDataset_0DRK8
   MODULE PROCEDURE readDataset_1DRK8
   MODULE PROCEDURE readDataset_2DRK8
   MODULE PROCEDURE readDataset_3DRK8
   MODULE PROCEDURE readDataset_4DRK8
   MODULE PROCEDURE readDataset_5DRK8
   MODULE PROCEDURE readDataset_6DRK8
   MODULE PROCEDURE readDataset_7DRK8

   MODULE PROCEDURE readDataset_0DCHAR
   MODULE PROCEDURE readDataset_1DCHAR
   MODULE PROCEDURE readDataset_2DCHAR
   MODULE PROCEDURE readDataset_3DCHAR
   MODULE PROCEDURE readDataset_4DCHAR
   MODULE PROCEDURE readDataset_5DCHAR
   MODULE PROCEDURE readDataset_6DCHAR
   MODULE PROCEDURE readDataset_7DCHAR
end interface


CONTAINS

!///////////////////////////////////////////////////////////////////////////////
!
!///////////////////////// P U B L I C   M E T H O D S /////////////////////////
!
!///////////////////////////////////////////////////////////////////////////////

function createDataset(parent_id, sds, &
                       L_createAttributes, dcpl) result(errorStatus)

use H5LT
implicit none

character(len=*), parameter :: routineName = "H5Util_createDataset"

integer(HID_T), intent(in) :: parent_id
type(H5SDS_T), intent(inout) :: sds
integer :: errorStatus
logical, intent(in), optional :: L_createAttributes
integer(HID_T), intent(in), optional :: dcpl

integer(HID_T) :: dataspace_id, datatype_id, property_id
integer :: hdferr, exists_status

character(len=256) :: msg


errorStatus = SUCCESS_STATE

exists_status = h5ltfind_dataset_f(parent_id, sds%name)
if(exists_status == 1) then
   errorStatus = WARNING_STATE
   write(msg,*) trim(sds%name), " already exists, cannot create it again" 
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

!
! verify I have enough information to create a dataset:
!
if(sds%rank < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide rank for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

if(sds%datatype_id < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide datatype_id for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif


!
! define dataspace, simple or scalar?
!
if(sds%rank == 0) then
   call h5screate_f(H5S_SCALAR_F, dataspace_id, hdferr)
else
   call H5Screate_simple_f(sds%rank, sds%dims, dataspace_id, hdferr)
endif
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "H5Screate_simple_f() failed for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
      goto 99999
endif

!
! create a copy of data type:
!
call h5tcopy_f(sds%datatype_id, datatype_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "h5tcopy_f() failed for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif


!
! create property
!
if(sds%rank > 0) then
   if (present(dcpl)) then
       call h5pcopy_f(dcpl, property_id, hdferr) 
   else
       call H5Pcreate_f(H5P_DATASET_CREATE_F, property_id, hdferr)
   endif

   if(hdferr < 0) then
      errorStatus = FAILURE_STATE
      write(msg,*) "H5Pcreate_f() failed for dataset: ", trim(sds%name)
      call Display_Message(routineName, trim(msg), errorStatus)
      goto 99999
   endif

   call H5Pset_fill_value_f(property_id,H5T_NATIVE_DOUBLE,sds%fill_value,hdferr)
   if(hdferr < 0) then
      errorStatus = FAILURE_STATE
      write(msg,*) "H5Pcreate_f() failed for dataset: ", trim(sds%name)
      call Display_Message(routineName, trim(msg), errorStatus)
      goto 99999
   endif
endif

!
! create dataset:
!
if(sds%rank > 0) then
   call H5Dcreate_f(parent_id, sds%name, &
                    datatype_id, dataspace_id, sds%dataset_id, &
                    hdferr, property_id)
else
   call H5Dcreate_f(parent_id, sds%name, &
                    datatype_id, dataspace_id, sds%dataset_id, &
                    hdferr)
endif

!
! write dataset attributes (optional):
!
if(present(L_createAttributes)) then
  if(L_createAttributes) call writeDatasetAttributes(sds, hdferr)
endif

!
! closing up shop:
!
call H5Tclose_f(datatype_id, hdferr)
if(sds%rank > 0) call H5Pclose_f(property_id, hdferr)
call H5Sclose_f(dataspace_id, hdferr)


99999 return

end function createDataset

!...............................................................................

function createDataset_char(parent_id, sds) result (errorStatus)

use H5LT
implicit none

character(len=*), parameter :: routineName = "createDataset_char"

integer(HID_T), intent(in) :: parent_id
type(H5SDS_CHAR_T), intent(inout) :: sds
integer :: errorStatus

integer(HID_T) :: dataspace_id, datatype_id, attr_id
character(len=256) :: msg
integer :: hdferr, exists_status
logical :: sameType

integer(SIZE_T) :: attr_len, characterLength
integer(HSIZE_T), dimension(1) :: dims 
integer :: rank


 
errorStatus = SUCCESS_STATE

exists_status = h5ltfind_dataset_f(parent_id, sds%name)
if(exists_status == 1) then
   errorStatus = WARNING_STATE
   write(msg,*) trim(sds%name), " already exists, cannot create it again"
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

!
! verify I have enough information to create a dataset:
!
if(sds%rank <= 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide rank for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

! if(sds%datatype_id /= H5T_NATIVE_CHARACTER) then
call h5tequal_f(sds%datatype_id, H5T_NATIVE_CHARACTER, sameType, hdferr)
if(.not. sameType) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide datatype_id for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

if(sds%character_length <= 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide character_length for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif


characterLength = sds%character_length
call h5tcopy_f(H5T_NATIVE_CHARACTER, datatype_id, hdferr)
call h5tset_size_f(datatype_id, characterLength, hdferr)
call H5screate_simple_f(sds%rank, sds%dims, dataspace_id, hdferr)
call h5dcreate_f(parent_id, sds%name, &
                 datatype_id, dataspace_id, sds%dataset_id, hdferr)
call h5sclose_f(dataspace_id, hdferr)
call h5tclose_f(datatype_id, hdferr)

! ACDD coverage content type (DISC):
! http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery
!     image, thematicClassification, physicalMeasurement, 
!     auxiliaryInformation, qualityInformation, referenceInformation, 
!     modelResult, or coordinate
!

attr_len = len_trim(sds%content_type)
if(attr_len > 1) then
!  dims(1) = 1
!  call H5Screate_f(H5S_SCALAR_F, dataspace_id, hdferr)
!  call H5Tcopy_f(H5T_NATIVE_CHARACTER, datatype_id, hdferr)
!  call H5Tset_strpad_f(datatype_id, H5T_STR_NULLTERM_F, hdferr)
!  call H5Tset_size_f(datatype_id, attr_len, hdferr)
!  call H5Acreate_f(sds%dataset_id, "coverage_content_type", &
!                   datatype_id, dataspace_id, attr_id, hdferr)
!  call H5Awrite_f(attr_id, datatype_id, sds%content_type, dims, hdferr)
!  call H5Aclose_f(attr_id, hdferr)
!  call H5Tclose_f(datatype_id, hdferr)
!  call H5Sclose_f(dataspace_id, hdferr)

   call H5Util_writeScalarAttribute(sds%dataset_id, 'coverage_content_type', &
                                    sds%content_type, hdferr)
endif

attr_len = len_trim(sds%comment)
if(attr_len > 1) then
   call H5Util_writeScalarAttribute(sds%dataset_id, 'comment', &
                                    sds%comment, hdferr)
endif


attr_len = len_trim(sds%title)
if(attr_len > 1) then
   call H5Util_writeScalarAttribute(sds%dataset_id,'long_name',sds%title,hdferr)
endif
99999 return

end function createDataset_char


!*******************************************************************************

function selectDataset(parent_id, sds) result(errorStatus)

implicit none

character(len=*), parameter :: routineName = "selectDataset"

integer(HID_T), intent(in) :: parent_id
type(H5SDS_T), intent(inout) :: sds
integer :: errorStatus

integer(HID_T) :: dataspace_id
integer(HSIZE_T), dimension(MAXRANK) :: maxdims
integer :: rank, hdferr


errorStatus = SUCCESS_STATE

if(sds%dataset_id > 0) then
   errorStatus = WARNING_STATE
   call Display_Message(routineName, &
               "dataset " // trim(sds%name) // " has already been selected", &
                errorStatus)
   goto 99999
endif

call h5dopen_f(parent_id, sds%name, sds%dataset_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "HDF5 error acquiring " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5dget_type_f(sds%dataset_id, sds%datatype_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5dget_type_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5dget_space_f(sds%dataset_id, dataspace_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5dget_space_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5sget_simple_extent_ndims_f(dataspace_id, sds%rank, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5sget_simple_extent_ndims_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif


call h5sget_simple_extent_dims_f(dataspace_id, sds%dims, maxdims, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5sget_simple_extent_dims_f error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

! read in attributes if they exist:

call readDatasetAttributes(sds, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "readDatasetAttributes error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif


99999 return

end function selectDataset

!...............................................................................

function selectDataset_char(parent_id, sds) result(errorStatus)

implicit none

character(len=*), parameter :: routineName = "selectDataset_char"

integer(HID_T), intent(in) :: parent_id
type(H5SDS_CHAR_T), intent(inout) :: sds
integer :: errorStatus

integer(HID_T) :: dataspace_id
integer(HSIZE_T), dimension(MAXRANK) :: maxdims
integer :: rank, hdferr


errorStatus = SUCCESS_STATE

if(sds%dataset_id > 0) then
   errorStatus = WARNING_STATE
   call Display_Message(routineName, &
               "dataset " // trim(sds%name) // " has already been selected", &
                errorStatus)
   goto 99999
endif

call h5dopen_f(parent_id, sds%name, sds%dataset_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "HDF5 error acquiring " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5dget_type_f(sds%dataset_id, sds%datatype_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5dget_type_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5dget_space_f(sds%dataset_id, dataspace_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5dget_space_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5sget_simple_extent_ndims_f(dataspace_id, sds%rank, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5sget_simple_extent_ndims_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5sget_simple_extent_dims_f(dataspace_id, sds%dims, maxdims, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5sget_simple_extent_dims_f error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

99999 return

end function selectDataset_char


!*******************************************************************************

function disposeDataset(sds) result(errorStatus)

implicit none

character(len=*), parameter :: routineName = "H5Util_disposeDataset"

type(H5SDS_T), intent(inout) :: sds
integer :: errorStatus, hdferr


errorStatus = SUCCESS_STATE

if(sds%dataset_id <= 0) then
   errorStatus = WARNING_STATE
   call Display_Message(routineName, &
               "invalid dataset_id for " // trim(sds%name), errorStatus)
   return
endif

call h5dclose_f(sds%dataset_id, hdferr)


sds%dataset_id = -1
sds%datatype_id = -1
sds%rank = 0
sds%dims = 0
sds%title = ""
sds%units = ""
sds%content_type = ""
sds%comment = ""
sds%valid_range = 0.0_DP
sds%fill_value = 0.0_DP

end function disposeDataset

!...............................................................................

function disposeDataset_char(sds) result(errorStatus)

implicit none

character(len=*), parameter :: routineName = "disposeDataset_char"

type(H5SDS_CHAR_T), intent(inout) :: sds
integer :: errorStatus, hdferr


errorStatus = SUCCESS_STATE

if(sds%dataset_id <= 0) then
   errorStatus = WARNING_STATE
   call Display_Message(routineName, &
               "invalid dataset_id for " // trim(sds%name), errorStatus)
   return
endif

call h5dclose_f(sds%dataset_id, hdferr)

sds%content_type = ""
sds%comment = ""
sds%dataset_id = -1
sds%datatype_id = -1
sds%rank = 0
sds%dims = 0
sds%character_length = 0

end function disposeDataset_char


!*******************************************************************************
!
! errorStatus = writeDataset_*(ds, data_in, start, count)
!
!*******************************************************************************

FUNCTION writeDataset_0DIK4( ds, data_in, start, count ) RESULT (status)

! scalar dataset, rank = 0 

IMPLICIT NONE
TYPE(H5SDS_T), INTENT( IN ) :: ds
INTEGER(I4B), INTENT(IN) :: data_in
INTEGER(HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER(I4B) :: status
CHARACTER(LEN=256) :: msg
INTEGER(HID_T) :: memspace , dataspace
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm


status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ",ds%name
   call Display_Message("writeDataset_0DIK4", msg, status)
   RETURN
ENDIF

dimsm = 1 
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call h5screate_f(H5S_SCALAR_F, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF


call h5dwrite_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_in, dimsm, error, &
                 memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_0DIK4", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT( start ) ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_0DIK4

!...............................................................................

FUNCTION writeDataset_1DIK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'writeDataset_1DIK4'

TYPE(H5SDS_T), INTENT( IN ) :: ds
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: data_in
INTEGER(HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER(I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HID_T) :: memtypeId, memspace, dataspace 
INTEGER :: rankm 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm

memtypeId = H5T_NATIVE_INTEGER

include 'writeDataset_template.inc'

END FUNCTION writeDataset_1DIK4

!...............................................................................

FUNCTION writeDataset_2DIK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'writeDataset_2DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HID_T) :: memtypeId, memspace, dataspace 
INTEGER :: rankm
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

memtypeId = H5T_NATIVE_INTEGER

include 'writeDataset_template.inc'

END FUNCTION writeDataset_2DIK4

!...............................................................................

FUNCTION writeDataset_3DIK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'writeDataset_3DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HID_T) :: memtypeId, memspace 
INTEGER(HID_T) :: dataspace 
INTEGER :: rankm
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

memtypeId = H5T_NATIVE_INTEGER

include 'writeDataset_template.inc'

END FUNCTION writeDataset_3DIK4

!...............................................................................

FUNCTION writeDataset_4DIK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'writeDataset_4DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

memtypeId = H5T_NATIVE_INTEGER

include 'writeDataset_template.inc'

END FUNCTION writeDataset_4DIK4

!...............................................................................

FUNCTION writeDataset_5DIK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'writeDataset_5DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER(I4B), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
INTEGER(HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER:: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memtypeId, dataspace 
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

memtypeId = H5T_NATIVE_INTEGER

include 'writeDataset_template.inc'

END FUNCTION writeDataset_5DIK4

!...............................................................................

FUNCTION writeDataset_6DIK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'writeDataset_6DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memtypeId, memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

memtypeId = H5T_NATIVE_INTEGER

include 'writeDataset_template.inc'

END FUNCTION writeDataset_6DIK4

!...............................................................................

FUNCTION writeDataset_7DIK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'writeDataset_7DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm

memtypeId = H5T_NATIVE_INTEGER

include 'writeDataset_template.inc'

END FUNCTION writeDataset_7DIK4

!...............................................................................

FUNCTION writeDataset_0DRK4( ds, data_in, start, count ) RESULT (status)

implicit none

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(SP), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HID_T) :: memspace, dataspace, native_datatype_id
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
logical :: sameType
status = SUCCESS_STATE


! IF( ds%datatype_id /= H5T_NATIVE_REAL) then
call h5tequal_f(ds%datatype_id, H5T_NATIVE_REAL, sameType, error) 
IF( .not. sameType) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("writeDataset_0DRK4", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_0DRK4", msg, status)
   RETURN
ENDIF

dimsm = 1
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call H5Screate_f(H5S_SCALAR_F, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dimsm, error, &
                 memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_0DRK4", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error )

END FUNCTION writeDataset_0DRK4

!...............................................................................

FUNCTION writeDataset_1DRK4( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_1DRK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(SP), DIMENSION(:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_REAL

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_1DRK4

!...............................................................................

FUNCTION writeDataset_2DRK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_2DRK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_REAL

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_2DRK4

!...............................................................................

FUNCTION writeDataset_3DRK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_3DRK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_REAL

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_3DRK4

!...............................................................................

FUNCTION writeDataset_4DRK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_4DRK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_REAL

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_4DRK4

!...............................................................................

FUNCTION writeDataset_5DRK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_5DRK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_REAL

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_5DRK4

!...............................................................................

FUNCTION writeDataset_6DRK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_6DRK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_REAL

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_6DRK4

!...............................................................................

FUNCTION writeDataset_7DRK4( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_7DRK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_REAL

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_7DRK4

!...............................................................................

FUNCTION writeDataset_0DRK8( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(DP), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HID_T) :: memspace,  dataspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
logical :: sameType

status = SUCCESS_STATE

! IF( ds%datatype_id /= H5T_NATIVE_DOUBLE) then
 call h5tequal_f(ds%datatype_id, H5T_NATIVE_DOUBLE, sameType, error)
IF(.not. sameType) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("writeDataset_0DRK8", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_0DRK8", msg, status)
   RETURN
ENDIF

dimsm = 1
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call H5Screate_f(H5S_SCALAR_F, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dimsm, error, &
                 memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_0DRK8", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error )

END FUNCTION writeDataset_0DRK8

!...............................................................................


FUNCTION writeDataset_1DRK8( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_1DRK8'

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(DP), DIMENSION(:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memtypeId, memspace,  dataspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
logical :: sameType


memtypeId = H5T_NATIVE_DOUBLE

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_1DRK8

!...............................................................................

FUNCTION writeDataset_2DRK8( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_2DRK8'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_DOUBLE

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_2DRK8

!...............................................................................

FUNCTION writeDataset_3DRK8( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_3DRK8'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_DOUBLE

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_3DRK8

!...............................................................................

FUNCTION writeDataset_4DRK8( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_4DRK8'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_DOUBLE

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_4DRK8

!...............................................................................

FUNCTION writeDataset_5DRK8( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_5DRK8'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
logical :: sameType

memtypeId = H5T_NATIVE_DOUBLE

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_5DRK8

!...............................................................................

FUNCTION writeDataset_6DRK8( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_6DRK8'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
logical :: sameType


memtypeId = H5T_NATIVE_DOUBLE

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_6DRK8

!...............................................................................

FUNCTION writeDataset_7DRK8( ds, data_in, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_7DRK8'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace, memtypeId 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
logical :: sameType


memtypeId = H5T_NATIVE_DOUBLE

include 'writeDataset_REAL_template.inc'

END FUNCTION writeDataset_7DRK8

!...............................................................................

FUNCTION writeDataset_0DCHAR( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_0DCHAR'

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
CHARACTER(LEN=*), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER(HID_T) :: dataspace, datatype_id, memtype_id 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
logical :: sameType
INTEGER(SIZE_T) :: dataInLen


! include 'writeDataset_CHAR_template.inc'

status = SUCCESS_STATE

! IF( ds%datatype_id /= H5T_NATIVE_CHARACTER) then
!call h5tequal_f(ds%datatype_id, H5T_NATIVE_CHARACTER, sameType, error)
!if(.not. sameType) then
!   status = FAILURE_STATE
!   write( msg,* ) "data type should be H5T_NATIVE_CHARACTER: ", ds%name
!   call Display_Message("writeDataset_0DCHAR", msg, status)
!   RETURN
!ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_0DCHAR", msg, status)
   RETURN
ENDIF

dimsm = 1
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call H5Screate_f(H5S_SCALAR_F, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF


call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
dataInLen = len(data_in) 
call h5tset_size_f(memtype_id, dataInLen, error)
call h5dwrite_f( ds%dataset_id, memtype_id, data_in, dimsm, error, &
                 memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_0DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_0DCHAR

!...............................................................................

FUNCTION writeDataset_1DCHAR( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_1DCHAR'

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER(HID_T) :: dataspace, datatype_id, memtype_id 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
logical :: sameType
INTEGER(SIZE_T) :: dataInLen

include 'writeDataset_CHAR_template.inc'

END FUNCTION writeDataset_1DCHAR

!...............................................................................

FUNCTION writeDataset_2DCHAR( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_2DCHAR'

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
CHARACTER(LEN=*), DIMENSION(:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER(HID_T) :: dataspace, datatype_id, memtype_id 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
logical :: sameType
INTEGER(SIZE_T) :: dataInLen


include 'writeDataset_CHAR_template.inc'

END FUNCTION writeDataset_2DCHAR

!...............................................................................

FUNCTION writeDataset_3DCHAR( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_3DCHAR'

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
CHARACTER(LEN=*), DIMENSION(:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER(HID_T) :: dataspace, datatype_id, memtype_id 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
logical :: sameType
INTEGER(SIZE_T) :: dataInLen

include 'writeDataset_CHAR_template.inc'


END FUNCTION writeDataset_3DCHAR

!...............................................................................

FUNCTION writeDataset_4DCHAR( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_4DCHAR'

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
CHARACTER(LEN=*), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER(HID_T) :: dataspace, datatype_id, memtype_id 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
logical :: sameType
INTEGER(SIZE_T) :: dataInLen

include 'writeDataset_CHAR_template.inc'

END FUNCTION writeDataset_4DCHAR

!...............................................................................

FUNCTION writeDataset_5DCHAR( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_5DCHAR'

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
CHARACTER(LEN=*), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER(HID_T) :: dataspace, datatype_id, memtype_id 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
logical :: sameType
INTEGER(SIZE_T) :: dataInLen


include 'writeDataset_CHAR_template.inc'

END FUNCTION writeDataset_5DCHAR

!...............................................................................

FUNCTION writeDataset_6DCHAR( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_6DCHAR'

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
CHARACTER(LEN=*), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER(HID_T) :: dataspace, datatype_id, memtype_id 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
logical :: sameType
INTEGER(SIZE_T) :: dataInLen


include 'writeDataset_CHAR_template.inc'

END FUNCTION writeDataset_6DCHAR

!...............................................................................

FUNCTION writeDataset_7DCHAR( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'writeDataset_7DCHAR'

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
CHARACTER(LEN=*), DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER :: error
INTEGER(HID_T) :: dataspace, datatype_id, memtype_id 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
logical :: sameType
INTEGER(SIZE_T) :: dataInLen

include 'writeDataset_CHAR_template.inc'

END FUNCTION writeDataset_7DCHAR

!*******************************************************************************
!
! errorStatus = readDataset_*(ds, data_out, start, count)
!
!*******************************************************************************

FUNCTION readDataset_0DIK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds 
INTEGER (I4B), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HID_T) :: memspace,  dataspace
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
integer :: typeClass


status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_0DIK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
IF(typeClass /=  H5T_INTEGER_F) then
   status = FAILURE_STATE
   write( msg,* ) "class type should be H5T_INTEGER_F: ", ds%name
   call Display_Message("readDataset_0DIK4", msg, status)
   RETURN
ENDIF

dimsm = 1
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
call H5Screate_f(H5S_SCALAR_F, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f(ds%dataset_id, H5T_NATIVE_INTEGER, data_out, dimsm, error, &
               memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_0DIK4", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error )

END FUNCTION readDataset_0DIK4

!...............................................................................


FUNCTION readDataset_1DIK4( ds, data_out, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'readDataset_1DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds 
INTEGER (I4B), DIMENSION(:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
integer :: typeClass

include 'readDataset_template.inc'

END FUNCTION readDataset_1DIK4

!...............................................................................

FUNCTION readDataset_2DIK4( ds, data_out, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'readDataset_2DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
integer :: typeClass

include 'readDataset_template.inc'

END FUNCTION readDataset_2DIK4

!...............................................................................

FUNCTION readDataset_3DIK4( ds, data_out, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'readDataset_3DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
integer :: typeClass

include 'readDataset_template.inc'

END FUNCTION readDataset_3DIK4

!...............................................................................

FUNCTION readDataset_4DIK4( ds, data_out, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'readDataset_4DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
integer :: typeClass

include 'readDataset_template.inc'

END FUNCTION readDataset_4DIK4

!...............................................................................

FUNCTION readDataset_5DIK4( ds, data_out, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'readDataset_5DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
integer :: typeClass

include 'readDataset_template.inc'

END FUNCTION readDataset_5DIK4

!...............................................................................

FUNCTION readDataset_6DIK4( ds, data_out, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'readDataset_6DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER(HID_T) :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
integer :: typeClass

include 'readDataset_template.inc'

END FUNCTION readDataset_6DIK4

!...............................................................................

FUNCTION readDataset_7DIK4( ds, data_out, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'readDataset_7DIK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER(HID_T) :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
integer :: typeClass


include 'readDataset_template.inc'

END FUNCTION readDataset_7DIK4

!...............................................................................

FUNCTION readDataset_0DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(SP), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
INTEGER(HID_T) :: memspace,  dataspace
integer :: typeClass
INTEGER(SIZE_T) :: typeSize



status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_0DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error) 
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_0DRK4", msg, status)
   RETURN
ENDIF

dimsm(1) = 1
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call H5Screate_f(H5S_SCALAR_F, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error , &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_0DRK4", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_0DRK4

!...............................................................................


FUNCTION readDataset_1DRK4( ds, data_out, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'readDataset_1DRK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(SP), DIMENSION(:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize

include 'readDataset_REAL_template.inc'

END FUNCTION readDataset_1DRK4

!...............................................................................

FUNCTION readDataset_2DRK4( ds, data_out, start, count ) RESULT (status)

CHARACTER(LEN=*), PARAMETER :: ROUTINE_NAME = 'readDataset_2DRK4'

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize


include 'readDataset_REAL_template.inc'

END FUNCTION readDataset_2DRK4

!...............................................................................

FUNCTION readDataset_3DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
integer :: typeClass
integer(SIZE_T) :: typeSize


status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_3DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_3DRK4", msg, status)
   RETURN
ENDIF

rankm = 3
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_3DRK4", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_3DRK4

!...............................................................................

FUNCTION readDataset_4DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
integer :: typeClass
integer(SIZE_T) :: typeSize



status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_4DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_4DRK4", msg, status)
   RETURN
ENDIF

rankm = 4
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_4DRK4", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_4DRK4

!...............................................................................

FUNCTION readDataset_5DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize


status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_5DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_5DRK4", msg, status)
   RETURN
ENDIF

rankm = 5
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_5DRK4", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_5DRK4

!...............................................................................

FUNCTION readDataset_6DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize


status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_6DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_6DRK4", msg, status)
   RETURN
ENDIF

rankm = 6
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_6DRK4", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_6DRK4

!...............................................................................

FUNCTION readDataset_7DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize


status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_7DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_7DRK4", msg, status)
   RETURN
ENDIF

rankm = 7
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_7DRK4", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_7DRK4

!...............................................................................

FUNCTION readDataset_0DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(DP), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count

INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
INTEGER(HID_T) :: memspace,  dataspace
integer :: typeClass
INTEGER(SIZE_T) :: typeSize



status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_0DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_0DRK8", msg, status)
   RETURN
ENDIF

dimsm = 1
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call H5Screate_f(H5S_SCALAR_F, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_0DRK8", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error )

END FUNCTION readDataset_0DRK8

!...............................................................................


FUNCTION readDataset_1DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(DP), DIMENSION(:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize



status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_1DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_1DRK8", msg, status)
   RETURN
ENDIF


rankm = 1
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_1DRK8", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_1DRK8

!...............................................................................

FUNCTION readDataset_2DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize



status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_2DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_2DRK8", msg, status)
   RETURN
ENDIF


rankm = 2
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_2DRK8", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_2DRK8

!...............................................................................

FUNCTION readDataset_3DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize


status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_3DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_3DRK8", msg, status)
   RETURN
ENDIF

rankm = 3
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_3DRK8", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_3DRK8

!...............................................................................

FUNCTION readDataset_4DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize



status = SUCCESS_STATE


IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_4DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_4DRK8", msg, status)
   RETURN
ENDIF

rankm = 4
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_4DRK8", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_4DRK8

!...............................................................................

FUNCTION readDataset_5DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
integer :: typeClass
INTEGER(SIZE_T) :: typeSize


status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_5DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_5DRK8", msg, status)
   RETURN
ENDIF

rankm = 5
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_5DRK8", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_5DRK8

!...............................................................................

FUNCTION readDataset_6DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
INTEGER :: typeClass
INTEGER(SIZE_T) :: typeSize



status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_6DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_6DRK8", msg, status)
   RETURN
ENDIF


rankm = 6
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_6DRK8", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_6DRK8

!...............................................................................

FUNCTION readDataset_7DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER :: rankm 
INTEGER(HID_T) :: memspace 
INTEGER(HID_T) :: dataspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
integer :: typeClass
integer(SIZE_T) :: typeSize



status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_7DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_7DRK8", msg, status)
   RETURN
ENDIF


rankm = 7
dimsm = SHAPE( data_out )
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dimsm, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_7DRK8", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_7DRK8

!...............................................................................

FUNCTION readDataset_0DCHAR( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds
CHARACTER(LEN=*), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error, status
INTEGER(HID_T) :: memtype_id, dataspace, memspace
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
INTEGER(SIZE_T) :: dataOutLen


status = SUCCESS_STATE


IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_0DCHAR", msg, status)
   RETURN
ENDIF

dimsm = 1
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call H5Screate_f(H5S_SCALAR_F, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
dataOutLen = len(data_out) 
call h5tset_size_f(memtype_id, dataOutLen, error)
CALL h5dread_f(ds%dataset_id, memtype_id, data_out, dimsm, error, &
               memspace, dataspace)
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_0DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)

END FUNCTION readDataset_0DCHAR

!...............................................................................

FUNCTION readDataset_1DCHAR( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds
CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status, rankm
INTEGER(HID_T) :: memtype_id, dataspace, memspace
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
INTEGER(SIZE_T) :: dataOutLen


status = SUCCESS_STATE


IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_1DCHAR", msg, status)
   RETURN
ENDIF

rankm = 1
dimsm = shape(data_out)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call h5screate_simple_f(rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
dataOutLen = len(data_out) 
call h5tset_size_f(memtype_id, dataOutLen, error)
CALL h5dread_f(ds%dataset_id, memtype_id, data_out, dimsm, &
               error, memspace, dataspace)
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_1DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)

END FUNCTION readDataset_1DCHAR

!...............................................................................

FUNCTION readDataset_2DCHAR( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds
CHARACTER(LEN=*), DIMENSION(:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status, rankm
INTEGER(HID_T) :: memtype_id, dataspace, memspace
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
INTEGER(SIZE_T) :: dataOutLen


status = SUCCESS_STATE


IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_2DCHAR", msg, status)
   RETURN
ENDIF

rankm = 2
dimsm = shape(data_out)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call h5screate_simple_f(rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
dataOutLen = len(data_out) 
call h5tset_size_f(memtype_id, dataOutLen, error)
call h5dread_f(ds%dataset_id, memtype_id, data_out, dimsm, &
               error, memspace, dataspace)
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_2DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)

END FUNCTION readDataset_2DCHAR


!...............................................................................

FUNCTION readDataset_3DCHAR( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds
CHARACTER(LEN=*), DIMENSION(:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status, rankm
INTEGER(HID_T) :: memtype_id, dataspace, memspace
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
INTEGER(SIZE_T) :: dataOutLen


status = SUCCESS_STATE


IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_3DCHAR", msg, status)
   RETURN
ENDIF

rankm = 3
dimsm = shape(data_out)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call h5screate_simple_f(rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
dataOutLen = len(data_out) 
call h5tset_size_f(memtype_id, dataOutLen, error)
CALL h5dread_f(ds%dataset_id, memtype_id, data_out, dimsm, &
               error, memspace, dataspace)
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_3DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)

END FUNCTION readDataset_3DCHAR


!...............................................................................

FUNCTION readDataset_4DCHAR( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds
CHARACTER(LEN=*), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status, rankm
INTEGER(HID_T) :: memtype_id, dataspace, memspace
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
INTEGER(SIZE_T) :: dataOutLen


status = SUCCESS_STATE


IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_4DCHAR", msg, status)
   RETURN
ENDIF

rankm = 4
dimsm = shape(data_out)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call h5screate_simple_f(rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
dataOutLen = len(data_out) 
call h5tset_size_f(memtype_id, dataOutLen, error)
CALL h5dread_f(ds%dataset_id, memtype_id, data_out, dimsm, &
               error, memspace, dataspace)
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_4DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)

END FUNCTION readDataset_4DCHAR


!...............................................................................

FUNCTION readDataset_5DCHAR( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds
CHARACTER(LEN=*), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status, rankm
INTEGER(HID_T) :: memtype_id, dataspace, memspace
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
INTEGER(SIZE_T) :: dataOutLen


status = SUCCESS_STATE


IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_5DCHAR", msg, status)
   RETURN
ENDIF

rankm = 5
dimsm = shape(data_out)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call h5screate_simple_f(rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

dataOutLen = len(data_out) 
call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
call h5tset_size_f(memtype_id, dataOutLen, error)
CALL h5dread_f(ds%dataset_id, memtype_id, data_out, dimsm, &
               error, memspace, dataspace)
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_5DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)

END FUNCTION readDataset_5DCHAR

!...............................................................................

FUNCTION readDataset_6DCHAR( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds
CHARACTER(LEN=*), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status, rankm
INTEGER(HID_T) :: memtype_id, dataspace, memspace
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
INTEGER(SIZE_T) :: dataOutLen


status = SUCCESS_STATE


IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_6DCHAR", msg, status)
   RETURN
ENDIF

rankm = 6
dimsm = shape(data_out)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call h5screate_simple_f(rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
dataOutLen = len(data_out)
call h5tset_size_f(memtype_id, dataOutLen, error)
CALL h5dread_f(ds%dataset_id, memtype_id, data_out, dimsm, &
               error, memspace, dataspace)
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_6DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)

END FUNCTION readDataset_6DCHAR


!...............................................................................

FUNCTION readDataset_7DCHAR( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds
CHARACTER(LEN=*), DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status, rankm
INTEGER(HID_T) :: memtype_id, dataspace, memspace
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
INTEGER(SIZE_T) :: dataOutLen

status = SUCCESS_STATE


IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_7DCHAR", msg, status)
   RETURN
ENDIF

rankm = 7
dimsm = shape(data_out)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   call h5screate_simple_f(rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
dataOutLen = len(data_out) 
call h5tset_size_f(memtype_id, dataOutLen, error)
CALL h5dread_f(ds%dataset_id, memtype_id, data_out, dimsm, &
               error, memspace, dataspace)
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_7DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT(start) ) call h5sclose_f(memspace, error)

END FUNCTION readDataset_7DCHAR


!///////////////////////////////////////////////////////////////////////////////
!
!//////////////////////// P R I V A T E   M E T H O D S ////////////////////////
!
!///////////////////////////////////////////////////////////////////////////////


subroutine readDatasetAttributes(sds, hdferr)

implicit none
type(H5SDS_T), intent(inout) :: sds
integer, intent(out) :: hdferr

integer(HID_T) :: memtype_id, attr_id
integer(SIZE_T) :: attr_len
integer(HSIZE_T), dimension(1) :: dims 

character(len=128) :: attr_name
logical :: attr_exists

!
! title:
!
attr_name = "long_name"
call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
if(attr_exists) then
   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
   dims(1) = 1
   attr_len = len(sds%title)
   call H5Tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
   call h5tset_size_f(memtype_id, attr_len, hdferr)
   call h5aread_f(attr_id, memtype_id,  sds%title, dims, hdferr) 
   call h5tclose_f(memtype_id, hdferr)
   call h5aclose_f(attr_id, hdferr)
endif

!
! units:
!
attr_name = "units"
call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
if(attr_exists) then
   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
   dims(1) = 1
   attr_len = len(sds%units)
   call H5Tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
   call h5tset_size_f(memtype_id, attr_len, hdferr)
   call h5aread_f(attr_id, memtype_id,  sds%units, dims, hdferr) 
   call h5tclose_f(memtype_id, hdferr)
   call h5aclose_f(attr_id, hdferr)
endif

!
! valid range:
!
attr_name = "valid_range"
call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
if(attr_exists) then
   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
   dims(1) = 2
   call h5aread_f(attr_id, H5T_NATIVE_DOUBLE,  sds%valid_range, dims, hdferr) 
   call h5aclose_f(attr_id, hdferr)
endif


!
! fill value attribute:
!
attr_name = "_FillValue"
call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
if(attr_exists) then
   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
   dims(1) = 1
   call h5aread_f(attr_id, H5T_NATIVE_DOUBLE,  sds%fill_value, dims, hdferr) 
   call h5aclose_f(attr_id, hdferr)
endif

end subroutine readDatasetAttributes


!*******************************************************************************

subroutine writeDatasetAttributes(sds, hdferr)

! Description: 
!    create dataset attributes: title, units, valid range and fill value. Very
!    similar to the set of attributes for HDFEOS datasets. 

implicit none
type(H5SDS_T), intent(in) :: sds
integer, intent(out) :: hdferr

integer(HID_T) :: datatype_id, dataspace_id, attr_id
integer(SIZE_T) :: attr_len
integer(HSIZE_T), dimension(1) :: dims 
integer :: rank

rank = 1

!
! title, aka long_name:
!
!dims(1) = 1
!attr_len = len_trim(sds%title)
!call H5Screate_f(H5S_SCALAR_F, dataspace_id, hdferr)
!call H5Tcopy_f(H5T_NATIVE_CHARACTER, datatype_id, hdferr)
!call H5Tset_strpad_f(datatype_id, H5T_STR_NULLTERM_F, hdferr)
!call H5Tset_size_f(datatype_id, attr_len, hdferr)
!call H5Acreate_f(sds%dataset_id, "long_name",&
!                 datatype_id, dataspace_id, attr_id, hdferr)
!call H5Awrite_f(attr_id, datatype_id, sds%title, dims, hdferr)
!call H5Aclose_f(attr_id, hdferr)
!call H5Tclose_f(datatype_id, hdferr)
!call H5Sclose_f(dataspace_id, hdferr)

call H5Util_writeScalarAttribute(sds%dataset_id, 'long_name', sds%title, hdferr)


!
! units:
!
!dims(1) = 1
!attr_len = len_trim(sds%units)
!call H5Screate_f(H5S_SCALAR_F, dataspace_id, hdferr)
!call H5Tcopy_f(H5T_NATIVE_CHARACTER, datatype_id, hdferr)
!call H5Tset_strpad_f(datatype_id, H5T_STR_NULLTERM_F, hdferr)
!call H5Tset_size_f(datatype_id, attr_len, hdferr)
!call H5Acreate_f(sds%dataset_id, "units", &
!                 datatype_id, dataspace_id, attr_id, hdferr)
!call H5Awrite_f(attr_id, datatype_id, sds%units, dims, hdferr)
!call H5Aclose_f(attr_id, hdferr)
!call H5Tclose_f(datatype_id, hdferr)
!call H5Sclose_f(dataspace_id, hdferr)

call H5Util_writeScalarAttribute(sds%dataset_id, 'units', sds%units, hdferr)




!
! ACDD coverage content type (DISC):
! http://wiki.esipfed.org/index.php/Attribute_Convention_for_Data_Discovery
!     image, thematicClassification, physicalMeasurement, 
!     auxiliaryInformation, qualityInformation, referenceInformation, 
!     modelResult, or coordinate
!

!dims(1) = 1
attr_len = len_trim(sds%content_type)
if(attr_len > 1) then
!  call H5Screate_f(H5S_SCALAR_F, dataspace_id, hdferr)
!  call H5Tcopy_f(H5T_NATIVE_CHARACTER, datatype_id, hdferr)
!  call H5Tset_strpad_f(datatype_id, H5T_STR_NULLTERM_F, hdferr)
!  call H5Tset_size_f(datatype_id, attr_len, hdferr)
!  call H5Acreate_f(sds%dataset_id, "coverage_content_type", &
!                   datatype_id, dataspace_id, attr_id, hdferr)
!  call H5Awrite_f(attr_id, datatype_id, sds%content_type, dims, hdferr)
!  call H5Aclose_f(attr_id, hdferr)
!  call H5Tclose_f(datatype_id, hdferr)
!  call H5Sclose_f(dataspace_id, hdferr)
   call H5Util_writeScalarAttribute(sds%dataset_id, 'coverage_content_type', &
                                    sds%content_type, hdferr)
endif

attr_len = len_trim(sds%comment)
if(attr_len > 1) then
   call H5Util_writeScalarAttribute(sds%dataset_id, 'comment', &
                                    sds%comment, hdferr)
endif


!
! valid range:
!
dims(1) = 2
call H5Screate_simple_f(rank, dims, dataspace_id, hdferr)
call h5Tcopy_f(sds%datatype_id, datatype_id, hdferr)
call H5Acreate_f(sds%dataset_id, "valid_range", &
                 datatype_id, dataspace_id, attr_id, hdferr)
call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, sds%valid_range, dims, hdferr)
call H5Aclose_f(attr_id, hdferr)
call H5Tclose_f(datatype_id, hdferr)
call H5Sclose_f(dataspace_id, hdferr)

!
! fill value attribute:
!

!dims(1) = 1
!call H5Screate_simple_f(rank, dims, dataspace_id, hdferr)
!call h5Tcopy_f(sds%datatype_id, datatype_id, hdferr)
!call H5Acreate_f(sds%dataset_id, "_FillValue", &
!                 datatype_id, dataspace_id, attr_id, hdferr)
!call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, sds%fill_value, dims, hdferr)
!call H5Aclose_f(attr_id, hdferr)
!call H5Tclose_f(datatype_id, hdferr)
!call H5Sclose_f(dataspace_id, hdferr)

call H5Util_writeScalarAttribute(sds%dataset_id, '_FillValue', &
                                 sds%fill_value, hdferr, sds%datatype_id)

end subroutine writeDatasetAttributes

!*******************************************************************************

subroutine writeScalarAttr_char(h5ObjectId, attrName, attrValue, errorFlag)

character(len=*), parameter :: ROUTINE_NAME = 'H5Util.writeScalarAttr_char' 

integer(HID_T),   intent(in ) :: h5ObjectId
character(len=*), intent(in ) :: attrName
character(len=*), intent(in ) :: attrValue
integer,          intent(out) :: errorFlag

integer(HID_T) :: dataSpaceId, dataTypeId, attrId
integer(SIZE_T) :: attrlen
integer(HSIZE_T), DIMENSION(1) :: dims = (/1/)

errorFlag = 0

call H5Screate_f(H5S_SCALAR_F, dataSpaceId, errorFlag)
call H5Tcopy_f(H5T_NATIVE_CHARACTER, dataTypeId, errorFlag)
call H5Tset_strpad_f(dataTypeId, H5T_STR_NULLTERM_F, errorFlag)

attrlen = len_trim(attrValue)
if(attrlen == 0) attrlen = 1
call H5Tset_size_f(dataTypeId, attrlen, errorFlag)

call H5Acreate_f(h5ObjectId,attrName,dataTypeId,dataSpaceId,attrId, errorFlag)

call H5Awrite_f(attrId, dataTypeId, attrValue, dims, errorFlag)
if(errorFlag /= 0) then
   errorFlag = FAILURE_STATE
   call Display_Message(ROUTINE_NAME, &
               'error writing attribute ' // trim(attrName), errorFlag)
   goto 9999
endif

call H5Aclose_f(attrId, errorFlag)
call H5Tclose_f(datatypeId, errorFlag)
call H5Sclose_f(dataSpaceId, errorFlag)

9999 return

end subroutine writeScalarAttr_char


!*******************************************************************************

subroutine writeScalarAttr_short(h5ObjectId, attrName, attrValue, errorFlag, &
                                 h5ObjectDataTypeId)

character(len=*), parameter :: ROUTINE_NAME = 'writeScalarAttr_short' 
 
integer(HID_T),   intent(in ) :: h5ObjectId
character(len=*), intent(in ) :: attrName
integer(I2B),     intent(in ) :: attrValue
integer(HID_T),   intent(in), optional :: h5ObjectDataTypeId 
integer,          intent(out) :: errorFlag


integer :: buffer
integer(HID_T) :: dataSpaceId, dataTypeIdOnDisk, dataTypeIdInMemory, attrId
integer(HID_T) :: defaultDataTypeId 
integer(SIZE_T) :: attrlen
integer(HSIZE_T), DIMENSION(1) :: dims = (/1/)
 
defaultDataTypeId = H5T_STD_I16LE

include 'writeAttribute_template.inc'

end subroutine writeScalarAttr_short

!*******************************************************************************

subroutine writeScalarAttr_long(h5ObjectId, attrName, attrValue, errorFlag, &
                                h5ObjectDataTypeId)


character(len=*), parameter :: ROUTINE_NAME = 'writeScalarAttr_long' 

integer(HID_T),   intent(in ) :: h5ObjectId
character(len=*), intent(in ) :: attrName
integer(I4B),     intent(in ) :: attrValue
integer(HID_T),   intent(in), optional :: h5ObjectDataTypeId 
integer,          intent(out) :: errorFlag

integer :: buffer
integer(HID_T) :: dataSpaceId, dataTypeIdOnDisk, dataTypeIdInMemory, attrId
integer(HID_T) :: defaultDataTypeId 
integer(SIZE_T) :: attrlen
integer(HSIZE_T), DIMENSION(1) :: dims = (/1/)
 

defaultDataTypeId = H5T_NATIVE_INTEGER

include 'writeAttribute_template.inc'

end subroutine writeScalarAttr_long

!*******************************************************************************

subroutine writeScalarAttr_float(h5ObjectId, attrName, attrValue, errorFlag, &
                                 h5ObjectDataTypeId)

character(len=*), parameter :: ROUTINE_NAME = 'writeScalarAttr_float' 

integer(HID_T),   intent(in ) :: h5ObjectId
character(len=*), intent(in ) :: attrName
real(SP),         intent(in ) :: attrValue
integer(HID_T),   intent(in), optional :: h5ObjectDataTypeId 
integer,          intent(out) :: errorFlag

real(SP) :: buffer
integer(HID_T) :: dataSpaceId, dataTypeIdOnDisk, dataTypeIdInMemory, attrId
integer(HID_T) :: defaultDataTypeId 
integer(SIZE_T) :: attrlen
integer(HSIZE_T), DIMENSION(1) :: dims = (/1/)
 
defaultDataTypeId = H5T_NATIVE_REAL

include 'writeAttribute_template.inc'

end subroutine writeScalarAttr_float

!*******************************************************************************

subroutine writeScalarAttr_double(h5ObjectId, attrName, attrValue, errorFlag, &
                                  h5ObjectDataTypeId)

character(len=*), parameter :: ROUTINE_NAME = 'writeScalarAttr_double' 

integer(HID_T),   intent(in ) :: h5ObjectId
character(len=*), intent(in ) :: attrName
real(DP),         intent(in ) :: attrValue
! desired output data type differs from attrValue, fill value for example. But
! very very rarely do we need such conversion.
integer(HID_T),   intent(in), optional :: h5ObjectDataTypeId 
integer,          intent(out) :: errorFlag

real(DP) :: buffer
integer(HID_T) :: dataSpaceId, dataTypeIdOnDisk, dataTypeIdInMemory, attrId
integer(HID_T) :: defaultDataTypeId  
integer(SIZE_T) :: attrlen
integer(HSIZE_T), DIMENSION(1) :: dims = (/1/)
 
defaultDataTypeId = H5T_NATIVE_DOUBLE

include 'writeAttribute_template.inc'

end subroutine writeScalarAttr_double

!*******************************************************************************


END MODULE H5Util_class

