MODULE ErrorHandler_class

implicit none

contains

!*******************************************************************************

subroutine Display_Message(routineName, messageText, errorStatus) 

implicit none

integer, parameter :: MAX_N_STATES = 4

character(11), parameter, dimension(MAX_N_STATES ) :: &
               State_Descriptor = (/ 'SUCCESS    ',   &
                                     'WARNING    ',   &
                                     'FAILURE    ',   &
                                     'UNDEFINED  '   /)

character(len=*), intent(in) :: routineName
character(len=*), intent(in) :: messageText
integer, intent(in), optional :: errorStatus

integer :: item
character(len=1024) :: label
character(len=2) :: punctuation


punctuation = ""

write(*,*)

!
! label:
!

if(present(errorStatus)) then
  select case(errorStatus)
         case(0)
               item = 1
               punctuation = " !"
         case(1) 
               item = 2
               punctuation = " !"
         case(-1) 
               item = 3
               punctuation = " !"
         case default
               item = 4
               punctuation = ". "
  end select
  write(*,*) "[" // trim(routineName) // "]: " // trim(State_Descriptor(item)) 
else
  write(*,*) "[" // trim(routineName) // "]: " 
endif

!
! content:
!
write(*,*) trim(messageText) // punctuation

write(*,*)

end subroutine Display_Message

!*******************************************************************************

END MODULE ErrorHandler_class
