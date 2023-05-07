      module compute_Boxes_and_edges
      implicit none
contains   
! ---------------------------------------------------------------------
!     
!!SUBROUTINE set_Boxes_and_edges
       Subroutine set_Boxes_and_edges(l1b_nLines,l1b_nXTrack,Ret_Xtrack,Ret_Lines)
      include 'mod04.inc' 
      INCLUDE 'read_Sat_MODIS.inc'   
      integer yy1,yy2,xx1,LinesPerRead,read_once_option 
      
              Ret_Xtrack=l1b_nXTrack/ILINE
              Ret_Lines= l1b_nLines/ILINE  
               return
              end  subroutine set_Boxes_and_edges
           end module  compute_Boxes_and_edges
           
               
                