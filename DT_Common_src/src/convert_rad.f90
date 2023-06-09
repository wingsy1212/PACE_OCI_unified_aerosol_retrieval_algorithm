 
module convert_rad
 
      implicit none

contains

! ---------------------------------------------------------------------
!  converting Rad into reflectances
!  
! ---------------------------------------------------------------------

       subroutine convert_rad_ref (l1b_nLines,l1b_nXTrack,OSOLA,&
       OW470,OW550,OW659,OW860,OW124,OW138,OW164,OW213,OW354,OW388,OW412)
       include 'read_Sat_MODIS.inc'  
       real, parameter :: pi = 3.1415927, dtr = pi/180.0
       real  cossza 
       integer i,j 
                           do j=1,l1b_nLines
                           do i=1,l1b_nXTrack
                           cossza = cos(dtr*OSOLA(i,j)) 
  ! for synthetic data                         
  !                          cossza = 1.
                           
  !                         OW354(i,j) = OW354(i,j)* PI
  !                         OW388(i,j) = OW388(i,j)* PI
                           OW412(i,j) = OW412(i,j)/cossza
                           OW470(i,j) = OW470(i,j)/cossza
                           OW550(i,j) = OW550(i,j)/cossza
                           OW659(i,j) = OW659(i,j)/cossza
                           OW860(i,j) = OW860(i,j)/cossza
                           OW124(i,j) = OW124(i,j)/cossza
                           OW138(i,j) = OW138(i,j)/cossza
                           OW164(i,j) = OW164(i,j)/cossza
                           OW213(i,j) = OW213(i,j)/cossza
                        enddo
                        enddo
                         
            return
            end  subroutine convert_rad_ref
            
                        
   end module  convert_rad
   
                           