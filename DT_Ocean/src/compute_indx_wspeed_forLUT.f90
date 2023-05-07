      module compute_indx_wspeed_forLUT
 
      implicit none

contains   
! ---------------------------------------------------------------------
!     
!!SUBROUTINE indx_wspeed
! ---------------------------------------------------------------------

      Subroutine indx_wspeed(ugrd,vgrd,Indx_wspeed1,Indx_wspeed2,&
              WSPEED,Wind,RTN_NCEP)
          IMPLICIT NONE
          SAVE

       INCLUDE 'mod04.inc'
        real ugrd,vgrd,WSPEED,Del,real,Wind(Lut_indx) 
         integer Indx_wspeed1,Indx_wspeed2,ij,inn,RTN_NCEP 
           
           If( RTN_NCEP .gt.0) then
              WSPEED =SQRT(ugrd*ugrd+vgrd*vgrd)  
!   Set indexes for wind speed IF <  )WIND(IJ) = 2            
              Indx_wspeed1= 1
              Indx_wspeed2= 2  
           do ij=1,lut_indx 
             if( ij .eq.1)WIND(IJ)=2
             if( ij .eq.2)WIND(IJ)=6
             if( ij .eq.3)WIND(IJ)=10
             if( ij .eq.4)WIND(IJ)=14
           enddo  
         do ij=1,lut_indx 
          if( WIND(ij) .LE. WSPEED)  then
          Indx_wspeed1=ij
          Indx_wspeed2=ij+1 
          endif
         enddo  
        if(Indx_wspeed2 .gt. lut_indx)then
            Indx_wspeed1=lut_indx-1
            Indx_wspeed2=lut_indx 
        endif   
! Default is 6 meters          
         Else
            WSPEED =6
            Indx_wspeed1= 1
            Indx_wspeed2= 2
         Endif 
         
        RETURN
         end  subroutine indx_wspeed
           end module compute_indx_wspeed_forLUT
      
