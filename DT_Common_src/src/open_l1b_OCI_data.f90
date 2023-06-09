 
module gather_l1b_OCI_data

      ! Declare use of Fortran 90 netCDF package
      ! If this part doesn't compile: module load sips/gcc
      use netcdf

      implicit none

contains
 
! ---------------------------------------------------------------------
!
! Begin subroutine open_l1b
!
!     Reads variables from VNP02MOD and VNP03MOD
!
! ---------------------------------------------------------------------

       subroutine open_l1b_OCI_data(Latitude_in,Longitude_in,SolarZenithAngle,SolarAzimuthAngle,&
                            ViewingZenithAngle,ViewingAzimuthAngle,TerrainHeight,&
                            UVtoSWIR_Reflectances,&
                            l1b_nLines,l1b_nXTrack,OLAT,&
                            OLON,OSOLA,OVIEW,OAZIM_View,OAZIM_Solar,&
                            OLand_ocean,OHeight,OW470,OW550,OW659,OW860,&
                            OW124,OW138,OW164,OW213,OW1100,OW659_Hres,OW354,OW388,OW412) 
 
      include 'read_Sat_MODIS.inc' 
      real, parameter :: pi = 3.1415927, dtr = pi/180.0
      integer ncid,grpid,dimid,varid,ifile,i,j,rad_flag 
      real sclfactor,offset,cossza,maxsolz,minsolz 
      integer ix,iy,jx,jy,Hqlines,Hqpixels 
      integer npixels,nlines,iobs,IXX,IYY
           npixels =  l1b_nXTrack
           nlines  =  l1b_nLines  
            
         ! Save datasets to output arrays
         OLAT(1:npixels,1:nlines)    = Latitude_in(1:npixels,1:nlines)
         OLON(1:npixels,1:nlines)     =  Longitude_in(1:npixels,1:nlines) 
         OHeight(1:npixels,1:nlines) = TerrainHeight(1:npixels,1:nlines) 
         OSOLA(1:npixels,1:nlines)   = SolarZenithAngle(1:npixels,1:nlines) 
         OVIEW(1:npixels,1:nlines)   = ViewingZenithAngle(1:npixels,1:nlines)
         OAZIM_View(1:npixels,1:nlines) = ViewingAzimuthAngle(1:npixels,1:nlines) 
         OAZIM_Solar(1:npixels,1:nlines) = SolarAzimuthAngle(1:npixels,1:nlines)   

           
! Save datasets to output arrays 
! 340, 354, 388, 412, 488, 55, 670, 680, 688, 860, 1250, 1378, 1615, 2260 nm. 
             
         OW354(1:npixels,1:nlines) = UVtoSWIR_Reflectances(2,1:npixels,1:nlines)
         OW388(1:npixels,1:nlines) = UVtoSWIR_Reflectances(3,1:npixels,1:nlines)
         OW412(1:npixels,1:nlines )= UVtoSWIR_Reflectances(4,1:npixels,1:nlines) 
         OW470(1:npixels,1:nlines) = UVtoSWIR_Reflectances(5,1:npixels,1:nlines)
         OW550(1:npixels,1:nlines) = UVtoSWIR_Reflectances(6,1:npixels,1:nlines)
         OW659(1:npixels,1:nlines) = UVtoSWIR_Reflectances(7,1:npixels,1:nlines)
         OW860(1:npixels,1:nlines) = UVtoSWIR_Reflectances(10,1:npixels,1:nlines)
         OW124(1:npixels,1:nlines) = UVtoSWIR_Reflectances(11,1:npixels,1:nlines)
         OW138(1:npixels,1:nlines) = UVtoSWIR_Reflectances(12,1:npixels,1:nlines) 
         OW164(1:npixels,1:nlines) = UVtoSWIR_Reflectances(13,1:npixels,1:nlines)
         OW213(1:npixels,1:nlines) = UVtoSWIR_Reflectances(14,1:npixels,1:nlines) 
          
!        OW1100(1:npixels,1:nlines) =Data_Viirs(10,1:npixels,1:nlines) 
         
          
                
                 
           Hqlines  = nlines*2 
           HQpixels = npixels*2
         
 ! making HRES 659....  Taking u0 out because it is already divided by U0        
              do j=1,Hqlines/2 
                  IY =  2*j-1
              do i=1,HQpixels/2
                 IX =  2*i-1 
!             cossza = cos(dtr*OSOLA(i,j))
                  do  jy = IY,2*j
                  do  jx = IX,2*i  
               OW659_Hres(jx,jy)= OW659(i,j)
                   enddo
                 enddo
              enddo
             enddo 
              
! save  Tropo data   
      end subroutine open_l1b_OCI_data
      
end module gather_l1b_OCI_data




