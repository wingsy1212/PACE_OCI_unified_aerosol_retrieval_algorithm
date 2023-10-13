       
module read_GMAO_anc
        use netcdf

        implicit none

contains

! ---------------------------------------------------------------------
!
!  
!  
!
! --------------------------------------------------------------------- 
       
       
        Subroutine Get_An_GMAO(lat,lon,met_ugrd,met_vgrd,&
            met_pwat,met_ozone,met_skinTemp,set_counter_for_anc,RTN_NCEP,anc_file)
       USE Linear_interpolation           
       include 'read_Sat_MODIS.inc' 
       CHARACTER(255) :: anc_file
       character(len=3000) :: filename_Anc
       integer ncid,grpid,dimid,varid,nlines,ifile,i,j,ct,npixels,tt 
!         Real , dimension(:,:), allocatable :: U,V,W,O
!          Real , dimension(:), allocatable ::LA,LO
           Real , dimension(:,:), allocatable ::ST
          Real    U(360*2,181*2),V(360*2,181*2),W(360*2,181*2),O(360*2,181*2)
          Real    LA(360*2),LO(360*2) 
        Real  Uwind(360*2,181*2),vwind(360*2,181*2),water_Vapor(360*2,181*2),Tozone(360*2,181*2) 
        Real  Latitude(360*2),Longitude(360*2) 
        integer dx,dy,x0,y0,IX
       real		met_pwat,met_skinTemp
       real		met_ugrd
       real		met_vgrd,x(2),y(2),xx1(2),yy1(2),y1
       real		met_ozone,lat,lon  
       integer set_counter_for_anc,RTN_NCEP,j1,j2,i1,i2,ilat,ilon,nn,mm,init,ik,ij
       integer KTHEM1,KTHEP1,istart, iend, jstart, jend,Imet_Var,Number_Lines,Number_pixels 
       SAVE
      ! Open the netCDF file
!     read only once  
    
      if( set_counter_for_anc .eq.1)then
!       open(ifile, FILE=trim(anc_file),status='old')
!             read(ifile,101)filename_Anc
!             rewind ifile
!         print*,'filename_Anc',trim(filename_Anc),set_counter_for_anc 
!  101     format(a3000)
        
      call check( nf90_open(trim(anc_file), NF90_NOWRITE, ncid) )
      ! Dimension number_of_lines varies; look it up so that variable
      ! dimensions can be allocated for cloud mask array 
      call check( nf90_inq_dimid(ncid,'lon',dimid) )
      call check( nf90_inquire_dimension(ncid,dimid,len=npixels) )
      call check( nf90_inq_dimid(ncid,'lat',dimid) )
      call check( nf90_inquire_dimension(ncid,dimid,len=nlines) )
      call check( nf90_inq_dimid(ncid,'time',dimid) )
      call check( nf90_inquire_dimension(ncid,dimid,len=tt) )
       LA =  scale_Double(ncid,'lat',nlines) 
       Lo =  scale_Double(ncid,'lon',npixels) 
       U = scale_float(ncid,'U10M',npixels,nlines,tt) 
       V = scale_float(ncid,'V10M',npixels,nlines,tt)   
       W = scale_float(ncid,'TQV',npixels,nlines,tt)
       O = scale_float(ncid,'TO3',npixels,nlines,tt) 
       ST =scale_float(ncid,'TS',npixels,nlines,tt)
       Number_pixels = npixels
       Number_Lines  = nlines 
        do Ik = 1,nlines
        Latitude(ik) = LA(ik) 
        do IJ = 1,npixels
        Longitude(ij) = LO(IJ)
       Uwind(ij,ik) = U(ij,ik)
       Vwind(ij,ik) = V(ij,IK)  
       water_Vapor(ij,ik) = W(ij,ik)
       Tozone(ij,ik) = O(ij,ik)
          Enddo
          ENDDO 
      CALL check(nf90_close(ncid))
          endif  
             
          call Find_index_for_inter(lat,lon,Number_pixels,Number_Lines,i1,i2,j1,j2)    
           do   Imet_Var = 1,5
             MM=0  
             DO  init = 1,2
             y1 =0
             x(init) = 0
             y(init)=0
             xx1(init)=0
             yy1(init)=0
             enddo  
            Do Ilat = j1,j2
             NN =0
            Do Ilon = I1,I2
            NN= NN +1
            X(NN)= Longitude(Ilon)
            if( Imet_Var .eq.1)Y(NN)= water_Vapor(Ilon,Ilat)
            if( Imet_Var .eq.2)Y(NN)= Uwind(Ilon,Ilat)
            if( Imet_Var .eq.3)Y(NN)= vwind(Ilon,Ilat)
            if( Imet_Var .eq.4)Y(NN)= O(Ilon,Ilat) 
            if( Imet_Var .eq.5)Y(NN)= ST(Ilon,Ilat) 
            ENDDO
            CALL INTERP(NN,lon,X,Y,Y1)  
             MM = MM +1
             XX1(MM)=Latitude(Ilat)
             YY1(MM)=Y1 
            ENDDO
             y1=0.0
             CALL INTERP(MM,lat,XX1,YY1,Y1)
              if( Imet_Var .eq.1) met_pwat = y1
              if( Imet_Var .eq.2) met_ugrd  =y1
              if( Imet_Var .eq.3) met_vgrd  =y1
              if( Imet_Var .eq.4) met_ozone  =y1 
              if( Imet_Var .eq.5) met_skinTemp  =y1 
         Enddo  
         RTN_NCEP=RTN_NCEP+1
             
      end subroutine  Get_An_GMAO 
       
! ---------------------------------------------------------------------
!
! Begin subroutine check(status)
!
!     Checks for errors in netCDF handling
!
! ---------------------------------------------------------------------

  
       subroutine check(status)

      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
         print *,TRIM(nf90_strerror(status))
         stop "Stopped"
      end if 
      end subroutine check
     
      
! ---------------------------------------------------------------------
!
! Begin subroutine scale_short
!
!     Read variable from netCDF to array,
!     applying scale factor and offset
!
! ---------------------------------------------------------------------

      function scale_float(fid,varname,rsx,rsy,tt)

      integer fid,rsx,rsy,vid,i,j,ex,tt
      character(len=*) varname
      real, dimension(rsx,rsy,tt) :: vdata 
      real, dimension(360*2,181*2) :: rdata,scale_float
      real  scl,ofs,fv 
       
      ! Read data array and necessary attributes
      call check( nf90_inq_varid(fid,varname,vid) ) 
      call check( nf90_get_var(fid,vid,vdata) ) 
       
      call check( nf90_get_att(fid,vid,'scale_factor',scl) ) 
       
      call check( nf90_get_att(fid,vid,'add_offset',ofs) )
      
      call check( nf90_get_att(fid,vid,'_FillValue',fv) )
      
      do i=1,rsx
         do j=1,rsy 
            ! Apply scale factor, offset and new fill value
            if (vdata(i,j,tt) .ne. fv) then
                rdata(i,j) = vdata(i,j,tt)*scl + ofs
            else
               rdata(i,j) = -9999.
            endif  
         enddo
      enddo 
      scale_float = rdata

end function scale_float
 ! ---------------------------------------------------------------------
!
! Begin subroutine scale_short
!
!     Read variable from netCDF to array,
!     applying scale factor and offset
!
! ---------------------------------------------------------------------

      function scale_Double(fid,varname,rsx)

      integer fid,rsx,rsy,vid,i,j,ex,tt
      character(len=*) varname
      real*8, dimension(rsx) :: vdata 
      real, dimension(360*2) :: rdata,scale_Double
      real  scl,ofs,fv 
       
      ! Read data array and necessary attributes
      call check( nf90_inq_varid(fid,varname,vid) ) 
      call check( nf90_get_var(fid,vid,vdata) )  
      do i=1,rsx 
                rdata(i) = vdata(i)  
         enddo
      scale_Double = rdata 
      end function scale_Double
      
      Subroutine Find_index_for_inter(lat,lon,npixels,nlines,i1,i2,j1,j2)
            integer nlines,npixels,i1,i2,j1,j2,y
            real lat,lon,Dellat,dellon 
            Dellat = 0.50
            Dellon = 0.625
          y = (npixels/2)+(lon/Dellon)  
           i1=y
           i2=y+1   
        IF(y.LE.0)THEN
         i1=1
         i2=2
         ENDIF 
        IF(i2.GE.npixels)THEN
         i1=y-1
         i2=Y    
       endif    
        
          y = (nlines/2)+(lat/Dellat)  
           j1=y
           j2=y+1
        IF(y.LE.0)THEN
         j1=1
         j2=2
         ENDIF
        IF(j2.GE.nlines)THEN
         j1=y-1
         j2=Y   
       endif   
         Return
         end   
  end module read_GMAO_anc 
       
