subroutine get_lut_igbp_land_cover( lut_file, status)
!
!f90
!description:
!   This subroutine reads core SDSs in the IGBP Land Cover
!   LUT hdf file.
!   file= 'MCD12C1.A2004001.005.Global_IGBP_Land_Cover_0.10deg.hdf'
!   includes:       IGBP_Land_Cover & Region_Index
!   array dimension [3600, 1800]

!   ver. 1.0  written by CES (jul 2010)
!   ver. 1.1  modified by MJ (aug 2011)
!
!   use GeneralAuxType
!   use lut_arrays
   
   implicit none
   
   include 'hdf.inc'
   include 'dffunc.inc'
  
   include 'newaottbl90.inc'
   
   character(*),           intent (in)      :: lut_file
   integer,                intent (out)     :: status
   integer, dimension (2)      		          :: start, edge, stride

!   integer(integer_fourbyte)		    :: number_type, nattrs
!   integer(integer_fourbyte)		    :: sds_id,sds_index,attr_index, hdfid
   integer		    :: number_type, nattrs
   integer		    :: sds_id,sds_index,attr_index, hdfid
   character(len=64)			          :: sds_name

   start  = (/ 0,0 /)
   edge   = (/ 3600,1800 /)
   stride = (/ 1,1 /)
    
   hdfid     = sfstart(lut_file,DFACC_READ)
   if (hdfid < 0) then
    print *, "ERROR: Unable to start veg MCD file: ", hdfid
    print *, "File: ", lut_file
    stop
   end if
   status    = sffinfo(hdfid, number_type, nattrs)
   if (status /= 0) then 
    print *, 'ERROR: Failed to get info on VEG MCD file: ', status
    stop
   end if
   sds_name  = 'IGBP_Land_Cover'
   sds_index = sfn2index(hdfid, sds_name)
   if (sds_index < 0) then 
    print *, "ERROR: Unable to get index of land cover SDS: ", status
    stop
   endif
   
   sds_id    = sfselect(hdfid , sds_index)
   if (sds_id < 0) then
    print *, "ERROR: Unable to select land cover SDS: ", sds_id
    stop
   end if
   
   status    = sfrdata(sds_id,start,stride,edge,xlcvr_2)
   if (status /= 0) then 
    print *, "ERROR: Unable to read veg land cover SDS: ", status
    stop
   end if
!   print *, xlcvr_2
!   print *, 'A:', status, hdfid, sds_index, sds_id
   status    = sfendacc(sds_id)
!  status    = sfend(hdfid) 
    
!  print *, 'AFTER'
!  print *, 'min ',minval(xlcvr_2)
!  print *, 'max ',maxval(xlcvr_2)
!  stop
  
!c... The lines below were added to read "Region_Index" in for selecting
!c... appropriate set of aerosol models in the over-vegetation retrievals.
!c... 08/05/2011 note by MJ.
!  hdfid     = sfstart(lut_file,DFACC_READ)
!  status    = sffinfo(hdfid, number_type, nattrs)
   sds_name  = 'Region_Index'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,regid_2)
!   print *, 'B:', status, hdfid, sds_index, sds_id
   status    = sfendacc(sds_id)

! close this HDF file
   status    = sfend(hdfid) 

!   print *, 'xlcvr_2(1800,900), xlcvr_2(1800,1000), xlcvr_2(1900,900)'
!   print *,  xlcvr_2(1800,900), xlcvr_2(1800,1000), xlcvr_2(1900,900)
!   print *, 'regid_2(1800,900), regid_2(1800,1000), regid_2(1900,900)'
!   print *,  regid_2(1800,900), regid_2(1800,1000), regid_2(1900,900) 

end subroutine get_lut_igbp_land_cover

