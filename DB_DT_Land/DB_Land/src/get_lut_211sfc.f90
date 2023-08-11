subroutine get_lut_211sfc( lut_file, status)
!
!f90
!description:
!   This subroutine reads LUTs for 2.1um sfc refl. calculation 
!   file= 'nvalx21um4sfc.hdf'
!   array dimension [10,46,30,4]

!   ver. 1.0  modified by MJ (aug 2011)
!
!   use GeneralAuxType
!   use lut_arrays

   implicit none
   include 'hdf.inc'
   include 'dffunc.inc'
   include 'sfc21tbl90.inc'
   
   character(*),           intent (in)      :: lut_file
   integer,                intent (out)     :: status
   integer, dimension (4)      	            :: start, edge, stride

   integer		    :: number_type, nattrs
   integer		    :: sds_id,sds_index,attr_index, hdfid
   character(len=64)			    :: sds_name

   start  = (/ 0,0,0,0 /)
   edge   = (/ 10,46,30,4 /)
   stride = (/ 1,1,1,1 /)
    
   hdfid     = sfstart(lut_file,DFACC_READ)
   status    = sffinfo(hdfid, number_type, nattrs)

   sds_name  = 'NVALX21_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,nvalx21)

!  print *, 'A:', status, hdfid, sds_index, sds_id
!  status    = sfendacc(sds_id)
  
   sds_name  = 'R0X21_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,r0x_21)

   sds_name  = 'SX21_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,sx_21)

   sds_name  = 'TX21_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,tx_21)

!jlee added

   sds_name  = 'NVALX672_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,nvalx672)

   sds_name  = 'R0X672_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,r0x_672)

   sds_name  = 'SX672_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,sx_672)

   sds_name  = 'TX672_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,tx_672)

   sds_name  = 'NVALX865_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,nvalx865)

   sds_name  = 'R0X865_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,r0x_865)

   sds_name  = 'SX865_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,sx_865)

   sds_name  = 'TX865_SFC'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status    = sfrdata(sds_id,start,stride,edge,tx_865)

!end jlee added

!  print *, 'B:', status, hdfid, sds_index, sds_id
!  status    = sfendacc(sds_id)

!c... 
!   print *, 'test reading nvalx21 ....'
!   print *, nvalx21(5,5,5,2),r0x_21(5,5,5,2),tx_21(5,5,5,2),sx_21(5,5,5,2)
!. print *, nvalx21(5,5,5,4),r0x_21(5,5,5,4),tx_21(5,5,5,4),sx_21(5,5,5,4)
!. print *, 'end of test reading nvalx21 ....'


! close this HDF file
   status    = sfend(hdfid) 

end subroutine get_lut_211sfc

