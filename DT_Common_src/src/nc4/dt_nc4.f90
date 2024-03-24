
! ---------------------------------------------------------------------
!  compute compute_Gascorrection



Subroutine compute_Gascorrection_nc4(Total_H2o,Total_O3,&
   SolZen,SatZen,set_counter_for_Gread,Multi_factor, RTN_NCEP)

   use netcdf
   USE OCIUAAER_Config_Module
   implicit None
   save
            include 'read_Sat_MODIS.inc'
   real Total_H2o,G_factor,SatZen,SolZen,LOGCON,LOGCON2,Total_O3
   integer iwave,wave_num,ifile,mband,vband,RTN_NCEP,Nfile
   parameter(wave_num=10)
   real EXPONENT,H2o_Coef(wave_num,3), Opt_H2O_Clim(wave_num)
   real O3_Coef(wave_num,2),Opt_O3_Clim(wave_num)
   real Opt_CO2_Clim(wave_num), RTrans_H2O(wave_num)
   real RTrans_O3(wave_num),RTrans_CO2(wave_num)
   real Multi_factor(wave_num),drygas
   real wave,mol(wave_num), DEGRAD
   real g_factor_flat,g_factor_h2o,g_factor_o3
   real g_factor_co2,c_vz,c_sz,hh,r,aa,ab
   integer set_counter_for_Gread,Viirs_Table,ik,ij
   character(len=10) :: Sat_Flag
   character * 1 line(132)

   integer, dimension (1) :: start1, edge1, stride1
   integer, dimension (2) :: start2, edge2, stride2

   integer               ::  status
   character(len=255)    ::  sds_name
   character(len=255)    ::  dset_name
   character(len=255)    ::  attr_name
   character(len=255)    ::  group_name

   integer               ::  nc_id
   integer               ::  dim_id
   integer               ::  dset_id
   integer               ::  grp_id
   integer               ::  cnt

   if (cnt == 0) then
      status = nf90_open(cfg%dt_nc4, nf90_nowrite, nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
         return
      end if

      group_name = 'gas_coeffs'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'MOL'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start1  = (/ 1 /)
      edge1   = (/ wave_num /)
      stride1 = (/ 1 /)
      status = nf90_get_var(grp_id, dset_id, mol, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'OPT_O3_CLIM'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      status = nf90_get_var(grp_id, dset_id, opt_o3_clim, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'OPT_H2O_CLIM'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      status = nf90_get_var(grp_id, dset_id, opt_h2o_clim, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'OPT_CO2_CLIM'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      status = nf90_get_var(grp_id, dset_id, opt_co2_clim, start=start1, &
         stride=stride1, count=edge1)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'O3_COEF'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      start2  = (/ 1,1 /)
      edge2   = (/ wave_num,2 /)
      stride2 = (/ 1,1 /)
      status = nf90_get_var(grp_id, dset_id, o3_coef, start=start2, &
         stride=stride2, count=edge2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
      dset_name = 'H2O_COEF'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      edge2   = (/ wave_num,3 /)
      status = nf90_get_var(grp_id, dset_id, h2o_coef, start=start2, &
         stride=stride2, count=edge2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if
       status = nf90_close(nc_id)
       if (status /= NF90_NOERR) then
           print *, "ERROR: Failed to close lut_nc4 file: ", status
           return
       end if
       cnt = 1
   endif

   Do Iwave =1,wave_num
      RTrans_H2O(Iwave) =1.
      RTrans_O3(Iwave)=1.
      RTrans_CO2(Iwave)=1.
   ENDDO
   ! Calculate gemoetric factor for 2-way transmission
   !
   DEGRAD=ACOS(-1.)/180.
   G_factor=-1.0
   G_factor_flat = -1.0
   G_factor_H2O = -1.0
   G_factor_O3 = -1.0
   G_factor_CO2 = -1.0


   IF(SatZen.GT.0.0.AND.SolZen.GT.0.0) THEN
      c_VZ = COS(DEGRAD*SatZen)
      c_SZ = COS(DEGRAD*SolZen)

      !  Calculate G_factors

      !  plane parallel g_factor
      G_factor_flat = (1./c_VZ) + (1./c_SZ)
      !  Spherical geometry g_factor
      hh = 9.  ! 9 km atmos scale height
      r = 6371./9.
      G_factor= (SQRT( (r*c_VZ)**2. + 2.*r + 1) - r*c_VZ)&
         +  (SQRT( (r*c_SZ)**2. + 2.*r + 1) - r*c_SZ)
      !  Kasten and Young g_factors (  1. / cosz + a1*z**a2 * (a3-z)**a4 )
      G_factor_H2O = &
         (1./(c_VZ + 0.0311*SatZen**(0.1) * (92.471-SatZen)**(-1.3814)))&
         + (1./(c_SZ + 0.0311*SolZen**(0.1) * (92.471-SolZen)**(-1.3814)))
      G_factor_O3 = &
         (1./(c_VZ + 268.45*SatZen**(0.5) * (115.42-SatZen)**(-3.2922)))&
         +    (1./(c_SZ + 268.45*SolZen**(0.5) * (115.42-SolZen)**(-3.2922)))
      G_factor_CO2 = &
         (1./(c_VZ + 0.4567*SatZen**(0.07) * (96.4836-SatZen)**(-1.6970)))&
         +    (1./(c_SZ + 0.4567*SolZen**(0.07) * (96.4836-SolZen)**(-1.6970)))

   ENDIF

   ! keep like old versio........
   !          G_factor_H2O =G_factor_flat
   !          G_factor_O3 =G_factor_flat
   !          G_factor_CO2=G_factor_flat

   IF(RTN_NCEP .gt. 0)then
      Total_H2o=Total_H2o/10.
   Else
      Total_H2o= 0
   endif

   ! If NCEP water is available compute Water transmission else use OptH20 from clim.
   IF(RTN_NCEP .gt. 0 .and. Total_H2O.GT.0.0.AND.G_factor.GT.0.0) THEN
      LOGCON=ALOG(Total_H2O*G_factor_H2O)
      LOGCON2=LOGCON*LOGCON
      Do Iwave =1,wave_num
         EXPONENT=H2o_Coef(Iwave,1)+H2o_Coef(Iwave,2)*LOGCON &
            +H2o_Coef(Iwave,3)*LOGCON2
         RTrans_H2O(Iwave)=EXP(EXP(EXPONENT))
      Enddo
   Else
      Do Iwave =1,wave_num
         RTrans_H2O(Iwave)=EXP(Opt_H2O_Clim(iwave)*G_factor_H2O)
      Enddo
   Endif




   ! If NCEP Ozone is available compute Ozonetransmission else use OptOzone from clim.
   IF(RTN_NCEP .gt.0 .and.Total_O3.GT.0.0.AND.G_factor.GT.0.0) THEN
      Do Iwave =1,wave_num
         EXPONENT=Total_O3*G_factor_O3
         RTrans_O3(Iwave)=EXP(O3_Coef(Iwave,1)+O3_Coef(Iwave,2)*EXPONENT)
      Enddo
   Else
      Do Iwave =1,wave_num
         RTrans_O3(Iwave)=EXP(Opt_O3_Clim(iwave)*G_factor_O3)
      Enddo
   Endif
   ! compute rest of gases from cli.
   Do Iwave =1,wave_num
      RTrans_CO2(iwave)=EXP(Opt_CO2_Clim(Iwave)*G_factor_CO2)
   Enddo
   ! compute total transmission

   Do Iwave =1,wave_num
      Multi_factor(Iwave)=RTrans_H2O(Iwave)*RTrans_O3(Iwave)* &
         RTrans_CO2(Iwave)

   Enddo

   close(Nfile)

   Return
end  subroutine compute_Gascorrection_nc4


!*********************************************************************

SUBROUTINE READ_LOOK_NC4(RGSS,SIGMAS,EXTSMALL,MOMENTSSMALL,&
   CCNSMALL,EXTNORSMALL,BACKSCTTSMALL,ASSYMSMALL,&
   RGSB,SIGMAB,EXTBIG,MOMENTSBIG,EXTNORBIG,BACKSCTTBIG,&
   ASSYMBIG,ALBEDOSMALL,ALBEDOBIG,&
   ALBEDO_R_SMALL,ALBEDO_R_BIG,ALBEDO_T_SMALL,ALBEDO_T_BIG,&
   PHC,THET,THET0,AINTS,TAUAS,WAVE,&
   AINTB,TAUAB,JPHI,ref_rayall,HANDLE_S,HANDLE_L,&
   HANDLE_Ext_554_O,Ext_554_small,Ext_554_large)

   !---------------------------------------------------------------------
   !!F90

   !
   !!INPUT PARAMETERS:
   !          HANDLE_S        Input read for small particles(lookup).
   !          HANDLE_L        Input read for large particles(lookup).
   !
   !!OUTPUT PARAMETERS:
   !      SMALL MODE.........
   !
   !                  RGSS        RADIUS
   !            SIGMAS        SIGMA
   ! EXTSMALL       EXTINCTION COEFF
   !       MOMENTSSMALL       MOMENTS SMALL MODE
   ! CCNSMALL       CLOUD CONDENSATION NUCLII
   !       EXTNORSMALL       EXTINCTION COEFF NORMALIZED
   !     BACKSCTTSMALL       BACKSCATTERING RATIO
   !         ASSYMSMALL       ASSYMETRY FACTOR
   !        ALBEDOSMALL       ALBEDO SINGLE SCATTERING
   !     ALBEDO_R_SMALL       REFLECTED ALBEDO
   !     ALBEDO_T_SMALL      TRANSMITTED ALBEDO
   !     LARGE MODE..........
   !
   !                   RGSB       RADIUS
   !   SIGMAB       SIGMA
   !   EXTBIG       EXTINCTION COEFF
   !         MOMENTSBIG       MOMENTS SMALL MODE
   !          EXTNORBIG       EXTINCTION COEFF NORMALIZED
   !        BACKSCTTBIG       BACKSCATTERING RATIO
   ! ASSYMBIG       ASSYMETRY FACTOR
   !          ALBEDOBIG       ALBEDO SINGLE SCATTERING
   !       ALBEDO_R_BIG       REFLECTED ALBEDO
   !       ALBEDO_T_BIG       TRANSMITTED ALBEDO
   !
   ! PHC       Azimuth angle.
   !
   !                   THET       view angle.
   !
   !                  THET0       Solar zenith angle.
   !
   !                  TFLUX       dowanward flux
   !
   !                   AINTS       radiance(l/fo),fo=1 small mode
   !
   !                   AINTB       radiance(l/fo),fo=1 large mode
   !
   !                   TAUAS       optical thickness SMALL MODE
   !
   !                   TAUAB       optical thickness LARGE MODE
   !
   !                   NWAV       number of wavelength
   !
   !                   NTAU       number of opticl thickness.
   !
   !                   NTH0       number of solar zenith angle.
   !
   !                  NTHET       number of view angle.
   !
   !                   NPHI       number of azimuth
   !
   !                   JPHI       azimuth
   !        ref_rayall      radiance(l/fo),fo=1 FOR RAYLEIGH TAUA=0.0


   use netcdf
   USE OCIUAAER_Config_Module
   IMPLICIT NONE
   SAVE

      INCLUDE 'mod04.inc'

   CHARACTER*132 LINE
   REAL EXTSMALL(Lut_indx,Numcases,NWAV),EXTbig(Lut_indx,Numcaseb,NWAV)
   REAL RGSS(NUMCASES),SIGMAS(NUMCASES)
   REAL RGSB(NUMCASEB),SIGMAB(NUMCASEB)
   REAL MOMENTSSMALL(Lut_indx,numcases,4),MOMENTSBIG(Lut_indx,NUMCASEB,4)
   REAL CCNSMALL(Lut_indx,nUMCASES),TR
   REAL EXTNORSMALL(NUMCASES,NWAV),EXTNORBIG(NUMCASEB,NWAV)
   REAL BACKSCTTSMALL(Lut_indx,NUMCASES,NWAV),BACKSCTTBIG(Lut_indx,NUMCASEB,NWAV)
   REAL ASSYMSMALL(Lut_indx,NUMCASES,NWAV),ASSYMBIG(Lut_indx,NUMCASEB,NWAV)
   REAL ALBEDOSMALL(Lut_indx,NUMCASES,NWAV),ALBEDOBIG(Lut_indx,NUMCASEB,NWAV)
   REAL ALBEDO_R_SMALL(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
   REAL ALBEDO_R_BIG(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)
   REAL ALBEDO_T_SMALL(Lut_indx,NTH0,NTAU,NWAV,NUMCASES)
   REAL ALBEDO_T_BIG(Lut_indx,NTH0,NTAU,NWAV,NUMCASEB)
   REAL ALBEDO_R_RAY(NWAV,NTH0)
   REAL ALBEDO_T_RAY(NWAV,NTH0)
   INTEGER ICASE,IFILE,IPHI,ITAU,ITH,ITH0,IWAV,IJ
   INTEGER JPHI(NPHI),Num_lut,HANDLE_Ext_554_O,Nfile
   REAL  PHC(NPHI),THET(NTHET),THET0(NTH0),WAVE(NWAV)
   REAL  AINTS(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASES)
   REAL  TAUAS(NUMCASES,NWAV,NTAU)
   REAL  AINTB(Lut_indx, NPHI, NTHET, NTH0, NTAU, NWAV, NUMCASEB)
   REAL  TAUAB(NUMCASEB,NWAV,NTAU),DUMMY(NTH0)
   REAL  CCNdummy,ref_rayall(Lut_indx,NPHI,NTHET,NTH0,NWAV)
   REAL EFFRADSMALL(NUMCASES),EFFRADBIG(NUMCASEB),QSCT
   Real Ext_554_small(Lut_indx,NUMCASES),Ext_554_large(Lut_indx,NUMCASEB)
   character * 1 cWS
   character ( len=10):: Sat_Flag

   integer, dimension (1)                   :: start1, edge1, stride1
   integer, dimension (2)                   :: start2, edge2, stride2
   integer, dimension (3)                   :: start3, edge3, stride3
   integer, dimension (4)                   :: start4, edge4, stride4
   integer, dimension (5)                   :: start5, edge5, stride5
   integer, dimension (7)                   :: start7, edge7, stride7

   integer               ::  status
   character(len=255)    ::  sds_name
   character(len=255)    ::  dset_name
   character(len=255)    ::  attr_name
   character(len=255)    ::  group_name

   integer               ::  nc_id
   integer               ::  dim_id
   integer               ::  dset_id
   integer               ::  grp_id

   status = nf90_open(cfg%dt_nc4, nf90_nowrite, nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
      return
   end if

   group_name = 'ext_coeffs'
   status = nf90_inq_ncid(nc_id, group_name, grp_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
      return
   end if

   dset_name = 'EXT_COEFFS'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start2  = (/ 1,1 /)
   edge2   = (/ 4,NUMCASES /)
   stride2 = (/ 1,1 /)
   status = nf90_get_var(grp_id, dset_id, Ext_554_small, start=start2, &
      stride=stride2, count=edge2)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start2  = (/ 1,NUMCASES+1 /)
   edge2   = (/ 4,NUMCASEB /)
   stride2 = (/ 1,1 /)
   status = nf90_get_var(grp_id, dset_id, Ext_554_large, start=start2, &
      stride=stride2, count=edge2)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   group_name = 'viirs_ocean_aerosol'
   status = nf90_inq_ncid(nc_id, group_name, grp_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
      return
   end if

   dset_name = 'THET0'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NTH0 /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, THET0, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'THET'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NTHET /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, THET, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'PHI'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NPHI /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, PHC, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   JPHI = PHC

   dset_name = 'WAVE'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NWAV /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, WAVE, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'RGS'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NUMCASES /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, RGSS, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ NUMCASES+1 /)
   edge1   = (/ NUMCASEB /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, RGSB, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'SIGMA'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NUMCASES /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, SIGMAS, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ NUMCASES+1 /)
   edge1   = (/ NUMCASEB /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, SIGMAB, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'EFFRAD'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NUMCASES /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, EFFRADSMALL, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ NUMCASES+1 /)
   edge1   = (/ NUMCASEB /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, EFFRADBIG, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'MOMENTS'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,1,1 /)
   edge3   = (/ LUT_INDX,NUMCASES,4 /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, MOMENTSSMALL, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,NUMCASES+1,1 /)
   edge3   = (/ LUT_INDX,NUMCASEB,4 /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, MOMENTSBIG, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'CCN'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start2  = (/ 1,1 /)
   edge2   = (/ LUT_INDX,NUMCASES /)
   stride2 = (/ 1,1 /)
   status = nf90_get_var(grp_id, dset_id, CCNSMALL, start=start2, &
      stride=stride2, count=edge2)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'ALBEDO'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,1,1 /)
   edge3   = (/ LUT_INDX,NUMCASES,NWAV /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ALBEDOSMALL, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,NUMCASES+1,1 /)
   edge3   = (/ LUT_INDX,NUMCASEB,NWAV /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ALBEDOBIG, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'ASSYM'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,1,1 /)
   edge3   = (/ LUT_INDX,NUMCASES,NWAV /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ASSYMSMALL, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,NUMCASES+1,1 /)
   edge3   = (/ LUT_INDX,NUMCASEB,NWAV /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ASSYMBIG, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'EXT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,1,1 /)
   edge3   = (/ LUT_INDX,NUMCASES,NWAV /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, EXTSMALL, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,NUMCASES+1,1 /)
   edge3   = (/ LUT_INDX,NUMCASEB,NWAV /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, EXTBIG, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'BACKSCTT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,1,1 /)
   edge3   = (/ LUT_INDX,NUMCASES,NWAV /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, BACKSCTTSMALL, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,NUMCASES+1,1 /)
   edge3   = (/ LUT_INDX,NUMCASEB,NWAV /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, BACKSCTTBIG, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'TAUA'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,1,1 /)
   edge3   = (/ NUMCASES,NWAV,NTAU /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, TAUAS, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ NUMCASES+1,1,1 /)
   edge3   = (/ NUMCASEB,NWAV,NTAU /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, TAUAB, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'ALBEDO_R'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start5  = (/ 1,1,1,1,1 /)
   edge5   = (/ LUT_INDX,NTH0,NTAU,NWAV,NUMCASES /)
   stride5 = (/ 1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ALBEDO_R_SMALL, start=start5, &
      stride=stride5, count=edge5)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start5  = (/ 1,1,1,1,NUMCASES+1 /)
   edge5   = (/ LUT_INDX,NTH0,NTAU,NWAV,NUMCASEB /)
   stride5 = (/ 1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ALBEDO_R_BIG, start=start5, &
      stride=stride5, count=edge5)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'ALBEDO_T'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start5  = (/ 1,1,1,1,1 /)
   edge5   = (/ LUT_INDX,NTH0,NTAU,NWAV,NUMCASES /)
   stride5 = (/ 1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ALBEDO_T_SMALL, start=start5, &
      stride=stride5, count=edge5)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start5  = (/ 1,1,1,1,NUMCASES+1 /)
   edge5   = (/ LUT_INDX,NTH0,NTAU,NWAV,NUMCASEB /)
   stride5 = (/ 1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ALBEDO_T_BIG, start=start5, &
      stride=stride5, count=edge5)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'REF_RAYALL'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start5  = (/ 1,1,1,1,1 /)
   edge5   = (/ LUT_INDX,NPHI,NTHET,NTH0,NWAV /)
   stride5 = (/ 1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, REF_RAYALL, start=start5, &
      stride=stride5, count=edge5)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'AINT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start7  = (/ 1,1,1,1,1,1,1 /)
   edge7   = (/ LUT_INDX,NPHI,NTHET,NTH0,NTAU,NWAV,NUMCASES /)
   stride7 = (/ 1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, AINTS, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start7  = (/ 1,1,1,1,1,1,NUMCASES+1 /)
   edge7   = (/ LUT_INDX,NPHI,NTHET,NTH0,NTAU,NWAV,NUMCASEB /)
   stride7 = (/ 1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, AINTB, start=start7, &
      stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   status = nf90_close(nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to close lut_nc4 file: ", status
      return
   end if

   RETURN
END



SUBROUTINE  READ_LOOK_UV_NC4(PHC,THET,THET0,&
   H1_AINTS_uv1,H1_AINTS_uv2,H2_AINTS_uv1,H2_AINTS_uv2,&
   H3_AINTS_uv1,H3_AINTS_uv2,H4_AINTS_uv1,H4_AINTS_uv2,TAUAS,WAVE,&
   H1_AINTB_uv1,H1_AINTB_uv2,H2_AINTB_uv1,H2_AINTB_uv2,&
   H3_AINTB_uv1,H3_AINTB_uv2,H4_AINTB_uv1,H4_AINTB_uv2,&
   TAUAB,JPHI,ref_rayall_uv,Iomega,&
   EXTSMALL,ALBEDOSMALL,EXTbig,ALBEDOBIG,Iheight,Num_iomega,&
   NUM_HEIGHT)

   use netcdf
   USE OCIUAAER_Config_Module
   IMPLICIT NONE
   SAVE
   INCLUDE 'mod04.inc'

   INTEGER Icase,IFILE,IPHI,ITAU,ITH,ITH0,&
      IWAV,IJ,Iopt_read,NNwave,Num_iomega,NUM_HEIGHT
   REAL  EXTSMALL(Num_iomega,NUM_HEIGHT,Numcases,NWAV_uv)
   REAL  ALBEDOSMALL(Num_iomega,NUM_HEIGHT,Numcases,NWAV_uv)
   REAL  EXTBIG(Num_iomega,NUM_HEIGHT,NUMCASEB,NWAV_uv)
   REAL  ALBEDOBIG(Num_iomega,Num_height,NUMCASEB,NWAV_uv)

   INTEGER JPHI(NPHI),IWS,Iomega,Mode_F,Mode_C,Iheight
   REAL  PHC(NPHI),THET(NTHET),THET0(NTH0),WAVE(NWAV)
   REAL  TAUAS(NUMCASES,NWAV,NTAU)
   REAL  TAUAB(NUMCASEB,NWAV,NTAU)
   real H1_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
   real H2_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
   real H3_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
   real H4_AINTS_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
   real H1_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
   real H2_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
   real H3_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
   real H4_AINTS_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASES)
   real H1_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
   real H2_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
   real H3_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
   real H4_AINTB_uv1(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
   real H1_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
   real H2_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
   real H3_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
   real H4_AINTB_uv2(Num_iomega,NWIND, NPHI, NTHET, NTH0, NTAU,NUMCASEB)
   Real  ref_rayall(Num_iomega,NUM_HEIGHT,NWIND,NPHI,NTHET,NTH0,NWAV),AA
   Real  ref_rayall_uv(Num_iomega,NWIND,NPHI,NTHET,NTH0,NWAV,NUM_HEIGHT)
   CHARACTER*1 cWS,kws,hhh
   integer  status
   integer, dimension (1)  :: start1, edge1, stride1
   integer, dimension (3)  :: start3, edge3, stride3
   integer, dimension (4)  :: start4, edge4, stride4
   integer, dimension (7)  :: start7, edge7, stride7
   integer, dimension (9)  :: start9, edge9, stride9

   character(len=255)     ::  sds_name
   character(len=255)    ::  dset_name
   character(len=255)    ::  attr_name
   character(len=255)    ::  group_name

   integer               ::  nc_id
   integer               ::  dim_id
   integer               ::  dset_id
   integer               ::  grp_id

   status = nf90_open(cfg%dt_nc4, nf90_nowrite, nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
      return
   end if

   group_name = 'pace_ocean_aerosol'
   status = nf90_inq_ncid(nc_id, group_name, grp_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
      return
   end if

   dset_name = 'THET0'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NTH0 /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, THET0, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'THET'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NTHET /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, THET, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'PHI'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NPHI /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, PHC, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   JPHI = PHC

   dset_name = 'WAVE'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NWAV_UV /)
   stride1 = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, WAVE, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   !  ----

   dset_name = 'ALBEDO'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start4  = (/ 1,1,1,1 /)
   edge4   = (/ NUM_IOMEGA,NUM_HEIGHT,NUMCASES,NWAV_UV /)
   stride4 = (/ 1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ALBEDOSMALL, start=start4, &
      stride=stride4, count=edge4)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start4  = (/ 1,1,NUMCASES+1,1 /)
   edge4   = (/ NUM_IOMEGA,NUM_HEIGHT,NUMCASEB,NWAV_UV /)
   stride4 = (/ 1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ALBEDOBIG, start=start4, &
      stride=stride4, count=edge4)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   !  ----

   dset_name = 'EXT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start4  = (/ 1,1,1,1 /)
   edge4   = (/ NUM_IOMEGA,NUM_HEIGHT,NUMCASES,NWAV_UV /)
   stride4 = (/ 1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, EXTSMALL, start=start4, &
      stride=stride4, count=edge4)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start4  = (/ 1,1,NUMCASES+1,1 /)
   edge4   = (/ NUM_IOMEGA,NUM_HEIGHT,NUMCASEB,NWAV_UV /)
   stride4 = (/ 1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, EXTBIG, start=start4, &
      stride=stride4, count=edge4)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   !  ----

   dset_name = 'TAUA'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,1,1 /)
   edge3   = (/ NUMCASES,NWAV_UV,NTAU /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, TAUAS(:,1:NWAV_UV,:), start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ NUMCASES+1,1,1 /)
   edge3   = (/ NUMCASEB,NWAV_UV,NTAU /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, TAUAB(:,1:NWAV_UV,:), start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   !  ----

   dset_name = 'REF_RAYALL'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if

   start7  = (/ 1,1,1,1,1,1,1 /)
   edge7   = (/ NUM_IOMEGA,NWIND,NPHI,NTHET,NTH0,NWAV_UV,NUM_HEIGHT /)
   stride7 = (/ 1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, ref_rayall_uv(:,:,:,:,:,1:NWAV_UV,:), &
      start=start7, stride=stride7, count=edge7)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   !  ----

   dset_name = 'AINT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if

   !   ---- H1

   start9  = (/ 1,1,1,1,1,1,1,1,1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASES /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H1_AINTS_UV1, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start9  = (/ 1,1,1,1,1,1,1,1,NUMCASES+1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASEB /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H1_AINTB_UV1, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   start9  = (/ 1,1,1,1,1,1,1,2,1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASES /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H1_AINTS_UV2, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start9  = (/ 1,1,1,1,1,1,1,2,NUMCASES+1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASEB /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H1_AINTB_UV2, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   !   ---- H2

   start9  = (/ 1,2,1,1,1,1,1,1,1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASES /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H2_AINTS_UV1, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start9  = (/ 1,2,1,1,1,1,1,1,NUMCASES+1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASEB /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H2_AINTB_UV1, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   start9  = (/ 1,2,1,1,1,1,1,2,1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASES /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H2_AINTS_UV2, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start9  = (/ 1,2,1,1,1,1,1,2,NUMCASES+1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASEB /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H2_AINTB_UV2, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if


   !   ---- H3

   start9  = (/ 1,3,1,1,1,1,1,1,1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASES /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H3_AINTS_UV1, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start9  = (/ 1,3,1,1,1,1,1,1,NUMCASES+1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASEB /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H3_AINTB_UV1, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   start9  = (/ 1,3,1,1,1,1,1,2,1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASES /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H3_AINTS_UV2, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start9  = (/ 1,3,1,1,1,1,1,2,NUMCASES+1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASEB /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H3_AINTB_UV2, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   !   ---- H4

   start9  = (/ 1,4,1,1,1,1,1,1,1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASES /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H4_AINTS_UV1, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start9  = (/ 1,4,1,1,1,1,1,1,NUMCASES+1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASEB /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H4_AINTB_UV1, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   start9  = (/ 1,4,1,1,1,1,1,2,1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASES /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H4_AINTS_UV2, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
   start9  = (/ 1,4,1,1,1,1,1,2,NUMCASES+1 /)
   edge9   = (/ NUM_IOMEGA,1,NWIND,NPHI,NTHET,NTH0,NTAU,1,NUMCASEB /)
   stride9 = (/ 1,1,1,1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, H4_AINTB_UV2, start=start9, &
      stride=stride9, count=edge9)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   status = nf90_close(nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to close lut_nc4 file: ", status
      return
   end if

   RETURN
END


!*************************************************
SUBROUTINE AEROSOL_MAP_NC4(HANDLE_LUTMAP, IMONTH,AEROSOL)
   !----------------------------------------------------
   !!F90

   !
   !!DESCRIPTION:  Subroutine AEROSOL_MAP reads
   !               the aerosol map (1 degree resolution)
   !               to determine expected aerosol type at
   !               a given location and season
   !
   !!INPUT PARAMETERS:
   !         HANDLE_LUTMAP   Logical Unit # for file
   !    IMONTH          Integer month 1-12
   !
   !!OUTPUT PARAMETERS
   !    AEROSOL     360x180 map of aerosol type for appropriate season

   use netcdf
   USE OCIUAAER_Config_Module
   IMPLICIT none

   INTEGER IMONTH, IMONTH1,iseason
   INTEGER HANDLE_LUTMAP
   INTEGER nlon, nlat, ilon, ilat
   PARAMETER (nlon = 360, nlat = 180)
   INTEGER AEROSOL(nlon,nlat)
   INTEGER AEROSOL_all(nlon,nlat,4)
   CHARACTER  (len=255) ::file_name ,Extension

   CHARACTER*3 MMMs(4), MMM
   DATA MMMs/'DJF','MAM','JJA','SON'/

   integer                :: status
   integer, dimension (3) :: start3, edge3, stride3

   character(len=255)    ::  sds_name
   character(len=255)    ::  dset_name
   character(len=255)    ::  attr_name
   character(len=255)    ::  group_name

   integer               ::  nc_id
   integer               ::  dim_id
   integer               ::  dset_id
   integer               ::  grp_id

   status = nf90_open(cfg%dt_nc4, nf90_nowrite, nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
      return
   end if

   group_name = 'aerosol_land_map'
   status = nf90_inq_ncid(nc_id, group_name, grp_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
      return
   end if

   dset_name = 'LAND_MAP'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,1,1 /)
   edge3   = (/ nlon,nlat,4 /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, AEROSOL_all, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if
    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to close lut_nc4 file: ", status
        return
    end if

   IMONTH1 = IMONTH
   IF (IMONTH .EQ. 12) THEN
      IMONTH1 = 0
   ENDIF
   iseason = IMONTH1/3 + 1
   MMM = MMMs(iseason)
   DO ilon = 1, nlon
      DO ilat = 1, nlat
         AEROSOL(ilon,ilat) = AEROSOL_all(ilon,ilat,iseason)
      ENDDO
   ENDDO

   RETURN
END

!***********************************************************************
SUBROUTINE RNLOOKUP_NC4(&
   HANDLE_LUT466,HANDLE_LUT553,HANDLE_LUT644,HANDLE_LUT213,&
   INT_NL0,Fd_NL0,T_NL0,OPTH_NL0,&
   SBAR_NL0,MASSCOEF_NL0,EXTNORM_NL0)

   !-----------------------------------------------------------------------
   !!F90
   !!DESCRIPTION:
   !               This subroutine reads 4 lookup tables
   !
   !!INPUT PARAMETERS:
   !
   !!OUTPUT PARAMETERS:
   !
   !          INT_NL0           Reflectance
   !          Fd_NL0           Flux down
   !          T_NL0           Transmission factor
   !          OPTH_NL0          Optical depth
   !          SBAR_NL0          Sbar
   !          MASSCOEF_NL0      Mass Concentration coeficient
   !          EXTNORM_NL0       Normalized extinction coeficient
   !
   !!OUTPUT PARAMETERS (INTO INCLUDE FILE)
   !          THET0           Array of LUT Solar Zenith Angles
   !          THE             Array of LUT Satellite View angles
   !          PHI             Array of LUT Azimuth angle differences
   !          WAV             Array of LUT wavelengths
   !
   !!DESIGN NOTES:
   !
   !  Read from look-up table.
   !   wav=wavelengths(2)(.466,.533,.644 and 2.13 micrometers)  4
   !   lthet0=number of theta0(solar zenith angle)              9
   !   thet0=(0,12,24,36,48,54,60,66,72   degrees)
   !   lthe =number of the(observation zenith angle)           16
   !   the=(0,6,.....72 degrees every 6 degrees)
   !   lphi =number of phi(observation azimuth angle)          31
   !   phi=(0,6,12,18,24,30,36,42,48,54,60,66,72,78,84,90,96,
   !        102,108,114,120,126,132,138,144,150,156,162,168,
   !       174,180)
   !   ltau =number of opth(optical thickness)                  7
   !   opth=(0.0,0.25,0.50,1.0,2.0,3.0,5.0)
   !   lhght = number of hght(observation heights)              1
   !   hght=(80.0 km)

   use netcdf
   USE OCIUAAER_Config_Module
   IMPLICIT NONE
   INCLUDE 'mod04_land.inc'
   INTEGER HANDLE_LUT466,HANDLE_LUT553
   INTEGER HANDLE_LUT644,HANDLE_LUT213
   INTEGER I,IPHI,ITAU,ITAB,NFILE
   INTEGER ITHE,ITHET0,IWAV
   REAL OPTH_NL0(NLTAU,NLWAV,NLTABLE)
   REAL MASSCOEF_NL0(NLTAU,NLWAV,NLTABLE)
   REAL EXTNORM_NL0(NLTAU,NLWAV,NLTABLE)
   REAL SSA_NL0(NLTAU,NLWAV,NLTABLE)
   REAL QEXT_NL0(NLTAU,NLWAV,NLTABLE)
   REAL BEXT_NL0(NLTAU,NLWAV,NLTABLE)
   REAL VEXT_NL0(NLTAU,NLWAV,NLTABLE)
   REAL MEXT_NL0(NLTAU,NLWAV,NLTABLE)
   REAL SBAR_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
   REAL INT_NL0(NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
   REAL Fd_NL0(NLTHET0,NLTAU,NLWAV,NLTABLE)
   REAL T_NL0(NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE)
   CHARACTER*1 LINE(132)
   character (len =10):: Sat_Flag
   REAL  Omega0(NLTAU,NLWAV,NLTABLE),ROD(NLWAV),GOD(NLWAV)
   SAVE

   integer                     :: status
   integer, dimension (1)      :: start1, edge1, stride1
   integer, dimension (3)      :: start3, edge3, stride3
   integer, dimension (4)      :: start4, edge4, stride4
   integer, dimension (5)      :: start5, edge5, stride5
   integer, dimension (6)      :: start6, edge6, stride6

   integer               ::  number_type, nattrs
   integer               ::  sds_id, sds_index, attr_index, hdfid
   character(len=255)    ::  sds_name
   character(len=255)    ::  dset_name
   character(len=255)    ::  attr_name
   character(len=255)    ::  group_name

   integer               ::  nc_id
   integer               ::  dim_id
   integer               ::  dset_id
   integer               ::  grp_id

   status = nf90_open(cfg%dt_nc4, nf90_nowrite, nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
      return
   end if

   group_name = 'land_aerosol'
   status = nf90_inq_ncid(nc_id, group_name, grp_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
      return
   end if

   dset_name = 'THET0'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NLTHET0 /)
   stride1  = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, THET0_NL, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'THET'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NLTHE /)
   stride1  = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, THE_NL, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'PHI'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NLPHI /)
   stride1  = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, PHI_NL, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'MU0'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NLTHET0 /)
   stride1  = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, MU0_NL, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'WAVE'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start1  = (/ 1 /)
   edge1   = (/ NLWAV /)
   stride1  = (/ 1 /)
   status = nf90_get_var(grp_id, dset_id, WAV_NL, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'ROD'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   status = nf90_get_var(grp_id, dset_id, ROD, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'GOD'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   status = nf90_get_var(grp_id, dset_id, GOD, start=start1, &
      stride=stride1, count=edge1)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'SSA'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start3  = (/ 1,1,1 /)
   edge3   = (/ NLTAU,NLWAV,NLTABLE /)
   stride3 = (/ 1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, SSA_NL0, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'QEXT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   status = nf90_get_var(grp_id, dset_id, QEXT_NL0, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'BEXT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   status = nf90_get_var(grp_id, dset_id, BEXT_NL0, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'VEXT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   status = nf90_get_var(grp_id, dset_id, VEXT_NL0, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'MEXT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   status = nf90_get_var(grp_id, dset_id, MEXT_NL0, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'MASSCOEF'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   status = nf90_get_var(grp_id, dset_id, MASSCOEF_NL0, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'OPTH'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   status = nf90_get_var(grp_id, dset_id, OPTH_NL0, start=start3, &
      stride=stride3, count=edge3)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'SBAR'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start4  = (/ 1,1,1,1 /)
   edge4   = (/ NLTHET0,NLTAU,NLWAV,NLTABLE /)
   stride4 = (/ 1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, SBAR_NL0, start=start4, &
      stride=stride4, count=edge4)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'FD'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   status = nf90_get_var(grp_id, dset_id, FD_NL0, start=start4, &
      stride=stride4, count=edge4)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'T'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start5  = (/ 1,1,1,1,1 /)
   edge5   = (/ NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE /)
   stride5 = (/ 1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, T_NL0, start=start5, &
      stride=stride5, count=edge5)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   dset_name = 'INT'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start6  = (/ 1,1,1,1,1,1 /)
   edge6   = (/ NLPHI,NLTHE,NLTHET0,NLTAU,NLWAV,NLTABLE /)
   stride6 = (/ 1,1,1,1,1,1 /)
   status = nf90_get_var(grp_id, dset_id, INT_NL0, start=start6, &
      stride=stride6, count=edge6)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   status = nf90_close(nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to close lut_nc4 file: ", status
      return
   end if

    QEXT_NL0(1,:,:) = QEXT_NL0(2,:,:)
    BEXT_NL0(1,:,:) = BEXT_NL0(2,:,:)
    VEXT_NL0(1,:,:) = VEXT_NL0(2,:,:)
    MEXT_NL0(1,:,:) = MEXT_NL0(2,:,:)
    MASSCOEF_NL0(1,:,:) = MASSCOEF_NL0(2,:,:)
    DO ITAB=1,NLTABLE
    DO ITAU=1,NLTAU
    DO IWAV=1,NLWAV
    EXTNORM_NL0(ITAU,IWAV,ITAB) = &
        QEXT_NL0(ITAU,IWAV,ITAB) / QEXT_NL0(ITAU,iwave_553,ITAB)
    ENDDO
    ENDDO
    ENDDO

   RETURN
END

subroutine  Read_urban_Table_nc4(Average_Urban, handle_Urban_Table_10km)
   use netcdf
   USE OCIUAAER_Config_Module
   IMPLICIT NONE
INCLUDE 'mod04.inc'
   Real Average_Urban(3601,1801)
   integer handle_Urban_Table_10km

   integer, dimension (2) :: start2, edge2, stride2

   integer               ::  status
   character(len=255)    ::  dset_name
   character(len=255)    ::  group_name

   integer               ::  nc_id
   integer               ::  dim_id
   integer               ::  dset_id
   integer               ::  grp_id
   integer               ::  cnt

   status = nf90_open(cfg%dt_nc4, nf90_nowrite, nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
      return
   end if

   group_name = 'urban_table'
   status = nf90_inq_ncid(nc_id, group_name, grp_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
      return
   end if

   dset_name = 'data'
   status = nf90_inq_varid(grp_id, dset_name, dset_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
   end if
   start2  = (/ 1, 1 /)
   edge2   = (/ 3601, 1801 /)
   stride2 = (/ 1, 1 /)
   status = nf90_get_var(grp_id, dset_id, Average_Urban, start=start2, &
      stride=stride2, count=edge2)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
   end if

   status = nf90_close(nc_id)
   if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to close lut_nc4 file: ", status
      return
   end if

   Return
end

