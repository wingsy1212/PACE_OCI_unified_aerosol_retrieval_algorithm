        module compute_Gascorrection_twoways
 
      implicit none

contains   
! ---------------------------------------------------------------------
!  compute compute_Gascorrection
           
           
           
           Subroutine compute_Gascorrection(Total_H2o,Total_O3,&
          SolZen,SatZen,set_counter_for_Gread,Multi_factor, RTN_NCEP)
          
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
      CHARACTER  (len=255) :: tbl,File_name
       
      
               
!       CALL get_environment_variable("Tables",tbl)  
            Nfile=30
! read the file  one time only
        if( set_counter_for_Gread .eq.1) then  
        
!        OPEN (Nfile, FILE = trim(tbl) // &
!              '/VIIRS_LUT/VIIRS_GAS_COEFS_LBL_v1.dat',& 
!                status = 'old',form='formatted') 

              file_name = cfg%dt_gas_coeff
              
              OPEN (Nfile, FILE = file_name,& 
                status = 'old',form='formatted')  
             
               read( Nfile,1)line 
               do Iwave=1,wave_num
        read(Nfile,*) mband,vband,wave,mol(iwave),&
        Opt_O3_Clim(iwave),Opt_H2O_Clim(iwave),&
        Opt_CO2_Clim(iwave),O3_Coef(iwave,1),O3_Coef(iwave,2),&
        H2o_Coef(iwave,1),H2o_Coef(iwave,2), H2o_Coef(iwave,3) 
               enddo   
          Endif
 
   1         format(132 a1)        
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
 end  subroutine compute_Gascorrection
           end module  compute_Gascorrection_twoways
           