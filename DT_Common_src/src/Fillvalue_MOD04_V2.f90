       SUBROUTINE FILLVALUE_Ocean(IDATA,SDS_ref,SDS_ref_STD,&
      SDSTAU_best,SDSTAUS_best,SDSTAUB_best,SDSTAU_average,&
      SDSTAUS_average,SDSTAUB_average,SDS_Least_error,&
      SDS_small_weighting,SDS_sol_INDX_small,SDS_sol_INDX_large,&
      SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,&
      SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,&
      SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,&
      SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,&
      SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_ocean,&
      SDS_Tau_Land_Ocean_img,Qcontrol_special,SDS_correc_small_weighting,&
      Dust_flag_10KM)

!-----------------------------------------------------------------------
! !F77
!
! !DESCRIPTION: This subroutine Fills the HDF array's with fill values.
!
! !INPUT PARAMETERS:
! IDATA        Index of box number
! !OUTPUT PARAMETERS:
! SDS_ref         Averaged Reflectance array
! SDS_ref_STD     Standard Deviation of Reflectance
! SDSTAU_best     Optical thickness for best solution
! SDSTAUS_best    Optical thickness contribution small particles for best solution
! SDSTAUB_best    Optical thickness contribution large particles for best solution
! SDSTAU_average  Optical thickness for best solution
! SDSTAUS_average Optical thickness contribution small particles for best solution
! SDSTAUB_average Optical thickness contribution large particles for best solution
! SDS_Least_error    Minimm Error function betwwen derived and computed radiances
! SDS_small_weighting Weight factor for large and small mode
! SDS_sol_INDX       Index for solution number
! SDSASSY_best       Assymetry Factor for best solution
! SDSASSY_average    Assymetry Factor for Average solution
! SDSBACK_best       Backscattering Ratio for best solution
! SDSBACK_average    Backscattering Ratio for Average solution
! SDS_effrad         Effective Radiance
! SDS_AOT_model Ration of optical Thickess small
! SDS_RefF_best      Reflected Flux Best solution
! SDS_RefF_average   Reflected Flux Average solution
! SDS_TranF_best     Transmitted Flux Best solution
! SDS_TranF_average  Transmitted Flux Average solution
! QCONTROL           Value for Quality Control
! SDS_SCAT_ANGLE_OCEAN Scattering angle ocean
! SDS_QCONTROL         Quality control SDS array
!SDS_NUMPIXELS        Number of Pixels used
! 

       IMPLICIT NONE
       SAVE

       INCLUDE 'mod04.inc'
       INCLUDE 'read_Sat_MODIS.inc'

       BYTE    SDS_QCONTROL_ocean(QA_LAND,NUMCELLS_B)
       INTEGER  IDATA,ICASE,IWAV,Dust_flag_10KM(NUMCELLS_B)
       Real       SDS_ref(NUMCELLS_B,NWAV+3),SDS_ref_STD(NUMCELLS_B,NWAV+3)
       Real        SDSTAU_best(NUMCELLS_B,NWAV),&
                  SDSTAU_average(NUMCELLS_B,NWAV),&
                  SDSTAUB_best(NUMCELLS_B,NWAV),&
                  SDSTAUB_average(NUMCELLS_B,NWAV),&
                  SDSTAUS_best(NUMCELLS_B,NWAV),&
                  SDSTAUS_average(NUMCELLS_B,NWAV),&
                  SDSASSY_best(NUMCELLS_B,NWAV),&
                  SDSASSY_average(NUMCELLS_B,NWAV),&
                  SDSBACK_best(NUMCELLS_B,NWAV),&
                  SDSBACK_average(NUMCELLS_B,NWAV),&
                  SDS_RefF_best(NUMCELLS_B,NWAV),&
                  SDS_RefF_average(NUMCELLS_B,NWAV),&
                  SDS_TranF_best(NUMCELLS_B,NWAV),&
                  SDS_TranF_average(NUMCELLS_B,NWAV)
                  

      Real        SDS_small_weighting(NUMCELLS_B,NUM_solutions),&
                  SDS_correc_small_weighting(NUMCELLS_B),&
                  SDS_Least_error(NUMCELLS_B,NUM_solutions),&
                  SDS_effrad(NUMCELLS_B,NUM_solutions),& 
                  SDS_angs_coeff1(NUMCELLS_B,NUM_solutions),&
                  SDS_angs_coeff2(NUMCELLS_B,NUM_solutions),&
                 SDS_AOT_model(NUMCELLS_B,num_model_index),&
                  SDS_CLDFRC_ocean(NUMCELLS_B),&
                  SDS_Tau_Land_Ocean_img(NUMCELLS_B)
       REAL SDS_mass_conc(NUMCELLS_B,NUM_solutions),&
            SDS_ccn(NUMCELLS_B,NUM_solutions)
        INTEGER  INT_Fill_value,Qcontrol_special
        integer    SDS_sol_INDX_small(NUMCELLS_B,NUM_solutions),&
                   SDS_sol_INDX_large(NUMCELLS_B,NUM_solutions) 
        REAL     FLOAT_Fill_value
       INTEGER   SDS_NUMPIXELS(NUMCELLS_B,NWAV+3),SDS_SCAT_ANGLE_OCEAN(NUMCELLS_B) 
!         Qcontrol_special .le.0 condition to fill everything with fill values 
          
          INT_Fill_value=-9999
          FLOAT_Fill_value=-999.0
          IF(Qcontrol_special .le.0 )then 
             SDS_SCAT_ANGLE_OCEAN(IDATA)=INT_Fill_value
            if(SDS_CLDFRC_ocean(IDATA) .lt.0  )SDS_CLDFRC_ocean(IDATA)=INT_Fill_value 
            Dust_flag_10KM(IDATA) = INT_Fill_value
            DO IWAV= 1,NWAV+3
            SDS_ref(IDATA,IWAV) =INT_Fill_value
            SDS_ref_STD(IDATA,IWAV) =INT_Fill_value
            SDS_NUMPIXELS(idata,IWAV)=INT_Fill_value
            Enddo
          DO IWAV= 1,NWAV 
            SDSASSY_best(IDATA,IWAV) =INT_Fill_value
            SDSASSY_average(IDATA,IWAV) =INT_Fill_value
            SDSBACK_BEST(IDATA,IWAV) =INT_Fill_value
            SDSBACK_AVERAGE(IDATA,IWAV) =INT_Fill_value 
            SDSTAU_BEST(IDATA,IWAV) = INT_Fill_value
            SDSTAUS_BEST(IDATA,IWAV) = INT_Fill_value
            SDSTAUB_BEST(IDATA,IWAV) = INT_Fill_value
            SDSTAU_AVERAGE(IDATA,IWAV) = INT_Fill_value
            SDSTAUS_AVERAGE(IDATA,IWAV) = INT_Fill_value
            SDSTAUB_AVERAGE(IDATA,IWAV) = INT_Fill_value
            ENDDO
            ENDIF
!        Qcontrol_special .eq.1 condition to fill optical thicknesses with zero for
!        ref at 0.865 is very small
             IF(Qcontrol_special .eq.1)then
!            SDS_SCAT_ANGLE_OCEAN(IDATA)=INT_Fill_value
!             if(SDS_CLDFRC_ocean(IDATA) .lt.0)SDS_CLDFRC_ocean(IDATA)=INT_Fill_value 
!            DO IWAV= 1,NWAV+3
!            SDS_ref(IDATA,IWAV) =INT_Fill_value
!            SDS_ref_STD(IDATA,IWAV) =INT_Fill_value
!            SDS_NUMPIXELS(idata,IWAV)=INT_Fill_value
!            Enddo
            Dust_flag_10KM(IDATA) = INT_Fill_value
            DO IWAV= 1,NWAV 
            SDSASSY_best(IDATA,IWAV) =INT_Fill_value
            SDSASSY_average(IDATA,IWAV) =INT_Fill_value
            SDSBACK_BEST(IDATA,IWAV) =INT_Fill_value
            SDSBACK_AVERAGE(IDATA,IWAV) =INT_Fill_value
            SDSTAU_BEST(IDATA,IWAV) = 0
            SDSTAUS_BEST(IDATA,IWAV) =0
            SDSTAUB_BEST(IDATA,IWAV) = 0
            SDSTAU_AVERAGE(IDATA,IWAV) = 0
            SDSTAUS_AVERAGE(IDATA,IWAV) = 0
            SDSTAUB_AVERAGE(IDATA,IWAV) = 0
          ENDDO
          ENDIF
          
           IF(Qcontrol_special .eq.2)then 
!            SDS_SCAT_ANGLE_OCEAN(IDATA)=INT_Fill_value 
!            SDS_CLDFRC_ocean(IDATA)=INT_Fill_value 
!            DO IWAV= 1,NWAV+3
!            SDS_ref(IDATA,IWAV) =INT_Fill_value
!            SDS_ref_STD(IDATA,IWAV) =INT_Fill_value
!            SDS_NUMPIXELS(idata,IWAV)=INT_Fill_value
!            Enddo
            Dust_flag_10KM(IDATA) = INT_Fill_value
            DO IWAV= 1,NWAV 
            SDSASSY_best(IDATA,IWAV) =INT_Fill_value
            SDSASSY_average(IDATA,IWAV) =INT_Fill_value
            SDSBACK_BEST(IDATA,IWAV) =INT_Fill_value
            SDSBACK_AVERAGE(IDATA,IWAV) =INT_Fill_value
            Enddo
          ENDIF
 
!         Qcontrol_special .eq.3 condition to fill Reflectances,sd and number of pixels
!          with real values over glint area 
!
            IF(Qcontrol_special .eq.3)then
            SDS_SCAT_ANGLE_OCEAN(IDATA)=INT_Fill_value
            SDS_CLDFRC_ocean(IDATA)=INT_Fill_value
            Dust_flag_10KM(IDATA) = INT_Fill_value
           DO IWAV= 1,NWAV
            SDSASSY_best(IDATA,IWAV) =INT_Fill_value
            SDSASSY_average(IDATA,IWAV) =INT_Fill_value
            SDSBACK_BEST(IDATA,IWAV) =INT_Fill_value
            SDSBACK_AVERAGE(IDATA,IWAV) =INT_Fill_value
            SDSTAU_BEST(IDATA,IWAV) = INT_Fill_value
            SDSTAUS_BEST(IDATA,IWAV) = INT_Fill_value
            SDSTAUB_BEST(IDATA,IWAV) = INT_Fill_value
            SDSTAU_AVERAGE(IDATA,IWAV) = INT_Fill_value
            SDSTAUS_AVERAGE(IDATA,IWAV) = INT_Fill_value
            SDSTAUB_AVERAGE(IDATA,IWAV) = INT_Fill_value
           ENDDO
           ENDIF
!  Fill with Fill values in all cases
           DO icase=1,NUM_solutions
            SDS_angs_coeff1(IDATA,ICASE) =INT_Fill_value
            SDS_angs_coeff2(IDATA,ICASE) =INT_Fill_value
            SDS_effrad(IDATA,ICASE) =INT_Fill_value
            SDS_Least_error(IDATA,ICASE) = INT_Fill_value
            SDS_small_weighting(IDATA,ICASE) = INT_Fill_value
            SDS_sol_INDX_small(IDATA,ICASE) = INT_Fill_value
            SDS_sol_INDX_large(IDATA,ICASE) = INT_Fill_value
            SDS_ccn(IDATA,ICASE) =FLOAT_Fill_value
            sds_mass_conc(IDATA,ICASE)= FLOAT_Fill_value
         ENDDO
!           SDS_correc_small_weighting(IDATA) = INT_Fill_value
           DO icase=1,num_model_index
            SDS_AOT_model(IDATA,ICASE) = INT_Fill_value
            enddo
            Qcontrol_special=0
           END


!***********************************************************************
      SUBROUTINE FILLVALUE_LAND(IDATA,SDS_Tau_Land_Ocean_img,&
      SDS_Aerosol_Type,SDSTAU_corrected_213,&
      SDS_SCAT_ANGLE_land,SDS_mass_conc_land,&
      SDS_angs_coeff_land,SDS_CLDFRC_land,SDS_dust_weighting,&
      SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,&
      SDS_NUMPIXELS_land,SDSTAU_corrected,SDS_ref_land,&
      SDS_ref_STD_land,SDS_QCONTROL_land,SDSTAU_small_land,&
      SDS_Surface_Reflectance_Land,&
      SDS_Fitting_Error_Land,Qcontrol_special_land)
!-----------------------------------------------------------------------
!!F77
!
!!DESCRIPTION:     This subroutine stores all the arrays to be written
!                  as output(HDF) file
!
!!INPUT PARAMETERS: all varaiables to be wrriten as output
!   LATM          Latitude of data cell
!   LONM          Longitude of data cell
!   IAER          Aerosol type
!   BARLBLUE      Average reflectance for blue channel
!   IDATA         Data cell number from 1 to NUMDATA
!   NUMDATA       Number of retrieval cells across orbit swath
!   TAUABLUE      Blue channel aerosol optical thickness (Continental)
!   BARLRED       Average reflectance for red channel
!   TAUARED       Red channel aerosol optical thickness (Continental)
!   SDBLUE        STD of blue channel reflectances
!   SDRED         STD of red channel reflectances
!   TAUABLUEC     Blue channel aerosol optical thickness (Corrected)
!   TAUAREDC      Red channel aerosol optical thickness (Corrected)
!   IBLUE         Number of blue channel observations
!   IRED          Number of red channel observations
!   IERROR        Error flag (0-4)
!   WEIGHT        Relative contribution of smoke/sulfate particles to
!                 dust in the computation of the aerosol optical depth
!   ISULP         Sulfate (1)/no sulfate flag (0)
!   ISMK          Smoke (1)/no smoke flag (0)
!   IPROCE        Aerosol retrieval procedure ID (0-4)
!   SAVERATIO     Aerosol path radiance ration (red to blue channel
!                 for continental model)
!   PIXELS        Number of cells across swath (same as NUMDATA)
!   LINES         Number of cells along swath (same as NUMSCAN)
!   NUMXBOX       Not used
!   NUMYBOX       Not used
!
!!OUTPUT PARAMETERS:ARRAY of 13 variables for output to HDF FILE
!  SDS1           HDF array of cell latitudes
!  SDS2           HDF array of cell longitudes
!  SDS3           HDF array of spectral reflectances
!  SDS4           HDF array of aerosol optical thicknesses for
!                     continental model
!  SDS5           HDF array of the STD of spectral reflectances
!  SDS6           HDF array of aerosol optical thicknesses for
!                     corrected model
!  SDS7           HDF array of STD for corrected optical thickness
!                     blue channel for continental model)
!  SDS8           HDF array of aerosol path radiance ratio (red to
!                     blue channel for continental model)
!  SDS9           HDF array of relative aerosol optical depth
!                     (smoke/sulfate particles to dust)
!  SDS10           HDF array of the number of blue channel clear pixels
!  SDS11          HDF array of the number of red channel clear pixels
!  SDS12          HDF array of retrieval procedure ID (0-4)
!  SDS13          HDF array of aerosol type (0-3)
!  SDS14          HDF array of Error flag (0-4)
 !--------------------------------------- 
      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'
      INCLUDE 'read_Sat_MODIS.inc'

      INTEGER IDATA,Num_Wave
      BYTE         SDS_QCONTROL_land(QA_LAND,NUMCELLS_B)
      REAL         SDS_mass_conc_land(NUMCELLS_B)
      Real         SDS_Tau_Land_Ocean_img(NUMCELLS_B),&
                   SDS_Aerosol_Type(NUMCELLS_B),&
                   SDS_SCAT_ANGLE_land(NUMCELLS_B),&
                   SDS_angs_coeff_land(NUMCELLS_B),&
                   SDS_CLDFRC_land(NUMCELLS_B),&
                   SDS_dust_weighting(NUMCELLS_B),& 
                   SDS_NUMPIXELS_land(NUMCELLS_B,Band_land+3),&
                   SDSTAU_corrected(NUMCELLS_B,Land_Sol3),&
                   SDS_ref_land(NUMCELLS_B,Band_land+3),&
                   SDS_ref_STD_land(NUMCELLS_B,Band_land+3),& 
                   SDS_Surface_Reflectance_Land(NUMCELLS_B,Land_Sol3),&
                   SDS_Fitting_Error_Land(NUMCELLS_B),&
                   SDSTAU_corrected_213(NUMCELLS_B),&
                   SDSTAU_small_land(NUMCELLS_B,Land_Sol4)  

!
! Obsolete (02/2006) Land SDS Arrays
!
      Real         &  
                 SDS_est_uncer(NUMCELLS_B,Land_Sol1),&
                  SDS_RefF_land(NUMCELLS_B,Land_Sol2),&
                  SDS_TranF_land(NUMCELLS_B,Land_Sol1)


      INTEGER  INT_Fill_value,Qcontrol_special_land
      REAL     FLOAT_Fill_value
 
! Set to Fill_values for integer and real variables
 

      INT_Fill_value=-9999
      FLOAT_Fill_value=-999.0


! Set Fill_Values to SDS arrays
 
      
       If ( Qcontrol_special_land .le.0) then
      if(SDS_CLDFRC_land(IDATA) .lt.0)SDS_CLDFRC_land(IDATA)=INT_Fill_value 
         SDS_mass_conc_land(IDATA)= FLOAT_Fill_value
        SDS_QCONTROL_land(1,IDATA)=0
        SDS_QCONTROL_land(2,IDATA)=0
        SDS_QCONTROL_land(3,IDATA)=0
        SDS_QCONTROL_land(4,IDATA)=0
        SDS_QCONTROL_land(5,IDATA)=0
         SDS_Aerosol_Type(IDATA)=INT_Fill_value
         SDS_SCAT_ANGLE_land(IDATA)=INT_Fill_value
         SDSTAU_corrected(IDATA,1)=INT_Fill_value
         SDSTAU_corrected(IDATA,2)=INT_Fill_value
         SDSTAU_corrected(IDATA,3)=INT_Fill_value
         SDSTAU_corrected_213(IDATA)=INT_Fill_value
         SDSTAU_small_land(IDATA,1)=INT_Fill_value
         SDSTAU_small_land(IDATA,2)=INT_Fill_value
         SDSTAU_small_land(IDATA,3)=INT_Fill_value
         SDSTAU_small_land(IDATA,4)=INT_Fill_value 
         SDS_angs_coeff_land(IDATA)=INT_Fill_value 
         SDS_dust_weighting(IDATA)=INT_Fill_value
         Do Num_Wave=1,Band_land+3
         SDS_NUMPIXELS_land(IDATA,Num_Wave)=INT_Fill_value 
         SDS_ref_land(IDATA,Num_Wave) =INT_Fill_value
         SDS_ref_STD_land(IDATA,Num_Wave)=INT_Fill_value
         Enddo 
         SDS_Surface_Reflectance_Land(IDATA,1) =INT_Fill_value
         SDS_Surface_Reflectance_Land(IDATA,2) =INT_Fill_value
         SDS_Surface_Reflectance_Land(IDATA,3) =INT_Fill_value 
         SDS_Fitting_Error_Land(IDATA)=INT_Fill_value 
          ENDIF
! if Path is B ie Iprocedure=2 then report only optical depth for 0.46 and 0.55 um
!
      IF(Qcontrol_special_land .eq.1)then
       if(SDS_CLDFRC_land(IDATA) .lt.0)SDS_CLDFRC_land(IDATA)=INT_Fill_value 
!        SDS_mass_conc_land(IDATA)= FLOAT_Fill_value
        SDS_QCONTROL_land(1,IDATA)=0
        SDS_QCONTROL_land(2,IDATA)=0
        SDS_QCONTROL_land(3,IDATA)=0
        SDS_QCONTROL_land(4,IDATA)=0
        SDS_QCONTROL_land(5,IDATA)=0 
         SDS_SCAT_ANGLE_land(IDATA)=INT_Fill_value 
          SDSTAU_corrected(IDATA,3)=INT_Fill_value
          SDSTAU_corrected_213(IDATA)=INT_Fill_value
          SDSTAU_small_land(IDATA,1)=INT_Fill_value
          SDSTAU_small_land(IDATA,2)=INT_Fill_value
          SDSTAU_small_land(IDATA,3)=INT_Fill_value
          SDSTAU_small_land(IDATA,4)=INT_Fill_value
          SDS_angs_coeff_land(IDATA)=INT_Fill_value 
          SDS_dust_weighting(IDATA)=INT_Fill_value
         Do Num_Wave=1,Band_land+3
         SDS_NUMPIXELS_land(IDATA,Num_Wave)=INT_Fill_value 
         SDS_ref_land(IDATA,Num_Wave) =INT_Fill_value
         SDS_ref_STD_land(IDATA,Num_Wave)=INT_Fill_value
         Enddo
         SDS_Surface_Reflectance_Land(IDATA,1) =INT_Fill_value
         SDS_Surface_Reflectance_Land(IDATA,2) =INT_Fill_value
         SDS_Surface_Reflectance_Land(IDATA,3) =INT_Fill_value 
         SDS_Fitting_Error_Land(IDATA)=INT_Fill_value
          ENDIF
         RETURN
         END

