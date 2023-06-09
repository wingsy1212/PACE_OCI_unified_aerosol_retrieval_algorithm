
!----------------------------------------------------------------------
!!F77-INC
!
!!DESCRIPTION:
!                This is an include file for the declaration of
!                variables used to write the Vdata and SDS arrays for
!                mod04 land algorithm.
!
!!Input Parameters:  N/A
!!Output Parameters: N/A
!
!!TEAM-UNIQUE HEADER:
!
!    This software is developed by the MODIS Science Data Support
!    Team for the National Aeronautics and Space Administration,
!    Goddard Space Flight Center, under contract NAS5-32373.
!
! 
!
!-----------------------------------------------------------------------

! Define parameters From Main

      INTEGER iwave_466, iwave_644, iwave_553, iwave_212
      PARAMETER(iwave_466=1, iwave_644=3, iwave_553=2, iwave_212=4)

	 
! Define parameters for new land (_NL) lookup tables
!

      INTEGER  NLWAV,NLTHET0,NLTHE,NLPHI,NLTAU,NLTABLE
      PARAMETER (NLWAV=4,NLTHET0=11,NLTHE=15,NLPHI=16,NLTAU=7,NLTABLE=5)
      INTEGER NLSIZE, NLETA
      PARAMETER (NLETA = 13, NLSIZE=2)
      INTEGER DTABLE_NL, MTABLE_NL(2)
      PARAMETER (DTABLE_NL = 5)
      Integer NUMCELLS_L,MaxPixels_left_L,Max_Pixels_L
      PARAMETER (NUMCELLS_L=1200,MaxPixels_left_L=1200,Max_Pixels_L=1200)

      REAL MINMTHET0_NL,MAXMTHET0_NL,MINMTHET_NL,&
         MAXMTHET_NL,MINMPHI_NL,MAXMPHI_NL
      PARAMETER(MINMTHET0_NL=0.0,MAXMTHET0_NL=84.0,MINMTHET_NL=0.0,&
         MAXMTHET_NL=72.0,MINMPHI_NL=0.0,MAXMPHI_NL=180.0)

      REAL REFW466L,REFW553L,REFW644L,REFW866L,REFW212L
      REAL REFW164L,REFW124L

!     Define parameters based on lookup tables
      REAL WAV_NL,PHI_NL,THE_NL,THET0_NL,MU0_NL

! 
	  
      COMMON/THREE/&
                    PHI_NL(NLPHI),&
                   THE_NL(NLTHE),&
                    THET0_NL(NLTHET0),&
                    MU0_NL(NLTHET0),&
                    WAV_NL(NLWAV),&
                    REFW466L,REFW644L,REFW212L,&
                    REFW553L,REFW866L,REFW164L,REFW124L




