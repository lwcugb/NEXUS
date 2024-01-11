!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_megan_mod.F90
!
! !DESCRIPTION: Module HCOX\_Megan\_Mod contains variables and routines
!  specifying the algorithms that control the MEGAN inventory of biogenic
!  emissions (as implemented into the GEOS-Chem model).
!\\
!\\
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!\\
!\\
! MEGAN calculates gamma activity factor based upon temperature and
! radiation information from the past. In the original GEOS-Chem
! code, the initial 10-d averages were explicitly calculated during
! initialization of MEGAN. This is not feasible in an ESMF environment,
! and the following restart variables can now be provided through the
! HEMCO configuration file:
! \begin{itemize}
! \item T\_DAVG: long-term historical temperature
! \item PARDR\_DAVG: long-term historical direct radiation
! \item PARDF\_DAVG: long-term historical diffuse radiation
! \item T\_PREVDAY: short-term historical temperature
! \end{itemize}
! These variables are automatically searched for on the first call of
! the run call. If not defined, default values will be used. The values
! of T\_DAVG, T\_PREVDAY, PARDR\_DAVG, and PARDF\_DAVG are continuously
! updated at the end of the run sequence, e.g. they represent the
! instantaneous running average. The e-folding times to be used when
! calculating the short=term and long-term running averages are defined
! as module parameter below (parameter TAU\_HOURS and TAU\_DAYS).
!\\
!\\
! A similar procedure is also applied to the leaf area index variables.
! The original GEOS-Chem MEGAN code used three LAI variables: current month
! LAI (LAI\_CM), previous month LAI (LAI\_PM), and instantaneous LAI (LAI),
! which was a daily interpolation of LAI\_CM and next month' LAI, (LAI\_NM).
! The HEMCO implementation uses only the instantaneous LAI, assuming it is
! updated every day. The short term historical LAI is kept in memory and
! used to determine the LAI change over time (used to calculate the gamma
! leaf age). It is also updated on every time step.
! For the first simulation day, the previous' day LAI is taken from the
! restart file (field LAI\_PREVDAY). If no restart variable is defined, a
! LAI change of zero is assumed (ckeller, 10/9/2014).
!\\
!\\
! !References:
!
!  \begin{itemize}
!  \item Guenther, A., et al., \emph{The Model of Emissions of Gases and
!        Aerosols from Nature version 2.1 (MEGAN2.1): an extended and updated
!        framework for modeling biogenic emissions}, \underline{Geosci. Model
!        Dev.}, \textbf{5}, 1471-1792, 2012.
!  \item Guenther, A., et al., \emph{A global model of natural volatile
!        organic compound emissions}, \underline{J.Geophys. Res.},
!        \textbf{100}, 8873-8892, 1995.
!  \item Wang, Y., D. J. Jacob, and J. A. Logan, \emph{Global simulation of
!        tropospheric O3-Nox-hydrocarbon chemistry: 1. Model formulation},
!        \underline{J. Geophys. Res.}, \textbf{103}, D9, 10713-10726, 1998.
!  \item Guenther, A., B. Baugh, G. Brasseur, J. Greenberg, P. Harley, L.
!        Klinger, D. Serca, and L. Vierling, \emph{Isoprene emission estimates
!        and uncertanties for the Central African EXPRESSO study domain},
!        \underline{J. Geophys. Res.}, \textbf{104}, 30,625-30,639, 1999.
!  \item Guenther, A. C., T. Pierce, B. Lamb, P. Harley, and R. Fall,
!        \emph{Natural emissions of non-methane volatile organic compounds,
!        carbon monoxide, and oxides of nitrogen from North America},
!        \underline{Atmos. Environ.}, \textbf{34}, 2205-2230, 2000.
!  \item Guenther, A., and C. Wiedinmyer, \emph{User's guide to Model of
!        Emissions of Gases and Aerosols from Nature}. http://cdp.ucar.edu.
!        (Nov. 3, 2004)
!  \item Guenther, A., \emph{AEF for methyl butenol}, personal commucation.
!        (Nov, 2004)
!  \item Sakulyanontvittaya, T., T. Duhl, C. Wiedinmyer, D. Helmig, S.
!        Matsunaga, M. Potosnak, J. Milford, and A. Guenther, \emph{Monoterpene
!        and sesquiterpene emission estimates for the United States},
!        \underline{Environ. Sci. Technol}, \textbf{42}, 1623-1629, 2008.
!  \end{itemize}
!
! !INTERFACE:
!
      MODULE HCOX_MEGAN_MOD
!
! !USES:
!
      !USE HCO_ERROR_MOD
      !USE HCO_DIAGN_MOD
      !USE HCOX_State_MOD,    ONLY : Ext_State
      !USE HCO_STATE_MOD,     ONLY : HCO_STATE

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
!      PUBLIC  :: HCOX_Megan_Init
!      PUBLIC  :: HCOX_Megan_Run
!      PUBLIC  :: HCOX_Megan_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PUBLIC :: GET_MEGAN_EMISSIONS  ! dbm, new MEGAN driver routine
                                      ! for all compounds (6/21/2012)
!      PRIVATE :: UPDATE_T_DAY
!      PRIVATE :: UPDATE_T_15_AVG
      PRIVATE :: GET_MEGAN_PARAMS
      !PRIVATE :: GET_MEGAN_AEF
      PRIVATE :: GET_GAMMA_PAR_PCEEA
      PRIVATE :: GET_GAMMA_T_LI
      PRIVATE :: GET_GAMMA_T_LD
      PRIVATE :: GET_GAMMA_LAI
      PRIVATE :: GET_GAMMA_AGE
      PRIVATE :: GET_GAMMA_SM
      PUBLIC :: CALC_NORM_FAC
      PRIVATE :: SOLAR_ANGLE
!      PRIVATE :: FILL_RESTART_VARS
!      PRIVATE :: CALC_AEF
      PRIVATE :: GET_GAMMA_CO2  ! (Tai, Jan 2013)
!
! !REVISION HISTORY:
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (6/2003).  Updated to
!        latest algorithm and modified for the standard code by May Fu
!        (11/2004).
!  (2 ) All emission are currently calculated using TS from DAO met field.
!        TS is the surface air temperature, which should be carefully
!        distinguished from TSKIN. (tmf, 11/20/2004)
!  (3 ) In GEOS4, the TS used here are the T2M in the A3 files, read in
!        'a3_read_mod.f'.
!  (4 ) Bug fix: change #if block to also cover GCAP met fields (bmy, 12/6/05)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Bug fix: Skip Feb 29th if GCAP in INIT_MEGAN (phs, 9/18/07)
!  (7 ) Added routine GET_AEF_05x0666 to read hi-res AEF data for the GEOS-5
!        0.5 x 0.666 nested grid simulations (yxw, dan, bmy, 11/6/08)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  09 Mar 2010 - R. Yantosca - Minor bug fix in GET_EMMONOT_MEGAN
!  17 Mar 2010 - H. Pye      - AEF_SPARE must be a scalar local variable
!                              in GET_EMMONOT_MEGAN for parallelization.
!  20 Aug 2010 - R. Yantosca - Move CMN_SIZE to top of module
!  20 Aug 2010 - R. Yantosca - Now set DAY_DIM = 24 for MERRA, since the
!                              surface temperature is now an hourly field.
!  01 Sep 2010 - R. Yantosca - Bug fix in INIT_MEGAN: now only read in
!                              NUM_DAYS (instead of 15) days of sfc temp data
!  22 Nov 2011 - R. Yantosca - Do not use erroneous AEF's for nested grids
!  06 Dec 2011 - E. Fischer  - Added Acetone emissions
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!  11 Apr 2012 - R. Yantosca - Replace lai_mod.F with modis_lai_mod.F90
!  13 Aug 2013 - M. Sulprizio- Modifications for updated SOA sim (H. Pye):
!                               Add sesquiterpenes to MEGAN group;
!                               Add plant functional types (PFT_xx);
!                               Rename GET_EMMONOG_MEGAN to GET_EMTERP_MEGAN;
!                               Add routines READ_PFT and GET_AEF_GEN
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  05 Oct 2013 - C. Keller   - Now a HEMCO extension
!  04 Aug 2014 - C. Keller   - Added 'manual' diagnostics for Acetone.
!  09 Oct 2014 - C. Keller   - Now use only GC_LAI (keep prev. LAI in memory)
!  22 Dec 2014 - C. Keller   - Now use flexible precision (hp) everywhere.
!                              Option to read temperature/irradiation from
!                              restart.
!  26 Jan 2015 - M. Sulprizio- Update from D. Millet (19 Jan 2013): Streamlined
!                              computations into a single driver routine
!                              and updated emissions according to MEGAN 2.1 as
!                              described in:
!                              Guenther et al., The Model of Emissions of Gases
!                              and Aerosols from Nature version 2.1 (MEGAN2.1):
!                              an extended and updated framework for modeling
!                              biogenic emissions, GMD, 5, 1471-1492, 2012.
!  12 Feb 2015 - M. Sulprizio- Remove GET_AEF_GEN routine. We now calculate AEFs
!                              for FARN, BCAR, and OSQT in CALC_AEF using
!                              parameters from Guenther et al., 2012.
!  18 Feb 2015 - M. Sulprizio- Remove LPECCA logical flag since we use this
!                              scheme exclusively now.
!                              Restore emissions of individual MEGAN species to
!                              diagnostics for consistency with pre-HEMCO code.
!  10 Jun 2015 - M. Sulprizio- Bug fix for SOA simulation: Now convert AEFs for
!                              sesquiterpenes to kg/m2/s.
!  15 Sep 2015 - M. Sulprizio- Add CO2 inhibition effect on isoprene emissions
!                              from Amos Tai (Jan 2013)
!  05 Nov 2015 - C. Keller   - Reorganize restart variables to running averages.
!  08 Dec 2015 - C. Keller   - Now treat previous' day LAI as running avg, too.
!  14 Oct 2016 - C. Keller   - Now use HCO_EvalFld instead of HCO_GetPtr.
!  05 Oct 2015 - M. Sulprizio- Activate MEGAN ethanol emissions for PAN updates
!                              from E. Fischer
!  17 Jul 2017 - C. Keller   - Now normalize LAI by PFTs.
!  05 Oct 2018 - R. Yantosca - For the standard (non-ESMF) environment, update
!                              variables from the restart file just once.  In
!                              the ESMF environment we have to do that on every
!                              call to get data from the External State object.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!     
      !------These are precision types in hco_error_mod.F90; 
      !I did not use the HCO_Error_Mod module because it may diff from model to model 
      ! Double and single precision definitions
      INTEGER, PARAMETER, PRIVATE  :: dp = KIND( REAL( 0.0, 8 ) ) ! Double (r8)
      INTEGER, PARAMETER, PRIVATE  :: sp = KIND( REAL( 0.0, 4 ) ) ! Single (r4)
!#if defined( USE_REAL8 )
      ! Use 8-byte floating point precision when asked.
      INTEGER, PARAMETER, PRIVATE :: fp = KIND( REAL( 0.0, 8 ) )
!#else
      ! Use 4-byte floating point by default.
      !INTEGER, PARAMETER, PRIVATE :: fp = KIND( REAL( 0.0, 4 ) )
!#endif
 
!#if defined( USE_REAL8 )
      INTEGER, PARAMETER, PRIVATE  :: hp = KIND( REAL( 0.0, 8 ) ) ! HEMCO prec r8
!#else
      !INTEGER, PARAMETER, PRIVATE  :: hp = KIND( REAL( 0.0, 4 ) ) ! HEMCO prec r4
!#endif
      TYPE :: MyMet
         !REAL(hp)  :: AEF
         ! Days between mid-months (updated by HEMCO clock)
         REAL(hp)  :: D_BTW_M
         REAL(hp)  :: TS
         REAL(hp)  :: SUNCOS
         !WINDSP seems not used anywhere, but inclided 
         REAL(hp)  :: WINDSP
         REAL(hp)  :: Q_DIR_2
         REAL(hp)  :: Q_DIFF_2
         REAL(hp)  :: ISOLAI, MISOLAI, PMISOLAI
         ! New restart variables (ckeller, 11/05/2015)
         REAL(sp)  :: T_LASTXDAYS            ! Avg. temperature of last NUM_DAYS
         REAL(sp)  :: T_LAST24H              ! Avg. temperature of last 24 hours
         REAL(sp)  :: PARDF_LASTXDAYS        ! Avg. PARDF of last NUM_DAYS
         REAL(sp)  :: PARDR_LASTXDAYS        ! Avg. PARDR of last NUM_DAYS
         LOGICAL   :: LISOPCO2               ! Include CO2 inhibition of ISOP?
         REAL(hp)  :: GLOBCO2                ! Global CO2 conc (ppmv)
         REAL(hp)  :: GWETROOT
         REAL(hp)  :: LAT, LocalHour
         INTEGER   :: DOY
         ! Physical parameter
         REAL(hp)  :: D2RAD  ! Degrees to radians
         REAL(hp)  :: RAD2D  ! Radians to degrees
      END TYPE MyMet

      ! Pointer to local MET used by MEGAN
      PUBLIC :: MyMet
      !TYPE(MyMet), POINTER      :: LocalMet => NULL()

      ! Scalars
!      INTEGER,  PARAMETER  :: DAY_DIM        = 24       ! # of 1-hr periods/day
!      INTEGER, PARAMETER  :: DAY_DIM        = 8        ! # of 3-hr periods/day
      !INTEGER,  PARAMETER  :: NUM_DAYS       = 10       ! # of days to avg
      REAL(hp), PARAMETER  :: TAU_DAYS       = 5.0_hp   ! e-folding time to be applied to
                                                        ! long-term past conditions (in days)
      REAL(hp), PARAMETER  :: TAU_HOURS      = 12.0_hp  ! e-folding time to be applied to
                                                        ! short-term past conditions (in hours)
      REAL(hp), PARAMETER  :: WM2_TO_UMOLM2S = 4.766_hp ! W/m2 -> umol/m2/s

      ! Maximum LAI value
      REAL(hp), PARAMETER  :: LAI_MAX        = 6.0_hp   ! cm2/cm2

      ! PI    : Double-Precision value of PI
      REAL(fp), PARAMETER :: PI       = 3.14159265358979323e+0_fp
      ! PI_180 : Number of radians per degree
      !REAL(fp), PARAMETER :: PI_180   = PI / 180e+0_fp


      !local RC values; always equal zero (Wei Li)
      INTEGER, PARAMETER, PUBLIC  :: RC_SUCCESS = 0

      ! testing only
      integer, parameter  :: ix = 20 !20 !25 !13 !20
      integer, parameter  :: iy = 43 !43 !22 !38 !31

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Megan_Emissions
!
! !DESCRIPTION: Subroutine Get\_Megan\_Emissions computes biogenic emissions in
!  units of [kgC/m2/s] or [kg/m2/s] using the MEGAN inventory. (dbm, 12/2012)
!\\
!\\
! !INTERFACE:
!
     ! SUBROUTINE GET_MEGAN_EMISSIONS( am_I_Root, HcoState, ExtState,
     !&                                Inst, I, J, CMPD, MEGAN_EMIS, RC )

      SUBROUTINE GET_MEGAN_EMISSIONS( LocalMet, I, J, CMPD, AEF, MEGAN_EMIS, RC )

! !INPUT PARAMETERS:
!
      !USE HCO_CLOCK_MOD, ONLY : HcoClock_Get, HcoClock_GetLocal

      !LOGICAL,          INTENT(IN)  :: am_I_Root
      !TYPE(HCO_STATE),  POINTER     :: HcoState
      !TYPE(Ext_State),  POINTER     :: ExtState
      TYPE(MyMet),     POINTER     :: LocalMet
      INTEGER,          INTENT(IN)  :: I, J      ! lon & lat indices (only need for logs at the end if necessary )
      CHARACTER(LEN=*), INTENT(IN)  :: CMPD      ! Compound name (dbm,6/21/2012)
      REAL(hp),         INTENT(IN)  :: AEF       ! directly add AEF here (not from GET_MEGAN_AEF) Wei Li
!
! !OUTPUT PARAMETERS:
!
      REAL(hp),         INTENT(OUT) :: MEGAN_EMIS ! VOC emission in kgC/m2/s
                                                  ! or kg/m2/s, depending on
                                                  ! units the compound is
                                                  ! carried in
!
! !INPUT/OUTPUT PARAMETERS:
!     Make it optional (Wei Li) 
      INTEGER,  INTENT(INOUT),  OPTIONAL :: RC
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2000, 2004, 2006
!  (2 ) Wang,    et al, 1998
!  (3 ) Guenther et al, 2007, MEGAN v2.1 User mannual
!  (4 ) Guenther et al, 2012 GMD MEGANv2.1 description and associated code at
!                                http://acd.ucar.edu/~guenther/MEGAN/
!
! !REVISION HISTORY:
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest
!        algorithm and modified for the standard code by May Fu (11/20/04)
!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO
!        met field. TS is the surface air temperature, which should be
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (3 ) Restructing of function & implementation of activity factors (mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  11 Apr 2012 - R. Yantosca - Now use data from modis_lai_mod.F90
!  11 Apr 2012 - R. Yantosca - Cosmetic changes
!  26 Jan 2015 - M. Sulprizio- Update from D. Millet (21 Jun 2012): New driver
!                              routine for all MEGAN compounds
!  15 Sep 2015 - M. Sulprizio- Add CO2 inhibition effect on isoprene emissions
!                              from Amos Tai (Jan 2013)
!  17 Jul 2017 - C. Keller   - Now normalize LAI by PFT's.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)  :: GAMMA_LAI
      REAL(hp)  :: GAMMA_AGE
      REAL(hp)  :: GAMMA_TP
      REAL(hp)  :: CDEA(5)
      REAL(hp)  :: VPGWT(5)
      REAL(hp)  :: GAMMA_PAR
      REAL(hp)  :: GAMMA_PAR_Sun, GAMMA_PAR_Shade
      REAL(hp)  :: GAMMA_T_LD
      REAL(hp)  :: GAMMA_T_LD_Sun, GAMMA_T_LD_Shade
      REAL(hp)  :: GAMMA_T_LI
      REAL(hp)  :: GAMMA_T_LI_Sun, GAMMA_T_LI_Shade
      REAL(hp)  :: GAMMA_SM
      REAL(hp)  :: GAMMA_CO2  ! (Tai, Jan 2013)
      !REAL(hp)  :: AEF
      !REAL(hp)  :: D_BTW_M
      !REAL(hp)  :: TS, SUNCOS, WINDSP
      !REAL(hp)  :: Q_DIR_2, Q_DIFF_2
      REAL(hp)  :: BETA, LDF, CT1, CEO
      REAL(hp)  :: ANEW, AGRO, AMAT, AOLD
      !REAL(hp)  :: MISOLAI,PMISOLAI
      !REAL(hp)  :: PFTSUM, LAT, LocalHour
      REAL(hp)  :: PSTD
      REAL(hp)  :: Ea1L, Ea2L, SINbeta, SunF
      LOGICAL   :: BIDIR
      INTEGER   :: K  !, DOY
      REAL(hp)  :: T_Leaf_Int_Sun(5)
      REAL(hp)  :: T_Leaf_Int_Shade(5)
      REAL(hp)  :: T_Leaf_Temp_Sun(5)
      REAL(hp)  :: T_Leaf_Temp_Shade(5)
      !REAL(hp)  :: T_Leaf_Wind_Sun(5)
      !REAL(hp)  :: T_Leaf_Wind_Shade(5)
      REAL(hp)  :: P_Leaf_Int_Sun(5)
      REAL(hp)  :: P_Leaf_Int_Shade(5)
      REAL(hp)  :: P_Leaf_LAI_Sun(5)
      REAL(hp)  :: P_Leaf_LAI_Shade(5)
      REAL(hp)  :: Distgauss(5)
      !=================================================================
      ! GET_MEGAN_EMISSIONS begins here!
      !=================================================================

      ! Initialize parameters, gamma values, and return value

      CDEA       = 0.0_hp
      MEGAN_EMIS = 0.0_hp
      GAMMA_LAI  = 0.0_hp
      GAMMA_TP  = 0.0_hp
      GAMMA_AGE  = 0.0_hp
      GAMMA_T_LD = 0.0_hp
      GAMMA_T_LI = 0.0_hp
      GAMMA_PAR  = 0.0_hp
      GAMMA_SM   = 0.0_hp
      GAMMA_CO2  = 0.0_hp  ! (Tai, Jan 2013)
      BETA       = 0.0_hp
      !AEF        = 0.0_hp
      LDF        = 0.0_hp
      CT1        = 0.0_hp
      CEO        = 0.0_hp
      ANEW       = 0.0_hp
      AGRO       = 0.0_hp
      AMAT       = 0.0_hp
      AOLD       = 0.0_hp
      BIDIR      = .FALSE.

      ! --------------------------------------------
      ! Get MEGAN parameters for this compound
      ! --------------------------------------------
      CALL GET_MEGAN_PARAMS ( CMPD, BETA, LDF,  CT1,  CEO,          &
                              ANEW, AGRO, AMAT, AOLD, BIDIR, RC )
      !IF ( RC /= HCO_SUCCESS ) RETURN
      IF ( RC /= RC_SUCCESS ) RETURN

      ! --------------------------------------------
      ! Get base emission factor for this compound and grid square
      ! Units: kgC/m2/s or kg/m2/s
      ! --------------------------------------------
      !CALL GET_MEGAN_AEF ( am_I_Root, HcoState, Inst,
      !&                     I, J, CMPD, AEF, RC )
      !IF ( RC /= HCO_SUCCESS ) RETURN

      !-----------------------------------------------------
      ! Only interested in terrestrial biosphere
      ! If (local LAI != 0 .AND. baseline emission !=0 )
      !-----------------------------------------------------
      IF ( LocalMet%ISOLAI * AEF > 0.0_hp ) THEN

         ! --------------------------------------------------
         ! GAMMA_par (light activity factor)
         ! --------------------------------------------------

         ! Calculate GAMMA PAR only during day
         IF ( LocalMet%SUNCOS > 0.0_hp ) THEN


            GAMMA_PAR = GET_GAMMA_PAR_PCEEA(LocalMet%Q_DIR_2,  LocalMet%Q_DIFF_2,               &
                                           LocalMet%PARDR_LASTXDAYS,LocalMet%PARDF_LASTXDAYS,   &
                                           LocalMet%LAT, LocalMet%DOY, LocalMet%LocalHour,      &
                                           LocalMet%D2RAD, LocalMet%RAD2D)

         ELSE

            ! If night
            GAMMA_PAR = 0.0_hp
         ENDIF

         ! --------------------------------------------------
         ! GAMMA_T_LI (temperature activity factor for
         ! light-independent fraction)
         ! --------------------------------------------------
!         GAMMA_T_LI = GET_GAMMA_T_LI( TS, BETA )

         ! --------------------------------------------------
         ! GAMMA_T_LD (temperature activity factor for
         ! light-dependent fraction)
         ! --------------------------------------------------
!         GAMMA_T_LD = GET_GAMMA_T_LD( TS, Inst%T_LASTXDAYS(I,J),
!     &        Inst%T_LAST24H(I,J), CT1, CEO )
!         GAMMA_T_LD = GET_GAMMA_T_LD( TS, T_15_AVG(I,J),
!     &        T_15(I,J,1), CT1, CEO )

         ! --------------------------------------------------
         ! GAMMA_LAI (leaf area index activity factor)
         ! --------------------------------------------------
         GAMMA_LAI = GET_GAMMA_LAI( LocalMet%MISOLAI, BIDIR )

         ! --------------------------------------------------
         ! GAMMA_AGE (leaf age activity factor)
         ! --------------------------------------------------
         GAMMA_AGE = GET_GAMMA_AGE( LocalMet%MISOLAI, LocalMet%PMISOLAI,              &
              LocalMet%D_BTW_M, LocalMet%T_LASTXDAYS, ANEW, AGRO, AMAT, AOLD )   
!     &        D_BTW_M, T_15_AVG(I,J), ANEW, AGRO, AMAT, AOLD )


         ! --------------------------------------------------
         ! GAMMA_SM (soil moisture activity factor)
         ! --------------------------------------------------
         GAMMA_SM = GET_GAMMA_SM( LocalMet%GWETROOT, CMPD )

         ! CO2 inhibition of isoprene (Tai, Jan 2013)
         IF ( LocalMet%LISOPCO2 ) THEN
            GAMMA_CO2 = GET_GAMMA_CO2( LocalMet%GLOBCO2 )
         ELSE
            GAMMA_CO2 = 1.0_hp
         ENDIF
            T_Leaf_Int_Sun  = (/-13.891_hp, -12.322_hp, -1.032_hp,     &
                               -5.172_hp, -5.589_hp/)
            T_Leaf_Int_Shade = (/-12.846_hp, -11.343_hp, -1.068_hp,    &
                                -5.551_hp, -5.955_hp/)

            T_Leaf_Temp_Sun = (/1.064_hp, 1.057_hp, 1.031_hp,          &
                               1.050_hp, 1.051_hp/)
            T_Leaf_Temp_Shade = (/1.060_hp, 1.053_hp, 1.031_hp,        &
                                 1.051_hp, 1.052_hp/)

            P_Leaf_Int_Sun  = (/1.0831_hp, 1.0964_hp, 1.1036_hp,       &
                               1.0985_hp, 1.0901_hp/)
            P_Leaf_Int_Shade = (/0.8706_hp, 0.8895_hp, 0.9160_hp,      &
                                0.9407_hp, 0.9564_hp/)

            P_Leaf_LAI_Sun = (/0.0018_hp, -0.1281_hp, -0.2977_hp,      &
                              -0.4448_hp, -0.5352_hp/)
            P_Leaf_LAI_Shade = (/0.0148_hp, -0.1414_hp, -0.3681_hp,    &
                               -0.5918_hp, -0.7425_hp/)
            

         VPGWT = (/0.1184635, 0.2393144, 0.284444444,                  &
                  0.2393144, 0.1184635/)

         Distgauss = (/0.0469101, 0.2307534, 0.5, 0.7692465,           &
                      0.9530899/)

         !LAT = HcoState%Grid%YMID%Val(I,J)

         ! Get day of year, local-time and latitude
         !CALL HcoClock_Get (am_I_Root, HcoState%Clock, cDOY = DOY,
     !&                      RC=RC )
     !    CALL HcoClock_GetLocal ( HcoState, I, J, cH = LocalHour,
     !&                            RC=RC )

         SINbeta =  SOLAR_ANGLE(LocalMet%DOY, LocalMet%LocalHour, LocalMet%LAT, LocalMet%D2RAD)

         CDEA = GET_CDEA( LocalMet%MISOLAI )
         GAMMA_TP = 0.0_hp

         DO K = 1, 5

           SunF = Calc_Sun_Frac(LocalMet%MISOLAI,SINbeta,Distgauss(K))

           PSTD = 200_hp
           GAMMA_PAR_Sun = GET_GAMMA_PAR_C(LocalMet%Q_DIR_2,  LocalMet%Q_DIFF_2,                     &
                                       LocalMet%PARDR_LASTXDAYS, LocalMet%PARDF_LASTXDAYS,           &
                                       P_Leaf_LAI_Sun(K), P_Leaf_Int_Sun(K), LocalMet%MISOLAI,       &
                                       PSTD)

           PSTD = 50_hp
           GAMMA_PAR_Shade = GET_GAMMA_PAR_C(LocalMet%Q_DIR_2,  LocalMet%Q_DIFF_2,                   &
                                       LocalMet%PARDR_LASTXDAYS, LocalMet%PARDF_LASTXDAYS,          &
                                       P_Leaf_LAI_Shade(K),P_Leaf_Int_Shade(K), LocalMet%MISOLAI,    &
                                       PSTD)

         GAMMA_T_LD_Sun = GET_GAMMA_T_LD_C( LocalMet%TS, LocalMet%T_LASTXDAYS, LocalMet%T_LAST24H,   &
                                           CT1, CEO,                                                 &
                                           T_Leaf_Int_Sun(K),                                        &
                                           T_Leaf_Temp_Sun(K) )                                      

         GAMMA_T_LD_Shade = GET_GAMMA_T_LD_C( LocalMet%TS, LocalMet%T_LASTXDAYS, LocalMet%T_LAST24H, &
                                           CT1, CEO,                                                 &
                                           T_Leaf_Int_Shade(K),                                      &
                                           T_Leaf_Temp_Shade(K) )                                    

         GAMMA_T_LI_Sun =  GET_GAMMA_T_LI( LocalMet%TS, BETA,                                        &
                                          T_Leaf_Int_Sun(K),                                         &
                                          T_Leaf_Temp_Sun(K) )

         GAMMA_T_LI_Shade =  GET_GAMMA_T_LI( LocalMet%TS, BETA,                                      &
                                          T_Leaf_Int_Shade(K),                                       &
                                          T_Leaf_Temp_Shade(K) )


           Ea1L  =  CDEA(K) * GAMMA_PAR_Sun * GAMMA_T_LD_Sun * SunF +                                &
                       GAMMA_PAR_Shade * GAMMA_T_LD_Shade * (1-SunF)

           Ea2L =  GAMMA_T_LI_Sun * SunF +                                                           &
                       GAMMA_T_LI_Shade * (1-SunF)

           GAMMA_TP  = GAMMA_TP +                                                                    &
                       (Ea1L*LDF + Ea2L*(1-LDF))* VPGWT(K)
!           IF ( i==ix .and. j==iy .and. K==1 .and. CMPD == 'ISOP') THEN
!         write(*,*) ' '
!         write(*,*) '--- GET_MEGAN_EMISSIONS --- '
!         write(*,*) 'HEMCO MEGAN @    ', i,j
!         write(*,*) 'Compound       : ', TRIM(CMPD)
!         write(*,*) 'TS        : ', TS
!         write(*,*) 'BETA : ', BETA
!         write(*,*) 'Inst%T_LAST24H(I,J)      : ', Inst%T_LAST24H(I,J)
!         write(*,*) 'CT1      : ', CT1
!         write(*,*) 'CEO     : ', CEO
!         write(*,*) 'Q_DIR_2     : ', Q_DIR_2
!         write(*,*) 'Q_DIFF_2      : ', Q_DIFF_2
!         write(*,*) 'PARDR_LAST       : ', Inst%PARDR_LASTXDAYS(I,J)
!         write(*,*) 'PARDF_L             : ', Inst%PARDF_LASTXDAYS(I,J)
!         write(*,*) 'MISOLAI        : ', MISOLAI
!         write(*,*) 'SunF        : ', SunF
!         write(*,*) 'SINbeta        : ', SINbeta
!         write(*,*) 'LDF        : ', LDF
!         write(*,*) 'CDEA       : ', CDEA(K)
!         write(*,*) 'VPGWT      : ', VPGWT(K)
!         write(*,*) 'Distgauss     : ', Distgauss(K)
!         write(*,*) 'GAMMA_TP        : ', GAMMA_TP
!         write(*,*) 'GAMMA_PAR_Sun        : ', GAMMA_PAR_Sun
!         write(*,*) 'GAMMA_PAR_Shade      : ', GAMMA_PAR_Shade
!         write(*,*) 'GAMMA_T_LD_Shade        : ', GAMMA_T_LD_Shade
!         write(*,*) 'GAMMA_T_LD_Sun        : ', GAMMA_T_LD_Sun
!         write(*,*) ' '
!            END IF

         ENDDO

      ELSE

         ! set activity factors to zero
         GAMMA_PAR  = 0.0_hp
         GAMMA_T_LI = 0.0_hp
         GAMMA_T_LD = 0.0_hp
         GAMMA_LAI  = 0.0_hp
         GAMMA_AGE  = 0.0_hp
         GAMMA_SM   = 0.0_hp
         GAMMA_CO2  = 0.0_hp  ! (Tai, Jan 2013)
         GAMMA_TP   = 0.0_hp
      END IF

      ! Emission is the product of all of these.
      ! Units here are kgC/m2/s or kg/m2/s as appropriate for the compound.
      ! Normalization factor ensures product of GAMMA values is 1.0 under
      !  standard conditions, Norm_FAC = 0.21.
      IF ( CMPD == 'ISOP' ) THEN
         ! Only apply CO2 inhibition to isoprene (mps, 9/15/15)
         ! Amos Tai wrote:
         !  In my opinion, we should not apply the CO2 inhibition factor on
         !  other monoterpene species yet, because the empirical data I've used
         !  are for isoprene emissions only. But we generally agree that CO2
         !  inhibition should affect monoterpenes too, so we'll leave room for
         !  future incorporation when new data arise.
!         MEGAN_EMIS = MISOLAI * AEF * GAMMA_AGE * GAMMA_SM *
!     &                GAMMA_LAI * Inst%NORM_FAC(1) * GAMMA_TP *GAMMA_CO2

         MEGAN_EMIS = LocalMet%MISOLAI * AEF * GAMMA_AGE * GAMMA_SM *     &
                      GAMMA_TP*GAMMA_CO2*GAMMA_LAI*0.21_hp
      ELSE
!         MEGAN_EMIS = MISOLAI * AEF * GAMMA_AGE * GAMMA_SM *
!     &                GAMMA_LAI * Inst%NORM_FAC(1) * GAMMA_TP
         MEGAN_EMIS = LocalMet%MISOLAI * AEF * GAMMA_AGE * GAMMA_SM *     &
                     GAMMA_TP * GAMMA_LAI * 0.21_hp
      ENDIF

!      ! testing only
!      if ( i==ix .and. j==iy ) then
!         write(*,*) ' '
!         write(*,*) '--- GET_MEGAN_EMISSIONS --- '
!         write(*,*) 'HEMCO MEGAN @    ', i,j
!         write(*,*) 'Compound       : ', TRIM(CMPD)
!         write(*,*) 'MEGAN_EMIS     : ', MEGAN_EMIS
!         write(*,*) 'SUNCOS         : ', SUNCOS
!         write(*,*) 'AEF [kgC/m2/s] : ', AEF
!         write(*,*) 'GAMMA_LAI      : ', GAMMA_LAI
!         write(*,*) 'GAMMA_AGE      : ', GAMMA_AGE
!         write(*,*) 'GAMMA_T_LI     : ', GAMMA_T_LI
!         write(*,*) 'GAMMA_T_LD     : ', GAMMA_T_LD
!         write(*,*) 'GAMMA_PAR      : ', GAMMA_PAR
!         write(*,*) 'GAMMA_SM       : ', GAMMA_SM
!         write(*,*) 'TS             : ', TS
!         write(*,*) 'T_15_AVG    : ', T_15_AVG(I,J)
!         write(*,*) 'HCOT_DAILY     : ', HCOT_DAILY(I,J)
!         write(*,*) 'PARDR_15_AVG: ', PARDR_15_AVG(I,J)
!         write(*,*) 'PARDF_15_AVG: ', PARDF_15_AVG(I,J)
!         write(*,*) 'ISOLAI         : ', ISOLAI
!         write(*,*) 'MISOLAI        : ', MISOLAI
!         write(*,*) 'PMISOLAI       : ', PMISOLAI
!         write(*,*) 'D_BTW_M        : ', D_BTW_M
!         write(*,*) ' '
!      endif

      ! Leave w/ success
      !RC = HCO_SUCCESS
      IF( PRESENT(RC) )  RC = RC_SUCCESS

      END SUBROUTINE GET_MEGAN_EMISSIONS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_NORM_FAC
!
! !DESCRIPTION: Function CALC\_NORM\_FAC calculates the normalization factor
!  needed to compute emissions. Called from GET\_MEGAN\_EMISSIONS.
!\\
!\\
! !INTERFACE:
!
      !SUBROUTINE CALC_NORM_FAC( am_I_Root, Inst, RC )
      SUBROUTINE CALC_NORM_FAC( D2RAD_FAC, NORM_FAC, RC )
!
! !INPUT PARAMETERS:
!
      !LOGICAL,         INTENT(IN)     :: am_I_Root
      !TYPE(MyInst),    POINTER        :: Inst
      REAL(hp),         INTENT(IN)     :: D2RAD_FAC
      REAL(hp),         INTENT(OUT)    :: NORM_FAC
!
! !INPUT/OUTPUT PARAMETERS
!
      INTEGER, INTENT(INOUT), optional  :: RC
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, (GMD 2012) and associated MEGANv2.1 source code
!
! !REVISION HISTORY:
!  (1 ) Created by dbm 11/2012. We calculate only 1 normalization factor for all
!       compounds based on the isoprene gamma values. Formally there should be a
!       different normalization factor for each compound, but we are following
!       Alex Guenther's approach here and the MEGAN source code.
!       "Hi Dylan, sorry for being so slow to get back to you.
!        Since the change is only a few percent or less, I didn't
!        bother to assign a different normalization factor to each
!        compound.  Since the MEGAN canopy environment model also
!        has 8 different canopy types (tropical broadleaf tree,
!        conifer tree, etc.) then to be correct we should have a
!        different CCE for each canopy type for each compound class
!        (which would be 160 slightly different values of CCE)."
!  07 Jan 2016 - E. Lundgren - Update ideal gas constant to NIST 2014 value
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp) :: PAC_DAILY, PHI, BBB, AAA, GAMMA_P_STANDARD
      REAL(hp) :: GAMMA_T_LI_STANDARD
      REAL(hp) :: GAMMA_SM_STANDARD
      REAL(hp) :: CMLAI, GAMMA_LAI_STANDARD
      REAL(hp) :: GAMMA_AGE_STANDARD
      REAL(hp) :: PT_15, T, R, CEO, CT1, E_OPT, T_OPT, CT2, X
      REAL(hp) :: GAMMA_T_LD_STANDARD
      REAL(hp) :: LDF, GAMMA_STANDARD
      REAL(hp) ::  SunF, GAMMA_TP_STANDARD
      REAL(hp)  :: T_Leaf_Int_Sun(5)
      REAL(hp)  :: T_Leaf_Int_Shade(5)
      REAL(hp)  :: T_Leaf_Temp_Sun(5)
      REAL(hp)  :: T_Leaf_Temp_Shade(5)
      REAL(hp)  :: T_Leaf_Wind_Sun(5)
      REAL(hp)  :: T_Leaf_Wind_Shade(5)
      REAL(hp)  :: P_Leaf_Int_Sun(5)
      REAL(hp)  :: P_Leaf_Int_Shade(5)
      REAL(hp)  :: P_Leaf_LAI_Sun(5)
      REAL(hp)  :: P_Leaf_LAI_Shade(5)
      REAL(hp)  :: Distgauss(5), CDEA(5), VPGWT(5)
      REAL(hp)  :: EA2L, EA1L, GAMMA_T_LD_SUN, GAMMA_T_LD_SHADE
      REAL(hp)  :: L_PT_T, L_T, C1, LAI, GAMMA_PAR_SUN
      REAL(hp)  :: GAMMA_PAR_SHADE, ALPHA, PAC_INSTANT
      INTEGER  :: Q

      !-----------------------------------------------------------------
      ! CALC_NORM_FAC
      !-----------------------------------------------------------------

      ! -----------------
      ! GAMMA_P for standard conditions
      ! -----------------
      ! Based on Eq. 11b from Guenther et al., 2006
      ! Using standard conditions of phi = 0.6, solar angle of 60 deg,
      ! and P_daily = 400
      ! Note corrigendum for Eq. 11b in that paper, should be a
      ! minus sign before the 0.9.
!      PAC_DAILY = 400.0_hp
!      PHI       = 0.6_hp
!      BBB       = 1.0_hp + 0.0005_hp *( PAC_DAILY - 400.0_hp )
!      AAA       = ( 2.46_hp * BBB * PHI ) - ( 0.9_hp * PHI**2 )
      ! sin(60) = 0.866
!      GAMMA_P_STANDARD = 0.866_hp * AAA

      ! -----------------
      ! GAMMA_T_LI for standard conditions
      ! -----------------
      ! gamma_t_li = EXP( Beta * ( T - T_Standard ) )
      ! This is 1.0 for T = T_Standard
      GAMMA_T_LI_STANDARD = 1.0_hp

      ! -----------------
      ! GAMMA_SM for standard conditions
      ! -----------------
      ! Standard condition is soil moisture = 0.3 m^3/m^3
      ! GAMMA_SM = 1.0 for all compounds under this condition
      GAMMA_SM_STANDARD = 1.0_hp

      ! -----------------
      ! GAMMA_TP for standard conditions
      ! -----------------
      ! gamma_t_li = EXP( Beta * ( T - T_Standard ) )
      ! This is 1.0 for T = T_Standard

      CDEA = GET_CDEA( 5.0_hp )
      Distgauss = (/0.0469101, 0.2307534, 0.5, 0.7692465,       &
                       0.9530899/)
      VPGWT = (/0.1184635, 0.2393144, 0.284444444,              &
                   0.2393144, 0.1184635/)
      P_Leaf_Int_Sun  = (/1.0831_hp, 1.0964_hp, 1.1036_hp,      &
                                1.0985_hp, 1.0901_hp/)
      P_Leaf_Int_Shade = (/0.8706_hp, 0.8895_hp, 0.9160_hp,     &
                                 0.9407_hp, 0.9564_hp/)

      P_Leaf_LAI_Sun = (/0.0018_hp, -0.1281_hp, -0.2977_hp,     &
                               -0.4448_hp, -0.5352_hp/)
      P_Leaf_LAI_Shade = (/0.0148_hp, -0.1414_hp, -0.3681_hp,   &
                                -0.5918_hp, -0.7425_hp/)

      T_Leaf_Int_Sun  = (/-13.891_hp, -12.322_hp, -1.032_hp,    &
                                -5.172_hp, -5.589_hp/)
      T_Leaf_Int_Shade = (/-12.846_hp, -11.343_hp, -1.068_hp,   &
                                 -5.551_hp, -5.955_hp/)

      T_Leaf_Temp_Sun = (/1.064_hp, 1.057_hp, 1.031_hp,         &
                                1.050_hp, 1.051_hp/)
      T_Leaf_Temp_Shade = (/1.060_hp, 1.053_hp, 1.031_hp,       &
                                  1.051_hp, 1.052_hp/)

      GAMMA_TP_STANDARD = 0.0_hp
      LDF = 1.0_hp
      LAI = 5.0_hp
      DO Q = 1, 5

            PAC_INSTANT  = 1500.0_hp/4.766_hp 
            PAC_DAILY = 740.0_hp/4.766_hp

            PAC_INSTANT = PAC_INSTANT * exp(P_Leaf_Int_Sun(Q) +      &
                               P_Leaf_LAI_Sun(Q) * LAI)
            PAC_DAILY = PAC_DAILY * exp(P_Leaf_Int_Sun(Q) +          &
                               P_Leaf_LAI_Sun(Q) * LAI)
            Alpha  = 0.004 - 0.0005*LOG(PAC_DAILY)
            C1 = 0.0468 * EXP(0.0005 * (PAC_DAILY - 200.0_hp)) *     &
                       (PAC_DAILY **  0.6)
            GAMMA_PAR_Sun = (Alpha * C1 * PAC_INSTANT) /             &
                       ((1 + Alpha**2. * PAC_INSTANT**2.)**0.5)

            PAC_INSTANT  = 1500.0_hp/4.766_hp
            PAC_DAILY = 740.0_hp/4.766_hp
            PAC_DAILY = PAC_DAILY * exp(P_Leaf_Int_Shade(Q) +        &
                               P_Leaf_LAI_Shade(Q) * LAI) 
            PAC_INSTANT = PAC_INSTANT * exp(P_Leaf_Int_Shade(Q) +    &
                               P_Leaf_LAI_Shade(Q) * LAI)
            Alpha  = 0.004 - 0.0005*LOG(PAC_DAILY)
            C1 = 0.0468 * EXP(0.0005 * (PAC_DAILY - 50.0_hp)) *      &
                       (PAC_DAILY **  0.6)
            GAMMA_PAR_Shade = (Alpha * C1 * PAC_INSTANT) /           &
                       ((1 + Alpha**2. * PAC_INSTANT**2.)**0.5)

            PT_15 = 298.5_hp
            T     = 303.0_hp
            R     = 8.3144598e-3_hp
            CEO = 2.0_hp
            CT1 = 95.0_hp
            CT2   = 230.0_hp

            L_T = T * T_Leaf_Temp_Sun(Q) + T_Leaf_Int_Sun(Q)
            L_PT_T = PT_15 * T_Leaf_Temp_Sun(Q) + T_Leaf_Int_Sun(Q)
            E_OPT = CEO * EXP( 0.1_hp * ( L_PT_T  - 2.97e2_hp ) )
            T_OPT = 3.125e2_hp + ( 6.0e-1_hp * ( L_PT_T - 2.97e2_hp ) )
            X     = ( 1.0_hp/T_OPT - 1.0_hp/L_T ) / R
            GAMMA_T_LD_Sun   = E_OPT * CT2 * EXP( CT1 * X ) /               &
              ( CT2 - CT1 * ( 1.0_hp - EXP( CT2 * X ) ) )

            L_T = T * T_Leaf_Temp_Shade(Q) + T_Leaf_Int_Shade(Q)
            L_PT_T = PT_15 * T_Leaf_Temp_Shade(Q) + T_Leaf_Int_Shade(Q)
            E_OPT = CEO * EXP( 0.08_hp * ( L_PT_T  - 2.97e2_hp ) )
            T_OPT = 3.125e2_hp + ( 6.0e-1_hp * ( L_PT_T - 2.97e2_hp ) )
            X     = ( 1.0_hp/T_OPT - 1.0_hp/L_T ) / R
            GAMMA_T_LD_Shade   = E_OPT * CT2 * EXP( CT1 * X ) /             &
              ( CT2 - CT1 * ( 1.0_hp - EXP( CT2 * X ) ) )

            SunF = Calc_Sun_Frac(5.0_hp,SIN(60.0_hp*D2RAD_FAC),            &
                             Distgauss(Q))
            Ea1L  =  CDEA(Q) * GAMMA_PAR_Sun * GAMMA_T_LD_Sun * SunF +      &
                        GAMMA_PAR_Shade * GAMMA_T_LD_Shade * (1-SunF)

            write(*,*) ' '
            write(*,*) '--- GET_MEGAN_NormFrac --- '
            write(*,*) 'GAMMA_TP_STANDARD      : ', GAMMA_TP_STANDARD
            write(*,*) 'Ea1L      : ', Ea1L
            write(*,*) 'SunF      : ', SunF
            write(*,*) 'CDEA      : ', CDEA(Q)
            write(*,*) 'VPGWT      : ', VPGWT(Q)
            write(*,*) 'Q      : ', Q
            write(*,*) 'Q      : ', Q
            write(*,*) 'GAMMA_PAR_Sun      : ', GAMMA_PAR_Sun
            write(*,*) 'GAMMA_PAR_Shade      : ', GAMMA_PAR_Shade
            write(*,*) 'GAMMA_T_LD_Sun      : ', GAMMA_T_LD_Sun
            write(*,*) 'GAMMA_T_LD_Shade      : ', GAMMA_T_LD_Shade

            GAMMA_TP_STANDARD  = GAMMA_TP_STANDARD +                       &
                        (Ea1L*LDF)* VPGWT(Q)
      ENDDO

            write(*,*) ' '
            write(*,*) '--- GET_MEGAN_NormFrac --- '
            write(*,*) 'GAMMA_TP_STANDARD      : ', GAMMA_TP_STANDARD
            write(*,*) 'Ea1L      : ', Ea1L
            write(*,*) 'SunF      : ', SunF
            write(*,*) 'CDEA      : ', CDEA(Q)
            write(*,*) 'VPGWT      : ', VPGWT(Q)
            write(*,*) 'Q      : ', Q

      ! -----------------
      ! GAMMA_LAI for standard conditions
      ! -----------------
      ! Standard condition is LAI = 5
      CMLAI = 5.0_hp
!      GAMMA_LAI_STANDARD = 0.49_hp *
!     &     CMLAI / SQRT( 1.0_hp + 0.2_hp * CMLAI*CMLAI )
      GAMMA_LAI_STANDARD = CMLAI
      ! -----------------
      ! GAMMA_AGE for standard conditions
      ! -----------------
      ! Standard condition is 0% new, 10% growing, 80% mature, 10% old foliage
      ! Isoprene uses A_NEW = 0.05d0, A_GRO = 0.6d0, A_MAT = 1.d0, A_OLD = 0.9d0
      GAMMA_AGE_STANDARD = 0.1_hp*0.6_hp + 0.8_hp*1.0_hp + 0.1_hp*0.9_hp

      ! -----------------
      ! GAMMA_T_LD for standard conditions
      ! -----------------
      ! Standard condition is
      ! PT_15 = average leaf temp over past 24-240 hours = 297K
      ! T = air temperature = 303K
!      PT_15 = 297.0_hp
!      T     = 303.0_hp
!      R     = 8.3144598e-3_hp
      ! parameters for isoprene
!      CEO = 2.0_hp
!      CT1 = 95.0_hp

!      E_OPT = CEO * EXP( 0.08_hp * ( PT_15  - 2.97e2_hp ) )
!      T_OPT = 3.13e2_hp + ( 6.0e-1_hp * ( PT_15 - 2.97e2_hp ) )
!      CT2   = 200.0_hp

      ! Variable related to temperature
!      X     = ( 1.0_hp/T_OPT - 1.0_hp/T ) / R

!      GAMMA_T_LD_STANDARD = E_OPT * CT2 * EXP( CT1 * X ) /
!     &     ( CT2 - CT1 * ( 1.0_hp - EXP( CT2 * X ) ) )

      ! -----------------
      ! Overall GAMMA_STANDARD
      ! -----------------
      ! LDF = 1d0 for isoprene
      GAMMA_STANDARD = GAMMA_AGE_STANDARD * GAMMA_SM_STANDARD *          &
           GAMMA_TP_STANDARD * GAMMA_LAI_STANDARD

      NORM_FAC = 1.0_hp / GAMMA_STANDARD

      write(*,*) ' '
      write(*,*) '--- GET_MEGAN_NormFrac --- '
      write(*,*) 'GAMMA_STANDARD      : ', GAMMA_STANDARD
      write(*,*) 'GAMMA_AGE_STANDARD      : ', GAMMA_AGE_STANDARD
      write(*,*) 'GAMMA_SM_STANDARD      : ', GAMMA_SM_STANDARD
      write(*,*) 'GAMMA_TP_STANDARD      : ', GAMMA_TP_STANDARD
      write(*,*) 'GAMMA_LAI_STANDARD      : ', GAMMA_LAI_STANDARD

      ! Return w/ success
      !RC = HCO_SUCCESS
      IF( PRESENT(RC) )  RC = RC_SUCCESS

      END SUBROUTINE CALC_NORM_FAC
!EOC


!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Sun_Frac
!
! !DESCRIPTION: Function Calc_Sun_Frac
!\\
!\\
! !INTERFACE:
!
      FUNCTION Calc_Sun_Frac( LAI, Sinbeta, Distgauss) RESULT( SunFrac )
!
! !INPUT PARAMETERS:
!
      ! Current leaf temperature, the surface air temperature field (TS)
      ! is assumed equivalent to the leaf temperature over forests.
      ! Temperature factor per species
      REAL(hp),  INTENT(IN) :: Distgauss
      REAL(hp),  INTENT(IN) :: SINbeta
      REAL(hp),  INTENT(IN) :: LAI
!
! !RETURN VALUE:
!
      ! Activity factor for the light-independent fraction of emissions
      REAL(hp)              :: SunFrac
!
! !DEFINED PARAMETERS:
!
      ! Standard reference temperature [K]
      REAL*8, PARAMETER   :: Cluster = 0.9
      REAL*8, PARAMETER   :: CANTRAN = 0.2
      REAL(hp)           :: Kb, LAIadj, LAIdepth

      Kb = Cluster * 0.5 / Sinbeta
      LAIadj = LAI / ( 1 - CANTRAN )
      LAIdepth   = LAIadj  * Distgauss

      IF ((Sinbeta  > 0.002) .AND. (LAIadj  > 0.001)) THEN
        SunFrac = EXP(-Kb * LAIdepth)
      ELSE
        SunFrac = 0.2
      ENDIF


      END FUNCTION Calc_Sun_Frac

!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Megan_Params
!
! !DESCRIPTION: Subroutine Get\_Megan\_Params returns the emission parameters
!  for each MEGAN compound needed to compute emissions. Called from
!  GET\_MEGAN\_EMISSIONS.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MEGAN_PARAMS( CPD,   BTA,   LIDF,  C_T1,  C_EO,   & 
                                  A_NEW, A_GRO, A_MAT, A_OLD, BI_DIR,   & 
                                  RC )
!
! !INPUT PARAMETERS:
!
      !LOGICAL,          INTENT(IN) :: am_I_Root ! Root CPU?
      !TYPE(HCO_State),  POINTER    :: HcoState
      CHARACTER(LEN=*), INTENT(IN) :: CPD       ! Compound name
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL(hp), INTENT(INOUT) :: BTA    ! Beta coefficient for temperature activity
                                        ! factor for light-independent fraction
      REAL(hp), INTENT(INOUT) :: LIDF   ! Light-dependent fraction of emissions
      REAL(hp), INTENT(INOUT) :: C_T1   ! CT1 parameter for temperature activity
                                        ! factor for light-dependent fraction
      REAL(hp), INTENT(INOUT) :: C_EO   ! Ceo parameter for temperature activity
                                        ! factor for light-dependent fraction
      REAL(hp), INTENT(INOUT) :: A_NEW  ! Relative emission factor (new leaves)
      REAL(hp), INTENT(INOUT) :: A_GRO  ! Relative emission factor (growing leaves)
      REAL(hp), INTENT(INOUT) :: A_MAT  ! Relative emission factor (mature leaves)
      REAL(hp), INTENT(INOUT) :: A_OLD  ! Relative emission factor (old leaves)
      LOGICAL,  INTENT(INOUT) :: BI_DIR ! Logical flag to indicate bidirectional exchange
      INTEGER,  INTENT(INOUT), OPTIONAL :: RC  !make it optional (Wei Li)
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, (GMD 2012) and associated MEGANv2.1 source code
!
! !REVISION HISTORY:
!  (1 ) Created by dbm 07/2012
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255):: MSG

      !=================================================================
      ! GET_MEGAN_PARAMS begins here!
      !=================================================================

      ! Initialize values
      BTA    = 0.0_hp
      LIDF   = 0.0_hp
      C_T1   = 0.0_hp
      C_EO   = 0.0_hp
      A_NEW  = 0.0_hp
      A_GRO  = 0.0_hp
      A_MAT  = 0.0_hp
      A_OLD  = 0.0_hp
      BI_DIR = .FALSE.

      ! ----------------------------------------------------------------
      ! Note that not all the above compounds are used in standard chemistry
      ! simulations, but they are provided here for future incorporation or
      ! specialized applications. More compounds can be added as needed
      ! by adding the corresponding CPD name and the appropriate paramaters.
      ! (dbm, 01/2013)
      !
      ! Values are from Table 4 in Guenther et al., 2012
      ! ----------------------------------------------------------------

      ! Isoprene, MBO
      IF ( TRIM(CPD) == 'ISOP' .OR.             &
           TRIM(CPD) == 'MBOX' ) THEN
         BTA    = 0.13_hp  ! Not actually used for ISOP, MBO
         LIDF   = 1.0_hp
         C_T1   = 95.0_hp
         C_EO   = 2.0_hp
         A_NEW  = 0.05_hp
         A_GRO  = 0.6_hp
         A_MAT  = 1.0_hp
         A_OLD  = 0.9_hp
         BI_DIR = .FALSE.

      ! Myrcene, sabinene, alpha-pinene
      ELSE IF ( TRIM(CPD) == 'MYRC' .OR.      &
                TRIM(CPD) == 'SABI' .OR.      &
                TRIM(CPD) == 'APIN' ) THEN
         BTA    = 0.10_hp
         LIDF   = 0.6_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 2.0_hp
         A_GRO  = 1.8_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.05_hp
         BI_DIR = .FALSE.

      ! Limonene, 3-carene, beta-pinene
      ELSE IF ( TRIM(CPD) == 'LIMO' .OR.     &
                TRIM(CPD) == 'CARE' .OR.     &
                TRIM(CPD) == 'BPIN' ) THEN
         BTA    = 0.10_hp
         LIDF   = 0.2_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 2.0_hp
         A_GRO  = 1.8_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.05_hp
         BI_DIR = .FALSE.

      ! t-beta-ocimene
      ELSE IF ( TRIM(CPD) == 'OCIM' ) THEN
         BTA    = 0.10_hp
         LIDF   = 0.8_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 2.0_hp
         A_GRO  = 1.8_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.05_hp
         BI_DIR = .FALSE.

      ! Other monoterpenes (lumped)
      ELSE IF ( TRIM(CPD) == 'OMON' ) THEN
         BTA    = 0.10_hp
         LIDF   = 0.4_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 2.0_hp
         A_GRO  = 1.8_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.05_hp
         BI_DIR = .FALSE.

      ! Methanol
      ELSE IF ( TRIM(CPD) == 'MOHX' ) THEN
         BTA    = 0.08_hp
         LIDF   = 0.8_hp
         C_T1   = 60.0_hp
         C_EO   = 1.6_hp
         A_NEW  = 3.5_hp
         A_GRO  = 3.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.2_hp
         BI_DIR = .FALSE.

      ! Acetone
      ELSE IF ( TRIM(CPD) == 'ACET' ) THEN
         BTA    = 0.1_hp
         LIDF   = 0.2_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 1.0_hp
         A_GRO  = 1.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.0_hp
         BI_DIR = .FALSE.

      ! Bidirectional VOC: Ethanol, formaldehyde, acetaldehyde, formic acid,
      ! acetic acid
      ELSE IF ( TRIM(CPD) == 'EOH'  .OR.     &
                TRIM(CPD) == 'CH2O' .OR.     &
                TRIM(CPD) == 'ALD2' .OR.     &
                TRIM(CPD) == 'FAXX' .OR.     &
                TRIM(CPD) == 'AAXX' ) THEN
         BTA    = 0.13_hp
         LIDF   = 0.8_hp
         C_T1   = 95.0_hp
         C_EO   = 2.0_hp
         A_NEW  = 1.0_hp
         A_GRO  = 1.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.0_hp
         BI_DIR = .TRUE.

      ! Stress VOCs: ethene, toluene, HCN
      ! There are others species in this category but none are currently
      ! used in GEOS-Chem
      ELSE IF ( TRIM(CPD) == 'C2H4' .OR.    &
                TRIM(CPD) == 'TOLU' .OR.    &
                TRIM(CPD) == 'HCNX' ) THEN
         BTA    = 0.1_hp
         LIDF   = 0.8_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 1.0_hp
         A_GRO  = 1.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.0_hp
         BI_DIR = .FALSE.

      ! Other VOCs: >C2 alkenes
      ! This includes propene, butene and very minor contribution from
      ! larger alkenes
      ELSE IF ( TRIM(CPD) == 'PRPE' ) THEN
         BTA    = 0.1_hp
         LIDF   = 0.2_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 1.0_hp
         A_GRO  = 1.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.0_hp
         BI_DIR = .FALSE.

      ! SOAupdate: Sesquiterpenes hotp 3/2/10
      ! alpha-Farnesene, beta-Caryophyllene, other sesquiterpenes
      ELSE IF ( TRIM(CPD) == 'FARN' .OR.     &
                TRIM(CPD) == 'BCAR' .OR.     &
                TRIM(CPD) == 'OSQT' ) THEN   
         BTA    = 0.17_hp
         LIDF   = 0.5_hp
         C_T1   = 130.0_hp
         C_EO   = 2.37_hp
         A_NEW  = 0.4_hp
         A_GRO  = 0.6_hp
         A_MAT  = 1.0_hp
         A_OLD  = 0.95_hp
         BI_DIR = .FALSE.

      ! Calls for any other MEGAN compounds (e.g. sesquiterpenes, CO, etc.) can
      ! be added following the above format based on the parameters in
      ! Guenther 2012 or the MEGAN source code (dbm, 6/21/2012).
      ELSE

         MSG = 'Invalid compound name'
         !CALL HCO_ERROR(HcoState%Config%Err,MSG, RC,
      !&                  THISLOC='GET_MEGAN_PARAMS' )
         write(*,*) MSG

         RETURN

      ENDIF

      ! Leave w/ success
      !RC = HCO_SUCCESS
      IF( PRESENT(RC) )  RC = RC_SUCCESS

      END SUBROUTINE GET_MEGAN_PARAMS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Megan_AEF (not used for now Wei Li)
!
! !DESCRIPTION: Function Get\_Megan\_AEF returns the appropriate AEF value
!  for a given compound and grid square.
!\\
!\\
! !INTERFACE:
!

!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_PAR_PCEEA
!
! !DESCRIPTION: Computes the PCEEA gamma activity factor with sensitivity
!  to LIGHT.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_PAR_PCEEA( Q_DIR_2, Q_DIFF_2, PARDR_AVG_SIM, PARDF_AVG_SIM,   &
                                    LAT, DOY, LocalHour, D2RAD, RAD2D )                &
      RESULT( GAMMA_P_PCEEA )
!
! !USES:
!
      !USE HCO_CLOCK_MOD, ONLY : HcoClock_Get, HcoClock_GetLocal
!
! !INPUT PARAMETERS:
!
      !LOGICAL,         INTENT(IN) :: am_I_Root
      !TYPE(HCO_State), POINTER    :: HcoState
      !TYPE(Ext_State), POINTER    :: ExtState
      !TYPE(MyInst),    POINTER    :: Inst
      !INTEGER,         INTENT(IN) :: I,  J             ! Lon & lat indices
      REAL(hp),        INTENT(IN) :: Q_DIR_2           ! Direct PAR [umol/m2/s]
      REAL(hp),        INTENT(IN) :: Q_DIFF_2          ! Diffuse PAR [umol/m2/s]
      REAL(sp),        INTENT(IN) :: PARDR_AVG_SIM     ! Avg direct PAR [W/m2]
      REAL(sp),        INTENT(IN) :: PARDF_AVG_SIM     ! Avg diffuse PAR [W/m2]
      REAL(hp),        INTENT(IN) :: LAT
      INTEGER ,        INTENT(IN) :: DOY
      REAL(hp),        INTENT(IN) :: LocalHour
      REAL(hp),        INTENT(IN) :: D2RAD, RAD2D
!
! !RETURN VALUE:
!
      REAL(hp)                    :: GAMMA_P_PCEEA     ! GAMMA factor for light
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 user guide
!
! !REVISION HISTORY:
!  (1 ) Here PAR*_AVG_SIM is the average light conditions over the simulation
!       period. I've set this = 10 days to be consistent with temperature & as
!       outlined in Guenther et al, 2006. (mpb,2009)
!  (2 ) Code was taken & adapted directly from the MEGAN v2.1 source code.
!       (mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)   :: mmPARDR_DAILY
      REAL(hp)   :: mmPARDF_DAILY
      REAL(hp)   :: PAC_DAILY, PAC_INSTANT, C_PPFD
      REAL(hp)   :: PTOA, PHI
      REAL(hp)   :: BETA,   SINbeta
      !INTEGER    :: RC
      REAL(hp)   :: AAA, BBB
      !REAL(hp)   :: LocalHour
      !REAL(hp)   :: LAT

      !-----------------------------------------------------------------
      ! Compute GAMMA_PAR_PCEEA
      !-----------------------------------------------------------------

      ! Initialize
      C_PPFD   = 0.0_hp
      PTOA     = 0.0_hp

      ! Convert past light conditions to micromol/m2/s
      mmPARDR_DAILY   = PARDR_AVG_SIM 
      mmPARDF_DAILY   = PARDF_AVG_SIM 

      ! Work out the light at the top of the canopy.
      PAC_DAILY    = mmPARDR_DAILY + mmPARDF_DAILY
      PAC_INSTANT  = Q_DIR_2       +  Q_DIFF_2

      ! Get latitude
      !LAT = HcoState%Grid%YMID%Val(I,J)

      ! Get day of year, local-time and latitude
      ! TODO: Evaluate RC?
      !CALL HcoClock_Get ( am_I_Root, HcoState%Clock, cDOY = DOY, RC=RC )
      !CALL HcoClock_GetLocal ( HcoState, I, J, cH = LocalHour, RC=RC )

      ! Get solar elevation angle
      SINbeta      =  SOLAR_ANGLE( DOY, LocalHour, LAT, D2RAD )
      BETA         =  ASIN( SINbeta ) * RAD2D

      IF ( SINbeta < 0.0_hp ) THEN

         GAMMA_P_PCEEA = 0.0_hp

      ELSEIF ( SINbeta > 0.0_hp ) THEN

         ! PPFD at top of atmosphere
         PTOA    = 3000.0_hp + 99.0_hp *                      &
                  COS( 2._hp * 3.14159265358979323_hp *       &
                  ( DOY - 10.0_hp ) / 365.0_hp )

         ! Above canopy transmission
         PHI     = PAC_INSTANT / ( SINbeta * PTOA )

         ! Work out gamma P
         BBB     = 1.0_hp + 0.0005_hp *( PAC_DAILY - 400.0_hp )
         AAA     = ( 2.46_hp * BBB * PHI ) - ( 0.9_hp * PHI**2 )

         GAMMA_P_PCEEA = SINbeta * AAA

      ENDIF

       ! Screen unforced errors. IF solar elevation angle is
       ! less than 1 THEN gamma_p can not be greater than 0.1.
       IF ( BETA < 1.0_hp .AND. GAMMA_P_PCEEA > 0.1_hp ) THEN
          GAMMA_P_PCEEA  = 0.0_hp
       ENDIF

      ! Prevent negative values
      GAMMA_P_PCEEA = MAX( GAMMA_P_PCEEA , 0.0_hp )

!      ! testing only
!      if ( i==ix .and. j==iy ) then
!         write(*,*) ' '
!         write(*,*) 'HEMCO GAMMA_PAR_PCEEA: ', GAMMA_P_PCEEA
!         write(*,*) 'DOY                  : ', DOY
!         write(*,*) 'LAT                  : ', LAT
!         write(*,*) 'SINbeta              : ', SINbeta
!         write(*,*) 'BETA                 : ', BETA
!         write(*,*) 'PTOA                 : ', PTOA
!         write(*,*) 'PHI                  : ', PHI
!         write(*,*) 'BBB                  : ', BBB
!         write(*,*) 'AAA                  : ', AAA
!         write(*,*) 'PAC_DAILY            : ', PAC_DAILY
!         write(*,*) 'PAC_INSTANT          : ', PAC_INSTANT
!      endif

      END FUNCTION GET_GAMMA_PAR_PCEEA
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Solar_Angle
!
! !DESCRIPTION: Function SOLAR\_ANGLE computes the local solar angle for a
!  given day of year, latitude and longitude (or local time).  Called from
!  routine Get\_Gamma\_P\_Pecca.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SOLAR_ANGLE( DOY, SHOUR, LAT, D2RAD ) RESULT(SINbeta)
!
! !INPUT PARAMETERS:
!
      ! Arguments
      !TYPE(HCO_State),   POINTER    :: HcoState
      !TYPE(MyInst),      POINTER    :: Inst
      INTEGER,           INTENT(IN) :: DOY       ! Day of year
      REAL(hp),          INTENT(IN) :: SHOUR     ! Local time
      REAL(hp),          INTENT(IN) :: LAT       ! Latitude
      REAL(hp),          INTENT(IN) :: D2RAD     ! Degree to radiance
! !RETURN VALUE:
!
      REAL(hp)                      :: SINbeta   ! Sin of the local solar angle
!
! !REMARKS:
!  References (see above for full citations):
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN v2.1 user mannual 2007-09
!
! !REVISION HISTORY:
!  (1 ) This code was taken directly from the MEGAN v2.1 source code.(mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !REAL(hp) :: BETA                        ! solar elevation angle
      REAL(hp) :: sindelta, cosdelta, A, B

      ! Calculation of sin beta
      sindelta = -SIN( 0.40907_hp ) *            &           
                 COS( 6.28_hp * ( DOY + 10_dp ) / 365_dp )

      cosdelta = (1-sindelta**2.0_hp)**0.5_hp

      A = SIN( LAT * D2RAD ) * sindelta
      B = COS( LAT * D2RAD ) * cosdelta

      SINbeta = A + B *      &
               COS( 2.0_hp * PI * ( SHOUR-12_dp )/24_dp )

      END FUNCTION SOLAR_ANGLE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_T_LI
!
! !DESCRIPTION: Function Get\_Gamma\_T\_LI computes the temperature activity
!  factor (GAMMA\_T\_LI) for the light-independent fraction of emissions
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_T_LI( T, BETA, T_Leaf_Int, T_Leaf_Temp) RESULT( GAMMA_T_LI )
!
! !INPUT PARAMETERS:
!
      ! Current leaf temperature, the surface air temperature field (TS)
      ! is assumed equivalent to the leaf temperature over forests.
      REAL(hp),  INTENT(IN) :: T
      REAL(hp), INTENT(IN) :: T_Leaf_Int
      REAL(hp), INTENT(IN) :: T_Leaf_Temp

      ! Temperature factor per species
      REAL(hp),  INTENT(IN) :: BETA
!
! !RETURN VALUE:
!
      ! Activity factor for the light-independent fraction of emissions
      REAL(hp)              :: GAMMA_T_LI
      REAL(hp)              :: L_T
!
! !REMARKS:
!  GAMMA_T =  exp[Beta*(T - T_Standard)]
!                                                                             .
!             where Beta   = temperature dependent parameter
!                   Ts     = standard temperature (normally 303K, 30C)
!                                                                             .
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
!
! !REVISION HISTORY:
!  (1 ) Original code by Michael Barkley (2009).
!       Note: If T = Ts  (i.e. standard conditions) then GAMMA_T = 1
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  (2 ) Modified to GAMMA_T_LI (dbm, 6/21/2012)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Standard reference temperature [K]
      REAL*8, PARAMETER   :: T_STANDARD = 303.d0
      L_T = T * T_Leaf_Temp + T_Leaf_Int

      !=================================================================
      ! GET_GAMMAT_T_LI begins here!
      !=================================================================

      GAMMA_T_LI = EXP( BETA * ( T - T_STANDARD ) )

      END FUNCTION GET_GAMMA_T_LI
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_T_LD
!
! !DESCRIPTION: Function Get\_Gamma\_T\_LD computes the temperature
!  sensitivity for the light-dependent fraction of emissions.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_T_LD( T, PT_15, PT_1, CT1, CEO )  RESULT( GAMMA_T_LD )
!
! !INPUT PARAMETERS:
!
      ! Current leaf temperature [K], the surface air temperature field (TS)
      ! is assumed equivalent to the leaf temperature over forests.
      REAL(hp), INTENT(IN) :: T

       ! Average leaf temperature over the past 15 days
      REAL(sp), INTENT(IN) :: PT_15

      ! Average leaf temperature over the past arbitray day(s).
      ! This is not used at present
      REAL(sp), INTENT(IN) :: PT_1

      ! Compound-specific parameters for light-dependent temperature activity
      ! factor (dbm, 6/21/2012)
      REAL(hp), INTENT(IN) :: CT1, CEO
!
! !RETURN VALUE:
!
      ! Temperature activity factor for the light-dependent fraction of
      ! emissions
      REAL(hp)             :: GAMMA_T_LD
!
! !REMARKS:
!  References (see above for full citations):
!  (1 ) Guenther et al, 1995
!  (2 ) Guenther et al, 2006
!  (3 ) Guenther et al, MEGAN v2.1 user mannual 2007-08
!  (4 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
!
! !REVISION HISTORY:
!  (1 ) Includes the latest MEGAN v2.1 temperature algorithm (mpb, 2009).
!       Note, this temp-dependence is the same for the PCEEA & hybrid models.
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  (2 ) Modified to gamma_t_ld and to permit compound specific parameters
!       CT1 and Ceo (dbm, 6/21/2012)
!  07 Jan 2016 - Update ideal gas constant to NIST 2014 value
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)              :: C_T,   CT2
      REAL(hp)              :: E_OPT, T_OPT, X
!
! !DEFINED PARAMETERS:
!
      ! Ideal gas constant [J/mol/K]
      REAL(hp), PARAMETER   :: R   = 8.3144598e-3_hp

      !=================================================================
      ! GET_GAMMA_T_LD begins here!
      !=================================================================
      E_OPT = CEO * EXP( 0.08_hp * ( PT_15  - 2.97e2_hp ) )
      T_OPT = 3.13e2_hp + ( 6.0e-1_hp * ( PT_15 - 2.97e2_hp ) )
      CT2   = 200.0_hp

      ! Variable related to temperature
      X     = ( 1.0_hp/T_OPT - 1.0_hp/T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including
      ! effect of average temperature over previous 15 days, based on
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) /              &
             ( CT2 - CT1 * ( 1.0_hp - EXP( CT2 * X ) ) )

      ! Hourly emission activity = C_T
      ! Prevent negative values
      GAMMA_T_LD = MAX( C_T, 0.0_hp )

      END FUNCTION GET_GAMMA_T_LD
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_Lai
!
! !DESCRIPTION: Function Get\_Gamma\_Lai computes the gamma exchange activity
!  factor which is sensitive to leaf area (= GAMMA\_LAI).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_LAI( CMLAI, BIDIREXCH )  RESULT( GAMMA_LAI )
!
! !INPUT PARAMETERS:
!
      REAL(hp),     INTENT(IN) :: CMLAI       ! Current month's LAI [cm2/cm2]
      LOGICAL,      INTENT(IN) :: BIDIREXCH   ! Logical flag indicating whether
                                              ! the compound undergoes bidirectional
                                              ! exchange
!
! !RETURN VALUE:
!
      REAL(hp)             :: GAMMA_LAI
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
!
! !REVISION HISTORY:
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard
!        code by May Fu (11/2004)
!  (2 ) Update to publically released (as of 11/2004) MEGAN algorithm and
!        modified for the standard code by May Fu (11/2004).
!  (3 ) Algorithm is based on the latest MEGAN v2.1 User's Guide (mpb,2009)
!  (4 ) Updated to treat bidirectional exchange compounds appropriately (dbm, 6/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !-----------------------
      ! Compute GAMMA_LAI
      !-----------------------

      ! Formulation for birectional compounds is as described for
      ! ALD2 in Millet et al., ACP 2010
      IF ( BIDIREXCH ) THEN

         IF ( CMLAI <= 6.0_hp) THEN

            ! if lai less than 2:
            IF ( CMLAI <= 2.0_hp ) THEN
               GAMMA_LAI = 0.5_hp * CMLAI

            ! if between 2 and 6:
            ELSE
               GAMMA_LAI = 1.0_hp - 0.0625_hp * ( CMLAI - 2.0_hp )
            END IF

         ELSE
            ! keep at 0.75 for LAI > 6
            GAMMA_LAI = 0.75_hp
         END IF

      ! For all other compounds use the standard gamma_lai formulation
      ELSE
!         GAMMA_LAI = 0.49_hp * CMLAI / SQRT( 1.0_hp + 0.2_hp *
!     &               CMLAI*CMLAI)
         GAMMA_LAI = 1.0_hp
      ENDIF

      END FUNCTION GET_GAMMA_LAI
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_T_LD_C
!
! !DESCRIPTION: Function Get\_Gamma\_T\_LD_C computes the temperature
!  sensitivity for the light-dependent fraction of emissions using the updated
!  Canopy Model (sjs).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_T_LD_C( T, PT_15, PT_24, CT1, CEO, T_Leaf_Int, T_Leaf_Temp)  &
               RESULT( GAMMA_T_LD_C )
!
! !INPUT PARAMETERS:
!
      ! Current leaf temperature [K], the surface air temperature field (TS)
      ! is assumed equivalent to the leaf temperature over forests.
      REAL(hp), INTENT(IN) :: T
      REAL(hp), INTENT(IN) :: T_Leaf_Int
      REAL(hp), INTENT(IN) :: T_Leaf_Temp

      ! Average leaf temperature over the past 15 days
      ! This is not used at present
      REAL(sp), INTENT(IN) :: PT_15

      ! Average leaf temperature over the past day.
      REAL(sp), INTENT(IN) :: PT_24

      ! Compound-specific parameters for light-dependent temperature activity
      ! factor (dbm, 6/21/2012)
      REAL(hp), INTENT(IN) :: CT1, CEO
!
! !RETURN VALUE:
!
      ! Temperature activity factor for the light-dependent fraction of
      ! emissions
      REAL(hp)             :: GAMMA_T_LD_C
      REAL(hp)              :: C_T,   CT2
      REAL(hp)              :: E_OPT, T_OPT, X
      REAL(hp)              :: L_T, L_PT_T

!
! !DEFINED PARAMETERS:
!
      ! Ideal gas constant [J/mol/K]
      REAL(hp), PARAMETER   :: R   = 8.3144598e-3_hp

      L_T = T * T_Leaf_Temp + T_Leaf_Int
      L_PT_T = PT_24 * T_Leaf_Temp + T_Leaf_Int

      !=================================================================
      ! GET_GAMMA_T_LD begins here!
      !=================================================================
      E_OPT = CEO * EXP( 0.1_hp * ( L_PT_T  - 2.97e2_hp ) )
      T_OPT = 3.125e2_hp + ( 6.0e-1_hp * ( L_PT_T - 2.97e2_hp ) )
      CT2   = 230.0_hp

      ! Variable related to temperature
      X     = ( 1.0_hp/T_OPT - 1.0_hp/L_T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including
      ! effect of average temperature over previous 15 days, based on
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) /      &
              ( CT2 - CT1 * ( 1.0_hp - EXP( CT2 * X ) ) )

      ! Hourly emission activity = C_T
      ! Prevent negative values
      IF (T < 260) THEN
        GAMMA_T_LD_C = 0.0_hp
      ELSE
        GAMMA_T_LD_C = MAX( C_T, 0.0_hp )
      ENDIF


      END FUNCTION GET_GAMMA_T_LD_C

!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_PAR_C
!
! !DESCRIPTION: Computes the PCEEA gamma activity factor with sensitivity
!  to LIGHT.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_PAR_C( Q_DIR_2, Q_DIFF_2,                     &
                               PARDR_AVG_SIM, PARDF_AVG_SIM,          &
                               P_Leaf_LAI, P_Leaf_Int, LAI, PSTD)     &
                   RESULT( GAMMA_P_C )
!
! !INPUT PARAMETERS:
!
      !LOGICAL,         INTENT(IN) :: am_I_Root
      !TYPE(HCO_State), POINTER    :: HcoState
      !TYPE(Ext_State), POINTER    :: ExtState
      !TYPE(MyInst),    POINTER    :: Inst
      !INTEGER,         INTENT(IN) :: I,  J             ! Lon & lat indices
      REAL(sp),        INTENT(IN) :: PARDR_AVG_SIM     ! Avg direct PAR [W/m2]
      REAL(sp),        INTENT(IN) :: PARDF_AVG_SIM     ! Avg diffuse PAR [W/m2]
      REAL(hp),        INTENT(IN) :: Q_DIR_2           ! Direct PAR [umol/m2/s]
      REAL(hp),        INTENT(IN) :: Q_DIFF_2          ! Diffuse PAR [umol/m2/s]
      REAL(hp),        INTENT(IN) :: P_Leaf_LAI
      REAL(hp),        INTENT(IN) :: P_Leaf_Int
      REAL(hp),        INTENT(IN) :: LAI
      REAL(hp),        INTENT(IN) :: PSTD
!
! !RETURN VALUE:
!
      REAL(hp)                    :: GAMMA_P_C     ! GAMMA factor for light
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 user guide
!
! !REVISION HISTORY:
!  (1 ) Here PAR*_AVG_SIM is the average light conditions over the simulation
!       period. I've set this = 10 days to be consistent with temperature & as
!       outlined in Guenther et al, 2006. (mpb,2009)
!  (2 ) Code was taken & adapted directly from the MEGAN v2.1 source code.
!       (mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)   :: mmPARDR_DAILY
      REAL(hp)   :: mmPARDF_DAILY
      REAL(hp)   :: PAC_DAILY, PAC_INSTANT, C_PPFD
      REAL(hp)   :: PTOA
      !REAL(hp)   :: LocalHour
      !REAL(hp)   :: LAT
      REAL(hp)   :: C1
      REAL(hp)   :: Alpha

      !-----------------------------------------------------------------
      ! Compute GAMMA_PAR_C
      !-----------------------------------------------------------------

      ! Initialize
      C_PPFD   = 0.0_hp
      PTOA     = 0.0_hp

      ! Convert past light conditions to micromol/m2/s
      mmPARDR_DAILY   = PARDR_AVG_SIM 
      mmPARDF_DAILY   = PARDF_AVG_SIM

      ! Work out the light at the top of the canopy.
      PAC_DAILY    = mmPARDR_DAILY + mmPARDF_DAILY
      PAC_INSTANT  = Q_DIR_2       +  Q_DIFF_2

      PAC_DAILY = PAC_DAILY * exp(P_Leaf_Int + P_Leaf_LAI * LAI)
      PAC_INSTANT = PAC_INSTANT * exp(P_Leaf_Int + P_Leaf_LAI * LAI)

      IF ( PAC_DAILY < 0.01_hp ) THEN

         GAMMA_P_C = 0.0_hp

      ELSE
        Alpha  = 0.004
        Alpha  = 0.004 - 0.0005*LOG(PAC_DAILY)
        C1 = 1.03
        C1 = 0.0468 * EXP(0.0005 * (PAC_DAILY - PSTD)) *           &
                      (PAC_DAILY **  0.6)
        GAMMA_P_C = (Alpha * C1 * PAC_INSTANT) /                   &
                        ((1 + Alpha**2. * PAC_INSTANT**2.)**0.5)
      ENDIF
      ! Prevent negative values
      GAMMA_P_C = MAX( GAMMA_P_C , 0.0_hp )

      END FUNCTION GET_GAMMA_PAR_C

!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_GET_CDEA
!
! !DESCRIPTION: Function Get\_GET_CDEA computes
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CDEA( CMLAI ) RESULT( CDEA )

!
! !INPUT PARAMETERS:
!
      REAL(hp), INTENT(IN) :: CMLAI       ! Current month's LAI [cm2/cm2]
      REAL(hp)             :: Cdepth(5)
      REAL(hp) :: LAIdepth
      INTEGER                 ::  K
!
! !RETURN VALUE:
!
      REAL(hp)             :: CDEA(5)

!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER   :: CCD1 = -0.2_hp
      REAL*8, PARAMETER   :: CCD2 = 1.3_hp
      !=================================================================
      ! GET_CDEA begins here!
      !=================================================================
      Cdepth (1)   = 0.0469101
      Cdepth (2)   = 0.2307534
      Cdepth (3)   = 0.5
      Cdepth (4)   = 0.7692465
      Cdepth (5)   = 0.9530899
      DO K = 1, 5
        LAIdepth = CMLAI * Cdepth(K)
        IF ( LAIdepth .GT. 3 ) THEN
           LAIdepth = 3.0
        ENDIF
        CDEA(K) = CCD1 * LAIdepth + CCD2
      ENDDO

      RETURN

      END FUNCTION GET_CDEA

!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_Age
!
! !DESCRIPTION: Function Get\_Gamma\_Age computes the gamma exchange
!  activity factor which is sensitive to leaf age (= Gamma\_Age).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_AGE( CMLAI, PMLAI, DBTWN, TT, AN, AG, AM, AO ) RESULT ( GAMMA_AGE )
!
! !INPUT PARAMETERS:
!
      REAL(hp), INTENT(IN) :: CMLAI     ! Current month's LAI [cm2/cm2]
      REAL(hp), INTENT(IN) :: PMLAI     ! Previous months LAI [cm2/cm2]
      REAL(hp), INTENT(IN) :: DBTWN     ! Number of days between
      REAL(sp), INTENT(IN) :: TT        ! Daily average temperature [K]
      REAL(hp), INTENT(IN) :: AN        ! Relative emiss factor (new leaves)
      REAL(hp), INTENT(IN) :: AG        ! Relative emiss factor (growing leaves)
      REAL(hp), INTENT(IN) :: AM        ! Relative emiss factor (mature leaves)
      REAL(hp), INTENT(IN) :: AO        ! Relative emiss factor (old leaves)
!
! !RETURN VALUE:
!
      REAL(hp)             :: GAMMA_AGE ! Activity factor
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
!
! !REVISION HISTORY:
!  (1 ) Original code by Dorian Abbot (9/2003). Modified for the standard
!        code by May Fu (11/2004)
!  (2 ) Update to publically released (as of 11/2004) MEGAN algorithm and
!        modified for the standard code by May Fu (11/2004).
!  (3 ) Algorithm is based on the latest User's Guide (tmf, 11/19/04)
!  (4 ) Renamed & now includes specific relative emission activity factors for
!       each BVOC based on MEGAN v2.1 algorithm (mpb,2008)
!  (5 ) Now calculate TI (number of days after budbreak required to induce
!       iso. em.) and TM (number of days after budbreak required to reach
!       peak iso. em. rates) using the daily average temperature, instead
!       of using fixed values (mpb,2008)
!       NOTE: Can create 20% increases in tropics (Guenther et al 2006)
!  (6 ) Implemented change for the calculation of FGRO if ( CMLAI > PMLAI ),
!       i.e. if LAI has increased with time, and used new values for
!       all foilage fractions if ( CMLAI = PMLAI ). Also removed TG variable
!       as not now needed. (mpb,2000)
!  (7 ) Updated to pass leaf age activity factors as arguments (dbm, 6/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  13 Aug 2013 - M. Sulprizio- Updated for sesquiterpenes (H. Pye)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)             :: FNEW  ! Fraction of new leaves in canopy
      REAL(hp)             :: FGRO  ! Fraction of growing leaves
      REAL(hp)             :: FMAT  ! Fraction of mature leaves
      REAL(hp)             :: FOLD  ! Fraction of old leaves

      ! TI: number of days after budbreak required to induce emissions
      REAL(hp)             :: TI

      ! TM: number of days after budbreak required to reach peak emissions
      REAL(hp)             :: TM

      !=================================================================
      ! GET_GAMMA_AGE begins here!
      !=================================================================

      !-----------------------
      ! Compute TI and TM
      ! (mpb,2009)
      !-----------------------

      IF ( TT <= 303.0_hp ) THEN
         TI = 5.0_hp + 0.7_hp * ( 300.0_hp - TT )
      ELSEIF ( TT >  303.0_hp ) THEN
         TI = 2.9_hp
      ENDIF
      TM = 2.3_hp * TI

      !-----------------------
      ! Compute GAMMA_AGE
      !-----------------------

      IF ( CMLAI == PMLAI ) THEN !(i.e. LAI stays the same)

         FNEW = 0.0_hp
         FGRO = 0.1_hp
         FMAT = 0.8_hp
         FOLD = 0.1_hp

      ELSE IF ( CMLAI > PMLAI ) THEN !(i.e. LAI has INcreased)

         ! Calculate Fnew
         IF ( DBTWN > TI ) THEN
            FNEW = ( TI / DBTWN ) * ( 1.0_hp -  PMLAI / CMLAI )
         ELSE
            FNEW = 1.0_hp - ( PMLAI / CMLAI )
         ENDIF

         ! Calculate FMAT
         IF ( DBTWN > TM ) THEN
            FMAT = ( PMLAI / CMLAI ) +              &
                 (( DBTWN - TM ) / DBTWN )*( 1.0_hp -  PMLAI / CMLAI )
         ELSE
            FMAT = ( PMLAI / CMLAI )
         ENDIF

         ! Calculate Fgro and Fold
         FGRO = 1.0_hp - FNEW - FMAT
         FOLD = 0.0_hp

      ELSE ! This is the case if  PMLAI > CMLAI (i.e. LAI has DEcreased)

         FNEW = 0.0_hp
         FGRO = 0.0_hp
         FOLD = ( PMLAI - CMLAI ) / PMLAI
         FMAT = 1.0_hp - FOLD

      ENDIF

      ! Age factor
      GAMMA_AGE = FNEW * AN + FGRO * AG +    &
                 FMAT * AM + FOLD * AO

      ! Prevent negative values
      GAMMA_AGE = MAX( GAMMA_AGE , 0.0_hp )

      END FUNCTION GET_GAMMA_AGE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gamma_sm
!
! !DESCRIPTION: Function GET\_GAMMA\_SM computes activity factor for soil
!  moisture
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_SM( GWETROOT, CMPD )   RESULT( GAMMA_SM )
!
! !INPUT PARAMETERS:
!
      !TYPE(Ext_State),  POINTER     :: ExtState
      !INTEGER,          INTENT(IN)  :: I, J     ! GEOS-Chem lon & lat indices
      CHARACTER(LEN=*), INTENT(IN)  :: CMPD     ! Compound name (dbm, 6/21/2012)
      REAL(hp),         INTENT(IN)  :: GWETROOT

! !RETURN VALUE:
!
      REAL(hp)                      :: GAMMA_SM ! Activity factor
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, ACP 2006
!  (2 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
!
! !REVISION HISTORY:
!  (1 ) Created by dbm (6/2012). We are not currently using a soil moisture
!       effect for isoprene. For all compounds other than acetaldehyde and
!       ethanol, gamma_sm =1 presently.
!  16 Apr 2015 - C. Keller   - Now restrict GWETROOT to values between 0.0 and
!                              1.0. This only seems to be a problem within the
!                              GEOS-5 ESM, where GWETROOT values over the ocean
!                              become 1e+15 (= missing value).
!  12 Aug 2015 - R. Yantosca - Extend #ifdef for MERRA2 meteorology

!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)  :: GWETROOT2

      !=================================================================
      ! GET_GAMMA_SM begins here!
      !=================================================================

      ! By default gamma_sm is 1.0
      GAMMA_SM = 1.0_hp

      ! Error trap: GWETROOT must be between 0.0 and 1.0 (ckeller, 4/16/15)
      GWETROOT2 = MIN(MAX(GWETROOT,0.0_hp),1.0_hp)

      IF ( TRIM( CMPD ) == 'ALD2' .OR. TRIM ( CMPD ) == 'EOH' ) THEN

         ! GWETROOT = degree of saturation or wetness in the root-zone
         ! (top meter of soil). This is defined as the ratio of the volumetric
         ! soil moisture to the porosity. We use a soil moisture activity factor
         ! for ALD2 to account for stimulation of emission by flooding.
         ! (Millet et al., ACP 2010)
         ! Constant value of 1.0 for GWETROOT = 0-0.9, increasing linearly to
         ! 3.0 at GWETROOT =1.
         GAMMA_SM = MAX( 20.0_hp * GWETROOT2 - 17.0_hp, 1.0_hp)

      ENDIF

      ! return to calling program
      END FUNCTION GET_GAMMA_SM

!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gamma_co2
!
! !DESCRIPTION: Function GET\_GAMMA\_CO2 computes the CO2 activity factor
!  associated with CO2 inhibition of isoprene emission. Called from
!  GET\_MEGAN\_EMISSIONS only.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_CO2( CO2a ) RESULT( GAMMA_CO2 )
!
! !INPUT PARAMETERS:
      REAL(hp), INTENT(IN) :: CO2a       ! Atmospheric CO2 conc [ppmv]
!
! !RETURN VALUE:
      REAL(hp)             :: GAMMA_CO2  ! CO2 activity factor [unitless]
!
! !LOCAL VARIABLES:
      REAL(hp)             :: CO2i       ! Intercellular CO2 conc [ppmv]
      REAL(hp)             :: ISMAXi     ! Asymptote for intercellular CO2
      REAL(hp)             :: HEXPi      ! Exponent for intercellular CO2
      REAL(hp)             :: CSTARi     ! Scaling coef for intercellular CO2
      REAL(hp)             :: ISMAXa     ! Asymptote for atmospheric CO2
      REAL(hp)             :: HEXPa      ! Exponent for atmospheric CO2
      REAL(hp)             :: CSTARa     ! Scaling coef for atmospheric CO2
      LOGICAL              :: LPOSSELL   ! Use Possell & Hewitt (2011)?
      LOGICAL              :: LWILKINSON ! Use Wilkinson et al. (2009)?

!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Heald, C. L., Wilkinson, M. J., Monson, R. K., Alo, C. A.,
!       Wang, G. L., and Guenther, A.: Response of isoprene emission
!       to ambient co(2) changes and implications for global budgets,
!       Global Change Biology, 15, 1127-1140, 2009.
!  (2 ) Wilkinson, M. J., Monson, R. K., Trahan, N., Lee, S., Brown, E.,
!       Jackson, R. B., Polley, H. W., Fay, P. A., and Fall, R.: Leaf
!       isoprene emission rate as a function of atmospheric CO2
!       concentration, Global Change Biology, 15, 1189-1200, 2009.
!  (3 ) Possell, M., and Hewitt, C. N.: Isoprene emissions from plants
!       are mediated by atmospheric co2 concentrations, Global Change
!       Biology, 17, 1595-1610, 2011.
!
! !REVISION HISTORY:
!  (1 ) Implemented in the standard code by A. Tai (Jun 2012).
!  15 Sep 2015 - M. Sulprizio- Implemented into hcox_megan_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC

      !----------------------------------------------------------
      ! Choose between two alternative CO2 inhibition schemes
      !----------------------------------------------------------

      ! Empirical relationship of Possell & Hewitt (2011) based on nine
      ! experimental studies including Wilkinson et al. (2009). This is
      ! especially recommended for sub-ambient CO2 concentrations:
      LPOSSELL    = .TRUE.   ! Default option

      ! Semi-process-based parameterization of Wilkinson et al. (2009),
      ! taking into account of sensitivity to intercellular CO2
      ! fluctuation, which is here set as a constant fraction of
      ! atmospheric CO2:
      LWILKINSON  = .FALSE.   ! Set .TRUE. only if LPOSSELL = .FALSE.

      !-----------------------
      ! Compute GAMMA_CO2
      !-----------------------

      IF ( LPOSSELL ) THEN

         ! Use empirical relationship of Possell & Hewitt (2011):
         GAMMA_CO2 = 8.9406_hp /       &
                    ( 1.0_hp + 8.9406_hp * 0.0024_hp * CO2a )

      ELSEIF ( LWILKINSON ) THEN

         ! Use parameterization of Wilkinson et al. (2009):

         ! Parameters for intercellular CO2 using linear interpolation:
         IF ( CO2a <= 600.0_hp ) THEN
            ISMAXi = 1.036_hp  - (1.036_hp - 1.072_hp) /            &
                     (600.0_hp - 400.0_hp) * (600.0_hp - CO2a)
            HEXPi  = 2.0125_hp - (2.0125_hp - 1.7000_hp) /          &
                     (600.0_hp - 400.0_hp) * (600.0_hp - CO2a)
            CSTARi = 1150.0_hp - (1150.0_hp - 1218.0_hp) /          &
                     (600.0_hp - 400.0_hp) * (600.0_hp - CO2a)
         ELSEIF ( CO2a > 600.0_hp .AND. CO2a < 800.0_hp ) THEN
            ISMAXi = 1.046_hp  - (1.046_hp - 1.036_hp) /            &
                     (800.0_hp - 600.0_hp) * (800.0_hp - CO2a)
            HEXPi  = 1.5380_hp - (1.5380_hp - 2.0125_hp) /          &
                     (800.0_hp - 600.0_hp) * (800.0_hp - CO2a)
            CSTARi = 2025.0_hp - (2025.0_hp - 1150.0_hp) /          &
                     (800.0_hp - 600.0_hp) * (800.0_hp - CO2a)
         ELSE
            ISMAXi = 1.014_hp - (1.014_hp - 1.046_hp) /             &
                     (1200.0_hp - 800.0_hp) * (1200.0_hp - CO2a)
            HEXPi  = 2.8610_hp - (2.8610_hp - 1.5380_hp) /          &
                     (1200.0_hp - 800.0_hp) * (1200.0_hp - CO2a)
            CSTARi = 1525.0_hp - (1525.0_hp - 2025.0_hp) /          &
                     (1200.0_hp - 800.0_hp) * (1200.0_hp - CO2a)
         ENDIF

         ! Parameters for atmospheric CO2:
         ISMAXa    = 1.344_hp
         HEXPa     = 1.4614_hp
         CSTARa    = 585.0_hp

         ! For now, set CO2_Ci = 0.7d0 * CO2_Ca as recommended by Heald
         ! et al. (2009):
         CO2i      = 0.7_hp * CO2a

         ! Compute GAMMA_CO2:
         GAMMA_CO2 = ( ISMAXi -  ISMAXi * CO2i**HEXPi /                 &
                     ( CSTARi**HEXPi + CO2i**HEXPi ) )                  &
                  *  ( ISMAXa - ISMAXa * ( 0.7_hp * CO2a )**HEXPa /     & 
                     ( CSTARa**HEXPa + ( 0.7_hp * CO2a )**HEXPa ) )

      ELSE

         ! No CO2 inhibition scheme is used; GAMMA_CO2 set to unity:
         GAMMA_CO2 = 1.0_hp

      ENDIF

      END FUNCTION GET_GAMMA_CO2
!EOC

      END MODULE HCOX_MEGAN_MOD

