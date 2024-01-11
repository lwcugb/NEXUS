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

!
! !INTERFACE:
!
      MODULE HCOX_MEGAN_WRAPPER
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_DIAGN_MOD
      USE HCOX_State_MOD,    ONLY : Ext_State
      USE HCO_STATE_MOD,     ONLY : HCO_STATE

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: HCOX_Megan_Init
      PUBLIC  :: HCOX_Megan_Run
      PUBLIC  :: HCOX_Megan_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
!      PRIVATE :: GET_MEGAN_EMISSIONS  ! dbm, new MEGAN driver routine
                                      ! for all compounds (6/21/2012)
!      PRIVATE :: UPDATE_T_DAY
!      PRIVATE :: UPDATE_T_15_AVG
!      PRIVATE :: GET_MEGAN_PARAMS
!      PRIVATE :: GET_MEGAN_AEF
!      PRIVATE :: GET_GAMMA_PAR_PCEEA
!      PRIVATE :: GET_GAMMA_T_LI
!      PRIVATE :: GET_GAMMA_T_LD
!      PRIVATE :: GET_GAMMA_LAI
!      PRIVATE :: GET_GAMMA_AGE
!      PRIVATE :: GET_GAMMA_SM
!      PRIVATE :: CALC_NORM_FAC
!      PRIVATE :: SOLAR_ANGLE
      PRIVATE :: FILL_RESTART_VARS
      PRIVATE :: CALC_AEF
!      PRIVATE :: GET_GAMMA_CO2  ! (Tai, Jan 2013)
!
! !REVISION HISTORY:
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (6/2003).  Updated to
!        latest algorithm and modified for the standard code by May Fu
!        (11/2004).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
      TYPE :: MyInst
         INTEGER                 :: Instance
         INTEGER                 :: ExtNr        ! Extension number for MEGAN
         INTEGER                 :: ExtNrMono    ! Extension number for MEGAN_Mono
         REAL(hp)                :: ISOP_SCALING ! Isoprene emissions scaling
         REAL(hp)                :: GLOBCO2      ! Global CO2 conc (ppmv)
         REAL(hp)                :: ISOPTOSOAP   ! SOAP Emission factor
         REAL(hp)                :: MONOTOSOAP   ! SOAP Emission factor
         REAL(hp)                :: OTHRTOSOAP   ! SOAP Emission factor
         REAL(hp)                :: ISOPTOSOAS   ! Direct SOA Emission factor
         REAL(hp)                :: MONOTOSOAS   ! Direct SOA Emission factor
         REAL(hp)                :: OTHRTOSOAS   ! Direct SOA Emission factor
         LOGICAL                 :: LISOPCO2     ! Include CO2 inhibition of ISOP?
         LOGICAL                 :: NORMLAI      ! Normalize LAI by PFT?

         ! Physical parameter
         REAL(hp)                :: D2RAD  ! Degrees to radians
         REAL(hp)                :: RAD2D  ! Radians to degrees

         ! HEMCO species IDs
         INTEGER                 :: IDTISOP
         INTEGER                 :: IDTACET
         INTEGER                 :: IDTCO
         INTEGER                 :: IDTPRPE
         INTEGER                 :: IDTC2H4
         INTEGER                 :: IDTALD2
         INTEGER                 :: IDTEOH
!         INTEGER                 :: IDTOCPI
         INTEGER                 :: IDTLIMO
         INTEGER                 :: IDTMTPA
         INTEGER                 :: IDTMTPO
         INTEGER                 :: IDTSESQ
         INTEGER                 :: IDTSOAP
         INTEGER                 :: IDTSOAS

         ! Pointers to annual emission factor arrays.
         ! These fields are obtained from ext. data (from config file)
         ! These arrays represent the base emission values in kg(C)/m2/s
         ! that will be scaled based upon meteorological conditions.
         REAL(hp), POINTER       :: AEF_ISOP(:,:)
         REAL(hp), POINTER       :: AEF_MBOX(:,:)
         REAL(hp), POINTER       :: AEF_BPIN(:,:)
         REAL(hp), POINTER       :: AEF_CARE(:,:)
         REAL(hp), POINTER       :: AEF_LIMO(:,:)
         REAL(hp), POINTER       :: AEF_OCIM(:,:)
         REAL(hp), POINTER       :: AEF_SABI(:,:)
         REAL(hp), POINTER       :: GEIA_ORVC(:,:)

         ! Annual emission factors arrays
         ! These fields are not read from file, but are computed in CALC_AEF
         REAL(hp), POINTER       :: AEF_APIN(:,:)              ! Alpha-pinene
         REAL(hp), POINTER       :: AEF_MYRC(:,:)              ! Myrcene
         REAL(hp), POINTER       :: AEF_OMON(:,:)              ! Other monoterpenes
         REAL(hp), POINTER       :: AEF_MOHX(:,:)              ! Methanol
         REAL(hp), POINTER       :: AEF_ACET(:,:)              ! Acetone
         REAL(hp), POINTER       :: AEF_EOH(:,:)               ! Ethanol
         REAL(hp), POINTER       :: AEF_CH2O(:,:)              ! Formaldehyde
         REAL(hp), POINTER       :: AEF_ALD2(:,:)              ! Acetaldehyde
         REAL(hp), POINTER       :: AEF_FAXX(:,:)              ! Formic acid
         REAL(hp), POINTER       :: AEF_AAXX(:,:)              ! Acetic acid
         REAL(hp), POINTER       :: AEF_C2H4(:,:)              ! Ethene
         REAL(hp), POINTER       :: AEF_TOLU(:,:)              ! Toluene
         REAL(hp), POINTER       :: AEF_HCNX(:,:)              ! HCN
         REAL(hp), POINTER       :: AEF_PRPE(:,:)              ! >= C3 alkenes
         REAL(hp), POINTER       :: AEF_FARN(:,:)              ! a-Farnesene
         REAL(hp), POINTER       :: AEF_BCAR(:,:)              ! b-Caryophyllene
         REAL(hp), POINTER       :: AEF_OSQT(:,:)              ! Other sesquiterp.

         ! Normalization factor
         REAL(hp), POINTER       :: NORM_FAC(:)

         ! New restart variables (ckeller, 11/05/2015)
         REAL(sp), POINTER       :: T_LASTXDAYS    (:,:)  ! Avg. temperature of last NUM_DAYS
         REAL(sp), POINTER       :: T_LAST24H      (:,:)  ! Avg. temperature of last 24 hours
         REAL(sp), POINTER       :: PARDF_LASTXDAYS(:,:)  ! Avg. PARDF of last NUM_DAYS
         REAL(sp), POINTER       :: PARDR_LASTXDAYS(:,:)  ! Avg. PARDR of last NUM_DAYS

         ! previous LAI values (ckeller, 10/09/2014)
         REAL(sp), POINTER       :: LAI_PREVDAY(:,:)      ! LAI of prev. day

         ! Array for PFT
         REAL(hp), POINTER       :: ARRAY_16(:,:,:)

         ! Days between mid-months (updated by HEMCO clock)
         INTEGER                 :: DAYS_BTW_M

         ! Input data suffix
         CHARACTER(LEN=255)      :: SUFFIX

         TYPE(MyInst), POINTER   :: NextInst => NULL()
      END TYPE MyInst

      ! Pointer to all instances
      TYPE(MyInst), POINTER      :: AllInst => NULL()

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
! !IROUTINE: HCOX_Megan_Run
!
! !DESCRIPTION: Subroutine HCOX\_MEGAN\_Run is the driver routine for
! the MEGAN model within the new emissions structure.
! Note that all diagnostics are commented since those are still written
! as part of the old emission structure.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_Megan_Run( ExtState, HcoState, RC )
!
! !USES:
!
      USE HCO_FLUXARR_MOD,     ONLY : HCO_EmisAdd
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_First
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_Rewind
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_NewHour
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_NewDay
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_Get    !Wei Li
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_GetLocal !Wei Li
      USE HCO_Restart_Mod,     ONLY : HCO_RestartWrite
      USE HCOX_MEGAN_MOD       !core megan module (Wei Li)
!
! !INPUT PARAMETERS:
!
      !LOGICAL,         INTENT(IN)    :: am_I_Root
      TYPE(Ext_State), POINTER       :: ExtState
      TYPE(HCO_State), POINTER       :: HcoState
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  05 Aug 2013 - C. Keller   - Initial version
!  08 Oct 2014 - C. Keller   - Now use HEMCO clock to get days between months
!  09 Oct 2014 - C. Keller   - Now use only GC_LAI. This makes days between months
!                              obsolete
!  16 Oct 2014 - C. Keller   - Initialize flux arrays to avoid float invalid
!   9 Mar 2015 - R. Yantosca - Bug fix: add EMIS_ALD2 to $OMP PRIVATE clause
!  08 May 2015 - C. Keller   - Now read/write restart variables from here to
!                              accomodate replay runs in GEOS-5.
!  30 Sep 2015 - C. Keller   - Now add OCPI flux to BIOGENIC_OCPI diagnostics.
!  12 Dec 2015 - C. Keller   - Apply e-folding times to past conditions.
!  07 Jan 2016 - E. Lundgren - Update Avogadro's # to NIST 2014
!  01 Jun 2016 - R. Yantosca - Bug fix: Now use ExtNrMono in the call to
!                              Hco_EmisAdd for the MONX (in order to make MONX
!                              show up as nonzero in HEMCO diagnostics).
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  02 Mar 2018 - M. Sulprizio- Remove species MONX and use MTPA/MTPO/LIMO;
!                              Consolidate MEGAN_SOA option into MEGAN_Mono
!                              since there is essentially no difference between
!                              the two now that SOA (MTPA/MTPO/LIMO) is standard
!  19 Mar 2019 - R. Yantosca - Bug fix: CO, MTPA, MTPO, LIMO, and SESQ should
!                              use Inst%ExtNrMono in calls to HCO_EmisAdd
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J
      REAL(hp)            :: EMIS_ISOP, EMIS_MBOX, EMIS_APIN, EMIS_BPIN
      REAL(hp)            :: EMIS_LIMO, EMIS_SABI, EMIS_MYRC, EMIS_CARE
      REAL(hp)            :: EMIS_OCIM, EMIS_OMON, EMIS_MONO, EMIS_ALD2
      REAL(hp)            :: EMIS_MOHX, EMIS_EOH,  EMIS_FAXX, EMIS_AAXX
      REAL(hp)            :: EMIS_ACET, EMIS_PRPE, EMIS_C2H4
      REAL(hp)            :: EMIS_FARN, EMIS_BCAR, EMIS_OSQT
      REAL(hp)            :: EMIS_OTHR
      REAL(hp)            :: BIOG_ALPH, BIOG_LIMO
      REAL(hp)            :: BIOG_ALCO
      REAL(hp)            :: MONO_MOL
      REAL(hp)            :: MBOX_MOL
      REAL(hp)            :: X, Y
      REAL(hp)            :: DIURFACT, DIURORVC
      REAL(hp)            :: TMP
      REAL(hp)            :: TS_EMIS
      REAL(hp)            :: DNEWFRAC
      REAL(hp)            :: DOLDFRAC
      REAL(hp)            :: HNEWFRAC
      REAL(hp)            :: HOLDFRAC
      LOGICAL             :: ERR
      LOGICAL             :: FIRST, IsNewHour, IsNewDay
      REAL(hp)            :: PFTSUM  !Wei Li
      REAL(hp)            :: D2RAD_FAC   !Wei Li
      REAL(hp)            :: NORM_FAC_OUT   !Wei Li
      REAL(hp)            :: LocalHour   !Wei Li
      Integer             :: DOY   !Wei Li
      TYPE(MyMet),POINTER :: LocalMet  !Wei Li

      ! Emission arrays (kgC/m2/s)
      REAL(hp), TARGET    :: FLUXISOP  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMONO  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXACETmo(HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXACETmb(HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXACETbg(HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXCO    (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXPRPE  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXC2H4  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXLIMO  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMTPA  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMTPO  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXSESQ  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXSOAP  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXSOAS  (HcoState%NX,HcoState%NY)

      ! Emission arrays for use in diagnostics only (kgC/m2/s)
      REAL(hp), TARGET    :: FLUXMBOX  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXAPIN  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXBPIN  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXSABI  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMYRC  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXCARE  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXOCIM  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXOMON  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXALD2  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMOHX  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXEOH   (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXFAXX  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXAAXX  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXFARN  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXBCAR  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXOSQT  (HcoState%NX,HcoState%NY)

      ! For diagnostics
      REAL(hp), POINTER   :: Arr2D(:,:)
      CHARACTER(LEN=63)   :: DiagnName
      CHARACTER(LEN=255)  :: MSG

      ! Molecules C / kg C
      REAL(hp), PARAMETER   :: XNUMOL_C = 6.022140857e+23_hp / 12e-3_hp
      REAL(hp), PARAMETER   :: UNITCONV = 1.0_hp

      ! Fraction of yield of OC (hydrophilic) from terpene emissions
      REAL(hp), PARAMETER   :: EMMT2OC = 1.0e-1_hp

      ! Conversion factors for OC/BC
      REAL(hp), PARAMETER   :: FC1 = 136.2364_hp / 120.11_hp
      REAL(hp), PARAMETER   :: FC2 = 154.2516_hp / 120.11_hp
      REAL(hp), PARAMETER   :: FC3 = 204.3546_hp / 180.165_hp
      REAL(hp), PARAMETER   :: FC4 = 152.0_hp    / 120.11_hp

      ! Conversion factors for acetone calculations
      REAL(hp), PARAMETER   :: YIELD_MO   = 0.116_hp
      REAL(hp), PARAMETER   :: MONO_SCALE = 0.89_hp
      REAL(hp), PARAMETER   :: MB_SCALE1  = 0.6_hp
      REAL(hp), PARAMETER   :: MB_SCALE2  = 0.76_hp

      TYPE(MyInst), POINTER :: Inst

      !=================================================================
      ! HCOX_Megan_Run begins here!
      !=================================================================

      ! Enter
      CALL HCO_ENTER( HcoState%Config%Err,
     &               'HCOX_Megan_Run (hcox_megan_mod.F)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      ERR = .FALSE.

      ! Nullify
      Arr2d => NULL()
      Inst  => NULL()

      ! Get instance
      CALL InstGet ( ExtState%Megan, Inst, RC )
      IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find Megan instance Nr. ', ExtState%Megan
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
      ENDIF

      ! Time information
      FIRST     = HcoClock_First  ( HcoState%Clock, .TRUE. )
      IsNewHour = HcoClock_NewHour( HcoState%Clock, .TRUE. )
      IsNewDay  = HcoClock_NewDay ( HcoState%Clock, .TRUE. )

      ! Initialize flux arrays
      FLUXISOP   = 0.0_hp
      FLUXMBOX   = 0.0_hp
      FLUXAPIN   = 0.0_hp
      FLUXBPIN   = 0.0_hp
      FLUXSABI   = 0.0_hp
      FLUXMYRC   = 0.0_hp
      FLUXCARE   = 0.0_hp
      FLUXOCIM   = 0.0_hp
      FLUXOMON   = 0.0_hp
      FLUXMONO   = 0.0_hp
      FLUXALD2   = 0.0_hp
      FLUXMOHX   = 0.0_hp
      FLUXEOH    = 0.0_hp
      FLUXFAXX   = 0.0_hp
      FLUXAAXX   = 0.0_hp
      FLUXACETmo = 0.0_hp
      FLUXACETmb = 0.0_hp
      FLUXACETbg = 0.0_hp
      FLUXCO     = 0.0_hp
      FLUXPRPE   = 0.0_hp
      FLUXC2H4   = 0.0_hp
      FLUXLIMO   = 0.0_hp
      FLUXMTPA   = 0.0_hp
      FLUXMTPO   = 0.0_hp
      FLUXFARN   = 0.0_hp
      FLUXBCAR   = 0.0_hp
      FLUXOSQT   = 0.0_hp
      FLUXSESQ   = 0.0_hp
      FLUXSOAP   = 0.0_hp
      FLUXSOAS   = 0.0_hp

#if defined( ESMF_ ) || defined( MODEL_GEOS )

      !----------------------------------------------------------------
      ! %%%%% ESMF environment: execute these on every call  %%%%%
      ! %%%%% because this will fill from the External State %%%%%
      !----------------------------------------------------------------

      ! Generate annual emission factors for MEGAN inventory
      CALL CALC_AEF( HcoState, ExtState, Inst, RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         MSG = 'Error encountered in MEGAN routine CALC_AEF!'
         CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
         RETURN
      ENDIF

      ! Calculate normalization factor (dbm, 11/2012)
      !CALL CALC_NORM_FAC( am_I_Root, Inst, RC )  !(Wei Li)
      D2RAD_FAC=Inst%D2RAD
      CALL CALC_NORM_FAC( D2RAD_FAC, NORM_FAC_OUT, RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         MSG = 'Error encountered in MEGAN routine CALC_NORM_FAC!'
         CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
         RETURN
      ENDIF
      Inst%NORM_FAC = NORM_FAC_OUT   !Wei Li

      ! Fill restart variables.
      CALL FILL_RESTART_VARS(  HcoState, ExtState, Inst, RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         MSG = 'Error encountered in MEGAN routine FILL_RESTART_VARS!'
         CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
         RETURN
      ENDIF

#else

      !----------------------------------------------------------------
      ! %%%%% Standard environment: Execute these only once     %%%%%
      ! %%%%% to avoid constantly overwriting restart variables %%%%%
      !----------------------------------------------------------------
      IF ( FIRST ) THEN

         ! Generate annual emission factors for MEGAN inventory
         CALL CALC_AEF(  HcoState, ExtState, Inst, RC )
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Error encountered in MEGAN routine CALC_AEF!'
            CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
            RETURN
         ENDIF

         ! Calculate normalization factor (dbm, 11/2012)
         !CALL CALC_NORM_FAC( am_I_Root, Inst, RC )   !Wei Li
         CALL CALC_NORM_FAC( D2RAD_FAC, NORM_FAC_OUT, RC )
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Error encountered in MEGAN routine CALC_NORM_FAC!'
            CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
            RETURN
         ENDIF
         Inst%NORM_FAC = NORM_FAC_OUT   !Wei Li

         ! Fill restart variables
         CALL FILL_RESTART_VARS( HcoState,
     &                           ExtState,  Inst,     RC )
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Error encountered in MEGAN routine ' //
     &            'FILL_RESTART_VARS!'
            CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
            RETURN
         ENDIF
      ENDIF

#endif

      ! LAI_PREVDAY is now continuously updated (ckeller, 12/8/15)
!      ! On a new day, copy previous LAI values to array that holds
!      ! previous day LAI.
!      ELSEIF ( IsNewDay ) THEN
!         LAI_PREVDAY = LAI_PREVSTEP
!      ENDIF

      ! Set to 1 day since we know that LAI is updated every day.
      Inst%DAYS_BTW_M = 1

      ! Calculate weights for running means of historic variables
      ! DNEWFRAC and DOLDFRAC are the weights given to the current
      ! and existing value, respectively, when updating running means
      ! over the last X days. HNEWFRAC and HOLDFRAC are the same but
      ! for the 24H means. (ckeller, 11/05/2015)
      TS_EMIS  = HcoState%TS_EMIS
      DNEWFRAC = TS_EMIS / ( TAU_DAYS * 24.0_hp * 3600.0_hp )
      DOLDFRAC = 1.0_hp - DNEWFRAC
      HNEWFRAC = TS_EMIS / ( TAU_HOURS * 3600.0_hp )
      HOLDFRAC = 1.0_hp - HNEWFRAC

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J,      EMIS_ISOP, EMIS_MBOX, EMIS_APIN, EMIS_BPIN   )
!$OMP+PRIVATE( EMIS_LIMO, EMIS_SABI, EMIS_MYRC, EMIS_CARE, EMIS_OCIM   )
!$OMP+PRIVATE( EMIS_OMON, EMIS_MONO, EMIS_MOHX, EMIS_EOH,  EMIS_FAXX   )
!$OMP+PRIVATE( EMIS_AAXX, EMIS_ACET, EMIS_PRPE, EMIS_C2H4              )
!$OMP+PRIVATE( EMIS_FARN, EMIS_BCAR, EMIS_OSQT, EMIS_ALD2, EMIS_OTHR   )
!$OMP+PRIVATE( MONO_MOL,  TMP,       MBOX_MOL                          )
!$OMP+PRIVATE( DIURFACT,  DIURORVC,  BIOG_ALPH                         )
!$OMP+PRIVATE( BIOG_LIMO, BIOG_ALCO, X, Y, RC                          )

      !-----------------------------------------------------------------
      ! Loop over all grid boxes
      !-----------------------------------------------------------------
      DO J = 1, HcoState%NY
      DO I = 1, HcoState%NX

         ! Zero biogenic species
         EMIS_ISOP  = 0.0_hp
         EMIS_MBOX  = 0.0_hp
         EMIS_APIN  = 0.0_hp
         EMIS_BPIN  = 0.0_hp
         EMIS_LIMO  = 0.0_hp
         EMIS_SABI  = 0.0_hp
         EMIS_MYRC  = 0.0_hp
         EMIS_CARE  = 0.0_hp
         EMIS_OCIM  = 0.0_hp
         EMIS_OMON  = 0.0_hp
         EMIS_MONO  = 0.0_hp
         EMIS_OTHR  = 0.0_hp
         EMIS_ALD2  = 0.0_hp
         EMIS_MOHX  = 0.0_hp
         EMIS_EOH   = 0.0_hp
         EMIS_FAXX  = 0.0_hp
         EMIS_AAXX  = 0.0_hp
         EMIS_ACET  = 0.0_hp
         EMIS_PRPE  = 0.0_hp
         EMIS_C2H4  = 0.0_hp
         EMIS_FARN  = 0.0_hp
         EMIS_BCAR  = 0.0_hp
         EMIS_OSQT  = 0.0_hp

         !-----------------------------------------------------------------
                        ! Get MET variables used in MEGAN core functions
                        !-----------------------------------------------------------------
         LocalMet => NULL()
         LocalMet%D_BTW_M = DBLE( Inst%DAYS_BTW_M )
         LocalMet%ISOLAI    = ExtState%LAI%Arr%Val(I,J)
         LocalMet%MISOLAI    = LocalMet%ISOLAI
         LocalMet%PMISOLAI   = Inst%LAI_PREVDAY(I,J)
         LocalMet%TS         = ExtState%T2M%Arr%Val(I,J)
         LocalMet%SUNCOS     = ExtState%SUNCOS%Arr%Val(I,J)
         LocalMet%WINDSP     = SQRT(ExtState%U10M%Arr%Val(I,J)**2 + 
     &                                ExtState%V10M%Arr%Val(I,J)**2)
         ! Convert Q_DIR and Q_DIFF from (W/m2) to (micromol/m2/s)
         LocalMet%Q_DIR_2  = ExtState%PARDR%Arr%Val(I,J)*WM2_TO_UMOLM2S
         LocalMet%Q_DIFF_2 = ExtState%PARDF%Arr%Val(I,J)*WM2_TO_UMOLM2S
         ! Eventually normalize LAI by PFT (if setting is set
         ! accordingly). ckeller, 7/17/17.
         IF ( Inst%NORMLAI ) THEN
             PFTSUM = SUM( Inst%ARRAY_16(I,J,2:16) )
             IF ( PFTSUM > 0.0_hp ) THEN
              LocalMet%MISOLAI = MIN( LocalMet%MISOLAI/PFTSUM,LAI_MAX )
              LocalMet%PMISOLAI= MIN( LocalMet%PMISOLAI/PFTSUM,LAI_MAX )
             ENDIF
         ENDIF
         LocalMet%PARDR_LASTXDAYS = Inst%PARDR_LASTXDAYS(I,J)
         LocalMet%PARDF_LASTXDAYS = Inst%PARDF_LASTXDAYS(I,J)
         LocalMet%T_LASTXDAYS     = Inst%T_LASTXDAYS(I,J)
         LocalMet%T_LAST24H       = Inst%T_LAST24H(I,J)
         LocalMet%LISOPCO2        = Inst%LISOPCO2
         LocalMet%GLOBCO2         = Inst%GLOBCO2
         ! Get day of year, local-time and latitude
         LocalMet%LAT             = HcoState%Grid%YMID%Val(I,J)
       
         CALL HcoClock_Get (HcoState%Clock,cDOY = DOY,RC=RC)
         CALL HcoClock_GetLocal( HcoState,I,J,cH = LocalHour, RC=RC )
         LocalMet%DOY             = DOY
         LocalMet%LocalHour       = LocalHour
         LocalMet%D2RAD           = Inst%D2RAD
         LocalMet%RAD2D           = Inst%RAD2D
         LocalMet%GWETROOT        = ExtState%GWETROOT%Arr%Val(I,J)

         !--------------------------------------------------------------------
         ! MEGAN Isoprene
         !--------------------------------------------------------------------

         ! Isoprene [kg C/m2/s]
         CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &      LocalMet,  I, J, 'ISOP', Inst%AEF_ISOP(I,J), EMIS_ISOP, RC )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                     'GET_MEGAN_EMISSIONS_ISOP', RC )
            ERR = .TRUE.
            EXIT
         ENDIF

         !FP_ISOP (6/2009)
         EMIS_ISOP = Inst%ISOP_SCALING * EMIS_ISOP

         ! Write isoprene emissions into tracer tendency array
         IF ( Inst%IDTISOP > 0 ) THEN
            FLUXISOP(I,J) = EMIS_ISOP

         ENDIF

         ! Biogenic emissions of SOA and SOA-Precursor from isopene
         IF ( (Inst%IDTISOP>0) .AND. (Inst%IDTSOAP>0) ) THEN
            FLUXSOAP(I,J) = FLUXSOAP(I,J) + EMIS_ISOP * Inst%ISOPTOSOAP
         ENDIF
         IF ( (Inst%IDTISOP>0) .AND. (Inst%IDTSOAS>0) ) THEN
            FLUXSOAS(I,J) = FLUXSOAS(I,J) + EMIS_ISOP * Inst%ISOPTOSOAS
         ENDIF

         !--------------------------------------------------------------------
         ! MEGAN monoterpenes
         !--------------------------------------------------------------------
         ! Updated the calls below for new MEGAN routine (dbm, 12/2012)
         IF ( Inst%ExtNrMono > 0 ) THEN

            ! ---------------------------------------------------
            ! Methyl butenol
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'MBOX',Inst%AEF_MBOX(I,J), EMIS_MBOX, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_MBOX', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXMBOX(I,J) = EMIS_MBOX

            ! ---------------------------------------------------
            ! Alpha Pinene emissions
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'APIN',Inst%AEF_APIN(I,J), EMIS_APIN, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_APIN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXAPIN(I,J) = EMIS_APIN

            ! ---------------------------------------------------
            ! Beta Pinene emissions
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'BPIN',Inst%AEF_BPIN(I,J), EMIS_BPIN, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_BPIN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXBPIN(I,J) = EMIS_BPIN

            ! ---------------------------------------------------
            ! Limonene emissions
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'LIMO',Inst%AEF_LIMO(I,J), EMIS_LIMO, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_LIMO', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXLIMO(I,J) = EMIS_LIMO

            ! ---------------------------------------------------
            ! Sabinene emissions
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'SABI',Inst%AEF_SABI(I,J), EMIS_SABI, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_SABI', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXSABI(I,J) = EMIS_SABI

            ! ---------------------------------------------------
            ! Mycrene emissions
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'MYRC',Inst%AEF_MYRC(I,J), EMIS_MYRC, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_MYRC', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXMYRC(I,J) = EMIS_MYRC

            ! ---------------------------------------------------
            ! 3-Carene emissions
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'CARE',Inst%AEF_CARE(I,J), EMIS_CARE, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_CARE', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXCARE(I,J) = EMIS_CARE

            ! ---------------------------------------------------
            ! Ocimene emissions
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'OCIM',Inst%AEF_OCIM(I,J), EMIS_OCIM, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_OCIM', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXOCIM(I,J) = EMIS_OCIM

            ! ---------------------------------------------------
            ! Other monoterpenes
            ! (added 12/2012; dbm)
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'OMON',Inst%AEF_OMON(I,J), EMIS_OMON, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_OMON', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXOMON(I,J) = EMIS_OMON

            ! ---------------------------------------------------
            ! Total monoterpenes = sum of individual
            ! dbm, now add other lumped monoterpenes (11/2012)
            EMIS_MONO = EMIS_APIN + EMIS_BPIN + EMIS_LIMO + EMIS_SABI +
     &                  EMIS_MYRC + EMIS_CARE + EMIS_OCIM + EMIS_OMON

            ! Add to tracer tendency array [kg C/m2/s]
            FLUXMONO(I,J) = EMIS_MONO
         ELSE
            EMIS_MONO = 0.0_hp
            EMIS_MBOX = 0.0_hp
         ENDIF


!         ! testing only
!         if ( i == ix .AND. j == iy .and. .FALSE. ) then
!                 write(*,*) ' '
!                 write(*,*) 'HEMCO MEGAN (kg/m2/s) @ ', i, j
!                 write(*,*) 'ISOP: ', EMIS_ISOP
!                 write(*,*) 'MBOX: ', EMIS_MBOX
!                 write(*,*) 'MONO: ', EMIS_MONO
!                 write(*,*) 'APIN: ', EMIS_APIN
!                 write(*,*) 'BPIN: ', EMIS_BPIN
!                 write(*,*) 'LIMO: ', EMIS_LIMO
!                 write(*,*) 'SABI: ', EMIS_SABI
!                 write(*,*) 'MYRC: ', EMIS_MYRC
!                 write(*,*) 'CARE: ', EMIS_CARE
!                 write(*,*) 'OCIM: ', EMIS_OCIM
!                 write(*,*) 'OMON: ', EMIS_OMON
!         endif

         !--------------------------------------------------------------------
         ! MEGAN Acetaldehyde
         !--------------------------------------------------------------------
         IF ( Inst%IDTALD2 > 0 ) THEN
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'ALD2',Inst%AEF_ALD2(I,J), EMIS_ALD2, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_ALD2', RC )
               ERR = .TRUE.
               EXIT
            ENDIF

            ! Add to tracer tendency array [kg C/m2/s]
            FLUXALD2(I,J) = EMIS_ALD2
         ENDIF

         ! ---------------------------------------------------
         ! MEGAN Ethanol
         IF ( Inst%IDTEOH > 0 ) THEN
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &           LocalMet, I, J, 'EOH',Inst%AEF_EOH(I,J), EMIS_EOH, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_EOH', RC )
               ERR = .TRUE.
               EXIT
            ENDIF

            ! Add to tracer tendency array [kg C/m2/s]
            FLUXEOH(I,J) = EMIS_EOH
         ENDIF

         !--------------------------------------------------------------------
         ! Other MEGAN biogenics
         ! Calls included here for future incorporation or
         ! specialized simulations. (dbm, 12/2012)
         !--------------------------------------------------------------------

         ! ---------------------------------------------------
         ! MEGAN Methanol
         ! Uncomment this to compute biogenic methanol emissions
!         CALL GET_MEGAN_EMISSIONS( am_I_Root, HcoState, ExtState,
!     &                             Inst, I, J, 'MOHX', EMIS_MOHX, RC )
!         IF ( RC /= HCO_SUCCESS ) THEN
!            CALL HCO_ERROR( HcoState%Config%Err,
!     &                        'GET_MEGAN_EMISSIONS_MOHX', RC )
!            ERR = .TRUE.
!            EXIT
!         ENDIF
!         FLUXMOHX(I,J) = EMIS_MOHX

         ! ---------------------------------------------------
         ! MEGAN Formic acid
         ! Uncomment this to compute biogenic formic acid emissions
!         CALL GET_MEGAN_EMISSIONS( am_I_Root, HcoState, ExtState,
!     &                             Inst, I, J, 'FAXX', EMIS_FAXX, RC )
!         IF ( RC /= HCO_SUCCESS ) THEN
!            CALL HCO_ERROR( HcoState%Config%Err,
!     &                        'GET_MEGAN_EMISSIONS_FAXX', RC )
!            ERR = .TRUE.
!            EXIT
!         ENDIF
!         FLUXFAXX(I,J) = EMIS_FAXX

         ! ---------------------------------------------------
         ! MEGAN Acetic acid
         ! Uncomment this to compute biogenic acetic acid emissions
!         CALL GET_MEGAN_EMISSIONS( am_I_Root, HcoState, ExtState,
!     &                             Inst, I, J, 'AAXX', EMIS_AAXX, RC )
!         IF ( RC /= HCO_SUCCESS ) THEN
!            CALL HCO_ERROR( HcoState%Config%Err,
!     &                        'GET_MEGAN_EMISSIONS_AAXX', RC )
!            ERR = .TRUE.
!            EXIT
!         ENDIF
!         FLUXAAXX(I,J) = EMIS_AAXX

         !--------------------------------------------------------------------
         ! Biogenic acetone emissions
         !--------------------------------------------------------------------
         IF ( Inst%IDTACET > 0 ) THEN

            !-----------------------------------------------------------------
            ! (1) BIOGENIC EMISSIONS OF ACETONE FROM MONOTERPENES
            ! Monoterpenes has same # molecules/kg of carbon as isoprene
            ! The yield for monoterpenes is .12 mol/mol from Reisell et.al.
            ! 1999 (this does not includes direct acetone emissions)

            ! Convert [kg C/m2/s] to [kg MONOTERPENE/m2/s]
            ! There are 10 C atoms per molecule of MONOTERPENE
            MONO_MOL = EMIS_MONO / 10.0_hp

            ! Apply yield from monoterpenes to get [kg ACET/m2/s]
            TMP      = MONO_MOL * YIELD_MO

            ! Convert acetone emissions back into [kg C/m2/s]
            TMP      = TMP * 3.0_hp

            ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
            TMP      = TMP * MONO_SCALE

            ! Add to total biogenic acetone emissions
            FLUXACETmo(I,J) = TMP

            !-----------------------------------------------------------------
            ! (2) BIOGENIC ACETONE FROM METHYL BUTENOL -- NORTH AMERICA
            !
            ! Methyl Butenol (a.k.a. MBO) produces acetone with a molar yield
            ! of 0.6 [Alvarado (1999)].  The biogenic source of MBO is thought
            ! to be restricted to North America.  According to Guenther (1999)
            ! North america emits 3.2Tg-C of MBO, producing 1.15 Tg-C of
            ! Acetone in North America.
            !=================================================================
            TMP = 0.0_hp

            ! Lon and lat of grid box (I,J) in degrees
            X = HcoState%Grid%XMID%Val( I, J )
            IF ( X >= 180.0_hp ) X = X - 360.0_hp
            Y = HcoState%Grid%YMID%Val( I, J )

            ! Methyl butenol is emitted only in North America, where
            ! ( -167.5 <= lon <= -52.5 ) and ( 16.0 <= lat <= 72.0 )
            IF ( ( X >= -167.5_hp .and. X <= -52.5_hp ) .AND.
     &           ( Y >=   16.0_hp .and. Y <=  72.0_hp ) ) THEN

               ! Convert from [kg C/m2/s] to [kg MBO/m2/s]
               ! There are 5 C atoms per molecule MBO
               MBOX_MOL = EMIS_MBOX / 5.0_hp

               ! Apply yield from MBO to get [kg ACET/m2/s]
               TMP      = MBOX_MOL * MB_SCALE1

               ! Convert from [kg ACET/m2/s] to [kg C/m2/s]
               ! There are 3 C atoms per acetone molecule
               TMP      = TMP * 3.0_hp

               ! Scale to a posteriori source from Jacob et al 2001 (bdf,
               ! 9/5/01)
               FLUXACETmb(I,J) = TMP * MB_SCALE2

            ENDIF

            !-----------------------------------------------------------------
            ! (3) BIOGENIC ACETONE -- DIRECT EMISSION
            ! evf, removed obsolete code, replaced with MEGAN acetone
            ! emissions (5/25/2011) Direct Emission now includes emission
            ! from grasses and emission from dry leaf matter
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'ACET',Inst%AEF_ACET(I,J), EMIS_ACET, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_ACET', RC )
               ERR = .TRUE.
               EXIT
            ENDIF

            ! Write to array
            FLUXACETbg(I,J) = EMIS_ACET
         ENDIF

         !--------------------------------------------------------------------
         !  Biogenic source of CO -- from oxidation of METHANOL and MONOTERPENES
         !
         !  CO from METHANOL oxidation -- scaled from ISOPRENE (bnd, 1/2/01)
         !
         !    We need to scale the Isoprene flux to get the CH3OH (methanol)
         !    flux. Currently, the annual isoprene flux in GEOS-CHEM is
         !    ~ 397 Tg C.
         !
         !    Daniel Jacob recommends a flux of 100 Tg/yr CO from CH3OH
         !    oxidation based on Singh et al. 2000 [JGR 105, 3795-3805] who
         !    estimate a global methanol source of 122 Tg yr-1, of which most
         !    (75 Tg yr-1) is "primary biogenic". He also recommends for now
         !    that the CO flux from CH3OH oxidation be scaled to monthly mean
         !    isoprene flux.
         !
         !    To get CO from METHANOL oxidation, we must therefore multiply
         !    the ISOPRENE flux by the following scale factor:
         !    ( 100 Tg CO from CH3OH Oxidation / 397 Tg C from Isoprene Flux ) *
         !    (  12 g C/mole                   / 28 g CO/mole                )
         !
         !  CO from MONOTERPENE oxidation (bnd, bmy, 1/2/01)
         !
         !    Assume the production of CO from monoterpenes is instantaneous
         !    even though the lifetime of intermediate species may be on the
         !    order of hours or days. This assumption will likely cause CO from
         !    monoterpene oxidation to be too high in the box in which the
         !    monoterpene is emitted.
         !
         !    The CO yield here is taken from:
         !
         !    Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
         !      "The ultimate yield of CO from the tropospheric oxidation of
         !       terpenes (including both O3 and OH reactions) was estimated
         !       to be 20% on the carbon number basis." They studied
         !       ALPHA- & BETA-pinene.
         !
         !    Vinckier et al. Fresenius Env. Bull., Vol. 7, p.361-368 (1998)
         !      "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
         !
         !    EMIS_MONO is now kgC, conversion factor then becomes:
         !    28g CO / 60g C = 46.666%
         !    (ckeller, 4/15/14).
         !--------------------------------------------------------------------

         ! NOTE: Don't save into EMISRR for the tagged CO simulation.
         ! Also for tagged CO we don't use monoterpenes  ??????
         ! (jaf, mak, bmy, 2/14/08)
         IF ( (Inst%IDTCO>0) .AND. (Inst%ExtNrMono>0) ) THEN
!            FLUX = EMMO * 0.2d0
            FLUXCO(I,J) = EMIS_MONO * 0.4666_hp

!            ! testing only
!            if ( i==ix .and. j==iy .and. .FALSE. ) then
!               write(*,*) ' '
!               write(*,*) 'HEMCO MEGAN @ ', ix, iy
!               write(*,*) 'CO (kg/m2/s): ', FLUXCO(I,J)
!            endif

         ENDIF

         !--------------------------------------------------------------------
         ! Biogenic emissions of SOA and SOA-Precursor from monoterpenes
         !--------------------------------------------------------------------
         IF ( (Inst%IDTSOAP>0) .AND. (Inst%ExtNrMono>0) ) THEN
            FLUXSOAP(I,J) = FLUXSOAP(I,J) + EMIS_MONO * Inst%MONOTOSOAP
         ENDIF
         IF ( (Inst%IDTSOAS>0) .AND. (Inst%ExtNrMono>0) ) THEN
            FLUXSOAS(I,J) = FLUXSOAS(I,J) + EMIS_MONO * Inst%MONOTOSOAS
         ENDIF

         !--------------------------------------------------------------------
         ! Biogenic emissions of PRPE
         !
         ! Now uses MEGAN2.1 (dbm, 12/2012)
         !--------------------------------------------------------------------
         IF ( Inst%IDTPRPE > 0 ) THEN
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'PRPE',Inst%AEF_PRPE(I,J), EMIS_PRPE, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_PRPE', RC )
               ERR = .TRUE.
               EXIT
            ENDIF

            ! Add to tracer tendency array [kg C/m2/s]
            FLUXPRPE(I,J) = EMIS_PRPE

!           ! testing only
!           if ( i==ix .and. j==iy .and. .FALSE. ) then
!              write(*,*) ' '
!              write(*,*) 'HEMCO MEGAN @ ', ix, iy
!              write(*,*) 'PRPE (kg/m2/s): ', FLUXPRPE(I,J)
!           endif

         ENDIF

         !--------------------------------------------------------------------
         ! Biogenic emissions of ethene (C2H4)
         !
         ! Now uses MEGAN2.1 (dbm, 12/2012)
         !--------------------------------------------------------------------
         IF ( Inst%IDTC2H4 > 0 ) THEN
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &         LocalMet, I, J, 'C2H4',Inst%AEF_C2H4(I,J), EMIS_C2H4, RC)
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_C2H4', RC )
               ERR = .TRUE.
               EXIT
            ENDIF

            ! Add to tracer tendency array [kg C/m2/s]
            FLUXC2H4(I,J) = EMIS_C2H4
         ENDIF

         ! Use CALTECH algorithm if MEGAN Monoterpene is turned on
         IF ( Inst%ExtNrMono > 0 ) THEN

            ! Error trap: check if SZAFACT is indeed defined
            IF ( .NOT. ASSOCIATED(ExtState%SZAFACT%Arr%Val) ) THEN
                  CALL HCO_ERROR( HcoState%Config%Err,
     &                           'SZAFACT not set - needed for SOA', RC)
                  ERR = .TRUE.
                  EXIT
            ENDIF

            ! Impose a diurnal variation on NVOC during the day
            DIURFACT = ExtState%SZAFACT%Arr%Val(I,J)
            DIURORVC = DIURFACT * Inst%GEIA_ORVC(I,J)

            ! ----------------------------------------------------------------
            ! The new MEGAN implementation has speciated information
            ! (hotp 3/7/10)
            !! For GCAP Meteorology year 2000 in Tg/yr:
            !! ------------------------------
            !! HC Class  New MEGAN  Old MEGAN
            !! --------  ---------  ---------
            !! ALPH        84         92
            !! LIMO        10         27
            !! TERP         3.2        3.5
            !! ALCO        47         38
            !! SESQ        15         15    (no change for SESQ)
            !!           -----      -----
            !! TOTAL      159        176
            ! ------------------------------
            ! updated again to include MEGAN sesquiterpenes
            ! and other monoterpenes (hotp 7/28/10)
            ! as of 7/28/10 for year 2000 GEOS4 2x2.5 in Tg/yr:
            ! ------------------------------
            ! HC Class  New MEGAN  Old MEGAN
            ! --------  ---------  ---------
            ! MTPA        73         86
            ! LIMO        10         25
            ! MTPO        33          3.2+38
            ! SESQ        13         15
            !           -----      --------
            ! TOTAL      129        167.2
            ! see Pye et al. 2010 ACP
            !
            ! MEGAN emission update for SESQ and other monoterps
            ! (hotp 3/10/10)
            ! terpenoid ketones no longer lumped into ALPH
            ! terpenoid alcohols no longer lumped in ALCO
            ! Other monoterpenes (OMON from MEGAN) all placed in TERP
            ! All sesq lumped together into SESQ
            !
            ! new mtp lumping (as of 5/20/10) hotp
            ! GEOS-Chem   MEGAN
            ! =========   ==========================================
            ! MTPA        a-pinene   (APIN), b-pinene  (BPIN),
            !             sabinene   (SABI), carene    (CARE)
            ! LIMO        limonene   (LIMO)
            ! MTPO        myrcene    (MYRC), ocimene   (OCIM),
            !             other mono (OMON)
            ! SESQ        farnesene  (FARN), b-caryoph (BCAR),
            !             other sesq (OSQT)
            ! =========   ==========================================

            !--------------------------------------------------------------
            ! MEGAN MTPA
            !--------------------------------------------------------------
            ! MTPA=a-,b-pinene,sabinene,carene (hotp 5/20/10)
            IF ( Inst%IDTMTPA > 0 ) THEN
               FLUXMTPA(I,J) = ( EMIS_APIN + EMIS_BPIN +
     &                           EMIS_SABI + EMIS_CARE ) * FC1
            ENDIF

            !--------------------------------------------------------------
            ! MEGAN Limonene
            !--------------------------------------------------------------
            ! [kg C/m2/s]
            IF ( Inst%IDTLIMO > 0 ) THEN
               FLUXLIMO(I,J) = EMIS_LIMO * FC1
            ENDIF

            !--------------------------------------------------------------
            ! MEGAN MTPO
            !--------------------------------------------------------------
            ! MTPO is all other monoterpenes (MEGAN categories:
            ! myrcene, ocimene, OMON) (hotp 5/20/10)
            ! All other monoterpenes (mostly camphene, linalool,
            ! terpinolene, terpinolene, phellandrene) (hotp 3/10/10)
            ! 14-18% of OMTP is terpinene and terpinolene
            IF ( Inst%IDTMTPO > 0 ) THEN
               FLUXMTPO(I,J) = ( EMIS_MYRC + EMIS_OCIM +
     &                           EMIS_OMON ) * FC1
            ENDIF

            !--------------------------------------------------------------
            ! MEGAN sesquiterpenes
            !--------------------------------------------------------------

            ! ---------------------------------------------------
            ! a-Farnesene
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &        LocalMet, I, J, 'FARN',Inst%AEF_FARN(I,J), EMIS_FARN, RC )
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_FARN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXFARN(I,J) = EMIS_FARN

            ! ---------------------------------------------------
            ! b_Caryophyllene
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &        LocalMet, I, J, 'BCAR',Inst%AEF_BCAR(I,J), EMIS_BCAR, RC )
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_BCAR', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXBCAR(I,J) = EMIS_BCAR

            ! ---------------------------------------------------
            ! Other sesquiterpene
            CALL GET_MEGAN_EMISSIONS( !am_I_Root, HcoState, ExtState,
     &        LocalMet, I, J, 'OSQT',Inst%AEF_OSQT(I,J), EMIS_OSQT, RC )
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( HcoState%Config%Err,
     &                        'GET_MEGAN_EMISSIONS_OSQT', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
            FLUXOSQT(I,J) = EMIS_OSQT

            ! ---------------------------------------------------
            ! Total sesquiterpenes from MEGAN (hotp 3/10/10)
            IF ( Inst%IDTSESQ > 0 ) THEN
               FLUXSESQ(I,J) = ( EMIS_FARN + EMIS_BCAR +
     &                           EMIS_OSQT ) * FC3
            ENDIF

            EMIS_OTHR = EMIS_FARN + EMIS_BCAR + EMIS_OSQT

         ELSE

            ! If MEGAN Monoterpene is turned off:
            ! NOTE: In the original GEOS-Chem code, many of these fluxes
            ! were derived from the monoterpenese fluxes. These fluxes
            ! are only calculated if the MEGAN_MONO option is turned on
            ! and their contributions are hence skipped in the following!
            ! (ckeller, 05/19/14)

            ! MTPA
            ! old ALPH is now MTPA (hotp 5/20/10)
            ! ALPHA-PINENE (0.35)
            ! BETA-PINENE lumped with ALPHA-PINENE (0.23)
            ! SABINENE    lumped with ALPHA-PINENE (0.05)
            ! D3-CARENE   lumped with ALPHA-PINENE (0.04)

            ! TERPENOID KETONE is lumped with SABINENE
            ! Then SABINENE    is lumped with ALPHA-PINENE
            IF ( Inst%IDTMTPA > 0 ) THEN
               FLUXMTPA(I,J) = DIURORVC * FC4 * 0.04_hp ! using campher
            ENDIF

            ! LIMONENE not defined (derived from TERP)
            IF ( Inst%IDTLIMO > 0 ) THEN
               FLUXLIMO(I,J) = 0.0_hp
            ENDIF

            ! MTPO
            IF ( Inst%IDTMTPO > 0 ) THEN
               FLUXMTPO(I,J) = DIURORVC * FC2 * 0.09_hp ! using LINALOOL
            ENDIF

            ! IDTSESQ
            IF ( Inst%IDTSESQ > 0 ) THEN
               FLUXSESQ(I,J) = DIURORVC * FC3 * 0.05_hp
            ENDIF

            EMIS_OTHR = 0.0_hp

         ENDIF

         !--------------------------------------------------------------------
         ! Biogenic emissions of SOA and SOA-Precursor from Other terpenes
         !--------------------------------------------------------------------
         IF ( (Inst%IDTSOAP>0) .AND. (Inst%ExtNrMono>0) ) THEN
            FLUXSOAP(I,J) = FLUXSOAP(I,J) + EMIS_OTHR * Inst%OTHRTOSOAP
         ENDIF
         IF ( (Inst%IDTSOAS>0) .AND. (Inst%ExtNrMono>0) ) THEN
            FLUXSOAS(I,J) = FLUXSOAS(I,J) + EMIS_OTHR * Inst%OTHRTOSOAS
         ENDIF

         !-----------------------------------------------------------------
         ! Update historical temperature / radiation values
         !-----------------------------------------------------------------
         ! Do this now on every time step. The arrays are simply the
         ! the running means over the intentend time window (24 hours
         ! and NUM_DAYS, respectively). This hugely faciliates warm
         ! model restarts, irrespective of simulation start/end dates
         ! and times. It also makes sure that all environmental
         ! variables are incorporated into the time averages (e.g. if
         ! emission time step is less than 60 minutes, all values will
         ! be used to calculate the daily mean).
         ! (ckeller, 11/05/2015)

         ! Updated LAI of last 24 hours
         Inst%LAI_PREVDAY(I,J) = ( HOLDFRAC * Inst%LAI_PREVDAY(I,J) )
     &                    + ( HNEWFRAC * ExtState%LAI%Arr%Val(I,J) )

         ! Updated temperature of last 24 hours
         Inst%T_LAST24H(I,J) = ( HOLDFRAC * Inst%T_LAST24H(I,J) )
     &                  + ( HNEWFRAC * ExtState%T2M%Arr%Val(I,J) )

         ! Updated temperature of last NUM_DAYS
         Inst%T_LASTXDAYS(I,J) = ( DOLDFRAC * Inst%T_LASTXDAYS(I,J) )
     &                    + ( DNEWFRAC * ExtState%T2M%Arr%Val(I,J) )

         ! Updated direct radiation of last NUM_DAYS
         Inst%PARDR_LASTXDAYS(I,J) =
     &                        ( DOLDFRAC * Inst%PARDR_LASTXDAYS(I,J)   )
     &                      + ( DNEWFRAC * ExtState%PARDR%Arr%Val(I,J) )

         ! Updated diffuse radiation of last NUM_DAYS
         Inst%PARDF_LASTXDAYS(I,J) =
     &                        ( DOLDFRAC * Inst%PARDF_LASTXDAYS(I,J)   )
     &                      + ( DNEWFRAC * ExtState%PARDF%Arr%Val(I,J) )

      ENDDO !I
      ENDDO !J
!$OMP END PARALLEL DO
      IF ( ERR ) THEN
         RC = HCO_FAIL
         RETURN
      ENDIF

      !=================================================================
      ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS
      !=================================================================

!      !!! testing only
!      if(am_I_Root) then
!       write(*,*) 'ISOP emissions instance ',ExtState%Megan
!       write(*,*) 'A MGN ',ExtState%Megan,Inst%IDTISOP,SUM(FLUXISOP)
!       write(*,*) 'B MGN ',ExtState%Megan,SUM(Inst%LAI_PREVDAY)
!       write(*,*) 'C MGN ',ExtState%Megan,SUM(Inst%T_LAST24H)
!       write(*,*) 'D MGN ',ExtState%Megan,SUM(Inst%T_LASTXDAYS)
!       write(*,*) 'E MGN ',ExtState%Megan,SUM(Inst%PARDR_LASTXDAYS)
!       write(*,*) 'F MGN ',ExtState%Megan,SUM(Inst%PARDF_LASTXDAYS)
!       write(*,*) 'G MGN ',ExtState%Megan,SUM(ExtState%T2M%Arr%Val)
!       write(*,*) 'H MGN ',ExtState%Megan,SUM(ExtState%SUNCOS%Arr%Val)
!       write(*,*) 'I MGN ',ExtState%Megan,SUM(ExtState%PARDR%Arr%Val)
!       write(*,*) 'J MGN ',ExtState%Megan,SUM(ExtState%PARDF%Arr%Val)
!       write(*,*) 'K MGN ',ExtState%Megan,SUM(ExtState%LAI%Arr%Val)
!       write(*,*) 'L MGN ',ExtState%Megan,SUM(ExtState%GWETROOT%Arr%Val)
!       write(*,*) 'M MGN ',ExtState%Megan,SUM(ExtState%SZAFACT%Arr%Val)
!      endif

      ! ----------------------------------------------------------------
      ! ISOPRENE
      IF ( Inst%IDTISOP > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXISOP,
     &                     Inst%IDTISOP, RC, ExtNr=Inst%ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                     'HCO_EmisAdd error: FLUXISOP', RC )
            RETURN
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ACETALDEHYDE
      IF ( Inst%IDTALD2 > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXALD2, Inst%IDTALD2,
     &                     RC,        ExtNr=Inst%ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXALD2', RC )
            RETURN
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ETHANOL
      IF ( Inst%IDTEOH > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXEOH, Inst%IDTEOH,
     &                     RC,        ExtNr=Inst%ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXEOH', RC )
            RETURN
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ACETONE
      IF ( Inst%IDTACET > 0 ) THEN

         ! Eventually add individual diagnostics. These are assumed to
         ! have names MEGAN_ACET_MONO, MEGAN_ACET_MBO, and
         ! MEGAN_ACET_DIRECT
         DiagnName =  'MEGAN_ACET_MONO'
         CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr,
     &      cName=TRIM(DiagnName), Array2D=FLUXACETmo, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         DiagnName =  'MEGAN_ACET_MBO'
         CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr,
     &      cName=TRIM(DiagnName), Array2D=FLUXACETmb, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         DiagnName =  'MEGAN_ACET_DIRECT'
         CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr,
     &      cName=TRIM(DiagnName), Array2D=FLUXACETbg, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! Add flux to emission array
         FLUXACETbg = FLUXACETbg + FLUXACETmo + FLUXACETmb
         CALL HCO_EmisAdd( HcoState, FLUXACETbg,
     &                     Inst%IDTACET, RC, ExtNr=Inst%ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXACETbg', RC )
            RETURN
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! SOA-Precursor (SOAP) and SOA-Simplified (SOAS)
      IF (  Inst%IDTSOAP > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXSOAP, Inst%IDTSOAP,
     &                     RC,        ExtNr=Inst%ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXSOAP', RC )
            RETURN
         ENDIF

      ENDIF

      IF (  Inst%IDTSOAS > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXSOAS, Inst%IDTSOAS,
     &                     RC,        ExtNr=Inst%ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXSOAS', RC )
            RETURN
         ENDIF

         ! Also archive the BIOGENIC_SOAS diagnostic
         DiagnName =  'BIOGENIC_SOAS'
         CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr,
     &      cName=TRIM(DiagnName), Array2D=FLUXSOAS, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN
         Arr2D => NULL()
      ENDIF

      ! ----------------------------------------------------------------
      ! CO
      IF ( ( Inst%ExtNrMono > 0 ) .AND.
     &     ( Inst%IDTCO     > 0 )      ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXCO, Inst%IDTCO,
     &                     RC,        ExtNr=Inst%ExtNrMono )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXCO', RC )
            RETURN
         ENDIF

         ! Also archive the BIOGENIC_CO diagnostic (bmy, 4/29/15)
         DiagnName =  'BIOGENIC_CO'
         CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXCO, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN
         Arr2D => NULL()
      ENDIF

      ! ----------------------------------------------------------------
      ! ALKENES
      IF ( Inst%IDTPRPE > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXPRPE, Inst%IDTPRPE,
     &                     RC,        ExtNr=Inst%ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXPRPE', RC )
            RETURN
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ETHENE
      IF ( Inst%IDTC2H4 > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXC2H4, INst%IDTC2H4,
     &                     RC,        ExtNr=Inst%ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXC2H4', RC )
            RETURN
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! MTPA
      IF ( ( Inst%ExtNrMono > 0 ) .AND.
     &     ( Inst%IDTMTPA   > 0 )      ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXMTPA, Inst%IDTMTPA,
     &                     RC,        ExtNr=Inst%ExtNrMono )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXMTPA', RC )
            RETURN
         ENDIF

      ENDIF

      ! ----------------------------------------------------------------
      ! MTPO
      IF ( ( Inst%ExtNrMono > 0 ) .AND.
     &     ( Inst%IDTMTPO   > 0 )      ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXMTPO, Inst%IDTMTPO,
     &                     RC,        ExtNr=Inst%ExtNrMono )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXMTPO', RC )
            RETURN
         ENDIF

      ENDIF

      ! ----------------------------------------------------------------
      ! LIMONENE
      IF ( ( Inst%ExtNrMono > 0 ) .AND.
     &     ( Inst%IDTLIMO   > 0 )      ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXLIMO, Inst%IDTLIMO,
     &                     RC,        ExtNr=Inst%ExtNrMono )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXLIMO', RC )
            RETURN
         ENDIF

      ENDIF

      ! ----------------------------------------------------------------
      ! SESQ
      IF ( ( Inst%ExtNrMono > 0 ) .AND.
     &     ( Inst%IDTSESQ   > 0 )      ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXSESQ, Inst%IDTSESQ,
     &                     RC,        ExtNr=Inst%ExtNrMono )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'HCO_EmisAdd error: FLUXSESQ', RC )
            RETURN
         ENDIF

      ENDIF

      ! Restore MEGAN monoterpenes and sesquiterpenes to diagnostics
      ! Add individual species to diagnostics rather than lumped species,
      ! following method used in pre-HEMCO diagnostics (mps, 2/18/15)
      IF ( Inst%ExtNrMono > 0 ) THEN

         ! -------------------------------------------------------------
         ! Methyl butenol
         DiagnName =  'BIOGENIC_MBOX'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &     cName=TRIM(DiagnName), Array2D=FLUXMBOX, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! Alpha Pinene
         DiagnName =  'BIOGENIC_APIN'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXAPIN, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! Beta Pinene
         DiagnName =  'BIOGENIC_BPIN'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXBPIN, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! Limonene
         DiagnName =  'BIOGENIC_LIMO'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXLIMO, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! Sabinene
         DiagnName =  'BIOGENIC_SABI'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXSABI, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! Mycrene
         DiagnName =  'BIOGENIC_MYRC'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXMYRC, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! 3-Carene
         DiagnName =  'BIOGENIC_CARE'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXCARE, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! Ocimene
         DiagnName =  'BIOGENIC_OCIM'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXOCIM, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! Other monoterpenes
         DiagnName =  'BIOGENIC_OMON'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXOMON, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! Total monoterpenes
         DiagnName =  'BIOGENIC_MONX'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXMONO, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! a-Farnesene
         DiagnName =  'BIOGENIC_FARN'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXFARN, RC=RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! b_Caryophyllene
         DiagnName =  'BIOGENIC_BCAR'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXBCAR, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! -------------------------------------------------------------
         ! Other sesquiterpenes
         DiagnName =  'BIOGENIC_OSQT'
         CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXOSQT, RC=RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

      ENDIF

      ! ----------------------------------------------------------------
      ! Other MEGAN biogenics
      !
      ! Calls included here for future incorporation or
      ! specialized simulations. (dbm, 12/2012)
      !
      ! These diagnostics will be zero because FLUXMOHX, FLUXFAXX,
      !  and FLUXAAXX are commented out above (mps, 2/19/15)
      ! ----------------------------------------------------------------
      ! Methanol
      DiagnName =  'BIOGENIC_MOHX'
      CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNr,
     &   cName=TRIM(DiagnName), Array2D=FLUXMOHX, RC=RC)
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! ----------------------------------------------------------------
      ! Formic acid
      DiagnName =  'BIOGENIC_FAXX'
      CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNr,
     &   cName=TRIM(DiagnName), Array2D=FLUXFAXX, RC=RC)
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! ----------------------------------------------------------------
      ! Acetic acid
      DiagnName =  'BIOGENIC_AAXX'
      CALL Diagn_Update(  HcoState, ExtNr=Inst%ExtNr,
     &   cName=TRIM(DiagnName), Array2D=FLUXAAXX, RC=RC)
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! ----------------------------------------------------------------
      ! Eventually copy internal values to ESMF internal state object
      ! ----------------------------------------------------------------

      ! LAI_PREVDAY
      CALL HCO_RestartWrite(  HcoState,
     &                       'LAI_PREVDAY', Inst%LAI_PREVDAY, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! T_LAST24H
      CALL HCO_RestartWrite(     HcoState,
     &                       'T_PREVDAY',  Inst%T_LAST24H, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! T_LASTXDAYS
      CALL HCO_RestartWrite(     HcoState,
     &                       'T_DAVG',     Inst%T_LASTXDAYS, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! PARDR_LASTXDAYS
      CALL HCO_RestartWrite(     HcoState,
     &                       'PARDR_DAVG', Inst%PARDR_LASTXDAYS, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! PARDF_LASTXDAYS
      CALL HCO_RestartWrite(     HcoState,
     &                       'PARDF_DAVG', Inst%PARDF_LASTXDAYS, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! ALL DONE!
      !=================================================================

!      ! testing only
!      write(*,*) ''
!      write(*,*) 'MEGAN done!'
!      write(*,*) 'total LAI : ', SUM(ExtState%GC_LAI%Arr%Val)
!      write(*,*) 'total ISOP: ', SUM(FLUXISOP)
!      write(*,*) 'total ACET: ', SUM(FLUXACETmo)+SUM(FLUXACETmb)
!     &                          +SUM(FLUXACETbg)
!      write(*,*) 'total CO  : ', SUM(FLUXCO  )
!      write(*,*) 'total PRPE: ', SUM(FLUXPRPE)
!      write(*,*) 'total OCPI: ', SUM(FLUXOCPI)
!      write(*,*) ''

      ! Cleanup
      Inst => NULL()

      ! Reset first-time flag
      FIRST = .FALSE.

      ! Return w/ success
      CALL HCO_LEAVE( HcoState%Config%Err,RC )

      END SUBROUTINE HCOX_Megan_Run
!EOC

!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fill_Restart_Vars
!
! !DESCRIPTION: Subroutine FILL\_RESTART\_VARS fills the megan restart
!  variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FILL_RESTART_VARS( HcoState,
     &                              ExtState, Inst, RC )
!
! !USES:
!
      USE HCO_Restart_Mod,    ONLY : HCO_RestartGet
!
! !INPUT PARAMETERS:
!
      !LOGICAL,         INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State), POINTER        :: HcoState
      TYPE(Ext_State), POINTER        :: ExtState
      TYPE(MyInst),    POINTER        :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  19 Dec 2014 - C. Keller    - Initial version
!  09 Mar 2015 - C. Keller    - Now use HCO_RestartGet to fill restarts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I
      LOGICAL           :: FND

      !=================================================================
      ! FILL_RESTART_VARS begins here!
      !=================================================================

      ! On first call, see if there are restart variables available for
      ! historic values of temperature and irradiation. If so, fill the
      ! daily arrays with those average values to start with. If not,
      ! simply use the current values.
      ! The same procedure is applied to the 15 day average values. For
      ! temperature, we also restart previous day temperature so that
      ! we can determine a trend in the daily average temperatures
      ! (between the long term daily averages and the most recent daily
      ! average).
      ! We now first call HCO_CopyFromInt_ESMF to check if the variable
      ! is found in the (ESMF) internal state object. In an non-ESMF
      ! environment, this call will always return FND=.FALSE.
      ! (ckeller, 3/9/15).

      !-----------------------------------------------------------------
      ! LAI
      !-----------------------------------------------------------------
      CALL HCO_RestartGet(      HcoState, 'LAI_PREVDAY',
     &                     Inst%LAI_PREVDAY, RC, FILLED=FND )

      ! Default value
      IF ( .NOT. FND ) THEN
         Inst%LAI_PREVDAY = ExtState%LAI%Arr%Val
      ENDIF

      !-----------------------------------------------------------------
      ! DIFFUSE RADIATION ( PARDF )
      !-----------------------------------------------------------------

      ! Temperature over last 24 hours
      CALL HCO_RestartGet(  HcoState, 'T_PREVDAY',
     &                     Inst%T_LAST24H, RC, DefVal = 288.15 )

      ! Temperature over last X days
      CALL HCO_RestartGet(    HcoState, 'T_DAVG',
     &                     Inst%T_LASTXDAYS, RC, DefVal = 288.15 )

      ! Direct radiation (PARDR) over last X days
      CALL HCO_RestartGet(        HcoState, 'PARDR_DAVG',
     &                     Inst%PARDR_LASTXDAYS, RC, DefVal = 30.0 )

      ! Diffuse radiation (PARDF) over last X days
      CALL HCO_RestartGet(        HcoState, 'PARDF_DAVG',
     &                     Inst%PARDF_LASTXDAYS, RC, DefVal = 48.0 )

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE FILL_RESTART_VARS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_AEF
!
! !DESCRIPTION: Subroutine CALC\_AEF reads Emission Factors for all
!  biogenic VOC species from disk.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_AEF( HcoState, ExtState, Inst, RC )
!
! !USES:
!
      ! References to F90 modules
      USE HCO_EMISLIST_MOD,    ONLY : HCO_GetPtr
      USE HCO_CALC_MOD,        ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
      !LOGICAL,         INTENT(IN)    :: am_I_Root
      TYPE(Ext_State), POINTER       :: ExtState
      TYPE(HCO_State), POINTER       :: HcoState
      TYPE(MyInst),    POINTER       :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,         INTENT(INOUT) :: RC
!
! !REMARKS:
!  Reference: (5 ) Guenther et al, 2004
!
! !REVISION HISTORY:
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard
!        code by May Fu (11/2004)
!  (2 ) AEF detailed in the latest MEGAN User's Guide (tmf, 11/19/04)
!  (3 ) Bug fix (tmf, 11/30/04)
!  (4 ) Now reads 1x1 files and regrids to current resolution (bmy, 10/24/05)
!  (5 ) Uses new v2.1 emission factors maps for isoprene, MBO and 7 monoterpene
!     species, download in 2009. (mpb,2009)
!  (6 ) Now use 2.1 emission factors for isoprene, MBO, and 7 monoterpenes. EFs for
!       other compounds are computed by reading in the PFT fractions and multiplying
!       the fractions by corresponding EF values. (dbm, 11/2012)
!  (7 ) Also, now read in the EF values already gridded to model resolution.
!       (dbm, 11/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  10 Jun 2015 - M. Sulprizio- Bug fix: Now convert sesquiterpenes (FARN, BCAR,
!                              OSQT) from ug compound/m2/h to kg C/m2/s
!  17 Jul 2017 - C. Keller   - PFT values (ARRAY_16) are now stored in Inst object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                 :: I, J, P, ARR_IND
      REAL(hp)                :: FACTOR, SPECIES2CARBON
      !REAL(hp)                :: ARRAY_16(HcoState%NX,HcoState%NY,16)
      REAL(hp)                :: PFT_EF_OMON(15), PFT_EF_MOHX(15)
      REAL(hp)                :: PFT_EF_ACET(15), PFT_EF_BIDR(15)
      REAL(hp)                :: PFT_EF_STRS(15), PFT_EF_OTHR(15)
! --->
! dbm, compute EF maps for a-pinene and myrcene as well since there seems to
! be an issue with the EF maps for these species provided on the MEGAN
! data portal
      REAL(hp)                :: PFT_EF_APIN(15), PFT_EF_MYRC(15)
! <---
      REAL(hp)                :: PFT_EF_FARN(15), PFT_EF_BCAR(15)
      REAL(hp)                :: PFT_EF_OSQT(15)
      REAL(hp)                :: EM_FRAC_ALD2(15), EM_FRAC_EOH(15)
      REAL(hp)                :: EM_FRAC_FAXX(15), EM_FRAC_AAXX(15)
      REAL(hp)                :: EM_FRAC_CH2O(15)

      ! Pointers to CLM4 plant functional types
      REAL(hp)  :: PFT_BARE(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_NDLF_EVGN_TMPT_TREE(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_NDLF_EVGN_BORL_TREE(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_NDLF_DECD_BORL_TREE(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_BDLF_EVGN_TROP_TREE(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_BDLF_EVGN_TMPT_TREE(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_BDLF_DECD_TROP_TREE(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_BDLF_DECD_TMPT_TREE(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_BDLF_DECD_BORL_TREE(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_BDLF_EVGN_SHRB(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_BDLF_DECD_TMPT_SHRB(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_BDLF_DECD_BORL_SHRB(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_C3_ARCT_GRSS(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_C3_NARC_GRSS(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_C4_GRSS(HcoState%NX,HcoState%NY)
      REAL(hp)  :: PFT_CROP(HcoState%NX,HcoState%NY)

      ! Suffix
      CHARACTER(LEN=255)      :: SFX

      !=================================================================
      ! CALC_AEF begins here!
      !=================================================================

      ! Suffix
      SFX = Inst%SUFFIX

      ! ----------------------------------------------------------------
      ! Note that not all these compounds are used in standard chemistry
      ! simulations, but they are provided here for future incorporation or
      ! specialized applications. More compounds can be added as needed
      ! by adding the corresponding PFT-specific emission factors and
      ! emission category fraction.
      ! (dbm, 01/2013)
      ! ----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Point to external data
      !-----------------------------------------------------------------
      CALL HCO_EvalFld( HcoState, 'MEGAN_AEF_ISOP'//TRIM(SFX),
     &                 Inst%AEF_ISOP, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

! ---> Now compute EF maps for acetone (dbm, 12/2012)
!     CALL HCO_GetPtr(  HcoState, 'MEGAN_AEF_ACET', AEF_ACET, RC )
!     IF ( RC /= HCO_SUCCESS ) RETURN
! <---

      ! Arrays used if monoterpene switch is turned on
      IF ( Inst%ExtNrMono > 0 ) THEN
         CALL HCO_EvalFld(  HcoState,
     &                    'MEGAN_AEF_MBOX'//TRIM(SFX),
     &                     Inst%AEF_MBOX, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

! ---> Now compute EF maps for a-pinene and myrcene (dbm, 12/2012)
!         CALL HCO_GetPtr(  HcoState, 'MEGAN_AEF_APIN',AEF_APIN, RC )
!         IF ( RC /= HCO_SUCCESS ) RETURN
! <---
         CALL HCO_EvalFld(  HcoState,
     &                    'MEGAN_AEF_BPIN'//TRIM(SFX),
     &                     Inst%AEF_BPIN, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld(  HcoState,
     &                    'MEGAN_AEF_CARE'//TRIM(SFX),
     &                     Inst%AEF_CARE, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld(  HcoState,
     &                    'MEGAN_AEF_LIMO'//TRIM(SFX),
     &                     Inst%AEF_LIMO, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

! ---> Now compute EF maps for a-pinene and myrcene (dbm, 12/2012)
!         CALL HCO_GetPtr(  HcoState, 'MEGAN_AEF_MYRC',AEF_MYRC, RC )
!         IF ( RC /= HCO_SUCCESS ) RETURN
! <---

         CALL HCO_EvalFld(  HcoState,
     &                    'MEGAN_AEF_OCIM'//TRIM(SFX),
     &                     Inst%AEF_OCIM, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld(  HcoState,
     &                    'MEGAN_AEF_SABI'//TRIM(SFX),
     &                     Inst%AEF_SABI, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Only used for MEGAN Monoterpene simulation:
      IF ( Inst%ExtNrMono > 0 ) THEN
         CALL HCO_EvalFld(  HcoState, 'MEGAN_ORVC'//TRIM(SFX),
     &                     Inst%GEIA_ORVC, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      !-----------------------------------------------------------------
      ! Point to PFT fractions
      !-----------------------------------------------------------------

      ! CLM4 PFT coverage (unitless)
      ! From Table 3 in Guenther et al., 2012
      ! PFT_BARE                : Bare
      ! PFT_NDLF_EVGN_TMPT_TREE : Needleleaf evergreen temperate tree
      ! PFT_NDLF_EVGN_BORL_TREE : Needleleaf evergreen boreal tree
      ! PFT_NDLF_DECD_BORL_TREE : Needleleaf deciduous boreal tree
      ! PFT_BDLF_EVGN_TROP_TREE : Broadleaf evergreen tropical tree
      ! PFT_BDLF_EVGN_TMPT_TREE : Broadleaf evergreen temperate tree
      ! PFT_BDLF_DECD_TROP_TREE : Broadleaf deciduous tropical tree
      ! PFT_BDLF_DECD_TMPT_TREE : Broadleaf deciduous temperate tree
      ! PFT_BDLF_DECD_BORL_TREE : Broadleaf deciduous boreal tree
      ! PFT_BDLF_EVGN_SHRB      : Broadleaf evergreen temperate shrub
      ! PFT_BDLF_DECD_TMPT_SHRB : Broadleaf deciduous temperate shrub
      ! PFT_BDLF_DECD_BORL_SHRB : Broadleaf deciduous boreal shrub
      ! PFT_C3_ARCT_GRSS        : Arctic C3 grass
      ! PFT_C3_NARC_GRSS        : Cool C3 grass
      ! PFT_C4_GRSS             : Warm C4 grass
      ! PFT_CROP                : Crop

      CALL HCO_EvalFld(  HcoState, 'CLM4_PFT_BARE'//TRIM(SFX),
     &                 PFT_BARE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_NDLF_EVGN_TMPT_TREE'//TRIM(SFX),
     &                 PFT_NDLF_EVGN_TMPT_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_NDLF_EVGN_BORL_TREE'//TRIM(SFX),
     &                 PFT_NDLF_EVGN_BORL_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_NDLF_DECD_BORL_TREE'//TRIM(SFX),
     &                 PFT_NDLF_DECD_BORL_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_BDLF_EVGN_TROP_TREE'//TRIM(SFX),
     &                 PFT_BDLF_EVGN_TROP_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_BDLF_EVGN_TMPT_TREE'//TRIM(SFX),
     &                 PFT_BDLF_EVGN_TMPT_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_BDLF_DECD_TROP_TREE'//TRIM(SFX),
     &                 PFT_BDLF_DECD_TROP_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_BDLF_DECD_TMPT_TREE'//TRIM(SFX),
     &                 PFT_BDLF_DECD_TMPT_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_BDLF_DECD_BORL_TREE'//TRIM(SFX),
     &                 PFT_BDLF_DECD_BORL_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_BDLF_EVGN_SHRB'//TRIM(SFX),
     &                 PFT_BDLF_EVGN_SHRB, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_BDLF_DECD_TMPT_SHRB'//TRIM(SFX),
     &                 PFT_BDLF_DECD_TMPT_SHRB, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                 'CLM4_PFT_BDLF_DECD_BORL_SHRB'//TRIM(SFX),
     &                 PFT_BDLF_DECD_BORL_SHRB, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                'CLM4_PFT_C3_ARCT_GRSS'//TRIM(SFX),
     &                 PFT_C3_ARCT_GRSS, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                'CLM4_PFT_C3_NARC_GRSS'//TRIM(SFX),
     &                 PFT_C3_NARC_GRSS, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState,
     &                'CLM4_PFT_C4_GRSS'//TRIM(SFX),
     &                 PFT_C4_GRSS, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_EvalFld(  HcoState, 'CLM4_PFT_CROP'//TRIM(SFX),
     &                 PFT_CROP, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Copy PFTs into ARRAY_16
      Inst%ARRAY_16(:,:, 1) = PFT_BARE
      Inst%ARRAY_16(:,:, 2) = PFT_NDLF_EVGN_TMPT_TREE
      Inst%ARRAY_16(:,:, 3) = PFT_NDLF_EVGN_BORL_TREE
      Inst%ARRAY_16(:,:, 4) = PFT_NDLF_DECD_BORL_TREE
      Inst%ARRAY_16(:,:, 5) = PFT_BDLF_EVGN_TROP_TREE
      Inst%ARRAY_16(:,:, 6) = PFT_BDLF_EVGN_TMPT_TREE
      Inst%ARRAY_16(:,:, 7) = PFT_BDLF_DECD_TROP_TREE
      Inst%ARRAY_16(:,:, 8) = PFT_BDLF_DECD_TMPT_TREE
      Inst%ARRAY_16(:,:, 9) = PFT_BDLF_DECD_BORL_TREE
      Inst%ARRAY_16(:,:,10) = PFT_BDLF_EVGN_SHRB
      Inst%ARRAY_16(:,:,11) = PFT_BDLF_DECD_TMPT_SHRB
      Inst%ARRAY_16(:,:,12) = PFT_BDLF_DECD_BORL_SHRB
      Inst%ARRAY_16(:,:,13) = PFT_C3_ARCT_GRSS
      Inst%ARRAY_16(:,:,14) = PFT_C3_NARC_GRSS
      Inst%ARRAY_16(:,:,15) = PFT_C4_GRSS
      Inst%ARRAY_16(:,:,16) = PFT_CROP

      ! --------------------------------------------------------------------------------
      ! PFT-specific EFs from Table 2 in Guenther et al., 2012
      ! in ug compound/m2/h
      ! PFTs 1-15 in the table correspond to #2-16
      ! (i.e., excluding bare ground #1) in the above array.
      ! --------------------------------------------------------------------------------
      ! Compound Class EF1 EF2 EF3 EF4 EF5 EF6 EF7 EF8 EF9 EF10 EF11 EF12 EF13 EF14 EF15
      ! --------------------------------------------------------------------------------
      ! Other Monoterp 180 180 170 150 150 150 150 150 110 200  110  5    5    5    5
      ! Methanol       900 900 900 500 900 500 900 900 900 900  900  500  500  500  900
      ! Acetone        240 240 240 240 240 240 240 240 240 240  240  80   80   80   80
      ! Bidirect VOC   500 500 500 500 500 500 500 500 500 500  500  80   80   80   80
      ! Stress VOC     300 300 300 300 300 300 300 300 300 300  300  300  300  300  300
      ! Other VOC      140 140 140 140 140 140 140 140 140 140  140  140  140  140  140
      ! a-Pinene       500 500 510 600 400 600 400 400 200 300  200    2    2    2    2
      ! Myrcene         70  70  60  80  30  80  30  30  30  50   30  0.3  0.3  0.3  0.3
      ! a-Farnesene     40  40  40  60  40  60  40  40  40  40   40    3    3    3    4
      ! b-Carophyllene  80  80  80  60  40  60  40  40  50  50   50    1    1    1    4
      ! Other sesqt.   120 120 120 120 100 120 100 100 100 100  100    2    2    2    2
      ! --------------------------------------------------------------------------------

      ! One thing to note is these are net emissions to the canopy atmosphere
      ! but not net emissions to the above canopy atmosphere since they don't
      !  account for within-canopy deposition. Only an issue for OVOCs.

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_OMON = (/180.0_hp, 180.0_hp, 170.0_hp, 150.0_hp, 150.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                150.0_hp, 150.0_hp, 150.0_hp, 110.0_hp, 200.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                110.0_hp, 5.0_hp  , 5.0_hp  , 5.0_hp  , 5.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_MOHX = (/900.0_hp, 900.0_hp, 900.0_hp, 500.0_hp, 900.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                500.0_hp, 900.0_hp, 900.0_hp, 900.0_hp, 900.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                900.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 900.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_ACET = (/240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                240.0_hp, 80.0_hp , 80.0_hp , 80.0_hp , 80.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_BIDR = (/500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                500.0_hp, 80.0_hp , 80.0_hp , 80.0_hp , 80.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_STRS = (/300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_OTHR = (/140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp/)

! ---> Now compute EFs for a-pinene and myrcene as well (dbm, 12/2012)
!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_APIN = (/500.0_hp, 500.0_hp, 510.0_hp, 600.0_hp, 400.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                600.0_hp, 400.0_hp, 400.0_hp, 200.0_hp, 300.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                200.0_hp, 2.0_hp,   2.0_hp,   2.0_hp,   2.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_MYRC = (/70.0_hp,  70.0_hp,  60.0_hp,  80.0_hp,  30.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                80.0_hp,  30.0_hp,  30.0_hp,  30.0_hp,  50.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                30.0_hp,  0.3_hp,   0.3_hp,   0.3_hp,   0.3_hp/)
! <---

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_FARN = (/40.0_hp,  40.0_hp,  40.0_hp,  60.0_hp,  40.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                60.0_hp,  40.0_hp,  40.0_hp,  40.0_hp,  40.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                40.0_hp,  3.0_hp,   3.0_hp,   3.0_hp,   4.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_BCAR = (/80.0_hp,  80.0_hp,  80.0_hp,  60.0_hp,  40.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                60.0_hp,  40.0_hp,  40.0_hp,  50.0_hp,  50.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                50.0_hp,  1.0_hp,   1.0_hp,   1.0_hp,   4.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5
      PFT_EF_OSQT = (/120.0_hp, 120.0_hp, 120.0_hp, 120.0_hp, 100.0_hp,
!                     EF6       EF7       EF8       EF9       EF10
     &                120.0_hp, 100.0_hp, 100.0_hp, 100.0_hp, 100.0_hp,
!                     EF11      EF12      EF13      EF14      EF15
     &                100.0_hp, 2.0_hp,   2.0_hp,   2.0_hp,   2.0_hp/)


      ! Other monoterpenes, methanol, acetone, MBO are each 100% of thier
      ! respective categories. The VOCs within the stress category each
      ! account for a specific fraction of emissions across all PFTs
      ! (ethene 58%, toluene 3%, HCN 1.5%). The VOCs within the
      ! other category also account for a given fraction of emissions
      ! across all PFTs (propene 48%, butene 24%, other alkenes 0.2%). But
      ! VOCs in the bidirectional category account for a different amount of
      ! the total flux for the different PFTs. So in this case we define a
      ! vector containing these fractions.

      ! Acetaldehyde: 40% of bidirectional category flux, except 25%
      ! for grasses and crops
      EM_FRAC_ALD2 = (/0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp,
     &                 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp,
     &                 0.40_hp, 0.25_hp, 0.25_hp, 0.25_hp, 0.25_hp/)

      ! Ethanol: 40% of bidirectional category flux, except 25%
      ! for grasses and crops
      EM_FRAC_EOH  = (/0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp,
     &                 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp,
     &                 0.40_hp, 0.25_hp, 0.25_hp, 0.25_hp, 0.25_hp/)

      ! Formic acid: 6% of bidirectional category flux, except 15%
      ! for grasses and crops
      EM_FRAC_FAXX = (/0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp,
     &                 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp,
     &                 0.06_hp, 0.15_hp, 0.15_hp, 0.15_hp, 0.15_hp/)

      ! Acetic acid: 6% of bidirectional category flux, except 15%
      ! for grasses and crops
      EM_FRAC_AAXX = (/0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp,
     &                 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp,
     &                 0.06_hp, 0.15_hp, 0.15_hp, 0.15_hp, 0.15_hp/)

      ! Formaldehyde: 8% of bidirectional category flux, except 20%
      ! for grasses and crops
      EM_FRAC_CH2O = (/0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp,
     &                 0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp,
     &                 0.08_hp, 0.20_hp, 0.20_hp, 0.20_hp, 0.20_hp/)

      ! Loop through plant types
      DO P = 1, 15

         ! Add 1 to Array_16 index to skip bare ground
         ARR_IND = P + 1

         ! Don't need to divide ARRAY_16 by 100 since data from netCDF
         ! file has already been converted to fraction (mps, 2/12/15)

         IF ( Inst%ExtNrMono > 0 ) THEN

! ---> Now compute EFs for a-pinene and myrcene as well (dbm, 12/2012)
            ! a-pinene: 100% of category
            Inst%AEF_APIN(:,:) = Inst%AEF_APIN(:,:) +
     &           Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_APIN(P)

            ! Myrcene: 100% of category
            Inst%AEF_MYRC(:,:) = Inst%AEF_MYRC(:,:) +
     &           Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_MYRC(P)
! <---
            ! Other monoterpenes: 100% of category
            Inst%AEF_OMON(:,:) = Inst%AEF_OMON(:,:) +
     &           Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_OMON(P)

            ! a-Farnesene: 100% of category
            Inst%AEF_FARN(:,:) = Inst%AEF_FARN(:,:) +
     &           Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_FARN(P)

            ! b-Caryophyllene: 100% of category
            Inst%AEF_BCAR(:,:) = Inst%AEF_BCAR(:,:) +
     &           Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_BCAR(P)

            ! Other sesquiterpenes: 100% of category
            Inst%AEF_OSQT(:,:) = Inst%AEF_OSQT(:,:) +
     &           Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_OSQT(P)

         ENDIF

         ! Methanol: 100% of category
         Inst%AEF_MOHX(:,:) = Inst%AEF_MOHX(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_MOHX(P)

         ! Acetone: 100% of category
         Inst%AEF_ACET(:,:) = Inst%AEF_ACET(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_ACET(P)

         ! Ethanol: variable fraction of category
         Inst%AEF_EOH(:,:) = Inst%AEF_EOH(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_EOH(P)

         ! Formaldehyde: variable fraction of category
         Inst%AEF_CH2O(:,:) = Inst%AEF_CH2O(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_CH2O(P)

         ! Acetaldehyde: variable fraction of category
         Inst%AEF_ALD2(:,:) = Inst%AEF_ALD2(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_ALD2(P)

         ! Formic acid: variable fraction of category
         Inst%AEF_FAXX(:,:) = Inst%AEF_FAXX(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_FAXX(P)

         ! Acetic acid: variable fraction of category
         Inst%AEF_AAXX(:,:) = Inst%AEF_AAXX(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_AAXX(P)

         ! Ethene: 58% of "stress" category for all PFTs
         Inst%AEF_C2H4(:,:) = Inst%AEF_C2H4(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_STRS(P) * 0.58_hp

         ! Toluene: 3% of "stress" category for all PFTs
         Inst%AEF_TOLU(:,:) = Inst%AEF_TOLU(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_STRS(P) * 0.03_hp

         ! HCN: 1.5% of "stress" category for all PFTs
         Inst%AEF_HCNX(:,:) = Inst%AEF_HCNX(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_STRS(P) * 0.015_hp

         ! Propene: 48% of "other" category for all PFTs
         ! Butene:  24% of "other" category for all PFTs
         ! Larger alkenes: 0.2% of "other" category for all PFTs
         ! Total: 72.2%
         Inst%AEF_PRPE(:,:) = Inst%AEF_PRPE(:,:) +
     &        Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_OTHR(P) * 0.722_hp

      ENDDO

!      ! Nullify pointers
!      PFT_BARE                => NULL()
!      PFT_NDLF_EVGN_TMPT_TREE => NULL()
!      PFT_NDLF_EVGN_BORL_TREE => NULL()
!      PFT_NDLF_DECD_BORL_TREE => NULL()
!      PFT_BDLF_EVGN_TROP_TREE => NULL()
!      PFT_BDLF_EVGN_TMPT_TREE => NULL()
!      PFT_BDLF_DECD_TROP_TREE => NULL()
!      PFT_BDLF_DECD_TMPT_TREE => NULL()
!      PFT_BDLF_DECD_BORL_TREE => NULL()
!      PFT_BDLF_EVGN_SHRB      => NULL()
!      PFT_BDLF_DECD_TMPT_SHRB => NULL()
!      PFT_BDLF_DECD_BORL_SHRB => NULL()
!      PFT_C3_ARCT_GRSS        => NULL()
!      PFT_C3_NARC_GRSS        => NULL()
!      PFT_C4_GRSS             => NULL()
!      PFT_CROP                => NULL()

      ! Conversion factor from [ug compound/m2/hr] to [kg compound/m2/s]
      FACTOR = 1.0e-9_hp / 3600.0_hp

      ! Loop over grid boxes
      DO J = 1, HcoState%NY
      DO I = 1, HcoState%NX

         ! Convert AEF arrays to [kgC/m2/s]
         ! Multiply arrays by FACTOR and ratio [g C/g compound]
         ! NOTE: AEFs for ISOP, MBOX, BPIN, CARE, LIMO, OCIM, SABI
         ! are read from file in [kgC/m2/s], so no need to convert here
         IF ( Inst%ExtNrMono > 0 ) THEN
            Inst%AEF_APIN(I,J) = Inst%AEF_APIN(I,J)
     &                         * FACTOR * 120.0_hp / 136.234_hp
            Inst%AEF_MYRC(I,J) = Inst%AEF_MYRC(I,J)
     &                         * FACTOR * 120.0_hp / 136.234_hp
            Inst%AEF_OMON(I,J) = Inst%AEF_OMON(I,J)
     &                         * FACTOR * 120.0_hp / 136.234_hp

            ! Sesquiterpenes
            SPECIES2CARBON = 15.d0 * 12.01d0 /
     &                     ( 15.d0 * 12.01d0 + 24.d0 * 1.01d0 )
            Inst%AEF_FARN(I,J) = Inst%AEF_FARN(I,J)
     &                         * FACTOR * SPECIES2CARBON
            Inst%AEF_BCAR(I,J) = Inst%AEF_BCAR(I,J)
     &                         * FACTOR * SPECIES2CARBON
            Inst%AEF_OSQT(I,J) = Inst%AEF_OSQT(I,J)
     &                         * FACTOR * SPECIES2CARBON

         ENDIF

         Inst%AEF_ACET(I,J) = Inst%AEF_ACET(I,J)
     &                      * FACTOR *  36.0_hp /  58.079_hp
         Inst%AEF_EOH(I,J)  = Inst%AEF_EOH(I,J)
     &                      * FACTOR *  24.0_hp /  46.068_hp
         Inst%AEF_ALD2(I,J) = Inst%AEF_ALD2(I,J)
     &                      * FACTOR *  24.0_hp /  44.053_hp
         Inst%AEF_C2H4(I,J) = Inst%AEF_C2H4(I,J)
     &                      * FACTOR *  24.0_hp /  28.053_hp
         Inst%AEF_TOLU(I,J) = Inst%AEF_TOLU(I,J)
     &                      * FACTOR *  84.0_hp /  92.138_hp
         Inst%AEF_PRPE(I,J) = Inst%AEF_PRPE(I,J)
     &                      * FACTOR *  36.0_hp /  42.080_hp

         ! Methanol, formaldehyde, formic acid, acetic acid, HCN are
         ! carried in kg, not kg C
         ! Convert AEF arrays to [kg/m2/s]
         Inst%AEF_MOHX(I,J) = Inst%AEF_MOHX(I,J) * FACTOR
         Inst%AEF_CH2O(I,J) = Inst%AEF_CH2O(I,J) * FACTOR
         Inst%AEF_FAXX(I,J) = Inst%AEF_FAXX(I,J) * FACTOR
         Inst%AEF_AAXX(I,J) = Inst%AEF_AAXX(I,J) * FACTOR
         Inst%AEF_HCNX(I,J) = Inst%AEF_HCNX(I,J) * FACTOR

      ENDDO
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE CALC_AEF
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Megan_Init
!
! !DESCRIPTION: Subroutine HCOX\_Megan\_Init allocates and initializes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_Megan_Init( HcoState, ExtName,
     &                            ExtState,    RC                )
!
! !USES:
!
      USE HCO_STATE_MOD,    ONLY : Hco_GetHcoID
      USE HCO_STATE_MOD,    ONLY : HCO_GetExtHcoID
      USE HCO_ExtList_Mod,  ONLY : GetExtNr, GetExtOpt
      USE HCO_Restart_Mod,  ONLY : HCO_RestartDefine
!
! !INPUT PARAMETERS:
!
      !LOGICAL,          INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State),  POINTER        :: HcoState
      CHARACTER(LEN=*), INTENT(IN   )  :: ExtName
      TYPE(Ext_State),  POINTER        :: ExtState
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER, INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  (1 ) Change the logic in the #if block for G4AHEAD. (bmy, 12/6/05)
!  (2 ) Bug fix: skip Feb 29th if GCAP (phs, 9/18/07)
!  (3 ) Now call GET_AEF_05x0666 for GEOS-5 nested grids (yxw,dan,bmy, 11/6/08)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  26 Aug 2010 - R. Yantosca - Now reference merra_a1_mod.f
!  01 Sep 2010 - R. Yantosca - Now read in NUM_DAYS of sfc temp data (this had
!                              been hardwired to 15 days previously)
!  07 Feb 2011 - R. Yantosca - Fix typos: make sure to zero out the proper
!                              PARDF_* and PARDR_* arrays after allocation
!  22 Nov 2011 - R. Yantosca - Do not use erroneous AEF's for nested grids
!  08 Feb 2012 - R. Yantosca - Now read surface temperature for GEOS-5.7.x
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  11 Apr 2012 - R. Yantosca - Now remove the call to INIT_LAI; we shall
!                              initialize the LAI arrays from main.F
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  11 Apr 2013 - R. Yantosca - Now pass directory info with Input_Opt
!  18 Feb 2015 - M. Sulprizio- Now allocate AEF arrays for species not read
!                              from file
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  17 Jul 2017 - C. Keller   - Added optional setting NORMLAI
!  30 Aug 2018 - C. Keller   - Added error trap for wrong config file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                        :: ExtNr, iExtNr
      INTEGER                        :: I, nSpc, AS, NX, NY
      INTEGER, ALLOCATABLE           :: HcoIDs(:)
      LOGICAL                        :: FOUND
      REAL*8                         :: PI_180
      REAL(hp), POINTER              :: Ptr2D(:,:)
      TYPE(MyInst), POINTER          :: Inst
      CHARACTER(LEN=255)             :: MSG
      CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
      LOGICAL                        :: Optfound

      !=================================================================
      ! HCOX_MEGAN_INIT begins here!
      !=================================================================

      ! Extension Nr.
      ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
      IF ( ExtNr <= 0 ) RETURN

      ! Enter
      CALL HCO_ENTER( HcoState%Config%Err,
     &                  'HCOX_Megan_Init (hcox_megan_mod.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Nullify
      Ptr2D => NULL()
      Inst  => NULL()

      ! Create AeroCom instance for this simulation
      CALL InstCreate ( ExtNr, ExtState%Megan, Inst, RC )
      IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err,
     &                 'Cannot create MEGAN instance', RC )
       RETURN
      ENDIF

      !-----------------------------------------------------------------
      ! Set settings and species IDs
      !-----------------------------------------------------------------

      ! Read settings specified in configuration file
      ! Note: the specified strings have to match those in
      !       the config. file!
      CALL GetExtOpt( HcoState%Config, ExtNr,
     &                 'Isoprene scaling',
     &                 OptValHp=Inst%ISOP_SCALING, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      CALL GetExtOpt( HcoState%Config, ExtNr,
     &                 'CO2 inhibition',
     &                 OptValBool=Inst%LISOPCO2, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      CALL GetExtOpt( HcoState%Config, ExtNr,
     &                 'CO2 conc (ppmv)',
     &                 OptValHp=Inst%GLOBCO2, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Normalize LAI by PFT? Default setting is 'yes'
      ! ckeller, 7/17/17.
      Inst%NORMLAI = .TRUE.
      CALL GetExtOpt( HcoState%Config, ExtNr,
     &                 'Normalize LAI', OptValBool=OptFound,
     &                 Found=FOUND, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      IF ( FOUND ) Inst%NORMLAI = OptFound

      ! Check GLOBCO2 if CO2 inhibition is turned on (LISOPCO2 = .TRUE.)
      ! GLOBCO2 should be between 150-1250 ppmv. Isoprene response to
      ! CO2 outside this range has no empirical basis. (Tai, Jan 2013)
      IF ( Inst%LISOPCO2 ) THEN
         IF ( Inst%GLOBCO2 <  150.0_hp .OR.
     &        Inst%GLOBCO2 > 1250.0_hp     ) THEN
            MSG = 'Global CO2 outside valid range of 150-1250 ppmv!'
            CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
            RETURN
         ENDIF
      ENDIF

      Optfound = .FALSE.
      !optional arguments for SOAP
      CALL GetExtOpt ( HcoState%Config, ExtNr, 'Isoprene to SOAP',
     &                 OptValHp=Inst%ISOPTOSOAP, FOUND=Optfound, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      IF ( Optfound ) THEN
         !convert from carbon basis to mass basis
         Inst%ISOPTOSOAP = Inst%ISOPTOSOAP * 1.134
      ELSE
         !set to zero if not specified
         Inst%ISOPTOSOAP = 0.0
      ENDIF

      Optfound = .FALSE.
      !optional arguments for SOAS
      CALL GetExtOpt ( HcoState%Config, ExtNr, 'Isoprene to SOAS',
     &                 OptValHp=Inst%ISOPTOSOAS, FOUND=Optfound, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      IF ( Optfound ) THEN
         !convert from carbon basis to mass basis
         Inst%ISOPTOSOAS = Inst%ISOPTOSOAS * 1.134
      ELSE
         !set to zero if not specified
         Inst%ISOPTOSOAS = 0.0
      ENDIF

      ! Verbose mode
      IF ( HcoState%amIRoot ) THEN
         MSG = 'Use MEGAN biogenic emissions (extension module)'
         CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
         WRITE(MSG,*) ' - This is instance ',ExtState%Megan
         CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
         MSG = ' - Use the following species:'
         CALL HCO_MSG(HcoState%Config%Err,MSG )
      ENDIF

      ! Get species IDs
      ! --> Assume that species are ordered ISOP, ACET, PRPE, C2H4 in
      !     config. file!
      CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC)
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Assign species IDs
      Inst%IDTISOP = -1
      Inst%IDTACET = -1
      Inst%IDTPRPE = -1
      Inst%IDTC2H4 = -1
      Inst%IDTALD2 = -1
      Inst%IDTEOH  = -1
      Inst%IDTSOAP = -1
      Inst%IDTSOAS = -1
      DO I = 1, nSpc
         SELECT CASE ( TRIM(SpcNames(I)) )
            CASE( 'ISOP' )
               Inst%IDTISOP = HcoIDs(I)
               WRITE(MSG,*) '   Isoprene   = ',
     &                      TRIM(SpcNames(I)),Inst%IDTISOP
            CASE( 'ACET' )
               Inst%IDTACET = HcoIDs(I)
               WRITE(MSG,*) '   Acetone    = ',
     &                      TRIM(SpcNames(I)),Inst%IDTACET
            CASE( 'PRPE' )
               Inst%IDTPRPE = HcoIDs(I)
               WRITE(MSG,*) '   C3 Alkenes = ',
     &                      TRIM(SpcNames(I)),Inst%IDTPRPE
            CASE( 'ALD2' )
               Inst%IDTALD2 = HcoIDs(I)
               WRITE(MSG,*) '   ALD2       = ',
     &                      TRIM(SpcNames(I)),Inst%IDTALD2
            CASE( 'EOH' )
               Inst%IDTEOH  = HcoIDs(I)
               WRITE(MSG,*) '   EOH        = ',
     &                      TRIM(SpcNames(I)),Inst%IDTEOH
            CASE( 'C2H4' )
               Inst%IDTC2H4 = HcoIDs(I)
               WRITE(MSG,*) '   Ethene     = ',
     &                      TRIM(SpcNames(I)),Inst%IDTC2H4
            CASE( 'SOAP' )
               Inst%IDTSOAP = HcoIDs(I)
               WRITE(MSG,*) 'SOA-Precursor = ',
     &                      TRIM(SpcNames(I)),Inst%IDTSOAP
            CASE( 'SOAS' )
               Inst%IDTSOAS = HcoIDs(I)
               WRITE(MSG,*) '  SOA-Simple  = ',
     &                      TRIM(SpcNames(I)),Inst%IDTSOAS
            CASE DEFAULT
               MSG = 'Invalid species names: ' // TRIM(SpcNames(I))
               CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
               RETURN
         END SELECT

         ! Verbose
         IF ( HcoState%amIRoot ) CALL HCO_MSG(HcoState%Config%Err,MSG )
      ENDDO

      ! Verbose mode
      IF ( HcoState%amIRoot ) THEN
         WRITE(MSG,*) ' --> Isoprene scale factor is ',Inst%ISOP_SCALING
         CALL HCO_MSG(HcoState%Config%Err,MSG )
         WRITE(MSG,*) ' --> Use CO2 inhibition on isoprene option ',
     &                Inst%LISOPCO2
         CALL HCO_MSG(HcoState%Config%Err,MSG )
         WRITE(MSG,*) ' --> Global atmospheric CO2 concentration : ',
     &                Inst%GLOBCO2, ' ppmv'
         CALL HCO_MSG(HcoState%Config%Err,MSG )
         WRITE(MSG,*) ' --> Normalize LAI by PFT: ',
     &                Inst%NORMLAI
         CALL HCO_MSG(HcoState%Config%Err,MSG )
      ENDIF

      IF ( HcoState%amIRoot ) THEN
         !convert back to direct mass basis just for show
         WRITE(MSG,*) ' --> Isoprene to SOA-Precursor',
     &       Inst%ISOPTOSOAP / 1.134
         CALL HCO_MSG( MSG )
         WRITE(MSG,*) ' --> Isoprene direct to SOA (Simple)',
     &       Inst%ISOPTOSOAS / 1.134
         CALL HCO_MSG( MSG )
      ENDIF

      ! Check for monoterpene option --> required for CO, MTPA, MTPO, LIMO, SESQ
      Inst%IDTCO     = -1
      Inst%IDTMTPA   = -1
      Inst%IDTMTPO   = -1
      Inst%IDTLIMO   = -1
      Inst%IDTSESQ   = -1
      Inst%ExtNrMono = GetExtNr( HcoState%Config%ExtList, 'MEGAN_Mono')

      ! Error trap: previous MEGAN versions had a 3rd option, SOA. This
      ! Is now merged with Mono. Check here if the HEMCO configuration file
      ! Still uses the old syntax. In that case stop with an error
      ! (ckeller, 8/30/18)
      iExtNr = GetExtNr( HcoState%Config%ExtList, 'MEGAN_SOA')
      IF ( iExtNr > 0 ) THEN
         IF ( HcoState%amIRoot ) THEN
            MSG = 'The HEMCO configuration file uses the old MEGAN' //
     &            'syntax - Please merge SOA into Mono, e.g.:'
            CALL HCO_MSG( MSG )
            MSG = '108     MEGAN             : on    '//
     &            'ISOP/ACET/PRPE/C2H4/ALD2/EOH/SOAP/SOAS'
            CALL HCO_MSG( MSG )
            MSG = '--> Isoprene scaling  :       1.0'
            CALL HCO_MSG( MSG )
            MSG = '--> CO2 inhibition    :       true'
            CALL HCO_MSG( MSG )
            MSG = '--> CO2 conc (ppmv)   :       390.0'
            CALL HCO_MSG( MSG )
            MSG = '--> Isoprene to SOAP  :       0.015'
            CALL HCO_MSG( MSG )
            MSG = '--> Isoprene to SOAS  :       0.015'
            CALL HCO_MSG( MSG )
            MSG = '--> Normalize LAI     :       true'
            CALL HCO_MSG( MSG )
            MSG='#3% mass yield from isoprene to SOA'
            CALL HCO_MSG( MSG )
            MSG='109     MEGAN_Mono      : on    CO/MTPA/MTPO/LIMO/SESQ'
            CALL HCO_MSG( MSG )
            MSG='--> Monoterp to SOAP  :       0.050'
            CALL HCO_MSG( MSG )
            MSG='--> Monoterp to SOAS  :       0.050'
            CALL HCO_MSG( MSG )
            MSG='#10% mass yield from monoterpenes to SOA'
            CALL HCO_MSG( MSG )
            MSG='--> Othrterp to SOAP  :       0.050'
            CALL HCO_MSG( MSG )
            MSG='--> Othrterp to SOAS  :       0.050'
            CALL HCO_MSG( MSG )
         ENDIF
         MSG = 'Your HEMCO configuration file uses old MEGAN syntax'
         CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
         RETURN
      ENDIF

      IF ( Inst%ExtNrMono > 0 ) THEN

         ! Verbose
         IF ( HcoState%amIRoot ) THEN
            MSG = ' - MEGAN monoterpene option enabled'
            CALL HCO_MSG(HcoState%Config%Err,MSG )
         ENDIF

         ! Get HEMCO species IDs
         CALL HCO_GetExtHcoID( HcoState, Inst%ExtNrMono, HcoIDs,
     &                         SpcNames, nSpc,           RC       )
         IF ( RC /= HCO_SUCCESS ) RETURN

         Optfound = .FALSE.
         CALL GetExtOpt ( HcoState%Config, Inst%ExtNrMono,
     &                 'Monoterp to SOAP',
     &                 OptValHp=Inst%MONOTOSOAP, FOUND=Optfound, RC=RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         IF ( Optfound ) THEN
            !convert from carbon basis to mass basis
            Inst%MONOTOSOAP = Inst%MONOTOSOAP! * 1.134
         ELSE
            !set to zero if not specified
            Inst%MONOTOSOAP = 0.0
         ENDIF

         Optfound = .FALSE.
         CALL GetExtOpt ( HcoState%Config, Inst%ExtNrMono,
     &                 'Monoterp to SOAS',
     &                 OptValHp=Inst%MONOTOSOAS, FOUND=Optfound, RC=RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         IF ( Optfound ) THEN
            !convert from carbon basis to mass basis
            Inst%MONOTOSOAS = Inst%MONOTOSOAS! * 1.134
         ELSE
            !set to zero if not specified
            Inst%MONOTOSOAS = 0.0
         ENDIF

         Optfound = .FALSE.
         CALL GetExtOpt ( HcoState%Config, Inst%ExtNrMono,
     &                 'Othrterp to SOAP',
     &                 OptValHp=Inst%OTHRTOSOAP, FOUND=Optfound, RC=RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         IF ( Optfound ) THEN
            !convert from carbon basis to mass basis
            Inst%OTHRTOSOAP = Inst%OTHRTOSOAP * 1.134
         ELSE
            !set to zero if not specified
            Inst%OTHRTOSOAP = 0.0
         ENDIF

         Optfound = .FALSE.
         CALL GetExtOpt ( HcoState%Config, Inst%ExtNrMono,
     &                 'Othrterp to SOAS',
     &                 OptValHp=Inst%OTHRTOSOAS, FOUND=Optfound, RC=RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         IF ( Optfound ) THEN
            !convert from carbon basis to mass basis
            Inst%OTHRTOSOAS = Inst%OTHRTOSOAS * 1.134
         ELSE
            !set to zero if not specified
            Inst%OTHRTOSOAS = 0.0
         ENDIF

         DO I = 1, nSpc
            SELECT CASE ( TRIM(SpcNames(I)) )
               CASE( 'CO' )
                  Inst%IDTCO = HcoIDs(I)
                  WRITE(MSG,*) '   CO         = ',
     &                         TRIM(SpcNames(I)),Inst%IDTCO
               CASE( 'MTPA' )
                  Inst%IDTMTPA = HcoIDs(I)
                  WRITE(MSG,*) '     a-,b-pinene    = ',
     &               TRIM(SpcNames(I)),Inst%IDTMTPA
               CASE( 'MTPO' )
                  Inst%IDTMTPO = HcoIDs(I)
                  WRITE(MSG,*) '     Other monoterp.= ',
     &               TRIM(SpcNames(I)),Inst%IDTMTPO
               CASE( 'LIMO' )
                  Inst%IDTLIMO = HcoIDs(I)
                  WRITE(MSG,*) '     Limonene       = ',
     &               TRIM(SpcNames(I)),Inst%IDTLIMO
               CASE( 'SESQ' )
                  Inst%IDTSESQ = HcoIDs(I)
                  WRITE(MSG,*) '     Sesquiterpenes = ',
     &               TRIM(SpcNames(I)),Inst%IDTSESQ
               CASE( 'SOAP' )
                  IF (Inst%IDTSOAP < 0) THEN
                     Inst%IDTSOAP = HcoIDs(I)
                  END IF
                  WRITE(MSG,*) '     SOA-Precursor = ',
     &               TRIM(SpcNames(I)),Inst%IDTSOAP
               CASE( 'SOAS' )
                  IF (Inst%IDTSOAS < 0) THEN
                     Inst%IDTSOAS = HcoIDs(I)
                  END IF
                  WRITE(MSG,*) '     SOA-Simple = ',
     &               TRIM(SpcNames(I)),Inst%IDTSOAS
               CASE DEFAULT
                  MSG = 'Invalid species names: ' // TRIM(SpcNames(I))
                  CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
                  RETURN
            END SELECT

            ! Verbose
            IF (HcoState%amIRoot) CALL HCO_MSG(HcoState%Config%Err,MSG)
         ENDDO

         ! Check if at least one of the Mono species is defined
         IF ( Inst%IDTCO   < 0 .AND.
     &        Inst%IDTMTPA < 0 .AND.
     &        Inst%IDTMTPO < 0 .AND.
     &        Inst%IDTLIMO < 0 .AND.
     &        Inst%IDTSESQ < 0 .AND.
     &        Inst%IDTSOAP < 0 .AND.
     &        Inst%IDTSOAS < 0  ) THEN
            MSG = 'MEGAN_MONO option is anabled but none of the '  //
     &            'species is found (CO, MTPA, MTPO, LIMO, SESQ, ' //
     &            'SOAP, SOAS)!'
            CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
            RETURN
         ENDIF

         IF ( HcoState%amIRoot ) THEN
            !convert back to direct mass basis just for show
            WRITE(MSG,*) ' --> Monoterpene to SOA-Precursor',
     &       Inst%MONOTOSOAP / 1.134
            CALL HCO_MSG( MSG )
            WRITE(MSG,*) ' --> Monoterpene direct to SOA (Simple)',
     &       Inst%MONOTOSOAS / 1.134
            CALL HCO_MSG( MSG )
            WRITE(MSG,*) ' --> Othrterpene to SOA-Precursor',
     &       Inst%OTHRTOSOAP / 1.134
            CALL HCO_MSG( MSG )
            WRITE(MSG,*) ' --> Othrterpene direct to SOA (Simple)',
     &       Inst%OTHRTOSOAS / 1.134
            CALL HCO_MSG( MSG )
         ENDIF

      ELSE

         ! CO emissions and acetone emissions will be zero/incomplete!
         IF ( Inst%IDTACET > 0 ) THEN
            MSG = 'Cannot fully calculate MEGAN acetone ' //
     &            'emissions because monoterpene switch is off!'
            CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
         ENDIF
         MSG = 'Will not calculate MEGAN CO emissions ' //
     &         'because monoterpene switch is off!'
         CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
      ENDIF

      CALL GetExtOpt( HcoState%Config, ExtNr, 'MEGAN_SUFFIX',
     &                OptValChar=Inst%SUFFIX, FOUND=FOUND, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      IF ( .NOT. FOUND ) Inst%SUFFIX = ''

      !-----------------------------------------------------------------
      ! Allocate module variables
      !-----------------------------------------------------------------

      ! Get horizontal grid extensions on this CPU
      NX = HcoState%NX
      NY = HcoState%NY

      ALLOCATE( Inst%T_LAST24H( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'T_LAST24H', RC )
         RETURN
      ENDIF
      Inst%T_LAST24H = 0.0_hp

      ALLOCATE( Inst%T_LASTXDAYS( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'T_LASTXDAYS', RC )
         RETURN
      ENDIF
      Inst%T_LASTXDAYS = 0.0_hp

      ALLOCATE( Inst%PARDR_LASTXDAYS( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'PARDR_LASTXDAYS', RC )
         RETURN
      ENDIF
      Inst%PARDR_LASTXDAYS = 0.0_hp

      ALLOCATE( Inst%PARDF_LASTXDAYS( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'PARDF_LASTXDAYS', RC )
         RETURN
      ENDIF
      Inst%PARDF_LASTXDAYS = 0.0_hp

      ALLOCATE( Inst%LAI_PREVDAY( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR( HcoState%Config%Err,
     &                  'LAI_PREVDAY', RC )
        RETURN
      ENDIF
      Inst%LAI_PREVDAY = 0.0_sp

      ALLOCATE( Inst%ARRAY_16( NX, NY, 16 ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR( HcoState%Config%Err,
     &                  'ARRAY_16', RC )
        RETURN
      ENDIF
      Inst%ARRAY_16 = 0.0_hp

      ! Normalization factor
      ! There should be a different normalization factor for each compound, but
      ! we calculate only 1 normalization factor for all compounds (dbm 11/2012)
      ALLOCATE( Inst%NORM_FAC( 1 ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR( HcoState%Config%Err,
     &                  'NORM_FAC', RC )
        RETURN
      ENDIF
      Inst%NORM_FAC = -99d0

      IF ( Inst%ExtNrMono > 0 ) THEN

         ALLOCATE( Inst%AEF_APIN( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_APIN', RC )
            RETURN
         ENDIF
         Inst%AEF_APIN = 0.0_hp

         ALLOCATE( Inst%AEF_MYRC( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_MYRC', RC )
            RETURN
         ENDIF
         Inst%AEF_MYRC = 0.0_hp

         ALLOCATE( Inst%AEF_OMON( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_OMON', RC )
            RETURN
         ENDIF
         Inst%AEF_OMON = 0.0_hp

         ALLOCATE( Inst%AEF_FARN( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_FARN', RC )
            RETURN
         ENDIF
         Inst%AEF_FARN = 0.0_hp

         ALLOCATE( Inst%AEF_BCAR( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_BCAR', RC )
            RETURN
         ENDIF
         Inst%AEF_BCAR = 0.0_hp

         ALLOCATE( Inst%AEF_OSQT( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN
            CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_OSQT', RC )
            RETURN
         ENDIF
         Inst%AEF_OSQT = 0.0_hp

      ENDIF

      ALLOCATE( Inst%AEF_MOHX( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_MOHX', RC )
         RETURN
      ENDIF
      Inst%AEF_MOHX = 0.0_hp

      ALLOCATE( Inst%AEF_ACET( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_ACET', RC )
         RETURN
      ENDIF
      Inst%AEF_ACET = 0.0_hp

      ALLOCATE( Inst%AEF_EOH( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_EOH', RC )
         RETURN
      ENDIF
      Inst%AEF_EOH = 0.0_hp

      ALLOCATE( Inst%AEF_CH2O( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_CH2O', RC )
         RETURN
      ENDIF
      Inst%AEF_CH2O = 0.0_hp

      ALLOCATE( Inst%AEF_ALD2( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_ALD2', RC )
         RETURN
      ENDIF
      Inst%AEF_ALD2 = 0.0_hp

      ALLOCATE( Inst%AEF_FAXX( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_FAXX', RC )
         RETURN
      ENDIF
      Inst%AEF_FAXX = 0.0_hp

      ALLOCATE( Inst%AEF_AAXX( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_AAXX', RC )
         RETURN
      ENDIF
      Inst%AEF_AAXX = 0.0_hp

      ALLOCATE( Inst%AEF_C2H4( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_C2H4', RC )
         RETURN
      ENDIF
      Inst%AEF_C2H4 = 0.0_hp

      ALLOCATE( Inst%AEF_TOLU( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_TOLU', RC )
         RETURN
      ENDIF
      Inst%AEF_TOLU = 0.0_hp

      ALLOCATE( Inst%AEF_HCNX( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_HCNX', RC )
         RETURN
      ENDIF
      Inst%AEF_HCNX = 0.0_hp

      ALLOCATE( Inst%AEF_PRPE( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF_PRPE', RC )
         RETURN
      ENDIF
      Inst%AEF_PRPE = 0.0_hp

      ALLOCATE( Inst%ARRAY_16( NX, NY, 16 ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'ARRAY_16', RC )
         RETURN
      ENDIF
      Inst%ARRAY_16 = 0.0_hp

      ALLOCATE ( Inst%AEF_ISOP ( NX, NY ),
     &           Inst%AEF_MBOX ( NX, NY ),
     &           Inst%AEF_BPIN ( NX, NY ),
     &           Inst%AEF_CARE ( NX, NY ),
     &           Inst%AEF_LIMO ( NX, NY ),
     &           Inst%AEF_OCIM ( NX, NY ),
     &           Inst%AEF_SABI ( NX, NY ),
     &           Inst%GEIA_ORVC( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN
         CALL HCO_ERROR( HcoState%Config%Err,
     &                  'AEF allocation error', RC )
         RETURN
      ENDIF
      Inst%AEF_ISOP  = 0.0_hp
      Inst%AEF_MBOX  = 0.0_hp
      Inst%AEF_BPIN  = 0.0_hp
      Inst%AEF_CARE  = 0.0_hp
      Inst%AEF_LIMO  = 0.0_hp
      Inst%AEF_OCIM  = 0.0_hp
      Inst%AEF_SABI  = 0.0_hp
      Inst%GEIA_ORVC = 0.0_hp

      !=================================================================
      ! Initialize internal diagnostics. These are the restart variables
      ! that can be used for a 'warm' start of MEGAN.
      !=================================================================
      CALL HCO_RestartDefine (      HcoState,
     &                         'LAI_PREVDAY', Inst%LAI_PREVDAY, '1', RC)
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_RestartDefine (        HcoState, 'T_PREVDAY',
     &                         Inst%T_LAST24H,  'K',      RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_RestartDefine (         HcoState, 'T_DAVG',
     &                         Inst%T_LASTXDAYS, 'K',      RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_RestartDefine (   HcoState, 'PARDR_DAVG',
     &                         Inst%PARDR_LASTXDAYS, 'W/m2',   RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_RestartDefine (  HcoState, 'PARDF_DAVG',
     &                         Inst%PARDF_LASTXDAYS, 'W/m2',   RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! The original MEGAN code used to read the emission factors here.
      ! We now get the emisson factors through HEMCO, hence no need to
      ! do this anymore (ckeller, 05/19/2014).
      !=================================================================

      ! Set physical constants
      PI_180     = HcoState%Phys%PI_180
      Inst%D2RAD = PI_180                ! Degrees to radians
      Inst%RAD2D = 1.0_hp / PI_180       ! Radians to degrees

      ! Enable met. fields
      ExtState%T2M%DoUse       = .TRUE.
      ExtState%SUNCOS%DoUse    = .TRUE.
      ExtState%PARDR%DoUse     = .TRUE.
      ExtState%PARDF%DoUse     = .TRUE.
      ExtState%LAI%DoUse       = .TRUE.
      ExtState%GWETROOT%DoUse  = .TRUE.

      ! SZAFACT field is only used for MEGAN Monoterpene simulation
      IF ( Inst%ExtNrMono > 0 ) THEN
         ExtState%SZAFACT%DoUse = .TRUE.
      ENDIF

      ! Cleanup
      Inst => NULL()

      ! Leave w/ success
      IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
      IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
      CALL HCO_LEAVE( HcoState%Config%Err,RC )

      END SUBROUTINE HCOX_Megan_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Megan_Final
!
! !DESCRIPTION: Subroutine HCOX\_Megan\_Final deallocates all allocated arrays
!  at the end of a GEOS-Chem model run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_MEGAN_FINAL ( HcoState, ExtState, RC )
!
! !USES
!
!
! !INPUT PARAMETERS:
!
      !LOGICAL,         INTENT(IN   )  :: am_I_Root     ! Root CPU?
      TYPE(HCO_State), POINTER        :: HcoState      ! HEMCO State obj
      TYPE(Ext_State), POINTER        :: ExtState      ! Extension State obj
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  17 Dec 2009 - R. Yantosca  - Added ProTeX headers
!  16 Aug 2013 - C. Keller    - Now a HEMCO extension
!  18 Feb 2015 - M. Sulprizio - Added AEF arrays
!  09 Mar 2015 - C. Keller    - Write variables to internal state in ESMF mode
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! HCOX_MEGAN_FINAL begins here!
      !=================================================================
      CALL InstRemove ( ExtState%Megan )

      END SUBROUTINE HCOX_MEGAN_FINAL
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a pointer to the desired instance.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
      INTEGER                             :: Instance
      TYPE(MyInst),     POINTER           :: Inst
      INTEGER                             :: RC
      TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      TYPE(MyInst),     POINTER    :: PrvInst

      !=================================================================
      ! InstGet begins here!
      !=================================================================

      ! Get instance. Also archive previous instance.
      PrvInst => NULL()
      Inst    => AllInst
      DO WHILE ( ASSOCIATED(Inst) )
         IF ( Inst%Instance == Instance ) EXIT
         PrvInst => Inst
         Inst    => Inst%NextInst
      END DO
      IF ( .NOT. ASSOCIATED( Inst ) ) THEN
         RC = HCO_FAIL
         RETURN
      ENDIF

      ! Pass output arguments
      IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

      ! Cleanup & Return
      PrvInst => NULL()
      RC = HCO_SUCCESS

      END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate adds a new instance to the list of
!  instances, assigns a unique instance number to this new instance, and
!  archives this instance number to output argument Instance.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
      INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
      INTEGER,       INTENT(  OUT)    :: Instance
      TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
      TYPE(MyInst), POINTER          :: TmpInst
      INTEGER                        :: nnInst

      !=================================================================
      ! InstCreate begins here!
      !=================================================================

      ! ----------------------------------------------------------------
      ! Generic instance initialization
      ! ----------------------------------------------------------------

      ! Initialize
      Inst => NULL()

      ! Get number of already existing instances
      TmpInst => AllInst
      nnInst = 0
      DO WHILE ( ASSOCIATED(TmpInst) )
         nnInst  =  nnInst + 1
         TmpInst => TmpInst%NextInst
      END DO

      ! Create new instance
      ALLOCATE(Inst)
      Inst%Instance = nnInst + 1
      Inst%ExtNr    = ExtNr

      ! Attach to instance list
      Inst%NextInst => AllInst
      AllInst       => Inst

      ! Update output instance
      Instance = Inst%Instance

      ! ----------------------------------------------------------------
      ! Type specific initialization statements follow below
      ! ----------------------------------------------------------------

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove removes an instance from the list of
! instances.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
      INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      INTEGER                     :: RC
      TYPE(MyInst), POINTER       :: PrevInst
      TYPE(MyInst), POINTER       :: Inst

      !=================================================================
      ! InstRemove begins here!
      !=================================================================

      ! Get instance. Also archive previous instance.
      PrevInst => NULL()
      Inst     => NULL()
      CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

      ! Instance-specific deallocation
      IF ( ASSOCIATED(Inst) ) THEN

       ! ----------------------------------------------------------------
       ! Type specific initialization statements follow below
       ! ----------------------------------------------------------------
       ! Nullify pointers
       IF( ASSOCIATED( Inst%GEIA_ORVC ) ) DEALLOCATE(Inst%GEIA_ORVC)
       IF( ASSOCIATED( Inst%AEF_ISOP  ) ) DEALLOCATE(Inst%AEF_ISOP)
       IF( ASSOCIATED( Inst%AEF_MBOX  ) ) DEALLOCATE(Inst%AEF_MBOX)
       IF( ASSOCIATED( Inst%AEF_BPIN  ) ) DEALLOCATE(Inst%AEF_BPIN)
       IF( ASSOCIATED( Inst%AEF_CARE  ) ) DEALLOCATE(Inst%AEF_CARE)
       IF( ASSOCIATED( Inst%AEF_LIMO  ) ) DEALLOCATE(Inst%AEF_LIMO)
       IF( ASSOCIATED( Inst%AEF_OCIM  ) ) DEALLOCATE(Inst%AEF_OCIM)
       IF( ASSOCIATED( Inst%AEF_SABI  ) ) DEALLOCATE(Inst%AEF_SABI)
       IF( ASSOCIATED( Inst%ARRAY_16  ) ) DEALLOCATE(Inst%ARRAY_16)

       ! Cleanup module arrays
       IF ( ASSOCIATED( Inst%LAI_PREVDAY     ) )
     &      DEALLOCATE( Inst%LAI_PREVDAY     )
       IF ( ASSOCIATED( Inst%T_LASTXDAYS     ) )
     &      DEALLOCATE( Inst%T_LASTXDAYS     )
       IF ( ASSOCIATED( Inst%T_LAST24H       ) )
     &      DEALLOCATE( Inst%T_LAST24H       )
       IF ( ASSOCIATED( Inst%PARDF_LASTXDAYS ) )
     &      DEALLOCATE( Inst%PARDF_LASTXDAYS )
       IF ( ASSOCIATED( Inst%PARDR_LASTXDAYS ) )
     &      DEALLOCATE( Inst%PARDR_LASTXDAYS )
       IF ( ASSOCIATED( Inst%NORM_FAC        ) )
     &      DEALLOCATE( Inst%NORM_FAC        )
       IF ( ASSOCIATED( Inst%AEF_APIN        ) )
     &      DEALLOCATE( Inst%AEF_APIN        )
       IF ( ASSOCIATED( Inst%AEF_MYRC        ) )
     &      DEALLOCATE( Inst%AEF_MYRC        )
       IF ( ASSOCIATED( Inst%AEF_OMON        ) )
     &      DEALLOCATE( Inst%AEF_OMON        )
       IF ( ASSOCIATED( Inst%AEF_MOHX        ) )
     &      DEALLOCATE( Inst%AEF_MOHX        )
       IF ( ASSOCIATED( Inst%AEF_ACET        ) )
     &      DEALLOCATE( Inst%AEF_ACET        )
       IF ( ASSOCIATED( Inst%AEF_EOH         ) )
     &      DEALLOCATE( Inst%AEF_EOH         )
       IF ( ASSOCIATED( Inst%AEF_CH2O        ) )
     &      DEALLOCATE( Inst%AEF_CH2O        )
       IF ( ASSOCIATED( Inst%AEF_ALD2        ) )
     &      DEALLOCATE( Inst%AEF_ALD2        )
       IF ( ASSOCIATED( Inst%AEF_FAXX        ) )
     &      DEALLOCATE( Inst%AEF_FAXX        )
       IF ( ASSOCIATED( Inst%AEF_AAXX        ) )
     &      DEALLOCATE( Inst%AEF_AAXX        )
       IF ( ASSOCIATED( Inst%AEF_C2H4        ) )
     &      DEALLOCATE( Inst%AEF_C2H4        )
       IF ( ASSOCIATED( Inst%AEF_TOLU        ) )
     &      DEALLOCATE( Inst%AEF_TOLU        )
       IF ( ASSOCIATED( Inst%AEF_HCNX        ) )
     &      DEALLOCATE( Inst%AEF_HCNX        )
       IF ( ASSOCIATED( Inst%AEF_PRPE        ) )
     &      DEALLOCATE( Inst%AEF_PRPE        )
       IF ( ASSOCIATED( Inst%AEF_FARN        ) )
     &      DEALLOCATE( Inst%AEF_FARN        )
       IF ( ASSOCIATED( Inst%AEF_BCAR        ) )
     &      DEALLOCATE( Inst%AEF_BCAR        )
       IF ( ASSOCIATED( Inst%AEF_OSQT        ) )
     &      DEALLOCATE( Inst%AEF_OSQT        )

       ! ----------------------------------------------------------------
       ! Pop off instance from list
       ! ----------------------------------------------------------------
       IF ( ASSOCIATED(PrevInst) ) THEN
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
      ENDIF

      END SUBROUTINE InstRemove
!EOC
      END MODULE HCOX_MEGAN_WRAPPER

