MODULE sbc_ice
   !!======================================================================
   !!                 ***  MODULE  sbc_ice  ***
   !! Surface module - LIM-3: parameters & variables defined in memory
   !!======================================================================
   !! History :  3.0  ! 2006-08  (G. Madec)  Surface module
   !!            3.2  ! 2009-06  (S. Masson) merge with ice_oce
   !!            3.3.1! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.4  ! 2011-11  (C. Harris) CICE added as an option
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option                      NO LIM 2.0 or 3.0 or CICE sea-ice model
   !!----------------------------------------------------------------------
   USE in_out_manager   ! I/O manager
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim2    = .FALSE.  !: no LIM-2 ice model
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim3    = .FALSE.  !: no LIM-3 ice model
   LOGICAL         , PUBLIC, PARAMETER ::   lk_cice    = .FALSE.  !: no CICE  ice model
   CHARACTER(len=1), PUBLIC, PARAMETER ::   cp_ice_msh = '-'      !: no grid ice-velocity
   REAL            , PUBLIC, PARAMETER ::   cldf_ice = 0.81       !: cloud fraction over sea ice, summer CLIO value   [-]
   INTEGER         , PUBLIC, PARAMETER ::   jpl = 1 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   u_ice, v_ice,fr1_i0,fr2_i0          ! jpi, jpj
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tn_ice, alb_ice, qns_ice, dqns_ice  ! (jpi,jpj,jpl)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   a_i
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qsr_ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ht_i, ht_s
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   topmelt, botmelt

   !!======================================================================
END MODULE sbc_ice
