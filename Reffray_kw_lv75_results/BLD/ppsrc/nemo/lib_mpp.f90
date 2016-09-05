












MODULE lib_mpp
   !!======================================================================
   !!                       ***  MODULE  lib_mpp  ***
   !! Ocean numerics:  massively parallel processing library
   !!=====================================================================
   !! History :  OPA  !  1994  (M. Guyon, J. Escobar, M. Imbard)  Original code
   !!            7.0  !  1997  (A.M. Treguier)  SHMEM additions
   !!            8.0  !  1998  (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI
   !!                 !  1998  (J.M. Molines) Open boundary conditions
   !!   NEMO     1.0  !  2003  (J.-M. Molines, G. Madec)  F90, free form
   !!                 !  2003  (J.M. Molines) add mpp_ini_north(_3d,_2d)
   !!             -   !  2004  (R. Bourdalle Badie)  isend option in mpi
   !!                 !  2004  (J.M. Molines) minloc, maxloc
   !!             -   !  2005  (G. Madec, S. Masson)  npolj=5,6 F-point & ice cases
   !!             -   !  2005  (R. Redler) Replacement of MPI_COMM_WORLD except for MPI_Abort
   !!             -   !  2005  (R. Benshila, G. Madec)  add extra halo case
   !!             -   !  2008  (R. Benshila) add mpp_ini_ice
   !!            3.2  !  2009  (R. Benshila) SHMEM suppression, north fold in lbc_nfd
   !!            3.2  !  2009  (O. Marti)    add mpp_ini_znl
   !!            4.0  !  2011  (G. Madec)  move ctl_ routines from in_out_manager
   !!            3.5  !  2012  (S.Mocavero, I. Epicoco) Add 'mpp_lnk_bdy_3d', 'mpp_lnk_obc_3d', 
   !!                          'mpp_lnk_bdy_2d' and 'mpp_lnk_obc_2d' routines and update
   !!                          the mppobc routine to optimize the BDY and OBC communications
   !!            3.5  !  2013  ( C. Ethe, G. Madec ) message passing arrays as local variables 
   !!            3.5  !  2013 (S.Mocavero, I.Epicoco - CMCC) north fold optimizations
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ctl_stop   : update momentum and tracer Kz from a tke scheme
   !!   ctl_warn   : initialization, namelist read, and parameters control
   !!   ctl_opn    : Open file and check if required file is available.
   !!   ctl_nam    : Prints informations when an error occurs while reading a namelist
   !!   get_unit   : give the index of an unused logical unit
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case:            Dummy module        share memory computing
   !!----------------------------------------------------------------------
   USE in_out_manager

   INTERFACE mpp_sum
      MODULE PROCEDURE mpp_sum_a2s, mpp_sum_as, mpp_sum_ai, mpp_sum_s, mpp_sum_i, mppsum_realdd, mppsum_a_realdd
   END INTERFACE
   INTERFACE mpp_max
      MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
   END INTERFACE
   INTERFACE mpp_min
      MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
   END INTERFACE
   INTERFACE mpp_minloc
      MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
   END INTERFACE
   INTERFACE mpp_maxloc
      MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
   END INTERFACE

   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .FALSE.      !: mpp flag
   LOGICAL, PUBLIC            ::   ln_nnogather          !: namelist control of northfold comms (needed here in case "key_mpp_mpi" is not used)
   INTEGER :: ncomm_ice
   INTEGER, PUBLIC            ::   mpi_comm_opa          ! opa local communicator
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION lib_mpp_alloc(kumout)          ! Dummy function
      INTEGER, INTENT(in) ::   kumout
      lib_mpp_alloc = 0
   END FUNCTION lib_mpp_alloc

   FUNCTION mynode( ldtxt, ldname, kumnam_ref, knumnam_cfg,  kumond , kstop, localComm ) RESULT (function_value)
      INTEGER, OPTIONAL            , INTENT(in   ) ::   localComm
      CHARACTER(len=*),DIMENSION(:) ::   ldtxt
      CHARACTER(len=*) ::   ldname
      INTEGER ::   kumnam_ref, knumnam_cfg , kumond , kstop
      IF( PRESENT( localComm ) ) mpi_comm_opa = localComm
      function_value = 0
      IF( .FALSE. )   ldtxt(:) = 'never done'
      CALL ctl_opn( kumond, TRIM(ldname), 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. , 1 )
   END FUNCTION mynode

   SUBROUTINE mppsync                       ! Dummy routine
   END SUBROUTINE mppsync

   SUBROUTINE mpp_sum_as( parr, kdim, kcom )      ! Dummy routine
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_as: You should not have seen this print! error?', kdim, parr(1), kcom
   END SUBROUTINE mpp_sum_as

   SUBROUTINE mpp_sum_a2s( parr, kdim, kcom )      ! Dummy routine
      REAL   , DIMENSION(:,:) :: parr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_a2s: You should not have seen this print! error?', kdim, parr(1,1), kcom
   END SUBROUTINE mpp_sum_a2s

   SUBROUTINE mpp_sum_ai( karr, kdim, kcom )      ! Dummy routine
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_ai: You should not have seen this print! error?', kdim, karr(1), kcom
   END SUBROUTINE mpp_sum_ai

   SUBROUTINE mpp_sum_s( psca, kcom )            ! Dummy routine
      REAL                  :: psca
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_s: You should not have seen this print! error?', psca, kcom
   END SUBROUTINE mpp_sum_s

   SUBROUTINE mpp_sum_i( kint, kcom )            ! Dummy routine
      integer               :: kint
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_i: You should not have seen this print! error?', kint, kcom
   END SUBROUTINE mpp_sum_i

   SUBROUTINE mppsum_realdd( ytab, kcom )
      COMPLEX(wp), INTENT(inout)         :: ytab    ! input scalar
      INTEGER , INTENT( in  ), OPTIONAL :: kcom
      WRITE(*,*) 'mppsum_realdd: You should not have seen this print! error?', ytab
   END SUBROUTINE mppsum_realdd

   SUBROUTINE mppsum_a_realdd( ytab, kdim, kcom )
      INTEGER , INTENT( in )                        ::   kdim      ! size of ytab
      COMPLEX(wp), DIMENSION(kdim), INTENT( inout ) ::   ytab      ! input array
      INTEGER , INTENT( in  ), OPTIONAL :: kcom
      WRITE(*,*) 'mppsum_a_realdd: You should not have seen this print! error?', kdim, ytab(1), kcom
   END SUBROUTINE mppsum_a_realdd

   SUBROUTINE mppmax_a_real( parr, kdim, kcom )
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmax_a_real: You should not have seen this print! error?', kdim, parr(1), kcom
   END SUBROUTINE mppmax_a_real

   SUBROUTINE mppmax_real( psca, kcom )
      REAL                  :: psca
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmax_real: You should not have seen this print! error?', psca, kcom
   END SUBROUTINE mppmax_real

   SUBROUTINE mppmin_a_real( parr, kdim, kcom )
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmin_a_real: You should not have seen this print! error?', kdim, parr(1), kcom
   END SUBROUTINE mppmin_a_real

   SUBROUTINE mppmin_real( psca, kcom )
      REAL                  :: psca
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmin_real: You should not have seen this print! error?', psca, kcom
   END SUBROUTINE mppmin_real

   SUBROUTINE mppmax_a_int( karr, kdim ,kcom)
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmax_a_int: You should not have seen this print! error?', kdim, karr(1), kcom
   END SUBROUTINE mppmax_a_int

   SUBROUTINE mppmax_int( kint, kcom)
      INTEGER               :: kint
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmax_int: You should not have seen this print! error?', kint, kcom
   END SUBROUTINE mppmax_int

   SUBROUTINE mppmin_a_int( karr, kdim, kcom )
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmin_a_int: You should not have seen this print! error?', kdim, karr(1), kcom
   END SUBROUTINE mppmin_a_int

   SUBROUTINE mppmin_int( kint, kcom )
      INTEGER               :: kint
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmin_int: You should not have seen this print! error?', kint, kcom
   END SUBROUTINE mppmin_int

   SUBROUTINE mpp_minloc2d( ptab, pmask, pmin, ki, kj )
      REAL                   :: pmin
      REAL , DIMENSION (:,:) :: ptab, pmask
      INTEGER :: ki, kj
      WRITE(*,*) 'mpp_minloc2d: You should not have seen this print! error?', pmin, ki, kj, ptab(1,1), pmask(1,1)
   END SUBROUTINE mpp_minloc2d

   SUBROUTINE mpp_minloc3d( ptab, pmask, pmin, ki, kj, kk )
      REAL                     :: pmin
      REAL , DIMENSION (:,:,:) :: ptab, pmask
      INTEGER :: ki, kj, kk
      WRITE(*,*) 'mpp_minloc3d: You should not have seen this print! error?', pmin, ki, kj, kk, ptab(1,1,1), pmask(1,1,1)
   END SUBROUTINE mpp_minloc3d

   SUBROUTINE mpp_maxloc2d( ptab, pmask, pmax, ki, kj )
      REAL                   :: pmax
      REAL , DIMENSION (:,:) :: ptab, pmask
      INTEGER :: ki, kj
      WRITE(*,*) 'mpp_maxloc2d: You should not have seen this print! error?', pmax, ki, kj, ptab(1,1), pmask(1,1)
   END SUBROUTINE mpp_maxloc2d

   SUBROUTINE mpp_maxloc3d( ptab, pmask, pmax, ki, kj, kk )
      REAL                     :: pmax
      REAL , DIMENSION (:,:,:) :: ptab, pmask
      INTEGER :: ki, kj, kk
      WRITE(*,*) 'mpp_maxloc3d: You should not have seen this print! error?', pmax, ki, kj, kk, ptab(1,1,1), pmask(1,1,1)
   END SUBROUTINE mpp_maxloc3d

   SUBROUTINE mppstop
      STOP      ! non MPP case, just stop the run
   END SUBROUTINE mppstop

   SUBROUTINE mpp_ini_ice( kcom, knum )
      INTEGER :: kcom, knum
      WRITE(*,*) 'mpp_ini_ice: You should not have seen this print! error?', kcom, knum
   END SUBROUTINE mpp_ini_ice

   SUBROUTINE mpp_ini_znl( knum )
      INTEGER :: knum
      WRITE(*,*) 'mpp_ini_znl: You should not have seen this print! error?', knum
   END SUBROUTINE mpp_ini_znl

   SUBROUTINE mpp_comm_free( kcom )
      INTEGER :: kcom
      WRITE(*,*) 'mpp_comm_free: You should not have seen this print! error?', kcom
   END SUBROUTINE mpp_comm_free

   !!----------------------------------------------------------------------
   !!   All cases:         ctl_stop, ctl_warn, get_unit, ctl_opn, ctl_nam   routines
   !!----------------------------------------------------------------------

   SUBROUTINE ctl_stop( cd1, cd2, cd3, cd4, cd5 ,   &
      &                 cd6, cd7, cd8, cd9, cd10 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE  stop_opa  ***
      !!
      !! ** Purpose :   print in ocean.outpput file a error message and
      !!                increment the error number (nstop) by one.
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
      !!----------------------------------------------------------------------
      !
      nstop = nstop + 1
      IF(lwp) THEN
         WRITE(numout,cform_err)
         IF( PRESENT(cd1 ) )   WRITE(numout,*) cd1
         IF( PRESENT(cd2 ) )   WRITE(numout,*) cd2
         IF( PRESENT(cd3 ) )   WRITE(numout,*) cd3
         IF( PRESENT(cd4 ) )   WRITE(numout,*) cd4
         IF( PRESENT(cd5 ) )   WRITE(numout,*) cd5
         IF( PRESENT(cd6 ) )   WRITE(numout,*) cd6
         IF( PRESENT(cd7 ) )   WRITE(numout,*) cd7
         IF( PRESENT(cd8 ) )   WRITE(numout,*) cd8
         IF( PRESENT(cd9 ) )   WRITE(numout,*) cd9
         IF( PRESENT(cd10) )   WRITE(numout,*) cd10
      ENDIF
                               CALL FLUSH(numout    )
      IF( numstp     /= -1 )   CALL FLUSH(numstp    )
      IF( numsol     /= -1 )   CALL FLUSH(numsol    )
      IF( numevo_ice /= -1 )   CALL FLUSH(numevo_ice)
      !
      IF( cd1 == 'STOP' ) THEN
         IF(lwp) WRITE(numout,*)  'huge E-R-R-O-R : immediate stop'
         CALL mppstop()
      ENDIF
      !
   END SUBROUTINE ctl_stop


   SUBROUTINE ctl_warn( cd1, cd2, cd3, cd4, cd5,   &
      &                 cd6, cd7, cd8, cd9, cd10 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE  stop_warn  ***
      !!
      !! ** Purpose :   print in ocean.outpput file a error message and
      !!                increment the warning number (nwarn) by one.
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
      !!----------------------------------------------------------------------
      !
      nwarn = nwarn + 1
      IF(lwp) THEN
         WRITE(numout,cform_war)
         IF( PRESENT(cd1 ) ) WRITE(numout,*) cd1
         IF( PRESENT(cd2 ) ) WRITE(numout,*) cd2
         IF( PRESENT(cd3 ) ) WRITE(numout,*) cd3
         IF( PRESENT(cd4 ) ) WRITE(numout,*) cd4
         IF( PRESENT(cd5 ) ) WRITE(numout,*) cd5
         IF( PRESENT(cd6 ) ) WRITE(numout,*) cd6
         IF( PRESENT(cd7 ) ) WRITE(numout,*) cd7
         IF( PRESENT(cd8 ) ) WRITE(numout,*) cd8
         IF( PRESENT(cd9 ) ) WRITE(numout,*) cd9
         IF( PRESENT(cd10) ) WRITE(numout,*) cd10
      ENDIF
      CALL FLUSH(numout)
      !
   END SUBROUTINE ctl_warn


   SUBROUTINE ctl_opn( knum, cdfile, cdstat, cdform, cdacce, klengh, kout, ldwp, karea )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_opn  ***
      !!
      !! ** Purpose :   Open file and check if required file is available.
      !!
      !! ** Method  :   Fortan open
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(  out) ::   knum      ! logical unit to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdfile    ! file name to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdstat    ! disposition specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdform    ! formatting specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdacce    ! access specifier
      INTEGER          , INTENT(in   ) ::   klengh    ! record length
      INTEGER          , INTENT(in   ) ::   kout      ! number of logical units for write
      LOGICAL          , INTENT(in   ) ::   ldwp      ! boolean term for print
      INTEGER, OPTIONAL, INTENT(in   ) ::   karea     ! proc number
      !!
      CHARACTER(len=80) ::   clfile
      INTEGER           ::   iost
      !!----------------------------------------------------------------------

      ! adapt filename
      ! ----------------
      clfile = TRIM(cdfile)
      IF( PRESENT( karea ) ) THEN
         IF( karea > 1 )   WRITE(clfile, "(a,'_',i4.4)") TRIM(clfile), karea-1
      ENDIF
      knum=get_unit()

      iost=0
      IF( cdacce(1:6) == 'DIRECT' )  THEN
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat, RECL=klengh, ERR=100, IOSTAT=iost )
      ELSE
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat             , ERR=100, IOSTAT=iost )
      ENDIF
      IF( iost == 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*) '     file   : ', clfile,' open ok'
            WRITE(kout,*) '     unit   = ', knum
            WRITE(kout,*) '     status = ', cdstat
            WRITE(kout,*) '     form   = ', cdform
            WRITE(kout,*) '     access = ', cdacce
            WRITE(kout,*)
         ENDIF
      ENDIF
100   CONTINUE
      IF( iost /= 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*)
            WRITE(kout,*) ' ===>>>> : bad opening file: ', clfile
            WRITE(kout,*) ' =======   ===  '
            WRITE(kout,*) '           unit   = ', knum
            WRITE(kout,*) '           status = ', cdstat
            WRITE(kout,*) '           form   = ', cdform
            WRITE(kout,*) '           access = ', cdacce
            WRITE(kout,*) '           iostat = ', iost
            WRITE(kout,*) '           we stop. verify the file '
            WRITE(kout,*)
         ENDIF
         STOP 'ctl_opn bad opening'
      ENDIF

   END SUBROUTINE ctl_opn

   SUBROUTINE ctl_nam ( kios, cdnam, ldwp )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_nam  ***
      !!
      !! ** Purpose :   Informations when error while reading a namelist
      !!
      !! ** Method  :   Fortan open
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(inout) ::   kios      ! IO status after reading the namelist
      CHARACTER(len=*) , INTENT(in   ) ::   cdnam     ! group name of namelist for which error occurs
      CHARACTER(len=4)                 ::   clios     ! string to convert iostat in character for print
      LOGICAL          , INTENT(in   ) ::   ldwp      ! boolean term for print
      !!----------------------------------------------------------------------

      ! 
      ! ----------------
      WRITE (clios, '(I4.0)') kios
      IF( kios < 0 ) THEN         
         CALL ctl_warn( 'W A R N I N G:  end of record or file while reading namelist ' &
 &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF

      IF( kios > 0 ) THEN
         CALL ctl_stop( 'E R R O R :   misspelled variable in namelist ' &
 &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF
      kios = 0
      RETURN
      
   END SUBROUTINE ctl_nam

   INTEGER FUNCTION get_unit()
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  get_unit  ***
      !!
      !! ** Purpose :   return the index of an unused logical unit
      !!----------------------------------------------------------------------
      LOGICAL :: llopn
      !!----------------------------------------------------------------------
      !
      get_unit = 15   ! choose a unit that is big enough then it is not already used in NEMO
      llopn = .TRUE.
      DO WHILE( (get_unit < 998) .AND. llopn )
         get_unit = get_unit + 1
         INQUIRE( unit = get_unit, opened = llopn )
      END DO
      IF( (get_unit == 999) .AND. llopn ) THEN
         CALL ctl_stop( 'get_unit: All logical units until 999 are used...' )
         get_unit = -1
      ENDIF
      !
   END FUNCTION get_unit

   !!----------------------------------------------------------------------
END MODULE lib_mpp
