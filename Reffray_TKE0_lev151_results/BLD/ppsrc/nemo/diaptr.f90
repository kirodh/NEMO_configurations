












MODULE diaptr
   !!======================================================================
   !!                       ***  MODULE  diaptr  ***
   !! Ocean physics:  Computes meridonal transports and zonal means
   !!=====================================================================
   !! History :  1.0  ! 2003-09  (C. Talandier, G. Madec)  Original code
   !!            2.0  ! 2006-01  (A. Biastoch)  Allow sub-basins computation
   !!            3.2  ! 2010-03  (O. Marti, S. Flavoni) Add fields
   !!            3.3  ! 2010-10  (G. Madec)  dynamical allocation
   !!            3.6  ! 2014-12  (C. Ethe) use of IOM
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_ptr      : Poleward Transport Diagnostics module
   !!   dia_ptr_init : Initialization, namelist read
   !!   ptr_sjk      : "zonal" mean computation of a field - tracer or flux array
   !!   ptr_sj       : "zonal" and vertical sum computation of a "meridional" flux array
   !!                   (Generic interface to ptr_sj_3d, ptr_sj_2d)
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE phycst           ! physical constants
   !
   USE iom              ! IOM library
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! MPP library
   USE timing           ! preformance summary

   IMPLICIT NONE
   PRIVATE

   INTERFACE ptr_sj
      MODULE PROCEDURE ptr_sj_3d, ptr_sj_2d
   END INTERFACE

   PUBLIC   ptr_sj         ! call by tra_ldf & tra_adv routines
   PUBLIC   ptr_sjk        ! 
   PUBLIC   dia_ptr_init   ! call in step module
   PUBLIC   dia_ptr        ! call in step module

   !                                  !!** namelist  namptr  **
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:) ::   htr_adv, htr_ldf   !: Heat TRansports (adv, diff, overturn.)
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:) ::   str_adv, str_ldf   !: Salt TRansports (adv, diff, overturn.)
   

   LOGICAL, PUBLIC ::   ln_diaptr   !  Poleward transport flag (T) or not (F)
   LOGICAL, PUBLIC ::   ln_subbas   !  Atlantic/Pacific/Indian basins calculation
   INTEGER         ::   nptr        ! = 1 (l_subbas=F) or = 5 (glo, atl, pac, ind, ipc) (l_subbas=T) 

   REAL(wp) ::   rc_sv    = 1.e-6_wp   ! conversion from m3/s to Sverdrup
   REAL(wp) ::   rc_pwatt = 1.e-15_wp  ! conversion from W    to PW (further x rau0 x Cp)
   REAL(wp) ::   rc_ggram = 1.e-6_wp   ! conversion from g    to Pg

   CHARACTER(len=3), ALLOCATABLE, SAVE, DIMENSION(:)     :: clsubb
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: btmsk   ! T-point basin interior masks
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:)   :: btm30   ! mask out Southern Ocean (=0 south of 30°S)

   REAL(wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:)     :: p_fval1d
   REAL(wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: p_fval2d


   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                    ***  domzgr_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
   !!      factors depending on the vertical coord. used, using CPP macro.
   !!----------------------------------------------------------------------
   !! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!----------------------------------------------------------------------

! z- or s-coordinate (1D or 3D + no time dependency) use reference in all cases







! This part should be removed one day ...
! ... In that case all occurence of the above statement functions
!     have to be replaced in the code by xxx_n

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 4488 2014-02-06 10:43:09Z rfurner $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop start/end indices with CPP macro
   !!                allow unrolling of do-loop (useful with vector processors)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: vectopt_loop_substitute.h90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: diaptr.F90 5147 2015-03-13 10:01:32Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_ptr( pvtr )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ptr  ***
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::   pvtr   ! j-effective transport
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   zv, zsfc               ! local scalar
      REAL(wp), DIMENSION(jpi,jpj)     ::  z2d   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  z3d   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zmask   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts) ::  zts   ! 3D workspace
      CHARACTER( len = 10 )  :: cl1
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dia_ptr')

      !
      IF( PRESENT( pvtr ) ) THEN
         IF( iom_use("zomsfglo") ) THEN    ! effective MSF
            z3d(1,:,:) = ptr_sjk( pvtr(:,:,:) )  ! zonal cumulative effective transport
            DO jk = 2, jpkm1 
              z3d(1,:,jk) = z3d(1,:,jk-1) + z3d(1,:,jk)   ! effective j-Stream-Function (MSF)
            END DO
            DO ji = 1, jpi
               z3d(ji,:,:) = z3d(1,:,:)
            ENDDO
            cl1 = TRIM('zomsf'//clsubb(1) )
            CALL iom_put( cl1, z3d * rc_sv )
            DO jn = 2, nptr                                    ! by sub-basins
               z3d(1,:,:) =  ptr_sjk( pvtr(:,:,:), btmsk(:,:,jn)*btm30(:,:) ) 
               DO jk = 2, jpkm1 
                  z3d(1,:,jk) = z3d(1,:,jk-1) + z3d(1,:,jk)    ! effective j-Stream-Function (MSF)
               END DO
               DO ji = 1, jpi
                  z3d(ji,:,:) = z3d(1,:,:)
               ENDDO
               cl1 = TRIM('zomsf'//clsubb(jn) )
               CALL iom_put( cl1, z3d * rc_sv )
            END DO
         ENDIF
         !
      ELSE
         !
         IF( iom_use("zotemglo") ) THEN    ! i-mean i-k-surface 
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zsfc = e1t(ji,jj) * e3t_0(ji,jj,jk)
                     zmask(ji,jj,jk)      = tmask(ji,jj,jk)      * zsfc
                     zts(ji,jj,jk,jp_tem) = tsn(ji,jj,jk,jp_tem) * zsfc
                     zts(ji,jj,jk,jp_sal) = tsn(ji,jj,jk,jp_sal) * zsfc
                  ENDDO
               ENDDO
            ENDDO
            DO jn = 1, nptr
               zmask(1,:,:) = ptr_sjk( zmask(:,:,:), btmsk(:,:,jn) )
               cl1 = TRIM('zosrf'//clsubb(jn) )
               CALL iom_put( cl1, zmask )
               !
               z3d(1,:,:) = ptr_sjk( zts(:,:,:,jp_tem), btmsk(:,:,jn) ) &
                  &            / MAX( zmask(1,:,:), 10.e-15 )
               DO ji = 1, jpi
                  z3d(ji,:,:) = z3d(1,:,:)
               ENDDO
               cl1 = TRIM('zotem'//clsubb(jn) )
               CALL iom_put( cl1, z3d )
               !
               z3d(1,:,:) = ptr_sjk( zts(:,:,:,jp_sal), btmsk(:,:,jn) ) &
                  &            / MAX( zmask(1,:,:), 10.e-15 )
               DO ji = 1, jpi
                  z3d(ji,:,:) = z3d(1,:,:)
               ENDDO
               cl1 = TRIM('zosal'//clsubb(jn) )
               CALL iom_put( cl1, z3d )
            END DO
         ENDIF
         !
         !                                ! Advective and diffusive heat and salt transport
         IF( iom_use("sophtadv") .OR. iom_use("sopstadv") ) THEN   
            z2d(1,:) = htr_adv(:) * rc_pwatt        !  (conversion in PW)
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sophtadv'                 
            CALL iom_put( TRIM(cl1), z2d )
            z2d(1,:) = str_adv(:) * rc_ggram        ! (conversion in Gg)
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sopstadv'
            CALL iom_put( TRIM(cl1), z2d )
         ENDIF
         !
         IF( iom_use("sophtldf") .OR. iom_use("sopstldf") ) THEN   
            z2d(1,:) = htr_ldf(:) * rc_pwatt        !  (conversion in PW) 
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sophtldf'
            CALL iom_put( TRIM(cl1), z2d )
            z2d(1,:) = str_ldf(:) * rc_ggram        !  (conversion in Gg)
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sopstldf'
            CALL iom_put( TRIM(cl1), z2d )
         ENDIF
         !
      ENDIF
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_ptr')
      !
   END SUBROUTINE dia_ptr


   SUBROUTINE dia_ptr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ptr_init  ***
      !!                   
      !! ** Purpose :   Initialization, namelist read
      !!----------------------------------------------------------------------
      INTEGER ::  jn           ! local integers
      INTEGER ::  inum, ierr   ! local integers
      INTEGER ::  ios          ! Local integer output status for namelist read
      !!
      NAMELIST/namptr/ ln_diaptr, ln_subbas
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namptr in reference namelist : Poleward transport
      READ  ( numnam_ref, namptr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namptr in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namptr in configuration namelist : Poleward transport
      READ  ( numnam_cfg, namptr, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namptr in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namptr )

      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_ptr_init : poleward transport and msf initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namptr : set ptr parameters'
         WRITE(numout,*) '      Poleward heat & salt transport (T) or not (F)      ln_diaptr  = ', ln_diaptr
         WRITE(numout,*) '      Global (F) or glo/Atl/Pac/Ind/Indo-Pac basins      ln_subbas  = ', ln_subbas
      ENDIF

      IF( ln_diaptr ) THEN  
         !
         IF( ln_subbas ) THEN 
            nptr = 5            ! Global, Atlantic, Pacific, Indian, Indo-Pacific
            ALLOCATE( clsubb(nptr) )
            clsubb(1) = 'glo' ;  clsubb(2) = 'atl'  ;  clsubb(3) = 'pac'  ;  clsubb(4) = 'ind'  ;  clsubb(5) = 'ipc'
         ELSE               
            nptr = 1       ! Global only
            ALLOCATE( clsubb(nptr) )
            clsubb(1) = 'glo' 
         ENDIF

         !                                      ! allocate dia_ptr arrays
         IF( dia_ptr_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dia_ptr_init : unable to allocate arrays' )

         rc_pwatt = rc_pwatt * rau0_rcp          ! conversion from K.s-1 to PetaWatt

         IF( lk_mpp )   CALL mpp_ini_znl( numout )     ! Define MPI communicator for zonal sum

         IF( ln_subbas ) THEN                ! load sub-basin mask
            CALL iom_open( 'subbasins', inum,  ldstop = .FALSE.  )
            CALL iom_get( inum, jpdom_data, 'atlmsk', btmsk(:,:,2) )   ! Atlantic basin
            CALL iom_get( inum, jpdom_data, 'pacmsk', btmsk(:,:,3) )   ! Pacific  basin
            CALL iom_get( inum, jpdom_data, 'indmsk', btmsk(:,:,4) )   ! Indian   basin
            CALL iom_close( inum )
            btmsk(:,:,5) = MAX ( btmsk(:,:,3), btmsk(:,:,4) )          ! Indo-Pacific basin
            WHERE( gphit(:,:) < -30._wp)   ;   btm30(:,:) = 0._wp      ! mask out Southern Ocean
            ELSE WHERE                     ;   btm30(:,:) = ssmask(:,:)
            END WHERE
         ENDIF
   
         btmsk(:,:,1) = tmask_i(:,:)                                   ! global ocean
      
         DO jn = 1, nptr
            btmsk(:,:,jn) = btmsk(:,:,jn) * tmask_i(:,:)               ! interior domain only
         END DO

         ! Initialise arrays to zero because diatpr is called before they are first calculated
         ! Note that this means diagnostics will not be exactly correct when model run is restarted.
         htr_adv(:) = 0._wp  ;  str_adv(:) =  0._wp  
         htr_ldf(:) = 0._wp  ;  str_ldf(:) =  0._wp 
         !
      ENDIF 
      ! 
   END SUBROUTINE dia_ptr_init


   FUNCTION dia_ptr_alloc()
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ptr_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER               ::   dia_ptr_alloc   ! return value
      INTEGER, DIMENSION(3) ::   ierr
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( btmsk(jpi,jpj,nptr) ,           &
         &      htr_adv(jpj) , str_adv(jpj) ,   &
         &      htr_ldf(jpj) , str_ldf(jpj) , STAT=ierr(1)  )
         !
      ALLOCATE( p_fval1d(jpj), p_fval2d(jpj,jpk), Stat=ierr(2))
      !
      ALLOCATE( btm30(jpi,jpj), STAT=ierr(3)  )

         !
      dia_ptr_alloc = MAXVAL( ierr )
      IF(lk_mpp)   CALL mpp_sum( dia_ptr_alloc )
      !
   END FUNCTION dia_ptr_alloc


   FUNCTION ptr_sj_3d( pva, pmsk )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_sj_3d  ***
      !!
      !! ** Purpose :   i-k sum computation of a j-flux array
      !!
      !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
      !!              pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk)       ::   pva   ! mask flux array at V-point
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj), OPTIONAL ::   pmsk   ! Optional 2D basin mask
      !
      INTEGER                  ::   ji, jj, jk   ! dummy loop arguments
      INTEGER                  ::   ijpj         ! ???
      REAL(wp), POINTER, DIMENSION(:) :: p_fval  ! function value
      !!--------------------------------------------------------------------
      !
      p_fval => p_fval1d

      ijpj = jpj
      p_fval(:) = 0._wp
      IF( PRESENT( pmsk ) ) THEN 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! Vector opt.
                  p_fval(jj) = p_fval(jj) + pva(ji,jj,jk) * tmask_i(ji,jj) * pmsk(ji,jj)
               END DO
            END DO
         END DO
      ELSE
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! Vector opt.
                  p_fval(jj) = p_fval(jj) + pva(ji,jj,jk) * tmask_i(ji,jj) 
               END DO
            END DO
         END DO
      ENDIF
      !
   END FUNCTION ptr_sj_3d


   FUNCTION ptr_sj_2d( pva, pmsk )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_sj_2d  ***
      !!
      !! ** Purpose :   "zonal" and vertical sum computation of a i-flux array
      !!
      !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
      !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!----------------------------------------------------------------------
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj)           ::   pva   ! mask flux array at V-point
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj), OPTIONAL ::   pmsk   ! Optional 2D basin mask
      !
      INTEGER                  ::   ji,jj       ! dummy loop arguments
      INTEGER                  ::   ijpj        ! ???
      REAL(wp), POINTER, DIMENSION(:) :: p_fval ! function value
      !!--------------------------------------------------------------------
      ! 
      p_fval => p_fval1d

      ijpj = jpj
      p_fval(:) = 0._wp
      IF( PRESENT( pmsk ) ) THEN 
         DO jj = 2, jpjm1
            DO ji = nldi, nlei   ! No vector optimisation here. Better use a mask ?
               p_fval(jj) = p_fval(jj) + pva(ji,jj) * tmask_i(ji,jj) * pmsk(ji,jj)
            END DO
         END DO
      ELSE
         DO jj = 2, jpjm1
            DO ji = nldi, nlei   ! No vector optimisation here. Better use a mask ?
               p_fval(jj) = p_fval(jj) + pva(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
      ENDIF
      ! 
   END FUNCTION ptr_sj_2d


   FUNCTION ptr_sjk( pta, pmsk )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_sjk  ***
      !!
      !! ** Purpose :   i-sum computation of an array
      !!
      !! ** Method  : - i-sum of pva using the interior 2D vmask (vmask_i).
      !!
      !! ** Action  : - p_fval: i-mean poleward flux of pva
      !!----------------------------------------------------------------------
      !!
      IMPLICIT none
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk)           ::   pta    ! mask flux array at V-point
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj)    , OPTIONAL ::   pmsk   ! Optional 2D basin mask
      !!
      INTEGER                           :: ji, jj, jk ! dummy loop arguments
      REAL(wp), POINTER, DIMENSION(:,:) :: p_fval     ! return function value
      !!--------------------------------------------------------------------
      !
      p_fval => p_fval2d

      p_fval(:,:) = 0._wp
      !
      IF( PRESENT( pmsk ) ) THEN 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
!!gm here, use of tmask_i  ==> no need of loop over nldi, nlei....
               DO ji =  nldi, nlei   ! No vector optimisation here. Better use a mask ?
                  p_fval(jj,jk) = p_fval(jj,jk) + pta(ji,jj,jk) * pmsk(ji,jj)
               END DO
            END DO
         END DO
      ELSE 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji =  nldi, nlei   ! No vector optimisation here. Better use a mask ?
                  p_fval(jj,jk) = p_fval(jj,jk) + pta(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
      END IF
      !
      !
   END FUNCTION ptr_sjk


   !!======================================================================
END MODULE diaptr
