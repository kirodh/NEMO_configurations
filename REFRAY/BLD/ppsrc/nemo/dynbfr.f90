












MODULE dynbfr
   !!==============================================================================
   !!                 ***  MODULE  dynbfr  ***
   !! Ocean dynamics :  bottom friction component of the momentum mixing trend
   !!==============================================================================
   !! History :  3.2  ! 2008-11  (A. C. Coward)  Original code
   !!            3.4  ! 2011-09  (H. Liu) Make it consistent with semi-implicit
   !!                            Bottom friction (ln_bfrimp = .true.) 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_bfr       : Update the momentum trend with the bottom friction contribution
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables 
   USE zdf_oce        ! ocean vertical physics variables
   USE zdfbfr         ! ocean bottom friction variables
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE timing         ! Timing
   USE wrk_nemo       ! Memory Allocation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_bfr   !  routine called by step.F90

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
   !!                    *** zdfddm_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
!   Defautl option :                     avs = avt
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfddm_substitute.h90 4152 2013-11-05 11:59:53Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
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
   !! $Id: dynbfr.F90 5120 2015-03-03 16:11:55Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE dyn_bfr( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_bfr  ***
      !!
      !! ** Purpose :   compute the bottom friction ocean dynamics physics.
      !!
      !! ** Action  :   (ua,va)   momentum trend increased by bottom friction trend
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !! 
      INTEGER  ::   ji, jj       ! dummy loop indexes
      INTEGER  ::   ikbu, ikbv   ! local integers
      REAL(wp) ::   zm1_2dt      ! local scalar
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_bfr')
      !
!!gm issue: better to put the logical in step to control the call of zdf_bfr
!!          ==> change the logical from ln_bfrimp to ln_bfr_exp !!
      IF( .NOT.ln_bfrimp) THEN     ! only for explicit bottom friction form
                                    ! implicit bfr is implemented in dynzdf_imp

!!gm bug : time step is only rdt (not 2 rdt if euler start !)
        zm1_2dt = - 1._wp / ( 2._wp * rdt )

        IF( l_trddyn )   THEN                      ! temporary save of ua and va trends
           CALL wrk_alloc( jpi,jpj,jpk, ztrdu, ztrdv )
           ztrdu(:,:,:) = ua(:,:,:)
           ztrdv(:,:,:) = va(:,:,:)
        ENDIF


        DO jj = 2, jpjm1
           DO ji = 2, jpim1
              ikbu = mbku(ji,jj)          ! deepest ocean u- & v-levels
              ikbv = mbkv(ji,jj)
              !
              ! Apply stability criteria on absolute value  : abs(bfr/e3) < 1/(2dt) => bfr/e3 > -1/(2dt)
              ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + MAX(  bfrua(ji,jj) / e3u_0(ji,jj,ikbu) , zm1_2dt  ) * ub(ji,jj,ikbu)
              va(ji,jj,ikbv) = va(ji,jj,ikbv) + MAX(  bfrva(ji,jj) / e3v_0(ji,jj,ikbv) , zm1_2dt  ) * vb(ji,jj,ikbv)
           END DO
        END DO
        
        IF ( ln_isfcav ) THEN
           DO jj = 2, jpjm1
              DO ji = 2, jpim1
                 ! (ISF) stability criteria for top friction
                 ikbu = miku(ji,jj)          ! first wet ocean u- & v-levels
                 ikbv = mikv(ji,jj)
                 !
                 ! Apply stability criteria on absolute value  : abs(bfr/e3) < 1/(2dt) => bfr/e3 > -1/(2dt)
                 ua(ji,jj,ikbu) = ua(ji,jj,ikbu) + MAX(  tfrua(ji,jj) / e3u_0(ji,jj,ikbu) , zm1_2dt  ) * ub(ji,jj,ikbu) &
                    &             * (1.-umask(ji,jj,1))
                 va(ji,jj,ikbv) = va(ji,jj,ikbv) + MAX(  tfrva(ji,jj) / e3v_0(ji,jj,ikbv) , zm1_2dt  ) * vb(ji,jj,ikbv) &
                    &             * (1.-vmask(ji,jj,1))
                 ! (ISF)
              END DO
           END DO
        END IF

        !
        IF( l_trddyn )   THEN                      ! save the vertical diffusive trends for further diagnostics
           ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
           ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
           CALL trd_dyn( ztrdu(:,:,:), ztrdv(:,:,:), jpdyn_bfr, kt )
           CALL wrk_dealloc( jpi,jpj,jpk, ztrdu, ztrdv )
        ENDIF
        !                                          ! print mean trends (used for debugging)
        IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' bfr  - Ua: ', mask1=umask,               &
           &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
        !
      ENDIF     ! end explicit bottom friction
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_bfr')
      !
   END SUBROUTINE dyn_bfr

   !!==============================================================================
END MODULE dynbfr
