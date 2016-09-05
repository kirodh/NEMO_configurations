MODULE dynnxt
   !!=========================================================================
   !!                       ***  MODULE  dynnxt  ***
   !! Ocean dynamics: time stepping
   !!=========================================================================
   !! History :  OPA  !  1987-02  (P. Andrich, D. L Hostis)  Original code
   !!                 !  1990-10  (C. Levy, G. Madec)
   !!            7.0  !  1993-03  (M. Guyon)  symetrical conditions
   !!            8.0  !  1997-02  (G. Madec & M. Imbard)  opa, release 8.0
   !!            8.2  !  1997-04  (A. Weaver)  Euler forward step
   !!             -   !  1997-06  (G. Madec)  lateral boudary cond., lbc routine
   !!    NEMO    1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-10  (C. Talandier, A-M. Treguier) Open boundary cond.
   !!            2.0  !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.3  !  2007-07  (D. Storkey) Calls to BDY routines. 
   !!            3.2  !  2009-06  (G. Madec, R.Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-09  (D. Storkey, E.O'Dea) Bug fix for BDY module
   !!            3.3  !  2011-03  (P. Oddo) Bug fix for time-splitting+(BDY-OBC) and not VVL
   !!            3.5  !  2013-07  (J. Chanut) Compliant with time splitting changes
   !!            3.7  !  2014-04  (G. Madec) add the diagnostic of the time filter trends
   !!-------------------------------------------------------------------------
  
   !!-------------------------------------------------------------------------
   !!   dyn_nxt      : obtain the next (after) horizontal velocity
   !!-------------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE phycst          ! physical constants
   USE dynspg_oce      ! type of surface pressure gradient
   USE dynadv          ! dynamics: vector invariant versus flux form
   USE domvvl          ! variable volume
   USE bdy_oce         ! ocean open boundary conditions
   USE bdydta          ! ocean open boundary conditions
   USE bdydyn          ! ocean open boundary conditions
   USE bdyvol          ! ocean open boundary condition (bdy_vol routines)
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   USE trdken          ! trend manager: kinetic energy
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O manager library
   USE lbclnk          ! lateral boundary condition (or mpp link)
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE prtctl          ! Print control
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC    dyn_nxt   ! routine called by step.F90

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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: dynnxt.F90 5628 2015-07-22 20:26:35Z mathiot $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_nxt ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_nxt  ***
      !!                   
      !! ** Purpose :   Compute the after horizontal velocity. Apply the boundary 
      !!             condition on the after velocity, achieved the time stepping 
      !!             by applying the Asselin filter on now fields and swapping 
      !!             the fields.
      !!
      !! ** Method  : * After velocity is compute using a leap-frog scheme:
      !!                       (ua,va) = (ub,vb) + 2 rdt (ua,va)
      !!             Note that with flux form advection and variable volume layer
      !!             (lk_vvl=T), the leap-frog is applied on thickness weighted
      !!             velocity.
      !!             Note also that in filtered free surface (lk_dynspg_flt=T),
      !!             the time stepping has already been done in dynspg module
      !!
      !!              * Apply lateral boundary conditions on after velocity 
      !!             at the local domain boundaries through lbc_lnk call,
      !!             at the one-way open boundaries (lk_bdy=T),
      !!             at the AGRIF zoom   boundaries (lk_agrif=T)
      !!
      !!              * Apply the time filter applied and swap of the dynamics
      !!             arrays to start the next time step:
      !!                (ub,vb) = (un,vn) + atfp [ (ub,vb) + (ua,va) - 2 (un,vn) ]
      !!                (un,vn) = (ua,va).
      !!             Note that with flux form advection and variable volume layer
      !!             (lk_vvl=T), the time filter is applied on thickness weighted
      !!             velocity.
      !!
      !! ** Action :   ub,vb   filtered before horizontal velocity of next time-step
      !!               un,vn   now horizontal velocity of next time-step
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   iku, ikv     ! local integers
      REAL(wp) ::   z2dt         ! temporary scalar
      REAL(wp) ::   zue3a, zue3n, zue3b, zuf, zec      ! local scalars
      REAL(wp) ::   zve3a, zve3n, zve3b, zvf, z1_2dt   !   -      -
      REAL(wp), POINTER, DIMENSION(:,:)   ::  zue, zve
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ze3u_f, ze3v_f, zua, zva 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dyn_nxt')
      !
      CALL wrk_alloc( jpi,jpj,jpk,  ze3u_f, ze3v_f, zua, zva )
      IF( lk_dynspg_ts )   CALL wrk_alloc( jpi,jpj, zue, zve )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_nxt : time stepping'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF




      ! Update after velocity on domain lateral boundaries
      ! --------------------------------------------------      
      CALL lbc_lnk( ua, 'U', -1. )     !* local domain boundaries
      CALL lbc_lnk( va, 'V', -1. ) 
      !
      !

      IF( l_trddyn ) THEN             ! prepare the atf trend computation + some diagnostics
         z1_2dt = 1._wp / (2. * rdt)        ! Euler or leap-frog time step 
         IF( neuler == 0 .AND. kt == nit000 )   z1_2dt = 1._wp / rdt
         !
         !                                  ! Kinetic energy and Conversion
         IF( ln_KE_trd  )   CALL trd_dyn( ua, va, jpdyn_ken, kt )
         !
         IF( ln_dyn_trd ) THEN              ! 3D output: total momentum trends
            zua(:,:,:) = ( ua(:,:,:) - ub(:,:,:) ) * z1_2dt
            zva(:,:,:) = ( va(:,:,:) - vb(:,:,:) ) * z1_2dt
            CALL iom_put( "utrd_tot", zua )        ! total momentum trends, except the asselin time filter
            CALL iom_put( "vtrd_tot", zva )
         ENDIF
         !
         zua(:,:,:) = un(:,:,:)             ! save the now velocity before the asselin filter
         zva(:,:,:) = vn(:,:,:)             ! (caution: there will be a shift by 1 timestep in the
         !                                  !  computation of the asselin filter trends)
      ENDIF

      ! Time filter and swap of dynamics arrays
      ! ------------------------------------------
      IF( neuler == 0 .AND. kt == nit000 ) THEN        !* Euler at first time-step: only swap
         DO jk = 1, jpkm1
            un(:,:,jk) = ua(:,:,jk)                          ! un <-- ua
            vn(:,:,jk) = va(:,:,jk)
         END DO
         IF (lk_vvl) THEN
            DO jk = 1, jpkm1
               e3t_0(:,:,jk) = e3t_0(:,:,jk)
               e3u_0(:,:,jk) = e3u_0(:,:,jk)
               e3v_0(:,:,jk) = e3v_0(:,:,jk)
            ENDDO
         ENDIF
      ELSE                                             !* Leap-Frog : Asselin filter and swap
         !                                ! =============!
         IF( .NOT. lk_vvl ) THEN          ! Fixed volume !
            !                             ! =============!
            DO jk = 1, jpkm1                              
               DO jj = 1, jpj
                  DO ji = 1, jpi    
                     zuf = un(ji,jj,jk) + atfp * ( ub(ji,jj,jk) - 2._wp * un(ji,jj,jk) + ua(ji,jj,jk) )
                     zvf = vn(ji,jj,jk) + atfp * ( vb(ji,jj,jk) - 2._wp * vn(ji,jj,jk) + va(ji,jj,jk) )
                     !
                     ub(ji,jj,jk) = zuf                      ! ub <-- filtered velocity
                     vb(ji,jj,jk) = zvf
                     un(ji,jj,jk) = ua(ji,jj,jk)             ! un <-- ua
                     vn(ji,jj,jk) = va(ji,jj,jk)
                  END DO
               END DO
            END DO
            !                             ! ================!
         ELSE                             ! Variable volume !
            !                             ! ================!
            ! Before scale factor at t-points
            ! (used as a now filtered scale factor until the swap)
            ! ----------------------------------------------------
            IF (lk_dynspg_ts.AND.ln_bt_fw) THEN
               ! No asselin filtering on thicknesses if forward time splitting
                  e3t_0(:,:,:) = e3t_0(:,:,:)
            ELSE
               e3t_0(:,:,:) = e3t_0(:,:,:) + atfp * ( e3t_0(:,:,:) - 2._wp * e3t_0(:,:,:) + e3t_0(:,:,:) )
               ! Add volume filter correction: compatibility with tracer advection scheme
               ! => time filter + conservation correction (only at the first level)
               IF ( nn_isf == 0) THEN   ! if no ice shelf melting
                  e3t_0(:,:,1) = e3t_0(:,:,1) - atfp * rdt * r1_rau0 * ( emp_b(:,:) - emp(:,:) &
                                 &                                          -rnf_b(:,:) + rnf(:,:) ) * tmask(:,:,1)
               ELSE                     ! if ice shelf melting
                  DO jj = 1,jpj
                     DO ji = 1,jpi
                        jk = mikt(ji,jj)
                        e3t_0(ji,jj,jk) = e3t_0(ji,jj,jk) - atfp * rdt * r1_rau0                       &
                                          &                          * ( (emp_b(ji,jj)    - emp(ji,jj)   ) &
                                          &                            - (rnf_b(ji,jj)    - rnf(ji,jj)   ) &
                                          &                            + (fwfisf_b(ji,jj) - fwfisf(ji,jj)) ) * tmask(ji,jj,jk)
                     END DO
                  END DO
               END IF
            ENDIF
            !
            IF( ln_dynadv_vec ) THEN
               ! Before scale factor at (u/v)-points
               ! -----------------------------------
               CALL dom_vvl_interpol( e3t_0(:,:,:), e3u_0(:,:,:), 'U' )
               CALL dom_vvl_interpol( e3t_0(:,:,:), e3v_0(:,:,:), 'V' )
               ! Leap-Frog - Asselin filter and swap: applied on velocity
               ! -----------------------------------
               DO jk = 1, jpkm1
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        zuf = un(ji,jj,jk) + atfp * ( ub(ji,jj,jk) - 2._wp * un(ji,jj,jk) + ua(ji,jj,jk) )
                        zvf = vn(ji,jj,jk) + atfp * ( vb(ji,jj,jk) - 2._wp * vn(ji,jj,jk) + va(ji,jj,jk) )
                        !
                        ub(ji,jj,jk) = zuf                      ! ub <-- filtered velocity
                        vb(ji,jj,jk) = zvf
                        un(ji,jj,jk) = ua(ji,jj,jk)             ! un <-- ua
                        vn(ji,jj,jk) = va(ji,jj,jk)
                     END DO
                  END DO
               END DO
               !
            ELSE
               ! Temporary filtered scale factor at (u/v)-points (will become before scale factor)
               !------------------------------------------------
               CALL dom_vvl_interpol( e3t_0(:,:,:), ze3u_f, 'U' )
               CALL dom_vvl_interpol( e3t_0(:,:,:), ze3v_f, 'V' )
               ! Leap-Frog - Asselin filter and swap: applied on thickness weighted velocity
               ! -----------------------------------             ===========================
               DO jk = 1, jpkm1
                  DO jj = 1, jpj
                     DO ji = 1, jpi                  
                        zue3a = ua(ji,jj,jk) * e3u_0(ji,jj,jk)
                        zve3a = va(ji,jj,jk) * e3v_0(ji,jj,jk)
                        zue3n = un(ji,jj,jk) * e3u_0(ji,jj,jk)
                        zve3n = vn(ji,jj,jk) * e3v_0(ji,jj,jk)
                        zue3b = ub(ji,jj,jk) * e3u_0(ji,jj,jk)
                        zve3b = vb(ji,jj,jk) * e3v_0(ji,jj,jk)
                        !
                        zuf = ( zue3n + atfp * ( zue3b - 2._wp * zue3n  + zue3a ) ) / ze3u_f(ji,jj,jk)
                        zvf = ( zve3n + atfp * ( zve3b - 2._wp * zve3n  + zve3a ) ) / ze3v_f(ji,jj,jk)
                        !
                        ub(ji,jj,jk) = zuf                     ! ub <-- filtered velocity
                        vb(ji,jj,jk) = zvf
                        un(ji,jj,jk) = ua(ji,jj,jk)            ! un <-- ua
                        vn(ji,jj,jk) = va(ji,jj,jk)
                     END DO
                  END DO
               END DO
               e3u_0(:,:,1:jpkm1) = ze3u_f(:,:,1:jpkm1)      ! e3u_b <-- filtered scale factor
               e3v_0(:,:,1:jpkm1) = ze3v_f(:,:,1:jpkm1)
            ENDIF
            !
         ENDIF
         !
         IF (lk_dynspg_ts.AND.ln_bt_fw) THEN
            ! Revert "before" velocities to time split estimate
            ! Doing it here also means that asselin filter contribution is removed  
            zue(:,:) = e3u_0(:,:,1) * ub(:,:,1) * umask(:,:,1)
            zve(:,:) = e3v_0(:,:,1) * vb(:,:,1) * vmask(:,:,1)    
            DO jk = 2, jpkm1
               zue(:,:) = zue(:,:) + e3u_0(:,:,jk) * ub(:,:,jk) * umask(:,:,jk)
               zve(:,:) = zve(:,:) + e3v_0(:,:,jk) * vb(:,:,jk) * vmask(:,:,jk)    
            END DO
            DO jk = 1, jpkm1
               ub(:,:,jk) = ub(:,:,jk) - (zue(:,:) * hur(:,:) - un_b(:,:)) * umask(:,:,jk)
               vb(:,:,jk) = vb(:,:,jk) - (zve(:,:) * hvr(:,:) - vn_b(:,:)) * vmask(:,:,jk)
            END DO
         ENDIF
         !
      ENDIF ! neuler =/0
      !
      ! Set "now" and "before" barotropic velocities for next time step:
      ! JC: Would be more clever to swap variables than to make a full vertical
      ! integration
      !
      !
      IF (lk_vvl) THEN
         hu(:,:) = 0.
         hv(:,:) = 0.
         DO jk = 1, jpkm1
            hu(:,:) = hu(:,:) + e3u_0(:,:,jk) * umask(:,:,jk)
            hv(:,:) = hv(:,:) + e3v_0(:,:,jk) * vmask(:,:,jk)
         END DO
         hur(:,:) = umask_i(:,:) / ( hu(:,:) + 1._wp - umask_i(:,:) )
         hvr(:,:) = vmask_i(:,:) / ( hv(:,:) + 1._wp - vmask_i(:,:) )
      ENDIF
      !
      un_b(:,:) = 0._wp ; vn_b(:,:) = 0._wp
      ub_b(:,:) = 0._wp ; vb_b(:,:) = 0._wp
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               un_b(ji,jj) = un_b(ji,jj) + e3u_0(ji,jj,jk) * un(ji,jj,jk) * umask(ji,jj,jk)
               vn_b(ji,jj) = vn_b(ji,jj) + e3v_0(ji,jj,jk) * vn(ji,jj,jk) * vmask(ji,jj,jk)
               !
               ub_b(ji,jj) = ub_b(ji,jj) + e3u_0(ji,jj,jk) * ub(ji,jj,jk) * umask(ji,jj,jk)
               vb_b(ji,jj) = vb_b(ji,jj) + e3v_0(ji,jj,jk) * vb(ji,jj,jk) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      !
      un_b(:,:) = un_b(:,:) * hur(:,:)
      vn_b(:,:) = vn_b(:,:) * hvr(:,:)
      ub_b(:,:) = ub_b(:,:) * hur(:,:)
      vb_b(:,:) = vb_b(:,:) * hvr(:,:)
      !
      !

      IF( l_trddyn ) THEN                ! 3D output: asselin filter trends on momentum
         zua(:,:,:) = ( ub(:,:,:) - zua(:,:,:) ) * z1_2dt
         zva(:,:,:) = ( vb(:,:,:) - zva(:,:,:) ) * z1_2dt
         CALL trd_dyn( zua, zva, jpdyn_atf, kt )
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=un, clinfo1=' nxt  - Un: ', mask1=umask,   &
         &                       tab3d_2=vn, clinfo2=' Vn: '       , mask2=vmask )
      ! 
      CALL wrk_dealloc( jpi,jpj,jpk,  ze3u_f, ze3v_f, zua, zva )
      IF( lk_dynspg_ts )   CALL wrk_dealloc( jpi,jpj, zue, zve )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_nxt')
      !
   END SUBROUTINE dyn_nxt

   !!=========================================================================
END MODULE dynnxt
