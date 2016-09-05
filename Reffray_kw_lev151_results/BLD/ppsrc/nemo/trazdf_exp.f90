












MODULE trazdf_exp
   !!==============================================================================
   !!                    ***  MODULE  trazdf_exp  ***
   !! Ocean  tracers:  vertical component of the tracer mixing trend using
   !!                  a split-explicit time-stepping 
   !!==============================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)
   !!                 !  1992-06  (M. Imbard)  correction on tracer trend loops
   !!                 !  1996-01  (G. Madec)  statement function for e3
   !!                 !  1997-05  (G. Madec)  vertical component of isopycnal
   !!                 !  1997-07  (G. Madec)  geopotential diffusion in s-coord
   !!                 !  2000-08  (G. Madec)  double diffusive mixing
   !!   NEMO     1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2004-08  (C. Talandier) New trends organisation
   !!             -   !  2005-11  (G. Madec)  New organisation
   !!            3.0  !  2008-04  (G. Madec)  leap-frog time stepping done in trazdf
   !!            3.3  !  2010-06  (C. Ethe, G. Madec) Merge TRA-TRC
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_zdf_exp  : compute the tracer the vertical diffusion trend using a
   !!                  split-explicit time stepping and provide the after tracer
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers 
   USE dom_oce         ! ocean space and time domain 
   USE domvvl          ! variable volume levels
   USE zdf_oce         ! ocean vertical physics
   USE zdfddm          ! ocean vertical physics: double diffusion
   USE trc_oce         ! share passive tracers/Ocean variables
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_zdf_exp   ! routine called by step.F90

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
   !! $Id: trazdf_exp.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_zdf_exp( kt, kit000, cdtype, p2dt, kn_zdfexp,   &
      &                                ptb , pta      , kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_exp  ***
      !!                   
      !! ** Purpose :   Compute the after tracer fields due to the vertical
      !!      tracer mixing alone, and then due to the whole tracer trend.
      !!
      !! ** Method  : - The after tracer fields due to the vertical diffusion
      !!      of tracers alone is given by:
      !!                zwx = ptb + p2dt difft
      !!      where difft = dz( avt dz(ptb) ) = 1/e3t dk+1( avt/e3w dk(ptb) )
      !!           (if lk_zdfddm=T use avs on salinity and passive tracers instead of avt)
      !!      difft is evaluated with an Euler split-explit scheme using a
      !!      no flux boundary condition at both surface and bottomi boundaries.
      !!      (N.B. bottom condition is applied through the masked field avt).
      !!              - the after tracer fields due to the whole trend is 
      !!      obtained in leap-frog environment by :
      !!          pta = zwx + p2dt pta
      !!              - in case of variable level thickness (lk_vvl=T) the 
      !!     the leap-frog is applied on thickness weighted tracer. That is:
      !!          pta = [ ptb*e3tb + e3tn*( zwx - ptb + p2dt pta ) ] / e3tn
      !!
      !! ** Action : - after tracer fields pta
      !!---------------------------------------------------------------------
      !
      INTEGER                              , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000      ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype      ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt        ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kn_zdfexp   ! number of sub-time step
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt        ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb         ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta         ! tracer trend 
      !
      INTEGER  ::  ji, jj, jk, jn, jl        ! dummy loop indices
      REAL(wp) ::  zlavmr, zave3r, ze3tr     ! local scalars
      REAL(wp) ::  ztra, ze3tb               !   -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwx, zwy
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf_exp')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zwx, zwy ) 
      !

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_zdf_exp : explicit vertical mixing on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      ! Initializations
      ! ---------------
      zlavmr = 1. / float( kn_zdfexp )         ! Local constant
      !
      !
      DO jn = 1, kjpt                          ! loop over tracers
         !
         zwy(:,:, 1 ) = 0.e0     ! surface boundary conditions: no flux
         zwy(:,:,jpk) = 0.e0     ! bottom  boundary conditions: no flux
         !
         zwx(:,:,:)   = ptb(:,:,:,jn)  ! zwx array set to before tracer values

         ! Split-explicit loop  (after tracer due to the vertical diffusion alone)
         ! -------------------
         !
         DO jl = 1, kn_zdfexp
            !                     ! first vertical derivative
            DO jk = 2, jpk
               DO jj = 2, jpjm1 
                  DO ji = 2, jpim1   ! vector opt.
                     zave3r = 1.e0 / e3w_0(ji,jj,jk) 
                     IF( cdtype == 'TRA' .AND. jn == jp_tem ) THEN  ! temperature : use of avt
                        zwy(ji,jj,jk) =   avt(ji,jj,jk) * ( zwx(ji,jj,jk-1) - zwx(ji,jj,jk) ) * zave3r
                     ELSE                                           ! salinity or pass. tracer : use of avs
                        zwy(ji,jj,jk) = avt(ji,jj,jk) * ( zwx(ji,jj,jk-1) - zwx(ji,jj,jk) ) * zave3r
                     END IF
                  END DO
               END DO
            END DO
            !
            DO jk = 1, jpkm1      ! second vertical derivative   ==> tracer at kt+l*2*rdt/nn_zdfexp
               DO jj = 2, jpjm1 
                  DO ji = 2, jpim1   ! vector opt.
                     ze3tr = zlavmr / e3t_0(ji,jj,jk)
                     zwx(ji,jj,jk) = zwx(ji,jj,jk) + p2dt(jk) * ( zwy(ji,jj,jk) - zwy(ji,jj,jk+1) ) * ze3tr
                  END DO
               END DO
            END DO
            !
         END DO

         ! After tracer due to all trends
         ! ------------------------------
         IF( lk_vvl ) THEN          ! variable level thickness : leap-frog on tracer*e3t
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1 
                  DO ji = 2, jpim1   ! vector opt.
                     ze3tb = e3t_0(ji,jj,jk) / e3t_0(ji,jj,jk)                          ! before e3t
                     ztra  = zwx(ji,jj,jk) - ptb(ji,jj,jk,jn) + p2dt(jk) * pta(ji,jj,jk,jn)       ! total trends * 2*rdt 
                     pta(ji,jj,jk,jn) = ( ze3tb * ptb(ji,jj,jk,jn) + ztra ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ELSE                       ! fixed level thickness : leap-frog on tracers
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1 
                  DO ji = 2, jpim1   ! vector opt.
                     pta(ji,jj,jk,jn) = ( zwx(ji,jj,jk) + p2dt(jk) * pta(ji,jj,jk,jn) ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF
         !
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zwx, zwy ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf_exp')
      !
   END SUBROUTINE tra_zdf_exp

   !!==============================================================================
END MODULE trazdf_exp
