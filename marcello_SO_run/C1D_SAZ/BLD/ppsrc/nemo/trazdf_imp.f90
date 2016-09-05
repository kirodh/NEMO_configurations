MODULE trazdf_imp
   !!======================================================================
   !!                 ***  MODULE  trazdf_imp  ***
   !! Ocean  tracers:  vertical component of the tracer mixing trend
   !!======================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)
   !!                 !  1992-06  (M. Imbard) correction on tracer trend loops
   !!                 !  1996-01  (G. Madec) statement function for e3
   !!                 !  1997-05  (G. Madec) vertical component of isopycnal
   !!                 !  1997-07  (G. Madec) geopotential diffusion in s-coord
   !!                 !  2000-08  (G. Madec) double diffusive mixing
   !!   NEMO     1.0  !  2002-08  (G. Madec) F90: Free form and module
   !!            2.0  !  2006-11  (G. Madec) New step reorganisation
   !!            3.2  !  2009-03  (G. Madec)  heat and salt content trends
   !!            3.3  !  2010-06  (C. Ethe, G. Madec) Merge TRA-TRC
   !!             -   !  2011-02  (A. Coward, C. Ethe, G. Madec) improvment of surface boundary condition
   !!----------------------------------------------------------------------
  
   !!----------------------------------------------------------------------
   !!   tra_zdf_imp : Update the tracer trend with the diagonal vertical  
   !!                 part of the mixing tensor.
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE trc_oce         ! share passive tracers/ocean variables
   USE domvvl          ! variable volume
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldftra          ! lateral mixing type
   USE ldfslp          ! lateral physics: slope of diffusion
   USE zdfddm          ! ocean vertical physics: double diffusion
   USE traldf_iso_grif ! active tracers: Griffies operator
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_zdf_imp   !  routine called by step.F90

   REAL(wp) ::  r_vvl     ! variable volume indicator, =1 if lk_vvl=T, =0 otherwise 

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
   !!                    *** ldftra_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_substitute.h90 6312 2016-02-15 11:43:52Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   Default option :             aht: Constant coefficient
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
   !! $Id: trazdf_imp.F90 5120 2015-03-03 16:11:55Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
 
   SUBROUTINE tra_zdf_imp( kt, kit000, cdtype, p2dt, ptb, pta, kjpt ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_imp  ***
      !!
      !! ** Purpose :   Compute the after tracer through a implicit computation
      !!     of the vertical tracer diffusion (including the vertical component 
      !!     of lateral mixing (only for 2nd order operator, for fourth order 
      !!     it is already computed and add to the general trend in traldf) 
      !!
      !! ** Method  :  The vertical diffusion of the tracer t  is given by:
      !!                  difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
      !!      It is computed using a backward time scheme (t=ta).
      !!      If lk_zdfddm=T, use avs for salinity or for passive tracers
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      If iso-neutral mixing, add to avt the contribution due to lateral mixing.
      !!
      !! ** Action  : - pta  becomes the after tracer
      !!---------------------------------------------------------------------
      USE oce     , ONLY:   zwd => ua       , zws => va         ! (ua,va) used as 3D workspace
      !
      INTEGER                              , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt     ! number of tracers
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt     ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb      ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta      ! tracer trend 
      !
      INTEGER  ::  ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::  zrhs, ze3tb, ze3tn, ze3ta   ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwi, zwt
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf_imp')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zwi, zwt ) 
      !
      IF( kt == kit000 )  THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'tra_zdf_imp : implicit vertical mixing on ', cdtype
         IF(lwp)WRITE(numout,*) '~~~~~~~~~~~ '
         !
         IF( lk_vvl ) THEN   ;    r_vvl = 1._wp       ! Variable volume indicator
         ELSE                ;    r_vvl = 0._wp       
         ENDIF
      ENDIF
      !
      !                                               ! ============= !
      DO jn = 1, kjpt                                 !  tracer loop  !
         !                                            ! ============= !
         !
         !  Matrix construction
         ! --------------------
         ! Build matrix if temperature or salinity (only in double diffusion case) or first passive tracer
         !
         IF(  ( cdtype == 'TRA' .AND. ( jn == jp_tem .OR. ( jn == jp_sal .AND. lk_zdfddm ) ) ) .OR.   &
            & ( cdtype == 'TRC' .AND. jn == 1 )  )  THEN
            !
            ! vertical mixing coef.: avt for temperature, avs for salinity and passive tracers
            IF( cdtype == 'TRA' .AND. jn == jp_tem ) THEN   ;   zwt(:,:,2:jpk) = avt  (:,:,2:jpk)
            ELSE                                            ;   zwt(:,:,2:jpk) = avt(:,:,2:jpk)
            ENDIF
            DO jj=1, jpj
               DO ji=1, jpi
                  zwt(ji,jj,1) = 0._wp
               END DO
            END DO
!
            ! Diagonal, lower (i), upper (s)  (including the bottom boundary condition since avt is masked)
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ze3ta =  ( 1. - r_vvl ) +        r_vvl   * e3t_0(ji,jj,jk)   ! after scale factor at T-point
                     ze3tn =         r_vvl   + ( 1. - r_vvl ) * e3t_0(ji,jj,jk)   ! now   scale factor at T-point
                     zwi(ji,jj,jk) = - p2dt(jk) * zwt(ji,jj,jk  ) / ( ze3tn * e3w_0(ji,jj,jk  ) )
                     zws(ji,jj,jk) = - p2dt(jk) * zwt(ji,jj,jk+1) / ( ze3tn * e3w_0(ji,jj,jk+1) )
                     zwd(ji,jj,jk) = ze3ta - zwi(ji,jj,jk) - zws(ji,jj,jk)
                 END DO
               END DO
            END DO
            !
            !! Matrix inversion from the first level
            !!----------------------------------------------------------------------
            !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
            !
            !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
            !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
            !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
            !        (        ...               )( ...  ) ( ...  )
            !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
            !
            !   m is decomposed in the product of an upper and lower triangular matrix.
            !   The 3 diagonal terms are in 3d arrays: zwd, zws, zwi.
            !   Suffices i,s and d indicate "inferior" (below diagonal), diagonal
            !   and "superior" (above diagonal) components of the tridiagonal system.
            !   The solution will be in the 4d array pta.
            !   The 3d array zwt is used as a work space array.
            !   En route to the solution pta is used a to evaluate the rhs and then 
            !   used as a work space array: its value is modified.
            !
            ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
            ! done once for all passive tracers (so included in the IF instruction)
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zwt(ji,jj,1) = zwd(ji,jj,1)
               END DO
            END DO
            DO jk = 2, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwt(ji,jj,jk-1)
                  END DO
               END DO
            END DO
            !
         END IF 
         !         
         ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ze3tb = ( 1. - r_vvl ) + r_vvl * e3t_0(ji,jj,1)
               ze3tn = ( 1. - r_vvl ) + r_vvl * e3t_0(ji,jj,1)
               pta(ji,jj,1,jn) = ze3tb * ptb(ji,jj,1,jn)                     &
                  &                      + p2dt(1) * ze3tn * pta(ji,jj,1,jn)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ze3tb = ( 1. - r_vvl ) + r_vvl * e3t_0(ji,jj,jk)
                  ze3tn = ( 1. - r_vvl ) + r_vvl * e3t_0(ji,jj,jk)
                  zrhs = ze3tb * ptb(ji,jj,jk,jn) + p2dt(jk) * ze3tn * pta(ji,jj,jk,jn)   ! zrhs=right hand side 
                  pta(ji,jj,jk,jn) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) * pta(ji,jj,jk-1,jn)
               END DO
            END DO
         END DO

         ! third recurrence:    Xk = (Zk - Sk Xk+1 ) / Tk   (result is the after tracer)
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               pta(ji,jj,jpkm1,jn) = pta(ji,jj,jpkm1,jn) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
            END DO
         END DO
         DO jk = jpk-2, 1, -1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  pta(ji,jj,jk,jn) = ( pta(ji,jj,jk,jn) - zws(ji,jj,jk) * pta(ji,jj,jk+1,jn) )   &
                     &             / zwt(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !                                            ! ================= !
      END DO                                          !  end tracer loop  !
      !                                               ! ================= !
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zwi, zwt ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf_imp')
      !
   END SUBROUTINE tra_zdf_imp

   !!==============================================================================
END MODULE trazdf_imp
