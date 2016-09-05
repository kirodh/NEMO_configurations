












MODULE icblbc

   !!======================================================================
   !!                       ***  MODULE  icblbc  ***
   !! Ocean physics:  routines to handle boundary exchanges for icebergs
   !!======================================================================
   !! History : 3.3.1 !  2010-01  (Martin&Adcroft) Original code
   !!            -    !  2011-03  (Madec)          Part conversion to NEMO form
   !!            -    !                            Removal of mapping from another grid
   !!            -    !  2011-04  (Alderson)       Split into separate modules
   !!            -    !  2011-05  (Alderson)       MPP exchanges written based on lib_mpp
   !!            -    !  2011-05  (Alderson)       MPP and single processor boundary
   !!            -    !                            conditions added
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   icb_lbc       : -  Pass icebergs across cyclic boundaries
   !!   icb_lbc_mpp   : -  In MPP pass icebergs from linked list between processors
   !!                      as they advect around
   !!                   -  Lagrangian processes cannot be handled by existing NEMO MPP
   !!                      routines because they do not lie on regular jpi,jpj grids
   !!                   -  Processor exchanges are handled as in lib_mpp whenever icebergs step 
   !!                      across boundary of interior domain (nicbdi-nicbei, nicbdj-nicbej)
   !!                      so that iceberg does not exist in more than one processor
   !!                   -  North fold exchanges controlled by three arrays:
   !!                         nicbflddest - unique processor numbers that current one exchanges with
   !!                         nicbfldproc - processor number that current grid point exchanges with
   !!                         nicbfldpts  - packed i,j point in exchanging processor
   !!----------------------------------------------------------------------

   USE par_oce                             ! ocean parameters
   USE dom_oce                             ! ocean domain
   USE in_out_manager                      ! IO parameters
   USE lib_mpp                             ! MPI code and lk_mpp in particular
   USE icb_oce                             ! define iceberg arrays
   USE icbutl                              ! iceberg utility routines

   IMPLICIT NONE
   PRIVATE


   PUBLIC   icb_lbc
   PUBLIC   icb_lbc_mpp

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2011)
   !! $Id: icblbc.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_lbc()
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE icb_lbc  ***
      !!
      !! ** Purpose :   in non-mpp case need to deal with cyclic conditions
      !!                including north-fold
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER ::   this
      TYPE(point)  , POINTER ::   pt
      INTEGER                ::   iine
      !!----------------------------------------------------------------------

      !! periodic east/west boundaries
      !! =============================

      IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN

         this => first_berg
         DO WHILE( ASSOCIATED(this) )
            pt => this%current_point
            iine = INT( pt%xi + 0.5 )
            IF( iine > mig(nicbei) ) THEN
               pt%xi = ricb_right + MOD(pt%xi, 1._wp ) - 1._wp
            ELSE IF( iine < mig(nicbdi) ) THEN
               pt%xi = ricb_left + MOD(pt%xi, 1._wp )
            ENDIF
            this => this%next
         END DO
         !
      ENDIF

      !! north/south boundaries
      !! ======================
      ! south symmetric
      IF( nperio == 2 )   CALL ctl_stop(' south symmetric condition not implemented for icebergs')
      ! north fold
      IF( nperio == 3 .OR. nperio == 4 .OR. nperio == 5 .OR. nperio == 6 )   CALL icb_lbc_nfld()
      !
   END SUBROUTINE icb_lbc


   SUBROUTINE icb_lbc_nfld()
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE icb_lbc_nfld  ***
      !!
      !! ** Purpose :   single processor north fold exchange
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER ::   this
      TYPE(point)  , POINTER ::   pt
      INTEGER                ::   iine, ijne, ipts
      INTEGER                ::   iiglo, ijglo
      !!----------------------------------------------------------------------
      !
      this => first_berg
      DO WHILE( ASSOCIATED(this) )
         pt => this%current_point
         ijne = INT( pt%yj + 0.5 )
         IF( ijne .GT. mjg(nicbej) ) THEN
            !
            iine = INT( pt%xi + 0.5 )
            ipts  = nicbfldpts (mi1(iine))
            !
            ! moving across the cut line means both position and
            ! velocity must change
            ijglo = INT( ipts/nicbpack )
            iiglo = ipts - nicbpack*ijglo
            pt%xi = iiglo - ( pt%xi - REAL(iine,wp) )
            pt%yj = ijglo - ( pt%yj - REAL(ijne,wp) )
            pt%uvel = -1._wp * pt%uvel
            pt%vvel = -1._wp * pt%vvel
         ENDIF
         this => this%next
      END DO
      !
   END SUBROUTINE icb_lbc_nfld

   !!----------------------------------------------------------------------
   !!   Default case:            Dummy module        share memory computing
   !!----------------------------------------------------------------------
   SUBROUTINE icb_lbc_mpp()
      WRITE(numout,*) 'icb_lbc_mpp: You should not have seen this message!!'
   END SUBROUTINE icb_lbc_mpp


   !!======================================================================
END MODULE icblbc
