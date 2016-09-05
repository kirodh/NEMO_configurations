MODULE ldfslp
   !!======================================================================
   !!                       ***  MODULE  ldfslp  ***
   !! Ocean physics: slopes of neutral surfaces
   !!======================================================================
   !! History :  OPA  ! 1994-12  (G. Madec, M. Imbard)  Original code
   !!            8.0  ! 1997-06  (G. Madec)  optimization, lbc
   !!            8.1  ! 1999-10  (A. Jouzeau)  NEW profile in the mixed layer
   !!   NEMO     1.0  ! 2002-10  (G. Madec)  Free form, F90
   !!             -   ! 2005-10  (A. Beckmann)  correction for s-coordinates
   !!            3.3  ! 2010-10  (G. Nurser, C. Harris, G. Madec)  add Griffies operator
   !!             -   ! 2010-11  (F. Dupond, G. Madec)  bug correction in slopes just below the ML
   !!----------------------------------------------------------------------
   !!------------------------------------------------------------------------
   !!   Dummy module :                 NO Rotation of lateral mixing tensor
   !!------------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_ldfslp = .FALSE.    !: slopes flag
CONTAINS
   SUBROUTINE ldf_slp( kt, prd, pn2 )   ! Dummy routine
      INTEGER, INTENT(in) :: kt
      REAL, DIMENSION(:,:,:), INTENT(in) :: prd, pn2
      WRITE(*,*) 'ldf_slp: You should not have seen this print! error?', kt, prd(1,1,1), pn2(1,1,1)
   END SUBROUTINE ldf_slp
   SUBROUTINE ldf_slp_grif( kt )        ! Dummy routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'ldf_slp_grif: You should not have seen this print! error?', kt
   END SUBROUTINE ldf_slp_grif
   SUBROUTINE ldf_slp_init              ! Dummy routine
   END SUBROUTINE ldf_slp_init

   !!======================================================================
END MODULE ldfslp
