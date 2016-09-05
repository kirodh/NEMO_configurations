












MODULE trabbl
   !!==============================================================================
   !!                       ***  MODULE  trabbl  ***
   !! Ocean physics :  advective and/or diffusive bottom boundary layer scheme
   !!==============================================================================
   !! History :  OPA  ! 1996-06  (L. Mortier)  Original code
   !!            8.0  ! 1997-11  (G. Madec)    Optimization
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  free form + modules
   !!             -   ! 2004-01  (A. de Miranda, G. Madec, J.M. Molines ) add advective bbl
   !!            3.3  ! 2009-11  (G. Madec)  merge trabbl and trabbl_adv + style + optimization
   !!             -   ! 2010-04  (G. Madec)  Campin & Goosse advective bbl
   !!             -   ! 2010-06  (C. Ethe, G. Madec)  merge TRA-TRC
   !!             -   ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!             -   ! 2013-04  (F. Roquet, G. Madec)  use of eosbn2 instead of local hard coded alpha and beta
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                      No bottom boundary layer scheme
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_trabbl = .FALSE.   !: bbl flag
CONTAINS
   SUBROUTINE tra_bbl_init               ! Dummy routine
   END SUBROUTINE tra_bbl_init
   SUBROUTINE tra_bbl( kt )              ! Dummy routine
      WRITE(*,*) 'tra_bbl: You should not have seen this print! error?', kt
   END SUBROUTINE tra_bbl

   !!======================================================================
END MODULE trabbl
