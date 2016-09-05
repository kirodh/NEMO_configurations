












MODULE dynldf_bilapg
   !!======================================================================
   !!                       ***  MODULE  dynldf_bilapg  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!======================================================================
   !! History :  OPA  !  1997-07  (G. Madec)  Original code
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!            2.0  !  2004-08  (C. Talandier) New trends organization
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                         NO rotation of mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE dyn_ldf_bilapg( kt )               ! Dummy routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'dyn_ldf_bilapg: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_ldf_bilapg

   !!======================================================================
END MODULE dynldf_bilapg
