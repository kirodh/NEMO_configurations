MODULE dynldf_iso
   !!======================================================================
   !!                     ***  MODULE  dynldf_iso  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!======================================================================
   !! History :  OPA  !  97-07  (G. Madec)  Original code
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2004-08  (C. Talandier) New trends organization
   !!            2.0  !  2005-11  (G. Madec)  s-coordinate: horizontal diffusion
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                           NO rotation of mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE dyn_ldf_iso( kt )               ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'dyn_ldf_iso: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_ldf_iso

   !!======================================================================
END MODULE dynldf_iso
