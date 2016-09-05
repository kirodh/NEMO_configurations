MODULE ldfeiv
   !!======================================================================
   !!                     ***  MODULE  ldfeiv  ***
   !! Ocean physics:  variable eddy induced velocity coefficients
   !!======================================================================
   !! History :  OPA  ! 1999-03  (G. Madec, A. Jouzeau)  Original code
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  Free form, F90
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option                                         Dummy module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ldf_eiv( kt )       ! Empty routine
      INTEGER :: kt
      WRITE(*,*) 'ldf_eiv: You should not have seen this print! error?', kt
   END SUBROUTINE ldf_eiv

   !!======================================================================
END MODULE ldfeiv
