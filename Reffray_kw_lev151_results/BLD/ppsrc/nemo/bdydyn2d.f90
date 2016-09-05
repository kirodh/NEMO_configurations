












MODULE bdydyn2d
   !!======================================================================
   !!                       ***  MODULE  bdydyn  ***
   !! Unstructured Open Boundary Cond. :   Apply boundary conditions to barotropic solution
   !!======================================================================
   !! History :  3.4  !  2011     (D. Storkey) new module as part of BDY rewrite
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!            3.5  !  2013-07  (J. Chanut) Compliant with time splitting changes
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dyn2d( kt )      ! Empty routine
      INTEGER, intent(in) :: kt
      WRITE(*,*) 'bdy_dyn2d: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn2d


   !!======================================================================
END MODULE bdydyn2d

