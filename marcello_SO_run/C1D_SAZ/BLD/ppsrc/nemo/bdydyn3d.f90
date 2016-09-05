MODULE bdydyn3d
   !!======================================================================
   !!                       ***  MODULE  bdydyn3d  ***
   !! Unstructured Open Boundary Cond. :   Flow relaxation scheme on baroclinic velocities
   !!======================================================================
   !! History :  3.4  !  2011     (D. Storkey) new module as part of BDY rewrite 
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dyn3d( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn3d: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn3d

   SUBROUTINE bdy_dyn3d_dmp( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn3d_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn3d_dmp


   !!======================================================================
END MODULE bdydyn3d
