MODULE bdytra
   !!======================================================================
   !!                       ***  MODULE  bdytra  ***
   !! Ocean tracers:   Apply boundary conditions for tracers
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_tra(kt)      ! Empty routine
      WRITE(*,*) 'bdy_tra: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_tra

   SUBROUTINE bdy_tra_dmp(kt)      ! Empty routine
      WRITE(*,*) 'bdy_tra_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_tra_dmp


   !!======================================================================
END MODULE bdytra
