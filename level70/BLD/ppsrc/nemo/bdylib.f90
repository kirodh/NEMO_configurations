












MODULE bdylib
   !!======================================================================
   !!                       ***  MODULE  bdylib  ***
   !! Unstructured Open Boundary Cond. :  Library module of generic boundary algorithms.
   !!======================================================================
   !! History :  3.6  !  2013     (D. Storkey) new module
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_orlanski_2d( idx, igrd, phib, phia, phi_ext  )      ! Empty routine
      WRITE(*,*) 'bdy_orlanski_2d: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_orlanski_2d
   SUBROUTINE bdy_orlanski_3d( idx, igrd, phib, phia, phi_ext  )      ! Empty routine
      WRITE(*,*) 'bdy_orlanski_3d: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_orlanski_3d

   !!======================================================================
END MODULE bdylib
