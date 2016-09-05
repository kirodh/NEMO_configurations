












MODULE bdyice_lim
   !!======================================================================
   !!                       ***  MODULE  bdyice_lim  ***
   !! Unstructured Open Boundary Cond. :  Open boundary conditions for sea-ice (LIM2 and LIM3)
   !!======================================================================
   !!  History :  3.3  !  2010-09 (D. Storkey)  Original code
   !!             3.4  !  2011    (D. Storkey)  rewrite in preparation for OBC-BDY merge
   !!              -   !  2012-01 (C. Rousset)  add lim3 and remove useless jk loop 
   !!----------------------------------------------------------------------
   !!---------------------------------------------------------------------------------
   !!   Default option                                                    Empty module
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_ice_lim( kt )      ! Empty routine
      WRITE(*,*) 'bdy_ice_lim: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_ice_lim

   !!=================================================================================
END MODULE bdyice_lim
