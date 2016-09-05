MODULE bdytides
   !!======================================================================
   !!                       ***  MODULE  bdytides  ***
   !! Ocean dynamics:   Tidal forcing at open boundaries
   !!======================================================================
   !! History :  2.0  !  2007-01  (D.Storkey)  Original code
   !!            2.3  !  2008-01  (J.Holt)  Add date correction. Origins POLCOMS v6.3 2007
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (D.Storkey and E.O'Dea)  bug fixes
   !!            3.4  !  2012-09  (G. Reffray and J. Chanut) New inputs + mods
   !!            3.5  !  2013-07  (J. Chanut) Compliant with time splitting changes
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module         NO Unstruct Open Boundary Conditions for tides
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdytide_init             ! Empty routine
      WRITE(*,*) 'bdytide_init: You should not have seen this print! error?'
   END SUBROUTINE bdytide_init
   SUBROUTINE bdytide_update( kt, jit )   ! Empty routine
      WRITE(*,*) 'bdytide_update: You should not have seen this print! error?', kt, jit
   END SUBROUTINE bdytide_update
   SUBROUTINE bdy_dta_tides( kt, kit, time_offset )     ! Empty routine
      INTEGER, INTENT( in )            ::   kt          ! Dummy argument empty routine      
      INTEGER, INTENT( in ),OPTIONAL   ::   kit         ! Dummy argument empty routine
      INTEGER, INTENT( in ),OPTIONAL   ::   time_offset ! Dummy argument empty routine
      WRITE(*,*) 'bdy_dta_tides: You should not have seen this print! error?', kt, jit
   END SUBROUTINE bdy_dta_tides

   !!======================================================================
END MODULE bdytides

