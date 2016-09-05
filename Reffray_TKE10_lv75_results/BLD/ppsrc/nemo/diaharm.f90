












MODULE diaharm 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Harmonic analysis of tidal constituents 
   !!======================================================================
   !! History :  3.1  !  2007  (O. Le Galloudec, J. Chanut)  Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_diaharm = .FALSE.
CONTAINS
   SUBROUTINE dia_harm ( kt )     ! Empty routine
      INTEGER, INTENT( IN ) :: kt  
      WRITE(*,*) 'dia_harm: you should not have seen this print'
   END SUBROUTINE dia_harm

   !!======================================================================
END MODULE diaharm
