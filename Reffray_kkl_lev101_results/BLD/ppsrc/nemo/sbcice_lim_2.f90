












MODULE sbcice_lim_2
   !!======================================================================
   !!                       ***  MODULE  sbcice_lim_2  ***
   !! Surface module :  update surface ocean boundary condition over ice covered area using LIM sea-ice model
   !! Sea-Ice model  :  LIM-2 Sea ice model time-stepping
   !!======================================================================
   !! History :  1.0   !  06-2006  (G. Madec)  from icestp_2.F90
   !!            3.0   !  08-2008  (S. Masson, E. .... ) coupled interface
   !!            3.3   !  05-2009  (G.Garric) addition of the lim2_evp case
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module      NO LIM 2.0 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE sbc_ice_lim_2 ( kt, ksbc )     ! Dummy routine
      INTEGER, INTENT(in) ::   kt, ksbc    
      WRITE(*,*) 'sbc_ice_lim_2: You should not have seen this print! error?', kt, ksbc
   END SUBROUTINE sbc_ice_lim_2

   !!======================================================================
END MODULE sbcice_lim_2
