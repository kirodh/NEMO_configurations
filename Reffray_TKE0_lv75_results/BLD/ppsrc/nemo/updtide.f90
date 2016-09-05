












MODULE updtide
   !!======================================================================
   !!                       ***  MODULE  updtide  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  9.0  !  07  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   Dummy module :                                        NO TIDE
  !!----------------------------------------------------------------------
CONTAINS
  SUBROUTINE upd_tide( kt, kit, time_offset )  ! Empty routine
    INTEGER, INTENT(in)           ::   kt      !  integer  arg, dummy routine
    INTEGER, INTENT(in), OPTIONAL ::   kit     !  optional arg, dummy routine
    INTEGER, INTENT(in), OPTIONAL ::   time_offset !  optional arg, dummy routine
    WRITE(*,*) 'upd_tide: You should not have seen this print! error?', kt
  END SUBROUTINE upd_tide


  !!======================================================================

END MODULE updtide
