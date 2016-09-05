












MODULE dynspg_flt
   !!======================================================================
   !!                   ***  MODULE  dynspg_flt  ***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
   !! History    OPA  !  1998-05  (G. Roullet)  free surface
   !!                 !  1998-10  (G. Madec, M. Imbard)  release 8.2
   !!   NEMO     O.1  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!            1.0  !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!             -   !  2006-08  (J.Chanut, A.Sellar) Calls to BDY routines. 
   !!            3.2  !  2009-03  (G. Madec, M. Leclair, R. Benshila) introduce sshwzv module
   !!            3.7  !  2014-04  (F. Roquet, G. Madec)  add some trends diag
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module   No standart free surface cst volume
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE dyn_spg_flt( kt, kindic )       ! Empty routine
      WRITE(*,*) 'dyn_spg_flt: You should not have seen this print! error?', kt, kindic
   END SUBROUTINE dyn_spg_flt
   SUBROUTINE flt_rst    ( kt, cdrw )         ! Empty routine
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      WRITE(*,*) 'flt_rst: You should not have seen this print! error?', kt, cdrw
   END SUBROUTINE flt_rst
   
   !!======================================================================
END MODULE dynspg_flt
