MODULE traldf_iso_grif
   !!======================================================================
   !!                   ***  MODULE  traldf_iso_grif  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!======================================================================
   !! History : 3.3  ! 2010-10  (G. Nurser, C. Harris, G. Madec)
   !!                !          Griffies operator version adapted from traldf_iso.F90
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   default option :   Dummy code   NO rotation of the diffusive tensor
   !!----------------------------------------------------------------------
   REAL, PUBLIC, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::   psix_eiv, psiy_eiv   !: eiv stream function (diag only)
CONTAINS
   SUBROUTINE tra_ldf_iso_grif( kt, kit000, cdtype, pgu, pgv,              &
       &                                   ptb, pta, kjpt, pahtb0 )
      CHARACTER(len=3) ::   cdtype
      INTEGER          ::   kit000     ! first time step index
      REAL, DIMENSION(:,:,:) ::   pgu, pgv   ! tracer gradient at pstep levels
      REAL, DIMENSION(:,:,:,:) ::   ptb, pta
      WRITE(*,*) 'tra_ldf_iso_grif: You should not have seen this print! error?', kt, cdtype,    &
         &                  pgu(1,1,1), pgv(1,1,1), ptb(1,1,1,1), pta(1,1,1,1), kjpt, pahtb0
   END SUBROUTINE tra_ldf_iso_grif

   !!==============================================================================
END MODULE traldf_iso_grif
