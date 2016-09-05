MODULE traldf_iso
   !!======================================================================
   !!                   ***  MODULE  traldf_iso  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!======================================================================
   !! History :  OPA  !  1994-08  (G. Madec, M. Imbard)
   !!            8.0  !  1997-05  (G. Madec)  split into traldf and trazdf
   !!            NEMO !  2002-08  (G. Madec)  Free form, F90
   !!            1.0  !  2005-11  (G. Madec)  merge traldf and trazdf :-)
   !!            3.3  !  2010-09  (C. Ethe, G. Madec) Merge TRA-TRC
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   default option :   Dummy code   NO rotation of the diffusive tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE tra_ldf_iso( kt, kit000,cdtype, pgu, pgv, pgui, pgvi, ptb, pta, kjpt, pahtb0 )      ! Empty routine
      INTEGER:: kt, kit000
      CHARACTER(len=3) ::   cdtype
      REAL, DIMENSION(:,:,:) ::   pgu, pgv, pgui, pgvi    ! tracer gradient at pstep levels
      REAL, DIMENSION(:,:,:,:) ::   ptb, pta
      WRITE(*,*) 'tra_ldf_iso: You should not have seen this print! error?', kt, kit000, cdtype,   &
         &                       pgu(1,1,1), pgv(1,1,1), ptb(1,1,1,1), pta(1,1,1,1), kjpt, pahtb0
   END SUBROUTINE tra_ldf_iso

   !!==============================================================================
END MODULE traldf_iso
