MODULE traldf_bilapg
   !!==============================================================================
   !!                       ***  MODULE  traldf_bilapg  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================
   !! History : 8.0   ! 1997-07  (G. Madec)  Original code
   !!   NEMO    1.0   ! 2002-08  (G. Madec)  F90: Free form and module
   !!           3.3   ! 2010-06  (C. Ethe, G. Madec) Merge TRA-TRC
   !!==============================================================================
   !!----------------------------------------------------------------------
   !!   Dummy module :             NO rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE tra_ldf_bilapg( kt, kit000, cdtype, ptb, pta, kjpt )      ! Empty routine
      INTEGER :: kt, kit000
      CHARACTER(len=3) ::   cdtype
      REAL, DIMENSION(:,:,:,:) ::   ptb, pta
      WRITE(*,*) 'tra_ldf_iso: You should not have seen this print! error?', &
        &         kt, kit000, cdtype, ptb(1,1,1,1), pta(1,1,1,1), kjpt
   END SUBROUTINE tra_ldf_bilapg

   !!==============================================================================
END MODULE traldf_bilapg
