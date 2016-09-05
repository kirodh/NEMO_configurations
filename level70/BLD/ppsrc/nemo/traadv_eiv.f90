












MODULE traadv_eiv
   !!======================================================================
   !!                    ***  MODULE  traadv_eiv  ***
   !! Ocean tracers:  advection trend - eddy induced velocity
   !!======================================================================
   !! History :  1.0  !  2005-11 (G. Madec)  Original code, from traldf and zdf _iso
   !!            3.3  !  2010-05 (C. Ethe, G. Madec)  merge TRC-TRA 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :             No rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE tra_adv_eiv( kt, kit000, pun, pvn, pwn, cdtype )              ! Empty routine
      INTEGER  ::   kt    
      INTEGER  ::   kit000    
      CHARACTER(len=3) ::   cdtype
      REAL, DIMENSION(:,:,:) ::   pun, pvn, pwn
      WRITE(*,*) 'tra_adv_eiv: You should not have seen this print! error?', &
          &  kt, cdtype, pun(1,1,1), pvn(1,1,1), pwn(1,1,1)
   END SUBROUTINE tra_adv_eiv

   !!==============================================================================
END MODULE traadv_eiv
