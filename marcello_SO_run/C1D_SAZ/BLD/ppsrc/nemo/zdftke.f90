MODULE zdftke
   !!======================================================================
   !!                       ***  MODULE  zdftke  ***
   !! Ocean physics:  vertical mixing coefficient computed from the tke 
   !!                 turbulent closure parameterization
   !!=====================================================================
   !! History :  OPA  !  1991-03  (b. blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)   bug fix
   !!            7.1  !  1992-10  (G. Madec)   new mixing length and eav
   !!            7.2  !  1993-03  (M. Guyon)   symetrical conditions
   !!            7.3  !  1994-08  (G. Madec, M. Imbard)  nn_pdl flag
   !!            7.5  !  1996-01  (G. Madec)   s-coordinates
   !!            8.0  !  1997-07  (G. Madec)   lbc
   !!            8.1  !  1999-01  (E. Stretta) new option for the mixing length
   !!  NEMO      1.0  !  2002-06  (G. Madec) add tke_init routine
   !!             -   !  2004-10  (C. Ethe )  1D configuration
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!            3.0  !  2008-05  (C. Ethe,  G.Madec) : update TKE physics:
   !!                 !           - tke penetration (wind steering)
   !!                 !           - suface condition for tke & mixing length
   !!                 !           - Langmuir cells
   !!             -   !  2008-05  (J.-M. Molines, G. Madec)  2D form of avtb
   !!             -   !  2008-06  (G. Madec)  style + DOCTOR name for namelist parameters
   !!             -   !  2008-12  (G. Reffray) stable discretization of the production term 
   !!            3.2  !  2009-06  (G. Madec, S. Masson) TKE restart compatible with key_cpl 
   !!                 !                                + cleaning of the parameters + bugs correction
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.6  !  2014-11  (P. Mathiot) add ice shelf capability
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        NO TKE scheme
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdftke = .FALSE.   !: TKE flag
CONTAINS
   SUBROUTINE zdf_tke_init           ! Dummy routine
   END SUBROUTINE zdf_tke_init
   SUBROUTINE zdf_tke( kt )          ! Dummy routine
      WRITE(*,*) 'zdf_tke: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_tke
   SUBROUTINE tke_rst( kt, cdrw )
     CHARACTER(len=*) ::   cdrw
     WRITE(*,*) 'tke_rst: You should not have seen this print! error?', kt, cdwr
   END SUBROUTINE tke_rst

   !!======================================================================
END MODULE zdftke
