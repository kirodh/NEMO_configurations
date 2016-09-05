












MODULE zdftmx
   !!========================================================================
   !!                       ***  MODULE  zdftmx  ***
   !! Ocean physics: vertical tidal mixing coefficient
   !!========================================================================
   !! History :  1.0  !  2004-04  (L. Bessieres, G. Madec)  Original code
   !!             -   !  2006-08  (A. Koch-Larrouy) Indonesian strait
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option          Dummy module                NO Tidal MiXing
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdftmx = .FALSE.   !: tidal mixing flag
CONTAINS
   SUBROUTINE zdf_tmx_init           ! Dummy routine
      WRITE(*,*) 'zdf_tmx: You should not have seen this print! error?'
   END SUBROUTINE zdf_tmx_init
   SUBROUTINE zdf_tmx( kt )          ! Dummy routine
      WRITE(*,*) 'zdf_tmx: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_tmx

   !!======================================================================
END MODULE zdftmx
