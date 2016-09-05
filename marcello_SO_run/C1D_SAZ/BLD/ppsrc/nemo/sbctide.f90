MODULE sbctide
   !!======================================================================
   !!                       ***  MODULE  sbctide  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  9.0  !  2007  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE phycst
   USE daymod
   USE dynspg_oce
   USE tideini
   !
   USE iom
   USE in_out_manager  ! I/O units
   USE ioipsl          ! NetCDF IPSL library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PUBLIC

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   pot_astro   !

  !!----------------------------------------------------------------------
  !!   Default case :   Empty module
  !!----------------------------------------------------------------------
  LOGICAL, PUBLIC, PARAMETER ::   lk_tide = .FALSE.
CONTAINS
  SUBROUTINE sbc_tide( kt )      ! Empty routine
    INTEGER         , INTENT(in) ::   kt         ! ocean time-step
    WRITE(*,*) 'sbc_tide: You should not have seen this print! error?', kt
  END SUBROUTINE sbc_tide

  !!======================================================================
END MODULE sbctide
