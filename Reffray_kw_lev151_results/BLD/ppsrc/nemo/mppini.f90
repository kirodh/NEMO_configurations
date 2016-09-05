












MODULE mppini
   !!==============================================================================
   !!                       ***  MODULE mppini   ***
   !! Ocean initialization : distributed memory computing initialization
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   mpp_init       : Lay out the global domain over processors
   !!   mpp_init2      : Lay out the global domain over processors 
   !!                    with land processor elimination
   !!   mpp_init_ioispl: IOIPSL initialization in mpp
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain 
   USE in_out_manager  ! I/O Manager
   USE lib_mpp         ! distribued memory computing library
   USE ioipsl

   IMPLICIT NONE
   PRIVATE

   PUBLIC mpp_init       ! called by opa.F90
   PUBLIC mpp_init2      ! called by opa.F90

   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                    ***  domzgr_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
   !!      factors depending on the vertical coord. used, using CPP macro.
   !!----------------------------------------------------------------------
   !! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!----------------------------------------------------------------------

! z- or s-coordinate (1D or 3D + no time dependency) use reference in all cases







! This part should be removed one day ...
! ... In that case all occurence of the above statement functions
!     have to be replaced in the code by xxx_n

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 4488 2014-02-06 10:43:09Z rfurner $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: mppini.F90 6413 2016-03-31 16:22:52Z lovato $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   Default option :                            shared memory computing
   !!----------------------------------------------------------------------

   SUBROUTINE mpp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init  ***
      !!
      !! ** Purpose :   Lay out the global domain over processors.
      !!
      !! ** Method  :   Shared memory computing, set the local processor
      !!      variables to the value of the global domain
      !!
      !! History :
      !!   9.0  !  04-01  (G. Madec, J.M. Molines)  F90 : free form, north fold jpni >1
      !!----------------------------------------------------------------------

      ! No mpp computation
      nimpp  = 1
      njmpp  = 1
      nlci   = jpi
      nlcj   = jpj
      nldi   = 1
      nldj   = 1
      nlei   = jpi
      nlej   = jpj
      nperio = jperio
      nbondi = 2
      nbondj = 2
      nidom  = FLIO_DOM_NONE
      npolj = jperio

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'mpp_init(2) : NO massively parallel processing'
         WRITE(numout,*) '~~~~~~~~~~~: '
         WRITE(numout,*) '         nperio = ', nperio
         WRITE(numout,*) '         npolj  = ', npolj
         WRITE(numout,*) '         nimpp  = ', nimpp
         WRITE(numout,*) '         njmpp  = ', njmpp
      ENDIF

      IF(  jpni /= 1 .OR. jpnj /= 1 .OR. jpnij /= 1 ) &
          CALL ctl_stop( 'equality  jpni = jpnj = jpnij = 1 is not satisfied',   &
          &              'the domain is lay out for distributed memory computing! ' )

   END SUBROUTINE mpp_init


   SUBROUTINE mpp_init2 
      CALL mpp_init                             ! same routine as mpp_init
   END SUBROUTINE mpp_init2


   !!======================================================================
END MODULE mppini
