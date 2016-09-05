












MODULE ldftra_oce
   !!=====================================================================
   !!                      ***  MODULE  ldftra_oce  ***
   !! Ocean physics :  lateral tracer mixing coefficient defined in memory 
   !!=====================================================================
   !! History :  9.0  !  2002-11  (G. Madec)  Original code
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC ldftra_oce_alloc ! called by nemo_init->nemo_alloc, nemogcm.F90

   !!----------------------------------------------------------------------
   !! Lateral eddy diffusivity coefficients (tracers)
   !!----------------------------------------------------------------------
   !                                     !!* Namelist namtra_ldf : lateral mixing *
   LOGICAL , PUBLIC ::   ln_traldf_lap    !: laplacian operator
   LOGICAL , PUBLIC ::   ln_traldf_bilap  !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_traldf_level  !: iso-level direction
   LOGICAL , PUBLIC ::   ln_traldf_hor    !: horizontal (geopotential) direction
   LOGICAL , PUBLIC ::   ln_traldf_iso    !: iso-neutral direction
   LOGICAL , PUBLIC ::   ln_traldf_grif   !: griffies skew flux
   LOGICAL , PUBLIC ::   ln_traldf_gdia   !: griffies skew flux streamfunction diagnostics
   REAL(wp), PUBLIC ::   rn_aht_0         !: lateral eddy diffusivity (m2/s)
   REAL(wp), PUBLIC ::   rn_ahtb_0        !: lateral background eddy diffusivity (m2/s)
   REAL(wp), PUBLIC ::   rn_aeiv_0        !: eddy induced velocity coefficient (m2/s)
   REAL(wp), PUBLIC ::   rn_slpmax        !: slope limit
   REAL(wp), PUBLIC ::   rn_chsmag        !:  multiplicative factor in Smagorinsky diffusivity
   REAL(wp), PUBLIC ::   rn_smsh          !:  Smagorinsky diffusivity: = 0 - use only sheer
   REAL(wp), PUBLIC ::   rn_aht_m         !:  upper limit or stability criteria for lateral eddy diffusivity (m2/s)

   REAL(wp), PUBLIC ::   aht0, ahtb0, aeiv0         !!: OLD namelist names

   LOGICAL , PUBLIC ::   ln_triad_iso    !: calculate triads twice
   LOGICAL , PUBLIC ::   ln_botmix_grif  !: mixing on bottom
   LOGICAL , PUBLIC ::   l_grad_zps      = .FALSE.   !: special treatment for Horz Tgradients w partial steps 

   REAL(wp), PUBLIC ::   rldf                        !: multiplicative factor of diffusive coefficient
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   r_fact_lap
                                                     !: Needed to define the ratio between passive and active tracer diffusion coef. 

   REAL(wp), PUBLIC                                      ::   ahtt, ahtu, ahtv, ahtw   !: ** 0D coefficients ** at T-,U-,V-,W-points

   !!----------------------------------------------------------------------
   !!   Default option :                           NO eddy induced velocity
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_traldf_eiv   = .FALSE.   !: eddy induced velocity flag
   REAL(wp), PUBLIC            ::   aeiu, aeiv, aeiw            !: eddy induced coef. (not used)

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_oce.F90 6312 2016-02-15 11:43:52Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ldftra_oce_alloc()
     !!----------------------------------------------------------------------
      !!                 ***  FUNCTION ldftra_oce_alloc  ***
     !!----------------------------------------------------------------------
     INTEGER, DIMENSION(4) :: ierr
     !!----------------------------------------------------------------------
     ierr(:) = 0

      !
      ALLOCATE( r_fact_lap(jpi,jpj,jpk), STAT=ierr(4) )
      ldftra_oce_alloc = MAXVAL( ierr )
      IF( ldftra_oce_alloc /= 0 )   CALL ctl_warn('ldftra_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION ldftra_oce_alloc

   !!=====================================================================
END MODULE ldftra_oce
