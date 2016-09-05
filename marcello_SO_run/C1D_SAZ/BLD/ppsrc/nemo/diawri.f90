MODULE diawri
   !!======================================================================
   !!                     ***  MODULE  diawri  ***
   !! Ocean diagnostics :  write ocean output files
   !!=====================================================================
   !! History :  OPA  ! 1991-03  (M.-A. Foujols)  Original code
   !!            4.0  ! 1991-11  (G. Madec)
   !!                 ! 1992-06  (M. Imbard)  correction restart file
   !!                 ! 1992-07  (M. Imbard)  split into diawri and rstwri
   !!                 ! 1993-03  (M. Imbard)  suppress writibm
   !!                 ! 1998-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
   !!                 ! 1999-02  (E. Guilyardi)  name of netCDF files + variables
   !!            8.2  ! 2000-06  (M. Imbard)  Original code (diabort.F)
   !!   NEMO     1.0  ! 2002-06  (A.Bozec, E. Durand)  Original code (diainit.F)
   !!             -   ! 2002-09  (G. Madec)  F90: Free form and module
   !!             -   ! 2002-12  (G. Madec)  merge of diabort and diainit, F90
   !!                 ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2008-11  (B. Lemaire) creation from old diawri
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_wri       : create the standart output files
   !!   dia_wri_state : create an output NetCDF file for a single instantaeous ocean state and forcing fields
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE dynadv, ONLY: ln_dynadv_vec
   USE zdf_oce         ! ocean vertical physics
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE traldf_iso_grif, ONLY : psix_eiv, psiy_eiv
   USE sol_oce         ! solver variables
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE icb_oce         ! Icebergs
   USE icbdia          ! Iceberg budgets
   USE sbcssr          ! restoring term toward SST/SSS climatology
   USE phycst          ! physical constants
   USE zdfmxl          ! mixed layer
   USE dianam          ! build name of file (routine)
   USE zdfddm          ! vertical  physics: double diffusion
   USE diahth          ! thermocline diagnostics
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE diadimg         ! dimg direct access file format output
   USE iom
   USE ioipsl
   USE dynspg_oce, ONLY: un_adv, vn_adv ! barotropic velocities     

   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE wrk_nemo        ! working array

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_wri                 ! routines called by step.F90
   PUBLIC   dia_wri_state
   PUBLIC   dia_wri_alloc           ! Called by nemogcm module

   INTEGER ::   nid_T, nz_T, nh_T, ndim_T, ndim_hT   ! grid_T file
   INTEGER ::          nb_T              , ndim_bT   ! grid_T file
   INTEGER ::   nid_U, nz_U, nh_U, ndim_U, ndim_hU   ! grid_U file
   INTEGER ::   nid_V, nz_V, nh_V, ndim_V, ndim_hV   ! grid_V file
   INTEGER ::   nid_W, nz_W, nh_W                    ! grid_W file
   INTEGER ::   ndex(1)                              ! ???
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT, ndex_hU, ndex_hV
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_T, ndex_U, ndex_V
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_bT

   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                    *** zdfddm_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
!   Defautl option :                     avs = avt
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfddm_substitute.h90 4152 2013-11-05 11:59:53Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
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
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop start/end indices with CPP macro
   !!                allow unrolling of do-loop (useful with vector processors)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.7 , NEMO Consortium (2014)
   !! $Id: vectopt_loop_substitute.h90 4990 2014-12-15 16:42:49Z timgraham $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: diawri.F90 6348 2016-02-24 10:44:07Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dia_wri_alloc()
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ierr
      !!----------------------------------------------------------------------
      ierr = 0
      ALLOCATE( ndex_hT(jpi*jpj) , ndex_T(jpi*jpj*jpk) ,     &
         &      ndex_hU(jpi*jpj) , ndex_U(jpi*jpj*jpk) ,     &
         &      ndex_hV(jpi*jpj) , ndex_V(jpi*jpj*jpk) , STAT=ierr(1) )
         !
      dia_wri_alloc = MAXVAL(ierr)
      IF( lk_mpp )   CALL mpp_sum( dia_wri_alloc )
      !
  END FUNCTION dia_wri_alloc

   !!----------------------------------------------------------------------
   !!   Default option                                   NetCDF output file
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_iomput'                                        use IOM library
   !!----------------------------------------------------------------------

   SUBROUTINE dia_wri( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!                   
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!
      !! ** Method  :  use iom_put
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER                      ::   ji, jj, jk              ! dummy loop indices
      INTEGER                      ::   jkbot                   !
      REAL(wp)                     ::   zztmp, zztmpx, zztmpy   ! 
      !!
      REAL(wp), POINTER, DIMENSION(:,:)   :: z2d      ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: z3d      ! 3D workspace
      !!----------------------------------------------------------------------
      ! 
      IF( nn_timing == 1 )   CALL timing_start('dia_wri')
      ! 
      CALL wrk_alloc( jpi , jpj      , z2d )
      CALL wrk_alloc( jpi , jpj, jpk , z3d )
      !
      ! Output the initial state and forcings
      IF( ninist == 1 ) THEN                       
         CALL dia_wri_state( 'output.init', kt )
         ninist = 0
      ENDIF

      ! Output of initial vertical scale factor
      CALL iom_put("e3t_0", e3t_0(:,:,:) )
      CALL iom_put("e3u_0", e3t_0(:,:,:) )
      CALL iom_put("e3v_0", e3t_0(:,:,:) )
      !
      CALL iom_put( "e3t" , e3t_0(:,:,:) )
      CALL iom_put( "e3u" , e3u_0(:,:,:) )
      CALL iom_put( "e3v" , e3v_0(:,:,:) )
      CALL iom_put( "e3w" , e3w_0(:,:,:) )
      IF( iom_use("e3tdef") )   &
         CALL iom_put( "e3tdef"  , ( ( e3t_0(:,:,:) - e3t_0(:,:,:) ) / e3t_0(:,:,:) * 100 * tmask(:,:,:) ) ** 2 )


      CALL iom_put( "ssh" , sshn )                 ! sea surface height
      
      CALL iom_put( "toce", tsn(:,:,:,jp_tem) )    ! 3D temperature
      CALL iom_put(  "sst", tsn(:,:,1,jp_tem) )    ! surface temperature
      IF ( iom_use("sbt") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               jkbot = mbkt(ji,jj)
               z2d(ji,jj) = tsn(ji,jj,jkbot,jp_tem)
            END DO
         END DO
         CALL iom_put( "sbt", z2d )                ! bottom temperature
      ENDIF
      
      CALL iom_put( "soce", tsn(:,:,:,jp_sal) )    ! 3D salinity
      CALL iom_put(  "sss", tsn(:,:,1,jp_sal) )    ! surface salinity
      IF ( iom_use("sbs") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               jkbot = mbkt(ji,jj)
               z2d(ji,jj) = tsn(ji,jj,jkbot,jp_sal)
            END DO
         END DO
         CALL iom_put( "sbs", z2d )                ! bottom salinity
      ENDIF

      IF ( iom_use("taubot") ) THEN                ! bottom stress
         z2d(:,:) = 0._wp
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zztmpx = (  bfrua(ji  ,jj) * un(ji  ,jj,mbku(ji  ,jj))  &
                      &  + bfrua(ji-1,jj) * un(ji-1,jj,mbku(ji-1,jj))  )      
               zztmpy = (  bfrva(ji,  jj) * vn(ji,jj  ,mbkv(ji,jj  ))  &
                      &  + bfrva(ji,jj-1) * vn(ji,jj-1,mbkv(ji,jj-1))  ) 
               z2d(ji,jj) = rau0 * SQRT( zztmpx * zztmpx + zztmpy * zztmpy ) * tmask(ji,jj,1) 
               !
            ENDDO
         ENDDO
         CALL lbc_lnk( z2d, 'T', 1. )
         CALL iom_put( "taubot", z2d )           
      ENDIF
         
      CALL iom_put( "uoce", un(:,:,:)         )    ! 3D i-current
      CALL iom_put(  "ssu", un(:,:,1)         )    ! surface i-current
      IF ( iom_use("sbu") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               jkbot = mbku(ji,jj)
               z2d(ji,jj) = un(ji,jj,jkbot)
            END DO
         END DO
         CALL iom_put( "sbu", z2d )                ! bottom i-current
      ENDIF
      CALL iom_put(  "ubar", un_b(:,:)        )    ! barotropic i-current
      
      CALL iom_put( "voce", vn(:,:,:)         )    ! 3D j-current
      CALL iom_put(  "ssv", vn(:,:,1)         )    ! surface j-current
      IF ( iom_use("sbv") ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               jkbot = mbkv(ji,jj)
               z2d(ji,jj) = vn(ji,jj,jkbot)
            END DO
         END DO
         CALL iom_put( "sbv", z2d )                ! bottom j-current
      ENDIF
      CALL iom_put(  "vbar", vn_b(:,:)        )    ! barotropic j-current

      CALL iom_put( "woce", wn )                   ! vertical velocity
      IF( iom_use('w_masstr') .OR. iom_use('w_masstr2') ) THEN   ! vertical mass transport & its square value
         ! Caution: in the VVL case, it only correponds to the baroclinic mass transport.
         z2d(:,:) = rau0 * e12t(:,:)
         DO jk = 1, jpk
            z3d(:,:,jk) = wn(:,:,jk) * z2d(:,:)
         END DO
         CALL iom_put( "w_masstr" , z3d )  
         IF( iom_use('w_masstr2') )   CALL iom_put( "w_masstr2", z3d(:,:,:) * z3d(:,:,:) )
      ENDIF

      CALL iom_put( "avt" , avt                        )    ! T vert. eddy diff. coef.
      CALL iom_put( "avm" , avmu                       )    ! T vert. eddy visc. coef.
      CALL iom_put( "avs" , avt(:,:,:)               )    ! S vert. eddy diff. coef. (useful only with key_zdfddm)
                                                            ! Log of eddy diff coef
      IF( iom_use('logavt') )   CALL iom_put( "logavt", LOG( MAX( 1.e-20_wp, avt  (:,:,:) ) ) )
      IF( iom_use('logavs') )   CALL iom_put( "logavs", LOG( MAX( 1.e-20_wp, avt(:,:,:) ) ) )

      IF ( iom_use("sstgrad") .OR. iom_use("sstgrad2") ) THEN
         DO jj = 2, jpjm1                                    ! sst gradient
            DO ji = 2, jpim1   ! vector opt.
               zztmp      = tsn(ji,jj,1,jp_tem)
               zztmpx     = ( tsn(ji+1,jj  ,1,jp_tem) - zztmp ) / e1u(ji,jj) + ( zztmp - tsn(ji-1,jj  ,1,jp_tem) ) / e1u(ji-1,jj  )
               zztmpy     = ( tsn(ji  ,jj+1,1,jp_tem) - zztmp ) / e2v(ji,jj) + ( zztmp - tsn(ji  ,jj-1,1,jp_tem) ) / e2v(ji  ,jj-1)
               z2d(ji,jj) = 0.25 * ( zztmpx * zztmpx + zztmpy * zztmpy )   &
                  &              * umask(ji,jj,1) * umask(ji-1,jj,1) * vmask(ji,jj,1) * umask(ji,jj-1,1)
            END DO
         END DO
         CALL lbc_lnk( z2d, 'T', 1. )
         CALL iom_put( "sstgrad2",  z2d               )    ! square of module of sst gradient
         z2d(:,:) = SQRT( z2d(:,:) )
         CALL iom_put( "sstgrad" ,  z2d               )    ! module of sst gradient
      ENDIF
         
      ! clem: heat and salt content
      IF( iom_use("heatc") ) THEN
         z2d(:,:)  = 0._wp 
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  z2d(ji,jj) = z2d(ji,jj) + e3t_0(ji,jj,jk) * tsn(ji,jj,jk,jp_tem) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL iom_put( "heatc", (rau0 * rcp) * z2d )    ! vertically integrated heat content (J/m2)
      ENDIF

      IF( iom_use("saltc") ) THEN
         z2d(:,:)  = 0._wp 
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  z2d(ji,jj) = z2d(ji,jj) + e3t_0(ji,jj,jk) * tsn(ji,jj,jk,jp_sal) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL iom_put( "saltc", rau0 * z2d )   ! vertically integrated salt content (PSU*kg/m2)
      ENDIF
      !
      IF ( iom_use("eken") ) THEN
         rke(:,:,jk) = 0._wp                               !      kinetic energy 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zztmp   = 1._wp / ( e1e2t(ji,jj) * e3t_0(ji,jj,jk) )
                  zztmpx  = 0.5 * (  un(ji-1,jj,jk) * un(ji-1,jj,jk) * e2u(ji-1,jj) * e3u_0(ji-1,jj,jk)    &
                     &             + un(ji  ,jj,jk) * un(ji  ,jj,jk) * e2u(ji  ,jj) * e3u_0(ji  ,jj,jk) )  &
                     &          *  zztmp 
                  !
                  zztmpy  = 0.5 * (  vn(ji,jj-1,jk) * vn(ji,jj-1,jk) * e1v(ji,jj-1) * e3v_0(ji,jj-1,jk)    &
                     &             + vn(ji,jj  ,jk) * vn(ji,jj  ,jk) * e1v(ji,jj  ) * e3v_0(ji,jj  ,jk) )  &
                     &          *  zztmp 
                  !
                  rke(ji,jj,jk) = 0.5_wp * ( zztmpx + zztmpy )
                  !
               ENDDO
            ENDDO
         ENDDO
         CALL lbc_lnk( rke, 'T', 1. )
         CALL iom_put( "eken", rke )           
      ENDIF
      !
      CALL iom_put( "hdiv", hdivn )                  ! Horizontal divergence
      !
      IF( iom_use("u_masstr") .OR. iom_use("u_heattr") .OR. iom_use("u_salttr") ) THEN
         z3d(:,:,jpk) = 0.e0
         DO jk = 1, jpkm1
            z3d(:,:,jk) = rau0 * un(:,:,jk) * e2u(:,:) * e3u_0(:,:,jk) * umask(:,:,jk)
         END DO
         CALL iom_put( "u_masstr", z3d )                  ! mass transport in i-direction
      ENDIF
      
      IF( iom_use("u_heattr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( z2d, 'U', -1. )
         CALL iom_put( "u_heattr", (0.5 * rcp) * z2d )    ! heat transport in i-direction
      ENDIF

      IF( iom_use("u_salttr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( z2d, 'U', -1. )
         CALL iom_put( "u_salttr", 0.5 * z2d )            ! heat transport in i-direction
      ENDIF

      
      IF( iom_use("v_masstr") .OR. iom_use("v_heattr") .OR. iom_use("v_salttr") ) THEN
         z3d(:,:,jpk) = 0.e0
         DO jk = 1, jpkm1
            z3d(:,:,jk) = rau0 * vn(:,:,jk) * e1v(:,:) * e3v_0(:,:,jk) * vmask(:,:,jk)
         END DO
         CALL iom_put( "v_masstr", z3d )                  ! mass transport in j-direction
      ENDIF
      
      IF( iom_use("v_heattr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_tem) + tsn(ji,jj+1,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( z2d, 'V', -1. )
         CALL iom_put( "v_heattr", (0.5 * rcp) * z2d )    !  heat transport in j-direction
      ENDIF

      IF( iom_use("v_salttr") ) THEN
         z2d(:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk) * ( tsn(ji,jj,jk,jp_sal) + tsn(ji,jj+1,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( z2d, 'V', -1. )
         CALL iom_put( "v_salttr", 0.5 * z2d )            !  heat transport in j-direction
      ENDIF
      !
      CALL wrk_dealloc( jpi , jpj      , z2d )
      CALL wrk_dealloc( jpi , jpj, jpk , z3d )
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_wri')
      !
   END SUBROUTINE dia_wri



   SUBROUTINE dia_wri_state( cdfile_name, kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ocean state and forcing fields.
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! ** Method  :   NetCDF files using ioipsl
      !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
      !!      File 'output.abort.nc' is created in case of abnormal job end
      !!----------------------------------------------------------------------
      CHARACTER (len=* ), INTENT( in ) ::   cdfile_name      ! name of the file created
      INTEGER           , INTENT( in ) ::   kt               ! ocean time-step index
      !! 
      CHARACTER (len=32) :: clname
      CHARACTER (len=40) :: clop
      INTEGER  ::   id_i , nz_i, nh_i       
      INTEGER, DIMENSION(1) ::   idex             ! local workspace
      REAL(wp) ::   zsto, zout, zmax, zjulian, zdt
      !!----------------------------------------------------------------------
      ! 
!     IF( nn_timing == 1 )   CALL timing_start('dia_wri_state') ! not sure this works for routines not called in first timestep

      ! 0. Initialisation
      ! -----------------

      ! Define name, frequency of output and means
      clname = cdfile_name
      IF( .NOT. Agrif_Root() ) clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      zdt  = rdt
      zsto = rdt
      clop = "inst(x)"           ! no use of the mask value (require less cpu time)
      zout = rdt
      zmax = ( nitend - nit000 + 1 ) * zdt

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
      IF(lwp) WRITE(numout,*) '                and named :', clname, '.nc'


      ! 1. Define NETCDF files and fields at beginning of first time step
      ! -----------------------------------------------------------------

      ! Compute julian date from starting date of the run
      CALL ymds2ju( nyear, nmonth, nday, rdt, zjulian )         ! time axis 
      zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment
      CALL histbeg( clname, jpi, glamt, jpj, gphit,   &
          1, jpi, 1, jpj, nit000-1, zjulian, zdt, nh_i, id_i, domain_id=nidom, snc4chunks=snc4set ) ! Horizontal grid : glamt and gphit
      CALL histvert( id_i, "deptht", "Vertical T levels",   &    ! Vertical grid : gdept
          "m", jpk, gdept_1d, nz_i, "down")

      ! Declare all the output fields as NetCDF variables

      CALL histdef( id_i, "vosaline", "Salinity"              , "PSU"    ,   &   ! salinity
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "votemper", "Temperature"           , "C"      ,   &   ! temperature
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "sossheig", "Sea Surface Height"    , "m"      ,   &  ! ssh
         &          jpi, jpj, nh_i, 1  , 1, 1  , nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "vozocrtx", "Zonal Current"         , "m/s"    ,   &   ! zonal current
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      CALL histdef( id_i, "vomecrty", "Meridional Current"    , "m/s"    ,   &   ! meridonal current
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
      CALL histdef( id_i, "vovecrtz", "Vertical Velocity"     , "m/s"    ,   &   ! vertical current
         &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout ) 
      CALL histdef( id_i, "sowaflup", "Net Upward Water Flux" , "Kg/m2/S",   &   ! net freshwater 
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sohefldo", "Net Downward Heat Flux", "W/m2"   ,   &   ! net heat flux
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "soshfldo", "Shortwave Radiation"   , "W/m2"   ,   &   ! solar flux
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "soicecov", "Ice fraction"          , "[0,1]"  ,   &   ! fr_i
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sozotaux", "Zonal Wind Stress"     , "N/m2"   ,   &   ! i-wind stress
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      CALL histdef( id_i, "sometauy", "Meridional Wind Stress", "N/m2"   ,   &   ! j-wind stress
         &          jpi, jpj, nh_i, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
      IF( lk_vvl ) THEN
         CALL histdef( id_i, "vovvldep", "T point depth"         , "m"      ,   &   ! t-point depth
            &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
         CALL histdef( id_i, "vovvle3t", "T point thickness"         , "m"      ,   &   ! t-point depth
            &          jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout )
      END IF

      CALL histend( id_i, snc4chunks=snc4set )

      ! 2. Start writing data
      ! ---------------------
      ! idex(1) est utilise ssi l'avant dernier argument est diffferent de 
      ! la taille du tableau en sortie. Dans ce cas , l'avant dernier argument
      ! donne le nombre d'elements, et idex la liste des indices a sortir
      idex(1) = 1   ! init to avoid compil warning

      ! Write all fields on T grid
      CALL histwrite( id_i, "votemper", kt, tsn(:,:,:,jp_tem), jpi*jpj*jpk, idex )    ! now temperature
      CALL histwrite( id_i, "vosaline", kt, tsn(:,:,:,jp_sal), jpi*jpj*jpk, idex )    ! now salinity
      CALL histwrite( id_i, "sossheig", kt, sshn             , jpi*jpj    , idex )    ! sea surface height
      CALL histwrite( id_i, "vozocrtx", kt, un               , jpi*jpj*jpk, idex )    ! now i-velocity
      CALL histwrite( id_i, "vomecrty", kt, vn               , jpi*jpj*jpk, idex )    ! now j-velocity
      CALL histwrite( id_i, "vovecrtz", kt, wn               , jpi*jpj*jpk, idex )    ! now k-velocity
      CALL histwrite( id_i, "sowaflup", kt, (emp-rnf )       , jpi*jpj    , idex )    ! freshwater budget
      CALL histwrite( id_i, "sohefldo", kt, qsr + qns        , jpi*jpj    , idex )    ! total heat flux
      CALL histwrite( id_i, "soshfldo", kt, qsr              , jpi*jpj    , idex )    ! solar heat flux
      CALL histwrite( id_i, "soicecov", kt, fr_i             , jpi*jpj    , idex )    ! ice fraction
      CALL histwrite( id_i, "sozotaux", kt, utau             , jpi*jpj    , idex )    ! i-wind stress
      CALL histwrite( id_i, "sometauy", kt, vtau             , jpi*jpj    , idex )    ! j-wind stress
      IF( lk_vvl ) THEN
         CALL histwrite( id_i, "vovvldep", kt, gdept_0(:,:,:), jpi*jpj*jpk, idex )!  T-cell depth       
         CALL histwrite( id_i, "vovvle3t", kt, e3t_0(:,:,:), jpi*jpj*jpk, idex )!  T-cell thickness  
      END IF

      ! 3. Close the file
      ! -----------------
      CALL histclo( id_i )
       
!     IF( nn_timing == 1 )   CALL timing_stop('dia_wri_state') ! not sure this works for routines not called in first timestep
      ! 

   END SUBROUTINE dia_wri_state
   !!======================================================================
END MODULE diawri
