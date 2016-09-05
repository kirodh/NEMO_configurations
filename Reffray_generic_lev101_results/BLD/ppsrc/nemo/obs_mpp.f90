












MODULE obs_mpp
   !!======================================================================
   !!                       ***  MODULE obs_mpp  ***
   !! Observation diagnostics: Various MPP support routines
   !!======================================================================
   !! History :  2.0  ! 2006-03  (K. Mogensen)  Original code
   !!             -   ! 2006-05  (K. Mogensen)  Reformatted
   !!             -   ! 2008-01  (K. Mogensen)  add mpp_global_max
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! obs_mpp_bcast_integer : Broadcast an integer array from a processor to all processors
   !! obs_mpp_max_integer   : Find maximum on all processors of each value in an integer on all processors
   !! obs_mpp_find_obs_proc : Find processors which should hold the observations
   !! obs_mpp_sum_integers  : Sum an integer array from all processors
   !! obs_mpp_sum_integer   : Sum an integer from all processors
   !!----------------------------------------------------------------------
   USE dom_oce, ONLY :   nproc, mig, mjg   ! Ocean space and time domain variables
   USE mpp_map, ONLY :   mppmap
   USE in_out_manager
   IMPLICIT NONE
   PRIVATE

   PUBLIC obs_mpp_bcast_integer, & !: Broadcast an integer array from a proc to all procs
      &   obs_mpp_max_integer,   & !: Find maximum across processors in an integer array
      &   obs_mpp_find_obs_proc, & !: Find processors which should hold the observations
      &   obs_mpp_sum_integers,  & !: Sum an integer array from all processors
      &   obs_mpp_sum_integer,   & !: Sum an integer from all processors
      &   mpp_alltoall_int,      &
      &   mpp_alltoallv_int,     &
      &   mpp_alltoallv_real,    &
      &   mpp_global_max

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_mpp.F90 2513 2010-12-23 16:01:47Z smasson $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE obs_mpp_bcast_integer( kvals, kno, kroot )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_bcast_integer ***
      !!          
      !! ** Purpose : Send array kvals to all processors
      !!
      !! ** Method  : MPI broadcast
      !!
      !! ** Action  : This does only work for MPI. 
      !!              MPI_COMM_OPA needs to be replace for OASIS4.!
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) ::   kno     ! Number of elements in array
      INTEGER                , INTENT(in   ) ::   kroot   ! Processor to send data
      INTEGER, DIMENSION(kno), INTENT(inout) ::   kvals   ! Array to send on kroot, receive for non-kroot
      !
      ! no MPI: empty routine
      !
   END SUBROUTINE obs_mpp_bcast_integer

  
   SUBROUTINE obs_mpp_max_integer( kvals, kno )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_bcast_integer ***
      !!          
      !! ** Purpose : Find maximum across processors in an integer array.
      !!
      !! ** Method  : MPI all reduce.
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!              MPI_COMM_OPA needs to be replace for OASIS4.!
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) ::   kno     ! Number of elements in array
      INTEGER, DIMENSION(kno), INTENT(inout) ::   kvals   ! Array to send on kroot, receive for non-kroot  
      !
      ! no MPI: empty routine
   END SUBROUTINE obs_mpp_max_integer


   SUBROUTINE obs_mpp_find_obs_proc( kobsp, kobsi, kobsj, kno )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_find_obs_proc ***
      !!          
      !! ** Purpose : From the array kobsp containing the results of the grid
      !!              grid search on each processor the processor return a
      !!              decision of which processors should hold the observation.
      !!
      !! ** Method  : A temporary 2D array holding all the decisions is
      !!              constructed using mpi_allgather on each processor.
      !!              If more than one processor has found the observation
      !!              with the observation in the inner domain gets it
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) ::   kno
      INTEGER, DIMENSION(kno), INTENT(in   ) ::   kobsi, kobsj
      INTEGER, DIMENSION(kno), INTENT(inout) ::   kobsp
      !
      ! no MPI: empty routine
      !
   END SUBROUTINE obs_mpp_find_obs_proc


   SUBROUTINE obs_mpp_sum_integers( kvalsin, kvalsout, kno )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_sum_integers ***
      !!          
      !! ** Purpose : Sum an integer array.
      !!
      !! ** Method  : MPI all reduce.
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) :: kno
      INTEGER, DIMENSION(kno), INTENT(in   ) ::   kvalsin
      INTEGER, DIMENSION(kno), INTENT(  out) ::   kvalsout
      !
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      kvalsout(:) = kvalsin(:)
      !
   END SUBROUTINE obs_mpp_sum_integers


   SUBROUTINE obs_mpp_sum_integer( kvalin, kvalout )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_sum_integers ***
      !!          
      !! ** Purpose : Sum a single integer
      !!
      !! ** Method  : MPI all reduce.
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kvalin
      INTEGER, INTENT(  out) ::   kvalout
      !
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      kvalout = kvalin
      !
   END SUBROUTINE obs_mpp_sum_integer


   SUBROUTINE mpp_global_max( pval )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_global_or ***
      !!          
      !! ** Purpose : Get the maximum value across processors for a global 
      !!              real array
      !!
      !! ** Method  : MPI allreduce
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      REAL(KIND=wp), DIMENSION(jpiglo,jpjglo), INTENT(inout) ::   pval
      !
      INTEGER :: ierr
      !
      ! no MPI: empty routine
      !
   END SUBROUTINE mpp_global_max


   SUBROUTINE mpp_alltoall_int( kno, kvalsin, kvalsout )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_allgatherv ***
      !!          
      !! ** Purpose : all to all.
      !!
      !! ** Method  : MPI alltoall
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                      , INTENT(in   ) ::   kno
      INTEGER, DIMENSION(kno*jpnij), INTENT(in   ) ::   kvalsin
      INTEGER, DIMENSION(kno*jpnij), INTENT(  out) ::   kvalsout
      !!
      INTEGER :: ierr
      !
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      kvalsout = kvalsin
      !
   END SUBROUTINE mpp_alltoall_int


   SUBROUTINE mpp_alltoallv_int( kvalsin, knoin , kinv , kvalsout,   &
      &                                   knoout, koutv )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_alltoallv_int ***
      !!          
      !! ** Purpose : all to all (integer version).
      !!
      !! ** Method  : MPI alltoall
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in) :: knoin
      INTEGER                   , INTENT(in) :: knoout
      INTEGER, DIMENSION(jpnij)                 ::   kinv, koutv
      INTEGER, DIMENSION(knoin) , INTENT(in   ) ::   kvalsin
      INTEGER, DIMENSION(knoout), INTENT(  out) ::   kvalsout
      !!
      INTEGER :: ierr
      INTEGER :: jproc
      !
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      kvalsout = kvalsin
      !
   END SUBROUTINE mpp_alltoallv_int


   SUBROUTINE mpp_alltoallv_real( pvalsin, knoin , kinv , pvalsout,   &
      &                                    knoout, koutv )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_alltoallv_real ***
      !!          
      !! ** Purpose : all to all (integer version).
      !!
      !! ** Method  : MPI alltoall
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) :: knoin
      INTEGER                    , INTENT(in   ) :: knoout
      INTEGER , DIMENSION(jpnij)                 ::   kinv, koutv
      REAL(wp), DIMENSION(knoin) , INTENT(in   ) ::   pvalsin
      REAL(wp), DIMENSION(knoout), INTENT(  out) ::   pvalsout
      !!
      INTEGER :: ierr
      INTEGER :: jproc
      !
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      pvalsout = pvalsin
      !
   END SUBROUTINE mpp_alltoallv_real

   !!======================================================================
END MODULE obs_mpp
