MODULE timing
   !!========================================================================
   !!                     ***  MODULE  timing  ***
   !!========================================================================
   !! History : 4.0  ! 2001-05  (R. Benshila)   
   !!------------------------------------------------------------------------

   !!------------------------------------------------------------------------
   !!   timming_init    : initialize timing process 
   !!   timing_start    : start Timer
   !!   timing_stop     : stop  Timer
   !!   timing_reset    : end timing variable creation
   !!   timing_finalize : compute stats and write output in calling w*_info 
   !!   timing_ini_var  : create timing variables 
   !!   timing_listing  : print instumented subroutines in ocean.output
   !!   wcurrent_info   : compute and print detailed stats on the current CPU
   !!   wave_info       : compute and print averaged statson all processors
   !!   wmpi_info       : compute and write global stats  
   !!   supress         : suppress an element of the timing linked list  
   !!   insert          : insert an element of the timing linked list  
   !!------------------------------------------------------------------------
   USE in_out_manager  ! I/O manager 
   USE dom_oce         ! ocean domain
   USE lib_mpp          
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   timing_init, timing_finalize   ! called in nemogcm module 
   PUBLIC   timing_reset                   ! called in step module 
   PUBLIC   timing_start, timing_stop      ! called in each routine to time 
   

   ! Variables for fine grain timing
   TYPE timer
      CHARACTER(LEN=20)  :: cname
   	  REAL(wp)  :: t_cpu, t_clock, tsum_cpu, tsum_clock, tmax_cpu, tmax_clock, tmin_cpu, tmin_clock, tsub_cpu, tsub_clock
      INTEGER :: ncount, ncount_max, ncount_rate  
      INTEGER :: niter
      LOGICAL :: l_tdone
      TYPE(timer), POINTER :: next => NULL()
      TYPE(timer), POINTER :: prev => NULL()
      TYPE(timer), POINTER :: parent_section => NULL()
   END TYPE timer
    
   TYPE alltimer
      CHARACTER(LEN=20), DIMENSION(:), POINTER :: cname => NULL()
   	  REAL(wp), DIMENSION(:), POINTER :: tsum_cpu   => NULL()
   	  REAL(wp), DIMENSION(:), POINTER :: tsum_clock => NULL()
   	  INTEGER, DIMENSION(:), POINTER :: niter => NULL()
      TYPE(alltimer), POINTER :: next => NULL()
      TYPE(alltimer), POINTER :: prev => NULL()
   END TYPE alltimer 
 
   TYPE(timer), POINTER :: s_timer_root => NULL()
   TYPE(timer), POINTER :: s_timer      => NULL()
   TYPE(timer), POINTER :: s_wrk        => NULL()
   REAL(wp) :: t_overclock, t_overcpu
   LOGICAL :: l_initdone = .FALSE.
   INTEGER :: nsize
   
   ! Variables for coarse grain timing
   REAL(wp) :: tot_etime, tot_ctime
   REAL(kind=wp), DIMENSION(2)     :: t_elaps, t_cpu
   REAL(wp), ALLOCATABLE, DIMENSION(:) :: all_etime, all_ctime
   INTEGER :: nfinal_count, ncount, ncount_rate, ncount_max
   INTEGER, DIMENSION(8)           :: nvalues
   CHARACTER(LEN=8), DIMENSION(2)  :: cdate
   CHARACTER(LEN=10), DIMENSION(2) :: ctime
   CHARACTER(LEN=5)                :: czone
    
   ! From of ouput file (1/proc or one global)   !RB to put in nammpp or namctl
   LOGICAL :: ln_onefile = .TRUE. 
   LOGICAL :: lwriter
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: timing.F90 5120 2015-03-03 16:11:55Z acc $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE timing_start(cdinfo)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_start  ***
      !! ** Purpose :   collect execution time
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdinfo
      !
       
      ! Create timing structure at first call
      IF( .NOT. l_initdone ) THEN
         CALL timing_ini_var(cdinfo)
      ELSE
         s_timer => s_timer_root
         DO WHILE( TRIM(s_timer%cname) /= TRIM(cdinfo) ) 
            IF( ASSOCIATED(s_timer%next) ) s_timer => s_timer%next
         END DO
      ENDIF         
      s_timer%l_tdone = .FALSE.
      s_timer%niter = s_timer%niter + 1
      s_timer%t_cpu = 0.
      s_timer%t_clock = 0.
                  
      ! CPU time collection
      CALL CPU_TIME( s_timer%t_cpu  )
      ! clock time collection
      CALL SYSTEM_CLOCK(COUNT_RATE=s_timer%ncount_rate, COUNT_MAX=s_timer%ncount_max)
      CALL SYSTEM_CLOCK(COUNT = s_timer%ncount)
      !
   END SUBROUTINE timing_start


   SUBROUTINE timing_stop(cdinfo, csection)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_stop  ***
      !! ** Purpose :   finalize timing and output
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdinfo
      CHARACTER(len=*), INTENT(in), OPTIONAL :: csection
      !
      INTEGER  :: ifinal_count, iperiods    
      REAL(wp) :: zcpu_end, zmpitime
      !
      s_wrk => NULL()

      ! clock time collection
      CALL SYSTEM_CLOCK(COUNT = ifinal_count)
      ! CPU time collection
      CALL CPU_TIME( zcpu_end )

      s_timer => s_timer_root
      DO WHILE( TRIM(s_timer%cname) /= TRIM(cdinfo) ) 
         IF( ASSOCIATED(s_timer%next) ) s_timer => s_timer%next
      END DO
 
      ! CPU time correction
      s_timer%t_cpu  = zcpu_end - s_timer%t_cpu - t_overcpu - s_timer%tsub_cpu
  
      ! clock time correction
      iperiods = ifinal_count - s_timer%ncount
      IF( ifinal_count < s_timer%ncount )  &
          iperiods = iperiods + s_timer%ncount_max 
      s_timer%t_clock  = REAL(iperiods) / s_timer%ncount_rate - t_overclock - s_timer%tsub_clock
      
      ! Correction of parent section
      IF( .NOT. PRESENT(csection) ) THEN
         s_wrk => s_timer
         DO WHILE ( ASSOCIATED(s_wrk%parent_section ) )
            s_wrk => s_wrk%parent_section
            s_wrk%tsub_cpu   = s_wrk%tsub_cpu   + s_timer%t_cpu 
            s_wrk%tsub_clock = s_wrk%tsub_clock + s_timer%t_clock              
         END DO
      ENDIF
            
      ! time diagnostics 
      s_timer%tsum_clock = s_timer%tsum_clock + s_timer%t_clock 
      s_timer%tsum_cpu   = s_timer%tsum_cpu   + s_timer%t_cpu
!RB to use to get min/max during a time integration
!      IF( .NOT. l_initdone ) THEN
!         s_timer%tmin_clock = s_timer%t_clock 
!         s_timer%tmin_cpu   = s_timer%t_cpu 
!      ELSE
!         s_timer%tmin_clock = MIN( s_timer%tmin_clock, s_timer%t_clock ) 
!         s_timer%tmin_cpu   = MIN( s_timer%tmin_cpu  , s_timer%t_cpu   ) 
!      ENDIF   
!      s_timer%tmax_clock = MAX( s_timer%tmax_clock, s_timer%t_clock ) 
!      s_timer%tmax_cpu   = MAX( s_timer%tmax_cpu  , s_timer%t_cpu   )  
      !
      s_timer%tsub_clock = 0.
      s_timer%tsub_cpu = 0.
      s_timer%l_tdone = .TRUE.
      !
   END SUBROUTINE timing_stop
 
 
   SUBROUTINE timing_init
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_init  ***
      !! ** Purpose :   open timing output file
      !!----------------------------------------------------------------------
      INTEGER :: iperiods, istart_count, ifinal_count
      REAL(wp) :: zdum
      LOGICAL :: ll_f
             
      IF( ln_onefile ) THEN
         IF( lwp) CALL ctl_opn( numtime, 'timing.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout,.TRUE., narea )
         lwriter = lwp
      ELSE
         CALL ctl_opn( numtime, 'timing.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout,.FALSE., narea )
         lwriter = .TRUE.
      ENDIF
      
      IF( lwriter) THEN      
         WRITE(numtime,*)
         WRITE(numtime,*) '      CNRS - NERC - Met OFFICE - MERCATOR-ocean - CMCC - INGV'
         WRITE(numtime,*) '                             NEMO team'
         WRITE(numtime,*) '                  Ocean General Circulation Model'
         WRITE(numtime,*) '                        version 3.6  (2015) '
         WRITE(numtime,*)
         WRITE(numtime,*) '                        Timing Informations '
         WRITE(numtime,*)
         WRITE(numtime,*)
      ENDIF   
      
      ! Compute clock function overhead
      CALL SYSTEM_CLOCK(COUNT_RATE=ncount_rate, COUNT_MAX=ncount_max)
      CALL SYSTEM_CLOCK(COUNT = istart_count)
      CALL SYSTEM_CLOCK(COUNT = ifinal_count)
      iperiods = ifinal_count - istart_count
      IF( ifinal_count < istart_count )  &
          iperiods = iperiods + ncount_max 
      t_overclock = REAL(iperiods) / ncount_rate

      ! Compute cpu_time function overhead
      CALL CPU_TIME(zdum)
      CALL CPU_TIME(t_overcpu)
      
      ! End overhead omputation  
      t_overcpu = t_overcpu - zdum        
      t_overclock = t_overcpu + t_overclock        

      ! Timing on date and time
      CALL DATE_AND_TIME(cdate(1),ctime(1),czone,nvalues)
    
      CALL CPU_TIME(t_cpu(1))      
      CALL SYSTEM_CLOCK(COUNT_RATE=ncount_rate, COUNT_MAX=ncount_max)
      CALL SYSTEM_CLOCK(COUNT = ncount)
      !
   END SUBROUTINE timing_init


   SUBROUTINE timing_finalize
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_finalize ***
      !! ** Purpose :  compute average time 
      !!               write timing output file
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER :: s_temp
      INTEGER :: idum, iperiods, icode
      LOGICAL :: ll_ord, ll_averep
      CHARACTER(len=120) :: clfmt            
      
      ll_averep = .TRUE.
    
      ! total CPU and elapse
      CALL CPU_TIME(t_cpu(2))
      t_cpu(2)   = t_cpu(2)    - t_cpu(1)   - t_overcpu
      CALL SYSTEM_CLOCK(COUNT = nfinal_count)
      iperiods = nfinal_count - ncount
      IF( nfinal_count < ncount )  &
          iperiods = iperiods + ncount_max 
      t_elaps(2) = REAL(iperiods) / ncount_rate - t_overclock

      ! End of timings on date & time
      CALL DATE_AND_TIME(cdate(2),ctime(2),czone,nvalues)
       
      ! Compute the numer of routines
      nsize = 0 
      s_timer => s_timer_root
      DO WHILE( ASSOCIATED(s_timer) )
         nsize = nsize + 1
         s_timer => s_timer%next
      END DO
      idum = nsize
      IF(lk_mpp) CALL mpp_sum(idum)
      IF( idum/jpnij /= nsize ) THEN
         IF( lwriter ) WRITE(numtime,*) '        ===> W A R N I N G: '
         IF( lwriter ) WRITE(numtime,*) ' Some CPU have different number of routines instrumented for timing'
         IF( lwriter ) WRITE(numtime,*) ' No detailed report on averaged timing can be provided'
         IF( lwriter ) WRITE(numtime,*) ' The following detailed report only deals with the current processor'
         IF( lwriter ) WRITE(numtime,*)
         ll_averep = .FALSE.
      ENDIF   

      tot_etime = t_elaps(2)
      tot_ctime = t_cpu  (2)           

      ! write output file
      IF( lwriter ) WRITE(numtime,*) 'Total timing (sum) :'
      IF( lwriter ) WRITE(numtime,*) '--------------------'
      IF( lwriter ) WRITE(numtime,"('Elapsed Time (s)  CPU Time (s)')")
      IF( lwriter ) WRITE(numtime,'(5x,f12.3,1x,f12.3)')  tot_etime, tot_ctime
      IF( lwriter ) WRITE(numtime,*) 
      IF( lwriter ) CALL wcurrent_info
      
      clfmt='(1X,"Timing started on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
      IF( lwriter ) WRITE(numtime, TRIM(clfmt)) &           
      &       cdate(1)(7:8), cdate(1)(5:6), cdate(1)(1:4),   &
      &       ctime(1)(1:2), ctime(1)(3:4), ctime(1)(5:6),   &
      &       czone(1:3),    czone(4:5)                     
      clfmt='(1X,  "Timing   ended on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
      IF( lwriter ) WRITE(numtime, TRIM(clfmt)) &           
      &       cdate(2)(7:8), cdate(2)(5:6), cdate(2)(1:4),   &
      &       ctime(2)(1:2), ctime(2)(3:4), ctime(2)(5:6),   &
      &       czone(1:3),    czone(4:5)

      IF( lwriter ) CLOSE(numtime) 
      !
   END SUBROUTINE timing_finalize
   

   SUBROUTINE wcurrent_info
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE wcurrent_info ***
      !! ** Purpose :  compute and write timing output file
      !!----------------------------------------------------------------------
      LOGICAL :: ll_ord
      CHARACTER(len=2048) :: clfmt            
   
      ! reorder the current list by elapse time      
      s_wrk => NULL()
      s_timer => s_timer_root
      DO
         ll_ord = .TRUE.
         s_timer => s_timer_root
         DO WHILE ( ASSOCIATED( s_timer%next ) )
         IF (.NOT. ASSOCIATED(s_timer%next)) EXIT
            IF ( s_timer%tsum_clock < s_timer%next%tsum_clock ) THEN 
               ALLOCATE(s_wrk)
               s_wrk = s_timer%next
               CALL insert  (s_timer, s_timer_root, s_wrk)
               CALL suppress(s_timer%next)            
               ll_ord = .FALSE.
               CYCLE            
            ENDIF           
         IF( ASSOCIATED(s_timer%next) ) s_timer => s_timer%next
         END DO         
         IF( ll_ord ) EXIT
      END DO
            
      ! write current info
      WRITE(numtime,*) 'Detailed timing for proc :', narea-1
      WRITE(numtime,*) '--------------------------'
      WRITE(numtime,*) 'Section             ',            &
      &   'Elapsed Time (s)  ','Elapsed Time (%)  ',   &
      &   'CPU Time(s)  ','CPU Time (%)  ','CPU/Elapsed  ','Frequency' 
      s_timer => s_timer_root  
      clfmt = '(1x,a,4x,f12.3,6x,f12.3,x,f12.3,2x,f12.3,6x,f7.3,2x,i9)'
      DO WHILE ( ASSOCIATED(s_timer) )
         WRITE(numtime,TRIM(clfmt))   s_timer%cname,   &
         &   s_timer%tsum_clock,s_timer%tsum_clock*100./t_elaps(2),            &
         &   s_timer%tsum_cpu  ,s_timer%tsum_cpu*100./t_cpu(2)    ,            &
         &   s_timer%tsum_cpu/s_timer%tsum_clock, s_timer%niter
         s_timer => s_timer%next
      END DO
      WRITE(numtime,*)
      !                  
   END SUBROUTINE wcurrent_info



   SUBROUTINE timing_ini_var(cdinfo)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_ini_var  ***
      !! ** Purpose :   create timing structure 
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdinfo
      LOGICAL :: ll_section
       
      !
      IF( .NOT. ASSOCIATED(s_timer_root) ) THEN
         ALLOCATE(s_timer_root)
         s_timer_root%cname       = cdinfo
         s_timer_root%t_cpu      = 0._wp
         s_timer_root%t_clock    = 0._wp
         s_timer_root%tsum_cpu   = 0._wp
         s_timer_root%tsum_clock = 0._wp
         s_timer_root%tmax_cpu   = 0._wp
         s_timer_root%tmax_clock = 0._wp
         s_timer_root%tmin_cpu   = 0._wp
         s_timer_root%tmin_clock = 0._wp
         s_timer_root%tsub_cpu   = 0._wp
         s_timer_root%tsub_clock = 0._wp
         s_timer_root%ncount      = 0
         s_timer_root%ncount_rate = 0
         s_timer_root%ncount_max  = 0
         s_timer_root%niter       = 0
         s_timer_root%l_tdone  = .FALSE.
         s_timer_root%next => NULL()
         s_timer_root%prev => NULL()
         s_timer => s_timer_root
         !
         ALLOCATE(s_wrk)
         s_wrk => NULL()
         
      ELSE
         s_timer => s_timer_root
         ! case of already existing area (typically inside a loop)
         DO WHILE( ASSOCIATED(s_timer) ) 
            IF( TRIM(s_timer%cname) .EQ. TRIM(cdinfo) ) RETURN
            s_timer => s_timer%next
         END DO
         
         ! end of the chain
         s_timer => s_timer_root
         DO WHILE( ASSOCIATED(s_timer%next) )
            s_timer => s_timer%next
         END DO
          
         ALLOCATE(s_timer%next)      
         s_timer%next%cname       = cdinfo
         s_timer%next%t_cpu      = 0._wp
         s_timer%next%t_clock    = 0._wp
         s_timer%next%tsum_cpu   = 0._wp
         s_timer%next%tsum_clock = 0._wp  
         s_timer%next%tmax_cpu   = 0._wp
         s_timer%next%tmax_clock = 0._wp
         s_timer%next%tmin_cpu   = 0._wp
         s_timer%next%tmin_clock = 0._wp
         s_timer%next%tsub_cpu   = 0._wp
         s_timer%next%tsub_clock = 0._wp
         s_timer%next%ncount      = 0
         s_timer%next%ncount_rate = 0
         s_timer%next%ncount_max  = 0
         s_timer%next%niter       = 0
         s_timer%next%l_tdone  = .FALSE.
         s_timer%next%parent_section => NULL()
         s_timer%next%prev => s_timer
         s_timer%next%next => NULL()
         s_timer => s_timer%next

         ! are we inside a section
         s_wrk => s_timer%prev
         ll_section = .FALSE.
         DO WHILE( ASSOCIATED(s_wrk) .AND. .NOT. ll_section )
            IF( .NOT. s_wrk%l_tdone ) THEN
               ll_section = .TRUE.
               s_timer%parent_section => s_wrk 
            ENDIF
            s_wrk => s_wrk%prev
         END DO 
      ENDIF         
      !
   END SUBROUTINE timing_ini_var


   SUBROUTINE timing_reset
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE timing_reset  ***
      !! ** Purpose :   go to root of timing tree 
      !!----------------------------------------------------------------------
      l_initdone = .TRUE. 
!      IF(lwp) WRITE(numout,*)
!      IF(lwp) WRITE(numout,*) 'timing_reset : instrumented routines for timing'
!      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
!      CALL timing_list(s_timer_root)
!      WRITE(numout,*)
      !
   END SUBROUTINE timing_reset


   RECURSIVE SUBROUTINE timing_list(ptr)
   
      TYPE(timer), POINTER, INTENT(inout) :: ptr
      !
      IF( ASSOCIATED(ptr%next) ) CALL timing_list(ptr%next)
      IF(lwp) WRITE(numout,*)'   ', ptr%cname   
      !
   END SUBROUTINE timing_list


   SUBROUTINE insert(sd_current, sd_root ,sd_ptr)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE insert  ***
      !! ** Purpose :   insert an element in timer structure
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(inout) :: sd_current, sd_root, sd_ptr
      !
     
      IF( ASSOCIATED( sd_current, sd_root ) ) THEN
         ! If our current element is the root element then
         ! replace it with the one being inserted
         sd_root => sd_ptr
      ELSE
         sd_current%prev%next => sd_ptr
      END IF
      sd_ptr%next     => sd_current
      sd_ptr%prev     => sd_current%prev
      sd_current%prev => sd_ptr
      ! Nullify the pointer to the new element now that it is held
      ! within the list. If we don't do this then a subsequent call
      ! to ALLOCATE memory to this pointer will fail.
      sd_ptr => NULL()
      !    
   END SUBROUTINE insert
  
  
   SUBROUTINE suppress(sd_ptr)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE suppress  ***
      !! ** Purpose :   supress an element in timer structure
      !!----------------------------------------------------------------------
      TYPE(timer), POINTER, INTENT(inout) :: sd_ptr
      !
      TYPE(timer), POINTER :: sl_temp
    
      sl_temp => sd_ptr
      sd_ptr => sd_ptr%next    
      IF ( ASSOCIATED(sl_temp%next) ) sl_temp%next%prev => sl_temp%prev
      DEALLOCATE(sl_temp)
      sl_temp => NULL()
      !
    END SUBROUTINE suppress

   !!=====================================================================
END MODULE timing
