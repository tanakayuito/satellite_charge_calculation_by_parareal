!>
!! Provides an MPI-based implementation of the parallel-in-time Parareal method introduced by Lions, Maday and Turinici in 2001.
!! Denote the fine propagator by \\( F \\) and the coarse propagator by \\( G \\). Then the Parareal iteration reads
!!
!! \\( y^{k+1}_{n+1} = G(y^{k+1}_{n}) + F(y^k_n) - G(y^k_n) \\)
!!
!! Here, communication of \\( y^{k+1}_{n+1} \\) to the next time slice is done via MPI.
MODULE parareal_mpi
USE mpi
use OML_currents
use global_variables
use flags


IMPLICIT NONE

!INCLUDE 'mpif.h'

PRIVATE
PUBLIC :: InitializePararealMPI, FinalizePararealMPI, PararealMPI

!> For the MPI version, each MPI process runs only a single OpenMP thread
INTEGER, PARAMETER :: Nthreads = 1

INTEGER, PARAMETER :: order_adv_c = 1, order_diff_c = 2, order_adv_f = 5, order_diff_f = 4

!> Three buffers for the solution used in Parareal
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: para_phi, Gphi, phi_end

!> Time step of the coarse propagator
DOUBLE PRECISION :: dt_coarse

!> Time step of the fine propagator
DOUBLE PRECISION :: dt_fine

!> Length of each time slice
DOUBLE PRECISION :: dt_slice

!> Start point of the time slice handled by this MPI process
DOUBLE PRECISION :: tstart_myslice

!> End time of the time slice handled by this MPI process
DOUBLE PRECISION :: tend_myslice

! Internal variables to be used for timers.
DOUBLE PRECISION :: timer_coarse, timer_fine, timer_comm, timer_all, T0, T1

! A number of MPI-related status variables
INTEGER :: mpi_thread_provided, dim, ierr,  Nproc, myrank, recv_status(MPI_STATUS_SIZE), send_status(MPI_STATUS_SIZE), k = 0

CHARACTER(len=64) :: filename


!TYPE parareal_parameter
    !INTEGER Nx
!END TYPE

!TYPE(parareal_parameter) :: param

CONTAINS

  !> Initialize the Parareal module and all used modules
  SUBROUTINE InitializePararealMPI()
    
    !INTEGER, INTENT(IN) :: iconductor

    !CALL InitializeTimestepper(iconductor)
    ALLOCATE(para_phi(1:10))
    ALLOCATE(phi_end(1:10))
    ALLOCATE(Gphi(1:10))

    dim = SIZE(para_phi)

    !param%Nx = Nx
    

  END SUBROUTINE InitializePararealMPI

  !> Finalize Parareal and used modules
  SUBROUTINE FinalizePararealMPI()
    DEALLOCATE(para_phi)
    DEALLOCATE(Gphi)
    DEALLOCATE(phi_end)
    !CALL FinalizeTimestepper()
  END SUBROUTINE FinalizePararealMPI

  !> Key routine to run Parareal with MPI.
  !> @param[in] Q_initial Initial value
  !> @param[in] Tend Final simulation time
  !> @param[in] N_fine Number of fine steps *per time slice*
  !> @param[in] N_coarse Number of coarse steps *per time slice*
  !> @param[in] Niter Number of Parareal iterations
  !> @param[in] do_io Whether to perform IO or not
  !> @param[in] be_verbose If true, Parareal gives several status messages that can aid in debugging

  double precision function PararealMPI(phi_initial, Tstart, Tend, N_fine, N_coarse, Niter, C_inv, rho1, rho2, Th1, Th2, iconductor,be_verbose)

    DOUBLE PRECISION, DIMENSION(1:10),  INTENT(IN) :: phi_initial
    DOUBLE PRECISION,                 INTENT(IN) :: Tstart, Tend
    INTEGER,                          INTENT(IN) :: N_fine, N_coarse, Niter
    DOUBLE PRECISION,                 INTENT(IN) :: C_inv, rho1, rho2, Th1, Th2
    INTEGER,                          INTENT(IN) :: iconductor
    LOGICAL,                          INTENT(IN) :: be_verbose
    
    
    INTEGER :: n
    DOUBLE PRECISION :: Tall

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nproc, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    
    ! Output
    IF ((myrank==0) .AND. (be_verbose)) THEN
      WRITE(*,'(A, I2)') '--- Running MPI parareal, no. of processes: ', Nproc
    END IF

    ! Divide time interval [0,T] into Nproc many timeslices
    Tall = Tend - Tstart

    dt_slice  = Tall/DBLE(Nproc)
    dt_fine   = dt_slice/DBLE(N_fine)
    dt_coarse = dt_slice/DBLE(N_coarse)

    tstart_myslice = Tstart + DBLE(myrank)*dt_slice
    tend_myslice   = Tstart + DBLE(myrank+1)*dt_slice

    IF (be_verbose) THEN
        PRINT *, 'Process', myrank,' from ', tstart_myslice, ' to ', tend_myslice
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
        IF (myrank==0) THEN
            print *, 'Time slice length:  ', dt_slice
            print *, 'Fine step length:   ', dt_fine
            print *, 'Coarse step length: ', dt_coarse
        END IF
    END IF

    ! --- START PARAREAL --- !

    timer_all    = MPI_WTIME()
    timer_fine   = 0.0
    timer_coarse = 0.0
    timer_comm   = 0.0

    T0 = MPI_WTIME()

    ! Q <- y^0_0
    para_phi = phi_initial
    IF (myrank>0) THEN
      ! coarse predictor to get initial guess at beginning of slice
      ! Q <- y^0_n = G(y^0_(n-1)) = G(G(y^0_(n-2))) = ...
      CALL RungeKutta_Method_Phi_For_Parareal(para_phi(iconductor),Tstart, tstart_myslice, myrank*N_coarse, C_inv, iconductor, rho1, rho2, Th1, Th2,myrank, Nproc, k, 1)
    END IF 

    
    if(just_para_flag == 0) then
      ! one more coarse step to get coarse value at end of slice
      ! GQ <- G(y^0_n)
      Gphi = para_phi
      CALL RungeKutta_Method_Phi_For_Parareal(Gphi(iconductor), tstart_myslice, tend_myslice, N_coarse, C_inv, iconductor, rho1, rho2, Th1, Th2,myrank, Nproc, k,2)

      ! Qend <- G(y^0_n) = y^0_(n+1)
      phi_end = Gphi

    else if(just_para_flag == 1) then
      !  何もしない
    end if    
    

    T1 = MPI_WTIME()
    timer_coarse = timer_coarse + (T1-T0)

    DO k=1,Niter   ! Parareal iteration
         
      ! Initial state: 
      ! Q     <- y^(k-1)_n 
      ! GQ    <- G(y^(k-1)_n) 
      ! Qend  <- y(k-1)_(n+1)
      !
      ! End state:
      ! Q     <- y^k_n
      ! GQ    <- G(y^k_n)
      ! Qend  <- y^k_(n+1) 
           

      T0 = MPI_WTIME()

      ! Run fine integrator:
      ! Q <- F(y^(k-1)_n)
      CALL RungeKutta_Method_Phi_For_Parareal(para_phi(iconductor), tstart_myslice, tend_myslice, N_fine, C_inv, iconductor, rho1, rho2, Th1, Th2, myrank,Nproc, k, 3)      

      ! Compute difference fine minus coarse
      ! Qend <- F(y^(k-1)_n) - G(y^(k-1)_n)
      if(just_para_flag == 0) then
        phi_end = para_phi - Gphi
      else if(just_para_flag == 1) then
        phi_end = para_phi
      end if
      
      T1 = MPI_WTIME()
      timer_fine = timer_fine + (T1-T0)

      ! Fetch updated value from previous process
      IF(myrank>0) THEN
        
        ! Q <- y^k_n 
        T0 = MPI_WTIME()
        CALL MPI_RECV(para_phi, dim, MPI_DOUBLE_PRECISION, myrank-1, 0, MPI_COMM_WORLD, recv_status, ierr)
        T1 = MPI_WTIME()
        timer_comm = timer_comm + (T1-T0)
           
      ELSE IF (myrank==0) THEN

        ! Fetch initial value again
        ! Q <- y^k_n with n=0
        T0 = MPI_WTIME()
        para_phi= phi_initial
        T1 = MPI_WTIME()
        timer_comm = timer_comm + (T1-T0)

      ELSE
        WRITE(*,*) 'Found negative value for myrank, now exiting.'
        STOP  
      END IF
      
      ! Run correction
      ! GQ <- G(y^k_n)
      T0 = MPI_WTIME()
      if(just_para_flag == 0 )then
        ! GQ <- y^k_n
        Gphi = para_phi
        CALL RungeKutta_Method_Phi_For_Parareal(Gphi(iconductor), tstart_myslice, tend_myslice, N_coarse, C_inv, iconductor, rho1, rho2, Th1, Th2,myrank, Nproc,k,4)  
      else if(just_para_flag == 1) then
        !何もしない
      end if
      
      
      T1 = MPI_WTIME()
       
      ! Correct
      ! Qend <- G(y^k_n) + F(y^(k-1))_n - G(y^(k-1)_n) = y^k_(n+1)
      if(just_para_flag == 0 )then
        phi_end = Gphi + phi_end
      else if(just_para_flag == 1) then
        !何もしない
      end if

      timer_coarse = timer_coarse + (T1-T0)
        
      ! Send forward updated value
      IF(myrank<Nproc-1) THEN
        T0 = MPI_WTIME()
        CALL MPI_SEND(phi_end, dim, MPI_DOUBLE_PRECISION, myrank+1, 0, MPI_COMM_WORLD, send_status, ierr)
        T1 = MPI_WTIME()
        timer_comm = timer_comm + (T1-T0)
      END IF 
      
      ! Final state:
      ! Q    <- y^k_n
      ! GQ   <- G(y^k_n)
      ! Qend <- y^k_(n+1)
      
    END DO ! End of parareal loop

    ! --- END PARAREAL --- !
    timer_all = MPI_WTIME() - timer_all

    ! Return final value in Q_initial    
    PararealMPI = phi_end(iconductor)

    if(myrank == Nproc -1)then
      WRITE(*,'(A, F35.25)') 'phi at Tend:  ', phi_end(iconductor)
    end if

    !timer_all = MPI_WTIME() - timer_all

    !WRITE(filename, '(A,I0.2,A,I0.2,A,I0.2,A)') 'timings_mpi', Niter, '_', myrank, '_', Nproc, '.dat'
    !OPEN(UNIT=myrank, FILE=filename, ACTION='write', STATUS='replace')
    !WRITE(myrank, '(F8.5)') timer_all
    !WRITE(myrank, '(F8.2)') timer_fine
    !WRITE(myrank, '(F8.2)') timer_coarse
    !WRITE(myrank, '(F8.2)') timer_comm
    !CLOSE(myrank)

  end function PararealMPI

END MODULE parareal_mpi
