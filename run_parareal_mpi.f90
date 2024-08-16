PROGRAM run_parareal_mpi

USE mpi
USE parareal_mpi, only: InitializePararealMPI, FinalizePararealMPI, PararealMPI
USE parareal_params, only : do_io_first_coarse, do_io_fine, do_io_cor_coarse, be_verbose, Tstart, Tend, N_fine, N_coarse, Niter
use global_variables
use flags
use OML_currents


IMPLICIT NONE

CHARACTER(LEN=20) :: inputfilename
INTEGER :: ierr,myrank

DOUBLE PRECISION :: Temp1, Temp2, rho1,rho2
DOUBLE PRECISION :: cal_time, T0, T1
double precision :: C_inv ,phi_ini

integer:: i

namelist /io/ do_io_first_coarse, do_io_fine, do_io_cor_coarse, be_verbose
namelist /sim_time/ Tstart,Tend
namelist /num_steps/ N_fine, N_coarse
namelist /num_iter/ Niter
namelist /flag/ just_para_flag, change_env_flag
namelist /env_data/ Temp1, Temp2, rho1, rho2, phi_ini

!> @todo docu

! Initialize
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

!コマンドライン引数から設定ファイル名を取得
call get_command_argument(1,inputfilename)

!設定ファイルから値取得
open(1,file=inputfilename)
read(1,io)
read(1,sim_time)
read(1,num_steps)
read(1,num_iter)
read(1,flag)
read(1,env_data)
close(1)



!電位の初期値を代入
do i = 1, 10
    phi(i) = phi_ini
end do
C_inv = 1.0d0/C(iconductor)

!parareal法の初期化
CALL InitializePararealMPI()


!ここから
T0 = MPI_WTIME()
phi(i) = PararealMPI(phi, DBLE(0.0), Tend, N_fine, N_coarse, Niter, C_inv, rho1, rho2, Temp1, Temp2, iconductor,be_verbose)
T1 = MPI_WTIME()

cal_time = T1 -T0

!データの標準出力
if(myrank == 0) then
    if(just_para_flag == 0) then
        print *, "Solution : parareal "
    else if(just_para_flag == 1) then
        print *, "Solution : just pararell in time "
    end if

    if(change_env_flag == 0) then
        print *, "Complementary Methods : Certain Environment "
        print *, "Initial value of phi : ",phi_ini
        print *, "Temperature : ", Temp1
        print *, "Density : ", rho1
        print *, "time constant : ", Time_Constant(iconductor,phi(i),rho1,Temp1)
    else if(change_env_flag == 1) then
        print *, "Complementary Methods : linear completion "
        print *, "Initial value of phi : ",phi_ini
        print *, "Temperature1 : ", Temp1
        print *, "Temperature2 : ", Temp2
        print *, "Density1 : ", rho1
        print *, "Density2 : ", rho2
        print *, "time constant : ", Time_Constant(iconductor,phi(i),(rho1+rho2)/2,(Temp1+Temp2)/2)
    else if(change_env_flag == 2) then
        print *, "Complementary Methods : sin completion "
        print *, "Initial value of phi : ",phi_ini
        print *, "Temperature1 : ", Temp1
        print *, "Temperature2 : ", Temp2
        print *, "Density1 : ", rho1
        print *, "Density2 : ", rho2
        print *, "time constant : ", Time_Constant(iconductor,phi(i),rho1,Temp1)
    end if
end if


! Finalize
CALL FinalizePararealMPI();
CALL MPI_FINALIZE(ierr)


END PROGRAM run_parareal_mpi
