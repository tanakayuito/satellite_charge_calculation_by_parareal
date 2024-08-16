PROGRAM main

USE mpi
use oml_currents
use global_variables
use flags


IMPLICIT NONE

CHARACTER(LEN=20) :: inputfilename, outputfilename

double precision, dimension(10) :: phi_tmp
double precision :: C_inv
double precision :: tstart_myslice,dt_slice,tend_myslice

real :: start_time, end_time, elapsed_time
DOUBLE PRECISION :: T0,T1, Tall

integer :: ierr,myrank
integer:: i

DOUBLE PRECISION :: Tstart,Tend
INTEGER :: N_fine, N_coarse
INTEGER :: Nprocs_for_seq
DOUBLE PRECISION :: Temp1, Temp2, rho1,rho2, phi_ini


namelist /sim_time/ Tstart,Tend
namelist /num_steps/ N_fine, N_coarse
namelist /num_proc/ Nprocs_for_seq
namelist /flag/ just_para_flag, change_env_flag
namelist /env_data/ Temp1, Temp2, rho1, rho2, phi_ini


!> @todo docu

! Initialize
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

!コマンドライン引数からファイル名を取得
call get_command_argument(1,inputfilename)
call get_command_argument(2,outputfilename)

!設定ファイルから値取得
open(1,file=inputfilename)
read(1,sim_time)
read(1,num_steps)
read(1,num_proc)
read(1,flag)
read(1,env_data)
close(1)
C_inv = 1.0d0/C(iconductor)


! Load initial data
do i = 1, 10    
    phi_tmp(i) = phi_ini
    phi(i) = phi_ini
end do

Tall = Tend - Tstart
dt_slice  = Tall/DBLE(Nprocs_for_seq)
    
OPEN(UNIT=4, FILE= outputfilename, ACTION='write', STATUS='replace')
WRITE(4, '(F35.25)') phi(iconductor)

T0 = MPI_WTIME()

do i = 0,Nprocs_for_seq-1
    tstart_myslice = Tstart + DBLE(i)*dt_slice
    tend_myslice   = Tstart + DBLE(i+1)*dt_slice
    CALL RungeKutta_Method_Phi_For_time_sequential_cal(phi(iconductor),tstart_myslice, tend_myslice, N_fine,C_inv,iconductor,rho1, rho2, Temp1, Temp2, i, Nprocs_for_seq, outputfilename) 
end do

T1 = MPI_WTIME()
close(4)

! 経過時間を計算
elapsed_time = T1 - T0

! 数値結果出力ファイルを閉じる
CALL MPI_FINALIZE(ierr)

END PROGRAM main
