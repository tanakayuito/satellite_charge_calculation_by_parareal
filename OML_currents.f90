module OML_currents
  use constants
  use global_variables
  use flags
  USE parareal_params, only : do_io_first_coarse, do_io_fine, do_io_cor_coarse
  implicit none

  CHARACTER(len=64) :: filename
  double precision :: x
  integer :: n

contains

  ! OMLで入射電子電流
  double precision function OML_I_e(phi, rho, T, iconductor)
    implicit none

    integer, intent(in) :: iconductor
    double precision, intent(in) :: phi, rho, T

    ! phiの正負で場合分け
    if (phi <= 0d0) then
      OML_I_e = -area_surface(iconductor)*rho*e*sqrt(e*T/2d0/pi/m_e)*exp(phi/T)
    else
      OML_I_e = -area_surface(iconductor)*rho*e*sqrt(e*T/2d0/pi/m_e)*((1d0 + phi/T)**beta(iconductor))
      ! OML_I_e = 0d0
    end if

  end function OML_I_e

  ! OMLで入射イオン電流
  double precision function OML_I_i(phi, rho, T, iconductor)
    implicit none

    integer, intent(in) :: iconductor
    double precision, intent(in) :: phi, rho, T

    ! phiの正負で場合分け
    if (phi <= 0d0) then
      OML_I_i = area_surface(iconductor)*rho*e*sqrt(e*T/2d0/pi/m_i)*((1d0 - phi/T)**beta(iconductor))
      
    else
      OML_I_i = area_surface(iconductor)*rho*e*sqrt(e*T/2d0/pi/m_i)*exp(-phi/T)
      ! OML_I_i = 0d0
    end if

  end function OML_I_i

  ! OMLで2次電子電流
  double precision function OML_I_se(phi, iconductor, I_e)
    implicit none

    integer, intent(in) :: iconductor
    double precision, intent(in) :: phi, I_e

    delta = 1.114d0*delta_max/cos(theta)*((E_max/T)**3.5d-1)&
           & *(1d0 - exp(-2.28d0*cos(theta)*((T/E_max)**1.35d0)))

    ! phiの正負で場合分け
    if (phi <= 0d0) then
      OML_I_se = -delta*I_e
      
    else
      OML_I_se = -delta*exp(-phi/T_sec)*I_e
    end if
    ! OML_I_se = 0.0d0

  end function OML_I_se

  ! OMLで後方散乱電子電流
  double precision function OML_I_back(phi, iconductor, I_e)
    implicit none

    integer, intent(in) :: iconductor
    double precision, intent(in) :: phi, I_e

    ! phiの正負で場合分け
    if (phi <= 0d0) then
      OML_I_back = -eta*I_e
     
    else
      OML_I_back = -eta*exp(-phi/T_back)*I_e
    end if
    ! OML_I_back = 0.0d0

  end function OML_I_back

  ! OMLで光電子電流
  double precision function OML_I_ph(phi, iconductor,light_flag)
    implicit none

    double precision, intent(in) :: phi
    integer, intent(in) :: iconductor
    integer, intent(in) :: light_flag

    !日光で場合分け
    !日光あり
    if (light_flag == 1) then
      ! phiの正負で場合分け
      if (phi <= 0d0) then
        OML_I_ph = light_area(iconductor)*J_ph
        
      else
        OML_I_ph = light_area(iconductor)*J_ph*exp(-phi/T_ph)
      end if
      ! OML_I_ph = 0.0d0
    !日光なし
    else if(light_flag == 0) then
      OML_I_ph = 0
    end if

  end function OML_I_ph

  !ルンゲクッタ法の計算
  double precision function RungeKutta_Method_Phi(dt,C_inv,I_e,I_i,I_se,I_back,I_ph,rho,T,iconductor,phi) 
    implicit none
    double precision, intent(in) :: phi,I_e,I_i,I_se,I_back,I_ph
    double precision, intent(in) :: dt, C_inv
    double precision, intent(in) :: rho, T
    integer, intent(in) :: iconductor
    double precision, dimension(4) :: k ! Runge-Kutta係数
    
    k(1) = dt*C_inv*(I_e + I_i + I_se + I_back + I_ph)
    k(2) = dt*C_inv*(OML_I_e(phi + k(1)/2.0d0, rho, T, iconductor)&
    & + OML_I_i(phi + k(1)/2.0d0, rho, T, iconductor)&
    & + OML_I_se(phi + k(1)/2.0d0, iconductor, OML_I_e(phi + k(1)/2.0d0, rho, T, iconductor))&
    & + OML_I_back(phi + k(1)/2.0d0, iconductor, OML_I_e(phi + k(1)/2.0d0, rho, T, iconductor))&
    & + OML_I_ph(phi + k(1)/2.0d0, iconductor, light_flag))
    k(3) = dt*C_inv*(OML_I_e(phi + k(2)/2.0d0, rho, T, iconductor)&
    & + OML_I_i(phi + k(2)/2.0d0, rho, T, iconductor)&
    & + OML_I_se(phi + k(2)/2.0d0, iconductor, OML_I_e(phi + k(2)/2.0d0, rho, T, iconductor))&
    & + OML_I_back(phi + k(2)/2.0d0, iconductor, OML_I_e(phi + k(2)/2.0d0, rho, T, iconductor))&
    & + OML_I_ph(phi + k(2)/2.0d0, iconductor, light_flag))
    k(4) = dt*C_inv*(OML_I_e(phi + k(3), rho, T, iconductor)&
    & + OML_I_i(phi + k(3), rho, T, iconductor)&
    & + OML_I_se(phi + k(3), iconductor, OML_I_e(phi + k(3), rho, T, iconductor))&
    & + OML_I_back(phi + k(3), iconductor, OML_I_e(phi + k(3), rho, T, iconductor))&
    & + OML_I_ph(phi + k(3), iconductor, light_flag))

    RungeKutta_Method_Phi = phi + (k(1) + k(2)*2.0d0 + k(3)*2.0d0 + k(4))/6.0d0
    
  end function RungeKutta_Method_Phi  


  subroutine RungeKutta_Method_Phi_For_Parareal(phi, Time0, Time1, Nsteps, C_inv, iconductor, rho_n, rho_n_pl1, T_n, T_n_pl1, myrank, Nproc, iter, num_integ) 
    double precision,                INTENT(INOUT):: phi 
    DOUBLE PRECISION,                INTENT(IN)    :: Time0, Time1
    INTEGER,                         INTENT(IN) :: Nsteps
    double precision ,               INTENT(IN)::  C_inv
    !double precision,                INTENT(IN) :: rho, T
    INTEGER,                         INTENT(IN) :: iconductor,myrank, Nproc,iter, num_integ !num_integ 1:first coarse 2:fine 3:correct coarse
    double precision,                 INTENT(IN):: rho_n, rho_n_pl1, T_n, T_n_pl1
    
    DOUBLE PRECISION :: dt_para
    INTEGER :: n,i
    !double precision :: I_e

    INTEGER :: sig_dig
    
    dt_para = (Time1-Time0)/DBLE(Nsteps)

    if(num_integ == 1 .AND. do_io_first_coarse == .TRUE.) then
      WRITE(filename, '(A,I0.2,A,I0.2,A,I0.2,A)') 'phi_first_coarse_', myrank, '_', Nproc, '_mpi.dat'
      OPEN(UNIT=myrank, FILE=filename, ACTION='write', STATUS='replace')
      WRITE(myrank, '(F35.25)') phi
    else if(num_integ == 3 .AND. do_io_fine == .TRUE.) then
      WRITE(filename, '(A,I0.2,A,I0.2,A,I0.2,A)') 'phi_fine_', iter, '_', myrank, '_', Nproc, '_mpi.dat'
      OPEN(UNIT=myrank, FILE=filename, ACTION='write', STATUS='replace')
      WRITE(myrank, '(F35.25)') phi
    else if(num_integ == 4 .AND. do_io_cor_coarse == .TRUE.) then
      WRITE(filename, '(A,I0.2,A,I0.2,A,I0.2,A)') 'phi_cor_coarse_', iter, '_', myrank, '_', Nproc, '_mpi.dat'
      OPEN(UNIT=myrank, FILE=filename, ACTION='write', STATUS='replace')
      WRITE(myrank, '(F35.25)') phi
    end if

    !一定の環境下の場合
    if(change_env_flag == 0) then
      T = T_n
      rho = rho_n
    end if

    DO n=1,Nsteps
      !環境データの補完        
      if(num_integ == 1) then
        !線形補完の場合
        if(change_env_flag == 1) then
          T = (T_n_pl1 - T_n)/((Nsteps/myrank)*Nproc) * n + T_n
          rho = (rho_n_pl1 - rho_n)/((Nsteps/myrank)*Nproc) * n + rho_n

        !sin関数で補完する場合
        else if(change_env_flag == 2) then
          !2が周期決める
          T = T_n + abs(T_n_pl1 - T_n) * sin(n * pi/(( (Nsteps/myrank) * Nproc)/2))
          rho = rho_n + abs(rho_n_pl1 - rho_n) * sin( n * pi/(( (Nsteps/myrank) *Nproc)/2))

        end if
      else
        !線形補完の場合
        if(change_env_flag == 1) then
          T = (T_n_pl1 - T_n)/(Nsteps*Nproc) * (myrank*Nsteps + n) + T_n
          rho = (rho_n_pl1 - rho_n)/(Nsteps*Nproc) * (myrank*Nsteps + n) + rho_n

        !sin関数で補完する場合
        else if(change_env_flag == 2) then
          !2が周期決める
          T = T_n + abs(T_n_pl1 - T_n) * sin( (myrank*Nsteps + n) * pi/((Nsteps*Nproc)/2) )
          rho = rho_n + abs(rho_n_pl1 - rho_n)* sin( (myrank*Nsteps + n) * pi/((Nsteps*Nproc)/2))
        end if
      end if

      ! OMLで電流を求める

      I_e = OML_I_e(phi, rho, T, iconductor)
      I_i = OML_I_i(phi, rho, T, iconductor)
      I_se = OML_I_se(phi, iconductor, I_e)
      I_back = OML_I_back(phi, iconductor, I_e)
      I_ph = OML_I_ph(phi, iconductor, light_flag)

  
      ! 微分方程式を解く
      ! Runge-Kutta 4次
      phi = RungeKutta_Method_Phi(dt_para,C_inv,I_e,I_i,I_se,I_back,I_ph,rho,T,iconductor,phi)

      !結果出力
      if((num_integ == 1 .OR. num_integ == 2).AND. do_io_first_coarse == .TRUE.) then
        WRITE(myrank, '(F35.25)') phi
      else if(num_integ == 3 .AND. do_io_fine == .TRUE.) then
        WRITE(myrank, '(F35.25)') phi
      else if(num_integ == 4 .AND. do_io_cor_coarse == .TRUE.) then
        WRITE(myrank, '(F35.25)') phi
      end if
      
    END DO

    if(.not. num_integ == 1) then
      close(myrank)
    end if

  end subroutine RungeKutta_Method_Phi_For_Parareal

  subroutine RungeKutta_Method_Phi_For_time_sequential_cal(phi,Time0,Time1,Nsteps,C_inv,iconductor,rho_n, rho_n_pl1, T_n, T_n_pl1, myrank, Nproc, output_filename) 
    double precision,   INTENT(INOUT) :: phi 
    DOUBLE PRECISION,   INTENT(IN)    :: Time0, Time1
    INTEGER,            INTENT(IN)    :: Nsteps
    double precision,   INTENT(IN)    :: C_inv
    INTEGER,            INTENT(IN)    :: iconductor
    double precision,   INTENT(IN)    :: rho_n, rho_n_pl1, T_n, T_n_pl1
    INTEGER,            INTENT(IN)    :: myrank, Nproc
    CHARACTER(LEN=20),  INTENT(IN)    :: output_filename
    
    DOUBLE PRECISION :: dt_seq
    INTEGER :: n,i  

    dt_seq = (Time1-Time0)/DBLE(Nsteps)

    if(change_env_flag == 0) then
      T = T_n
      rho = rho_n
    end if

    DO n=1,Nsteps
      if(change_env_flag == 1) then
        T = (T_n_pl1 - T_n)/(Nsteps*Nproc) * (myrank*Nsteps + n) + T_n
        rho = (rho_n_pl1 - rho_n)/(Nsteps*Nproc) * (myrank*Nsteps + n) + rho_n
      else if(change_env_flag == 2) then
        T = T_n + abs(T_n_pl1 - T_n) * sin( (myrank*Nsteps + n) * pi/((Nsteps*Nproc)/2) )
        rho = rho_n + abs(rho_n_pl1 - rho_n)* sin( (myrank*Nsteps + n) * pi/((Nsteps*Nproc)/2))
      end if

      I_e = OML_I_e(phi, rho, T, iconductor)
      I_i = OML_I_i(phi, rho, T, iconductor)
      I_se = OML_I_se(phi, iconductor, I_e)
      I_back = OML_I_back(phi, iconductor, I_e)
      I_ph = OML_I_ph(phi, iconductor, light_flag)
  
      ! 微分方程式を解く
      ! Runge-Kutta 4次
      phi = RungeKutta_Method_Phi(dt_seq,C_inv,I_e,I_i,I_se,I_back,I_ph,rho,T,iconductor,phi)
      WRITE(4, '(F35.25)') phi
      
    END DO

  end subroutine RungeKutta_Method_Phi_For_time_sequential_cal
  

  !ニュートン法用の計算式１ Y_sunada ver210
  double precision function F_Phi(delta,iconductor,phi,rho,T)
    implicit none 
    double precision,intent(in) :: phi,rho,T,delta
    integer, intent(in):: iconductor
    double precision::Ie0,Ii0,Iph

    Ie0 = area_surface(iconductor) * rho * e * sqrt(e*T/2d0/pi/m_e)
    Ii0 = area_surface(iconductor) * rho * e * sqrt(e*T/2d0/pi/m_i)
    Iph = light_area(iconductor) * J_ph
    if(phi >= 0d0) then
      F_Phi = (delta*exp(-phi/T_sec) + eta*exp(-phi/T_back)-1d0) * Ie0 * (1d0+phi/T) + Ii0*exp(-phi/T) + Iph*exp(-phi/T_ph)
    else
      F_Phi = (delta+eta-1d0) * Ie0 * exp(phi/T) + (1d0-phi/T)*Ii0 + Iph
    end if
  end function F_Phi

  !ニュートン法用の計算式２ Y_sunada ver210
  double precision function F_Phi_d(delta,iconductor,phi,rho,T)
    implicit none 
    double precision,intent(in) :: phi,rho,T,delta
    integer, intent(in):: iconductor
    double precision::Ie0,Ii0,Iph

    Ie0 = area_surface(iconductor) * rho * e * sqrt(e*T/2d0/pi/m_e)
    Ii0 = area_surface(iconductor) * rho * e * sqrt(e*T/2d0/pi/m_i)
    Iph = light_area(iconductor) * J_ph
    if(phi >= 0d0) then
      F_Phi_d = -Ie0/T - Ii0/T*exp(-phi/T) - Iph/T_ph*exp(-phi/T_ph)  + delta * Ie0 *exp(-phi/T_sec) * (-1/T_sec + 1/T - phi/T/T_sec) + eta * Ie0 * exp(-phi/T_back)* (-1/T_back+ 1/T - phi/T/T_back)
    else
      F_Phi_d = ((delta+eta-1d0)*Ie0/T*exp(phi/T)-Ii0/T)
    end if
  end function F_Phi_d

  !ニュートン法による定常値の計算 Y_sunada ver210
  double precision function Newton_Method_Phi(iconductor,phi,rho,T)
    implicit none
    double precision,intent(in) :: rho,T,phi
    integer, intent(in) :: iconductor
    integer :: count = 0
    integer :: count_end = 1000000
    double precision :: phi_new,phi_old
    double precision :: decF,sig,acc
    
    decF = 1d0 !減速定数　初期値
    sig  = 4d0   !減速定数　規定値
    acc  = 1d-10 !収束精度
    delta = 1.114d0*delta_max/cos(theta)*((E_max/T)**3.5d-1)*(1d0 - exp(-2.28d0*cos(theta)*((T/E_max)**1.35d0)))
    phi_old = phi
    !メインループ
    !ニュートン法減速法　F(phi) = Ii - Ie + Iback + Iph + Ise0
    do while(.true.)
      count = count + 1
      !念の為の非常口
      if(count > count_end) then
        Newton_Method_Phi = phi_new
        exit
      end if
      phi_new = phi_old -F_Phi(delta,iconductor,phi_old,rho,T)/F_Phi_d(delta,iconductor,phi_old,rho,T) * decF
      !減速判定
      if(abs(F_Phi(delta,iconductor,phi_new,rho,T)) > (1d0-decF/sig)*abs(F_Phi(delta,iconductor,phi_old,rho,T)) )then
        decF = decF/2d0
        cycle
      end if
      !収束判定
      if(abs(F_Phi(delta,iconductor,phi_new,rho,T)) < acc) then
        Newton_Method_Phi = phi_new
        exit
      end if
      phi_old = phi_new
      decF = 1
    end do !メインループ終了
  end function Newton_Method_Phi

  !時定数の割り出し（定常値利用） Y_sunada ver210
  double precision function Time_Constant(iconductor,phi,rho,T)
    implicit none
    double precision,intent(in) :: rho,T,phi
    integer, intent(in) :: iconductor
    double precision :: Ie0,Ii0,Iph,delta
   
    delta = 1.114d0*delta_max/cos(theta)*((E_max/T)**3.5d-1)*(1d0 - exp(-2.28d0*cos(theta)*((T/E_max)**1.35d0)))
    Ie0 = area_surface(iconductor) * rho * e * sqrt(e*T/2d0/pi/m_e)
    Ii0 = area_surface(iconductor) * rho * e * sqrt(e*T/2d0/pi/m_i)
    Iph = light_area(iconductor) * J_ph
    if(phi>=0) then
      Time_Constant =  C(iconductor)/abs(Ii0/T/exp(phi/T) + Iph/T_ph/exp(phi/T_ph) + delta*Ie0/T_sec/exp(phi/T_sec) + eta*Ie0/T_back/exp(phi/T_back) - Ie0/T - delta*Ie0/T/exp(phi/T_sec) + delta*Ie0/T_sec*phi/T/exp(phi/T_sec) - eta*Ie0/T/exp(phi/T_back) + eta*Ie0/T_back*phi/T/exp(phi/T_back) )
    else
      Time_Constant =  C(iconductor)/abs(Ii0/T + (1-delta-eta)*Ie0/T * exp(phi/T))
    end if
  end function Time_Constant
end module OML_currents
