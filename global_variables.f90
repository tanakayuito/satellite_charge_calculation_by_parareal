! いろんなところで使う変数と配列はここに書く

module global_variables
  implicit none

  ! MHDのデータを扱うための配列(密度[1/m^3]と温度[V]と時刻[sec.])
  !double precision, dimension(1000) :: mhd_rho, mhd_T
  double precision, dimension(1400) :: mhd_rho, mhd_T

  ! 電位計算のためのやつ
  double precision :: rho, T ! 密度[1/m^3]と温度[eV]
  double precision :: I_e, I_i, I_se, I_back, I_ph, I_sum ! 電流 [A]
  double precision :: delta ! 2次電子放出係数
  integer :: iconductor = 1! メインループ用

  ! conductors
  integer :: num_conductors = 1 ! 電位を計算する導体の数(最大10)
  double precision, dimension(10) :: C = 4.0d-12 ! 静電容量 [F]
  !double precision, dimension(10) :: phi = -9.99953669475895d3 ! 衛星電位 [V]
  double precision, dimension(10) :: phi ! 衛星電位 [V]
  double precision, dimension(10) :: area_surface = 9.6d1  ! 衛星の表面積[m^2]
  double precision, dimension(10) :: light_area = 1.6d1  ! 光が当たる面積[m^2]
  double precision, dimension(10) :: beta = 1d0 ! 衛星の形状で決まる係数(0<beta<1)

  ! ニュートン法と時定数
  double precision :: T_constant, Phi_constant

  !日光当たってるか当たってないか判断する変数。当たってたら1,当たってなかったら0
  integer :: light_flag = 1

end module global_variables
