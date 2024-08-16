! まず変えないだろう定数はここに書く

module constants
  implicit none

  ! 基本の定数
  double precision, parameter :: pi = 3.14159265358979328d0 ! 円周率
  double precision, parameter :: e = 1.602176634d-19 ! 電気素量 [C]
  double precision, parameter :: m_e = 9.1093837015d-31 ! 電子の質量 [kg]
  double precision, parameter :: m_i = 1.67262192369d-27 ! 水素イオンの質量 [kg]

  ! 2次電子に関する定数
  double precision, parameter :: delta_max = 2.1d0
  double precision, parameter :: E_max = 1.5d2
  double precision, parameter :: theta = 5.235d-1   ! [rad] (30度)
  double precision, parameter :: T_sec = 2.5d0

  ! 後方散乱電子のパラメータ
  double precision, parameter :: eta = 1.0d-1
  double precision, parameter :: T_back = 2.5d0

  ! 光電子のパラメータ
  !double precision, parameter :: J_ph = 1d-20  ! 光電子電流密度[A/m^2]
  double precision, parameter :: J_ph = 4d-6  ! 光電子電流密度[A/m^2]
  double precision, parameter :: T_ph = 2.5d0

end module constants
