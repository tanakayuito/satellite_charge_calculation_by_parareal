constants.f90                   ：定数
flags.f90                       ：解法やデータ作成法を選択するフラッグ
global_variables.f90            ：global変数
parareal_params.f90             ：反復数やステップ数
OML_currents.f90                ：OML理論の電流の計算式、数値積分手法
parareal_mpi.f90                ：mpiを用いたparareal法のアルゴリズム
run_parareal_mpi.f90            ：parareal法による計算を行うメインプログラム
run_seq1.f90                    ：dt_fineと同じ時間刻み幅を用いて時間逐次計算を行うメインプログラム
run_seq2.f90                    ：dt_coarseと同じ時間刻み幅を用いて時間逐次計算を行うメインプログラム

create_fig_seq_results.py       ：run_seq1の帯電数値計算結果、run_seq2の帯電数値計算結果、run_seq1の結果に対するrun_seq2の結果の相対誤差のグラフ作成
create_fig_parareal_results.py  ：parareal法の帯電計算結果、run_seq1の結果に対するparareal法の結果の相対誤差のグラフ作成

params.inp(test_params.inp)     ：頻繁に変更するパラメータはここで設定する。入力データ。
内容
&io
  do_io_first_coarse (logical)  ：1回目の粗視化積分の結果を出力
  do_io_fine (logical)          ：厳密積分の結果出力
  do_io_cor_coarse (logical)    ：修正のため(2回目以降)の粗視化積分の結果を出力
  be_verbose (logical)          ：刻み幅や環境データなどの情報を出力
/

&sim_time
  Tstart (DOUBLE PRECISION)     ：シミュレーションの開始時刻
  Tend (DOUBLE PRECISION)       ：シミュレーションの終了時刻
/

&num_steps
  N_fine (int)                  ：1つの部分時間領域における厳密積分のステップ数
  N_coarse (int)                ：1つの部分時間領域における粗視化積分のステップ数
/

&num_iter
  Niter (int)                   ：反復数
/

&num_proc
  Nprocs_for_seq (int)          ：プロセス数。（Parareal法のプロセス数をここで決められるわけではない。ただし、ここにも正しく設定しないと時間逐次計算や可視化が上手くいかない）
/

&flag
  just_para_flag (int)          ：  0：parareal法の更新式(Qend <- G(y^k_n) + F(y^(k-1))_n - G(y^(k-1)_n) = y^k_(n+1))、1：細かい積分の数値結果のみで更新(Qend <- F(y^(k-1))_n = y^k_(n+1))
  change_env_flag (int)         ：  0：一定の環境下、1：線形な環境データ、2:sin関数の環境データ
/

&env_data
  Temp1 (DOUBLE PRECISION)      ：プラズマの温度の値。一定の環境下の場合はこの値を使用する。線形の場合、さらにsin関数の場合は開始時刻の値。 
  Temp2 (DOUBLE PRECISION)      ：プラズマの温度の値。一定の環境下の場合はこの値は無関係。線形の場合は最終時刻の値。つまり、Temp1とTemp2の間で線形補完する。sin関数の場合はTemp1+(Temp2-Temp1)*sin(x)
  rho1 (DOUBLE PRECISION)       ：プラズマの密度の値。一定の環境下の場合はこの値を使用する。線形の場合、さらにsin関数の場合は開始時刻の値。
  rho2 (DOUBLE PRECISION)       ：プラズマの温度の値。一定の環境下の場合はこの値は無関係。線形の場合は最終時刻の値。つまり、rho1とrho2の間で線形補完する。sin関数の場合はrho1+(rho2-rho1)*sin(x)
  phi_ini (DOUBLE PRECISION)    ：電位の初期値
/

&visualization
  match_y_axis_range (logical)  ：parareal法の結果のy軸の範囲を反復数1回目の結果に合わせる。
/