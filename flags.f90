module flags
  implicit none

    !0:parareal法の更新式、1:細かい積分の数値結果のみで更新
    INTEGER :: just_para_flag
    
    !環境データの設定
    !0:一定の環境下、1:線形の環境データ、2:sin関数で環境データ作成
    INTEGER :: change_env_flag
    
end module flags
