#!/usr/bin/env python3
import numpy as np
import parameters as P
# from set_allocate import allocate_arrays

def calc_constant():
    """
    物理定数やシミュレーション定数を計算し、parameters.py の変数を更新する。
    """
    # arrays = allocate_arrays(P.Natom, P.Nbead, P.Nnhc, P.Nys, P.Ncent, P.Isimulation)
    freq1 = 10.0  # [THz]

    # タイムステップ変換：fs → 原子単位
    P.dt *= P.fs2AU

    # 逆温度 β（原子単位）
    P.beta = 1.0 / (P.KtoAU * P.temperature)

    # 系の固有振動数（[fs^-1] → AU^-1）
    P.omega_system = 2.0 * P.pi / (freq1 * P.fs2AU)

    # 核運動の平均運動エネルギー（古典）
    P.gnkt = 3.0 * float(P.Natom) / P.beta
    P.gkt  = 1.0 / P.beta

    # 一時ファイルの出力先のパス（文字列処理）
    P.address0    = P.dir_scr.strip() + "/"
    P.laddress    = len(P.dir_scr.strip()) + 1
    P.addresstmp  = P.address0  # 同じ意味

    # パスインテグラル関連の角振動数
    dp = float(P.Nbead)
    P.dp     = dp
    P.dp_inv = 1.0 / dp
    P.omega_p2 = dp / (P.beta * P.beta)
    P.omega2   = P.omega_system ** 2
