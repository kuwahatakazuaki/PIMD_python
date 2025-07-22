#!/usr/bin/env python3
import parameters as P
import numpy as np
from mod_md_subroutine import get_kinetic_ene, norm_seq
# from ase import Atoms
# from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

def setup_time_mass():
    if P.Isimulation == 10:
        P.dt_ref = P.dt
    else:
        P.dt_ref = P.dt / float(P.Nref)

    # セントロイド熱浴質量の設定
    if P.Ncent == 3:
        P.qmcent31[0] = 3.0 * float(P.Natom) / P.beta / P.omega2
        for inhc in range(1, P.Nnhc):
            P.qmcent31[inhc] = 1.0 / P.beta / P.omega2

    # 非セントロイドモードの熱浴質量
    if P.Isimulation != 10:
        P.qmass[0] = 0.0
        for imode in range(1, P.Nbead):
            P.qmass[imode] = 1.0 / P.beta / P.omega_p2

    # セントロイドMDの場合は質量スケーリング
    P.gamma2 = P.gamma1 ** 2
    if P.Isimulation == 2:
        for imode in range(1, P.Nbead):
            P.qmass[imode] *= P.gamma2

    # 高次分解（Yoshida–Suzuki）重みの設定
    if P.Nys == 1:
        P.ysweight[0] = 1.0
    elif P.Nys == 3:
        w1 = 1.0 / (2.0 - 2.0**(1.0/3.0))
        P.ysweight[0] = w1
        P.ysweight[1] = 1.0 - 2.0 * w1
        P.ysweight[2] = w1
    elif P.Nys == 5:
        w1 = 1.0 / (4.0 - 4.0**(1.0/3.0))
        P.ysweight[0] = w1
        P.ysweight[1] = w1
        P.ysweight[2] = 1.0 - 4.0 * w1
        P.ysweight[3] = w1
        P.ysweight[4] = w1


def normal_mode():
    dp = float(P.Nbead)
    sqp = np.sqrt(dp)
    sqpinv = 1.0 / sqp
    dnorm = np.sqrt(2.0 / dp)

    # u(imode,1) = 1/√N, u(imode,N) = ±1/√N の設定
    dum = -1.0
    for imode in range(P.Nbead):
        P.u[imode, 0] = sqpinv
        P.u[imode, P.Nbead - 1] = dum * sqpinv
        dum *= -1.0

    # コサイン・サインの成分設定（Fourier-like行列）
    for imode in range(1, (P.Nbead - 2) // 2 + 1):
        di = float(imode)
        for jmode in range(P.Nbead):
            dj = float(jmode + 1)  # Fortranは1始まり、Pythonは0始まり
            P.u[jmode, 2 * imode - 1] = dnorm * np.cos(2.0 * np.pi * di * dj / dp)
            P.u[jmode, 2 * imode]     = dnorm * np.sin(2.0 * np.pi * di * dj / dp)

    # uinv の転置版として初期化
    P.uinv = P.u.T.copy()

    # tnm と tnminv の初期化
    P.tnm = sqp * P.u
    P.tnminv = sqpinv * P.uinv


def init_mass():
    twopi = 2.0 * np.pi
    dp = float(P.Nbead)

    # dnmmass[:, 0] = 0.0（Pythonの0-indexに合わせる）
    P.dnmmass[:, 0] = 0.0

    for iatom in range(P.Natom):
        # 最後のモード（nbead番目）用：PythonではNbead-1
        P.dnmmass[iatom, P.Nbead - 1] = 4.0 * dp * P.physmass[iatom]
        for imode in range(1, (P.Nbead - 2) // 2 + 1):
            di = float(imode)
            val = 2.0 * (1.0 - np.cos(twopi * di / dp)) * dp * P.physmass[iatom]
            P.dnmmass[iatom, 2 * imode - 1] = val
            P.dnmmass[iatom, 2 * imode]     = val

    # fictmass の設定
    if P.Isimulation in [1, 10]:  # RPMD or classical MD
        for iatom in range(P.Natom):
            P.fictmass[iatom, :] = P.physmass[iatom]
    else:  # PIMD or CMD
        for iatom in range(P.Natom):
            P.fictmass[iatom, 0] = P.physmass[iatom]
            for imode in range(1, P.Nbead):
                P.fictmass[iatom, imode] = P.gamma2 * P.dnmmass[iatom, imode]


def nm_position():
    """
    非セントロイドモード（j ≥ 2）の初期座標 ur にガウス揺らぎを設定。
    Lrandom_coor = False の場合はゼロに初期化。
    """
    if P.Lrandom_coor:
        for j in range(1, P.Nbead):  # Fortran: 2〜Nbead → Python: 1〜Nbead-1
            for i in range(P.Natom):
                usigma = np.sqrt(1.0 / P.beta / P.omega_p2 / P.dnmmass[i, j])
                for k in range(3):  # x, y, z
                    P.ur[k, i, j] = usigma * np.random.normal(0, 1)
    else:
        P.ur[:, :, 1:P.Nbead] = 0.0  # Python: 2〜Nbead に対応


def remove_translation_rotation():
    delta = 1.0e-5

    for j in range(P.Nbead):
        comr = np.zeros(3)
        sumvr = np.zeros(3)
        ang_vel = np.zeros(3)

        # --- 中心座標 & 並進速度の計算 ---
        for i in range(P.Natom):
            sumvr += P.vur[:, i, j] * P.fictmass[i, j]
            comr += P.ur[:, i, j] * P.fictmass[i, j]
        totmas = np.sum(P.fictmass[:, j])
        sumvr /= totmas
        comr /= totmas

        # --- 慣性テンソルの計算 ---
        inertia = np.zeros((3, 3))
        for i in range(P.Natom):
            r = P.ur[:, i, j] - comr
            m = P.fictmass[i, j]
            inertia += m * (np.dot(r, r) * np.identity(3) - np.outer(r, r))

        # --- 慣性テンソルの逆行列計算 ---
        detmat = np.linalg.det(inertia)
        if detmat <= delta:
            dinvmat = np.zeros((3, 3))
            for k in range(3):
                if abs(inertia[k, k]) > delta:
                    dinvmat[k, k] = 1.0 / inertia[k, k]
        else:
            dinvmat = np.linalg.inv(inertia)

        # --- 角運動量ベクトルの計算 ---
        moment = np.zeros(3)
        for i in range(P.Natom):
            r = P.ur[:, i, j] - comr
            v = P.vur[:, i, j]
            m = P.fictmass[i, j]
            moment += m * np.cross(r, v)

        # --- 角速度ベクトル ---
        ang_vel = dinvmat @ moment

        # --- 並進と回転の速度成分を除去 ---
        for i in range(P.Natom):
            r = P.ur[:, i, j] - comr
            P.vur[:, i, j] -= sumvr
            P.vur[:, i, j] -= np.cross(ang_vel, r)


def init_velocity():
    """
    各原子・各ビーズに対して Maxwell 分布に基づいた初速度を設定する。
    結果は P.vur[3, Natom, Nbead] に格納される。
    """
    np.random.seed(P.Iseed)
    for imode in range(P.Nbead):
        for iatom in range(P.Natom):
            vsigma = np.sqrt(1.0 / P.beta / P.fictmass[iatom, imode])
            for k in range(3):  # x, y, z
                P.vur[k, iatom, imode] = vsigma * np.random.normal(0, 1)

    remove_translation_rotation()
    

def init_bath():
    """
    ノーズ・フーバー熱浴の初期化。
    初期速度を Maxell 分布（平均0, 分散1）からサンプリングし、vrbath および vrbc31 に代入。
    """
    # セントロイドの初期値ゼロ化（imode=1）
    P.rbath[:, :, :, 0] = 0.0
    P.vrbath[:, :, :, 0] = 0.0

    # 非セントロイドモードの熱浴速度初期化
    for imode in range(1, P.Nbead):  # Fortran: 2〜 → Python: 1〜
        vsigma = np.sqrt(1.0 / P.beta / P.qmass[imode])
        for inhc in range(P.Nnhc):
            for iatom in range(P.Natom):
                P.vrbath[0, iatom, inhc, imode] = vsigma * np.random.normal(0, 1)
                P.vrbath[1, iatom, inhc, imode] = vsigma * np.random.normal(0, 1)
                P.vrbath[2, iatom, inhc, imode] = vsigma * np.random.normal(0, 1)
                P.rbath[:, iatom, inhc, imode] = 0.0

    # セントロイドモード用熱浴（Ncent == 3 のとき）
    if P.Ncent == 3:
        for inhc in range(P.Nnhc):
            for iatom in range(P.Natom):
                vsigma = np.sqrt(1.0 / P.beta / P.qmcent31[inhc])
                P.vrbc31[0, iatom, inhc] = vsigma * np.random.normal(0, 1)
                P.vrbc31[1, iatom, inhc] = vsigma * np.random.normal(0, 1)
                P.vrbc31[2, iatom, inhc] = vsigma * np.random.normal(0, 1)
                P.rbc31[:, iatom, inhc] = 0.0


def init_bath_cl():
    if P.Ncent == 3:
        for inhc in range(P.Nnhc):
            for iatom in range(P.Natom):
                vsigma = np.sqrt(1.0 / P.beta / P.qmcent31[inhc])
                P.vrbc31[0, iatom, inhc] = vsigma * np.random.normal(0, 1)
                P.vrbc31[1, iatom, inhc] = vsigma * np.random.normal(0, 1)
                P.vrbc31[2, iatom, inhc] = vsigma * np.random.normal(0, 1)
                P.rbc31[:, iatom, inhc] = 0.0



def ham_temp_cl():
    """
    Classical Hamiltonian and thermostat bath energy computation.
    Corresponds to Fortran subroutine Ham_Temp_Classical.
    """
    # Kinetic energy
    dkinetic = get_kinetic_ene()  # assumes get_kinetic_ene() available
    P.dkinetic = dkinetic

    # Temperature (not often used in classical MD, but computed for logging)
    P.temp = 2.0 * dkinetic / (P.Natom * P.KtoAU * 3.0)
    P.temp /= P.Nbead

    # Hamiltonian (classical part)
    P.hamiltonian = dkinetic + P.potential

    # Centroid bath energy (if thermostat on centroid)
    ebath_cent = 0.0
    if P.Ncent == 3:
        for inhc in range(P.Nnhc):
            for iatom in range(P.Natom):
                ebath_cent += (
                    0.5 * P.qmcent31[inhc] * norm_seq(P.vrbc31[:, iatom, inhc])
                    + P.gkt * np.sum(P.rbc31[:, iatom, inhc])
                )

    # Add bath contribution
    P.ebath_cent = ebath_cent
    P.hamiltonian += ebath_cent
