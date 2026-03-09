#!/usr/bin/env python3
import numpy as np
import parameters as P
from mod_md_subroutine import get_kinetic_ene, norm_seq
# from virial_estimator import virial_estimator  # call Virial_Estimator に対応

def ham_temp():
    """
    ハミルトニアン、温度、量子運動エネルギー、バスエネルギーの計算。
    Fortran サブルーチン `Ham_Temp` に相当。
    """

    # --- 運動エネルギー（classical beads） ---
    dkinetic = get_kinetic_ene()
    P.dkinetic = dkinetic

    # --- 温度計算（K） ---
    P.temp = 2.0 * dkinetic / (P.Natom * P.KtoAU * 3.0)
    P.temp /= P.Nbead

    # --- 量子運動エネルギー（qkinetic） ---
    ur_noncent = P.ur[:, :, 1:]
    ur_sq = np.sum(ur_noncent**2, axis=0)
    factqk = 0.5 * P.dnmmass[:, 1:] * P.omega_p2
    P.qkinetic = np.sum(factqk * ur_sq)

    # --- ハミルトニアン（暫定） ---
    P.hamiltonian = dkinetic + P.potential + P.qkinetic

    # --- バスエネルギー（非中心） ---
    if P.Ncent > 0:
        vrbath_sq = np.sum(P.vrbath[:, :, :, 1:]**2, axis=0)
        rbath_sum = np.sum(P.rbath[:, :, :, 1:], axis=0)
        P.ebath = (
            0.5 * np.sum(P.qmass[1:][None, None, :] * vrbath_sq)
            + P.gkt * np.sum(rbath_sum)
        )
    else:
        P.ebath = 0.0

    # --- バスエネルギー（中心モード用） ---
    if P.Ncent == 3:
        vrbc31_sq = np.sum(P.vrbc31**2, axis=0)
        rbc31_sum = np.sum(P.rbc31, axis=0)
        P.ebath_cent = (
            0.5 * np.sum(P.qmcent31[None, :] * vrbc31_sq)
            + P.gkt * np.sum(rbc31_sum)
        )
    else:
        P.ebath_cent = 0.0

    # --- 全体ハミルトニアン更新 ---
    P.hamiltonian += P.ebath + P.ebath_cent

    # --- ビリアル推定子計算 ---
    virial_estimator()


import numpy as np
import parameters as P

def virial_estimator():
    """
    Fortran の subroutine Virial_Estimator に対応。
    ビリアル推定を行い、P.E_Virial を更新。
    """

    diff = P.r - P.ur[:, :, [0]]
    E_Virial = 0.5 * np.sum(P.fr * diff)

    e_virial1 = 1.5 * float(P.Natom) / P.beta
    P.E_Virial = e_virial1 - E_Virial
