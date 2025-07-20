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
    qkinetic = 0.0
    for imode in range(1, P.Nbead):  # 2 to Nbead in Fortran → 1 to Nbead-1 in Python
        for iatom in range(P.Natom):
            factqk = 0.5 * P.dnmmass[iatom, imode] * P.omega_p2
            qkinetic += factqk * norm_seq(P.ur[:, iatom, imode])
    P.qkinetic = qkinetic

    # --- ハミルトニアン（暫定） ---
    P.hamiltonian = dkinetic + P.potential + qkinetic

    # --- バスエネルギー（非中心） ---
    ebath = 0.0
    if P.Ncent > 0:
        for imode in range(1, P.Nbead):
            qdummy = P.qmass[imode]
            for inhc in range(P.Nnhc):
                for iatom in range(P.Natom):
                    ebath += (
                        0.5 * qdummy * norm_seq(P.vrbath[:, iatom, inhc, imode])
                        + P.gkt * np.sum(P.rbath[:, iatom, inhc, imode])
                    )
    P.ebath = ebath

    # --- バスエネルギー（中心モード用） ---
    ebath_cent = 0.0
    if P.Ncent == 3:
        for inhc in range(P.Nnhc):
            for iatom in range(P.Natom):
                ebath_cent += (
                    0.5 * P.qmcent31[inhc] * norm_seq(P.vrbc31[:, iatom, inhc])
                    + P.gkt * np.sum(P.rbc31[:, iatom, inhc])
                )
    P.ebath_cent = ebath_cent

    # --- 全体ハミルトニアン更新 ---
    P.hamiltonian += ebath + ebath_cent

    # --- ビリアル推定子計算 ---
    virial_estimator()


import numpy as np
import parameters as P

def virial_estimator():
    """
    Fortran の subroutine Virial_Estimator に対応。
    ビリアル推定を行い、P.E_Virial を更新。
    """

    E_Virial = 0.0

    for imode in range(P.Nbead):
        for iatom in range(P.Natom):
            diff = P.r[:, iatom, imode] - P.ur[:, iatom, 0]  # ur(:,:,1) → ur[:, :, 0]
            E_Virial += np.dot(P.fr[:, iatom, imode], diff)

    E_Virial /= 2.0

    e_virial1 = 1.5 * float(P.Natom) / P.beta
    P.E_Virial = e_virial1 - E_Virial
