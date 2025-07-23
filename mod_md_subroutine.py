#!/usr/bin/env python3
# utility.py
import numpy as np
import parameters as P

def norm_seq(vec):
    """
    ベクトル vec のノルム二乗（x^2 + y^2 + z^2）を返す
    """
    return np.sum(vec ** 2)

def get_kinetic_ene():
    """
    系の全運動エネルギーを返す。
    kine = 0.5 * Σ_ij [fictmass(i,j) * |v(i,j)|^2]
    """
    kine = 0.0
    for i in range(P.Natom):
        for j in range(P.Nbead):
            kine += P.fictmass[i, j] * norm_seq(P.vur[:, i, j])
    return 0.5 * kine


def temp_ctr():
    """
    初期速度 vur を運動エネルギーから計算された温度でスケーリングして、
    設定温度 P.temperature に合わせる。
    """
    tempi = get_kinetic_ene()

    # 系の実温度を算出
    tempi = 2.0 * tempi / (3.0 * float(P.Natom)) / P.KtoAU
    tempi = tempi / float(P.Nbead)
    temp_scale = np.sqrt(P.temperature / tempi)

    # 初速度をスケーリング
    P.vur *= temp_scale


def nmtrans_ur2r():
    """
    通常座標 r に正準モード ur を変換する（r(xyz,i,j) = Σ_k tnm(j,k) * ur(xyz,i,k)）。
    """
    P.r[:, :, :] = np.einsum("jk, xik -> xij", P.tnm, P.ur)
    # P.r[:, :, :] = 0.0
    # for i in range(P.Natom):
    #     for j in range(P.Nbead):
    #         for xyz in range(3):
    #             P.r[xyz, i, j] = np.dot(P.tnm[j, :], P.ur[xyz, i, :])

def nmtrans_r2ur():
    """
    正準モード ur に通常座標 r を変換・加算する（ur(xyz,i,j) += Σ_k tnminv(j,k) * r(xyz,i,k)）。
    """
    P.ur = np.einsum("jk, xik -> xij", P.tnminv, P.r)
    # P.ur[:,:,:] = 0.0
    # for i in range(P.Natom):
    #     for j in range(P.Nbead):
    #         for xyz in range(3):
    #             P.ur[xyz, i, j] += np.dot(P.tnminv[j, :], P.r[xyz, i, :])


def nmtrans_fr2fur():
    """
    fur(iatom, imode) = sum_j fr(iatom, jmode) * tnm(jmode, imode)
    フーリエ変換的に fr を正準モード座標 fur に変換する。
    """
    P.fur[:, :, :] = np.einsum("kij,jm->kim", P.fr, P.tnm)
    # P.fur[:, :, :] = 0.0
    # for iatom in range(P.Natom):
    #     for imode in range(P.Nbead):
    #         for jmode in range(P.Nbead):
    #             P.fur[:, iatom, imode] += P.fr[:, iatom, jmode] * P.tnm[jmode, imode]


def get_force_ref():
    """
    Python equivalent of the Fortran subroutine Getforce_Ref.
    Computes reference force `fur_ref` from normal-mode coordinates `ur`.
    """
    P.fur_ref[:, :, 0] = 0.0
    # for j in range(1, P.Nbead):  # Fortran index 2 to Nbead → Python 1 to Nbead-1
    #     for i in range(P.Natom):
    #         P.fur_ref[:, i, j] = -P.dnmmass[i, j] * P.omega_p2 * P.ur[:, i, j]
    j_slice = slice(1, P.Nbead)
    P.fur_ref[:, :, j_slice] = -P.dnmmass[:, j_slice] * P.omega_p2 * P.ur[:, :, j_slice]


def Uupdate():
    """
    Update position ur using time step dt_ref and velocity vur
    ur[:, iatom, imode] += dt_ref * vur[:, iatom, imode]
    """
    # for imode in range(P.Nbead):
    #     for iatom in range(P.Natom):
    #         P.ur[:, iatom, imode] += P.dt_ref * P.vur[:, iatom, imode]
    P.ur[:,:,:] += P.dt_ref * P.vur[:,:,:]

def Vupdate():
    """
    Update velocity vur using time step dt and force fur
    vur[:, iatom, imode] += 0.5 * dt * fur[:, iatom, imode] / fictmass[iatom, imode]
    """
    # for imode in range(P.Nbead):
    #     for iatom in range(P.Natom):
    #         P.vur[:, iatom, imode] += 0.5 * P.dt * P.fur[:, iatom, imode] / P.fictmass[iatom, imode]
    P.vur += 0.5 * P.dt * P.fur / P.fictmass


def Vupdate_Ref():
    """
    Reference velocity update:
    vur[:, i, j] += 0.5 * dt_ref * fur_ref[:, i, j] / fictmass[i, j]
    for j = 2 to Nbead
    """
    # for j in range(1, P.Nbead):  # Fortranの2からnbead → Pythonの1からNbead-1
    #     for i in range(P.Natom):
    #         P.vur[:, i, j] += 0.5 * P.dt_ref * P.fur_ref[:, i, j] / P.fictmass[i, j]
    j_slice = slice(1, P.Nbead)
    P.vur[:, :, j_slice] += 0.5 * P.dt_ref * P.fur_ref[:, :, j_slice] / P.fictmass[:, j_slice]
