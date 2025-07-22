#!/usr/bin/env python3
# force_nnp_matlantis.py
import numpy as np
import os
import subprocess
import parameters as P
import run_matlantis
from utility import program_abort

def input_periodic():
    Fout = np.empty(P.Nbead,dtype=object)
    for Imode in range(P.Nbead):
        fname = f"str{Imode+1:05d}.vasp"
        Fout[Imode] = fname  # Fout[] は Nbead 要素のリストと仮定

    for Imode in range(P.Nbead):
        src = "LATTICE"
        dst = os.path.join(P.addresstmp, Fout[Imode])
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        os.system(f"cp {src} {dst}")

        with open(dst, "a") as f:
            for i in range(P.Natom):
                pos = P.r[:, i, Imode] * P.AUtoAng
                f.write(f"{pos[0]:20.12f} {pos[1]:20.12f} {pos[2]:20.12f}\n")
            f.write("\n")

def input_nonperio():
    Fout = np.empty(P.Nbead,dtype=object)
    for Imode in range(P.Nbead):
        fname = f"str{Imode+1:05d}.xyz"
        Fout[Imode] = fname

    for Imode in range(P.Nbead):
        dst = os.path.join(P.addresstmp, Fout[Imode])
        os.makedirs(os.path.dirname(dst), exist_ok=True)

        with open(dst, "w") as f:
            f.write(f"{P.Natom}\n")
            f.write("Properties=species:S:1:pos:R:3 pbc=\"F F F\"\n")
            for i in range(P.Natom):
                pos = P.r[:, i, Imode] * P.AUtoAng
                label = P.alabel[i]
                f.write(f"{label:2s} {pos[0]:20.12f} {pos[1]:20.12f} {pos[2]:20.12f}\n")


def force_nnp_matlantis():
    Fforce = "forces.out"
    Fenergy = "energy.out"
    Nfile = P.Nbead

    if P.Lperiodic:
        input_periodic()
    else:
        input_nonperio()

    # 古い出力ファイルがあれば削除
    if os.path.exists(Fforce):
        os.remove(Fforce)
    if os.path.exists(Fenergy):
        os.remove(Fenergy)

    run_matlantis.run_cal()

    # # 力の読み込み
    # if not os.path.exists(Fforce):
    #     program_abort(f'ERROR!!! There is no "{Fforce}"')

    # with open(Fforce, "r") as f:
    #     for j in range(P.Nbead):
    #         for i in range(P.Natom):
    #             line = f.readline()
    #             P.fr[:, i, j] = np.fromstring(line, sep=" ")

    # # エネルギーの読み込み
    # if not os.path.exists(Fenergy):
    #     program_abort(f'There is no "{Fenergy}"')

    # with open(Fenergy, "r") as f:
    #     for j in range(P.Nbead):
    #         P.Eenergy[j] = float(f.readline())

    # P.fr *= P.eVAng2AU * P.dp_inv
    # P.Eenergy *= P.eVtoAU
    # P.potential = np.sum(P.Eenergy) * P.dp_inv


