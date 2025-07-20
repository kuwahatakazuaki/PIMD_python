#!/usr/bin/env python3
import parameters as P
import datetime
from pathlib import Path
from datetime import datetime

def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def print_start(output_file="output.log"):
    """
    実行開始メッセージをログファイルに出力
    Parameters:
        output_file (str): 出力ファイル名
        restart (bool): 再スタートかどうか（True: 追記モード, False: 上書き）
    """
    mode = 'a' if P.Lrestart else 'w'
    with open(output_file, mode) as f:
        f.write("***********************\n")
        f.write("   Simulation Start!   \n")
        f.write("***********************\n")
        f.write(f" Simulation Started at {get_time()}\n\n")


def print_end(output_file="output.log"):
    """
    実行終了メッセージをログファイルに追記する
    Parameters:
        output_file (str): 出力ファイル名
    """
    with open(output_file, 'a') as f:
        f.write("***********************\n")
        f.write("    Simulation End!    \n")
        f.write("***********************\n")
        f.write(f"   Simulation Ended at {get_time()}\n\n")



def print_ini():
    ham_path = Path(P.dir_result)/ "ham.dat"
    fout_path = Path(P.Fout)

    if P.Nbead == 1:
        # ham.dat (replace mode: overwrite)
        with open(ham_path, "w") as f:
            f.write("# " + "*" * 87 + "\n")
            f.write("#  1Step 2 Hamiltonian   3 Temperature   4 Potential     5 DKinetic      6 EBath_Cent \n")
            f.write("# " + "*" * 87 + "\n")

        # Fout (append mode)
        with open(fout_path, "a") as f:
            f.write(" " + "*" * 95 + "\n")
            f.write("   Step  Hamiltonian  Temperature  Potential    DKinetic     EBath_Cent    Date & Time\n")
            f.write(" " + "*" * 95 + "\n")

    elif P.Nbead > 1:
        # ham.dat (replace mode: overwrite)
        with open(ham_path, "w") as f:
            f.write("# " + "*" * 132 + "\n")
            f.write(
                "# 1Step  2 Hamiltonian   3 Temperature   4 Potential     5 DKinetic      "
                "6 QKinetic      7 EBath         8 EBath_Cent    9 E_Virial\n"
            )
            f.write("# " + "*" * 132 + "\n")

        # Fout (append mode)
        with open(fout_path, "a") as f:
            f.write(" " + "*" * 121 + "\n")
            f.write(
                "   Step  Hamiltonian  Temperature  Potential    DKinetic     QKinetic     "
                "EBath        EBath_Cent     Date & Time\n"
            )
            f.write(" " + "*" * 121 + "\n")



def print_result_qm():
    """
    Fortran サブルーチン `print_result_qm` に対応。
    座標（coor.xyz）と、必要であれば力（force.dat）を出力。
    """
    # 出力ファイルパス
    coor_path = f"{P.dir_result}/coor.xyz"
    force_path = f"{P.dir_result}/force.dat"

    # === 座標出力 (coor.xyz) ===
    with open(coor_path, "a") as fcoor:
        fcoor.write(f"{P.Natom * P.Nbead:5d}\n")
        fcoor.write(f"{P.istepsv:10d}\n")
        for imode in range(P.Nbead):
            for iatom in range(P.Natom):
                x, y, z = P.r[:, iatom, imode] * P.AUtoAng
                fcoor.write(f"{P.alabel[iatom]:s} {x:15.9E} {y:15.9E} {z:15.9E}\n")

    # === 力出力 (force.dat) ===
    if P.Lsave_force:
        with open(force_path, "a") as ffor:
            ffor.write(f"# {P.istepsv:10d}\n")
            for imode in range(P.Nbead):
                for iatom in range(P.natom):
                    fx, fy, fz = P.fr[:, iatom, imode]
                    ffor.write(f"{fx:15.10f}{fy:15.10f}{fz:15.10f}\n")



def print_ham(Istep):
    """
    Print Hamiltonian and thermodynamic quantities.
    Corresponds to Fortran: print_ham
    """

    def print_ham_qm():
        if Istep % P.out_step == 0:
            with open(f"{P.dir_result}/ham.dat", "a") as f:
                f.write(f"{Istep:7d} {P.hamiltonian:16.8e} {P.temp:16.8e} {P.potential:16.8e} "
                        f"{P.dkinetic:16.8e} {P.qkinetic:16.8e} {P.ebath:16.8e} "
                        f"{P.ebath_cent:16.8e} {P.E_Virial:16.8e}\n")

        if Istep % 100 == 0:
            with open(P.Fout, "a") as f:
                f.write(f"{Istep:7d} {P.hamiltonian:13.5e} {P.temp:13.5e} {P.potential:13.5e} "
                        f"{P.dkinetic:13.5e} {P.qkinetic:13.5e} {P.ebath:13.5e} "
                        f"{P.ebath_cent:13.5e} {get_time():>20s}\n")

    def print_ham_cl():
        if Istep % P.out_step == 0:
            with open(f"{P.dir_result}/ham.dat", "a") as f:
                f.write(f"{Istep:7d} {P.hamiltonian:16.8e} {P.temp:16.8e} {P.potential:16.8e} "
                        f"{P.dkinetic:16.8e} {P.ebath_cent:16.8e}\n")

        if Istep % 100 == 0:
            with open(P.Fout, "a") as f:
                f.write(f"{Istep:7d} {P.hamiltonian:13.5e} {P.temp:13.5e} {P.potential:13.5e} "
                        f"{P.dkinetic:13.5e} {P.ebath_cent:13.5e} {get_time():>20s}\n")

    if P.Isimulation in [0, 1, 2]:
        print_ham_qm()
    elif P.Isimulation == 10:
        print_ham_cl()

