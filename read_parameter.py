#!/usr/bin/env python3
import os
import numpy as np
import parameters as P

def read_parameter(filename="input.inp"):
    """
    入力ファイル input.inp を読み込み、各種パラメータを個別変数に代入する。
    $Coords ～ $end Coords の行数から Natom を自動的に設定する。
    """

    # global Natom, Nbead, Nstep, temperature, dt
    # global Isimulation, Nref, out_step, Nys, Nnhc, Ncent
    # 初期化
    # Natom = Nbead = Nstep = Nref = out_step = Nys = Nnhc = Ncent = None
    # temperature = dt = Isimulation = None

    if not os.path.exists(filename):
        raise FileNotFoundError(f'ERROR!!: There is no input file named "{filename}"')

    with open(filename, 'r') as f:
        lines = f.readlines()

    # パラメータと Coords 検出用
    i = 0
    coords_start = coords_end = -1

    def get_value():
        nonlocal i
        i += 1
        val = lines[i].strip()
        return float(val) if '.' in val or 'e' in val.lower() else int(val)
    
    def get_bool():
        nonlocal i
        i += 1
        val = lines[i].strip().lower()
        return val in ['true', '1', 'yes']


    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("!"):
            i += 1
            continue

        # print(line, P.Lrestart, P.Lperiodic)
        # if line.startswith("$Lrestart"):
        #     print(line,"Find!!")
        #     # exit()

        if line.startswith("$Natom"):   # Plrase cut this
            P.Natom = get_value()
        elif line.startswith("$Nbead"):
            P.Nbead = get_value()
        elif line.startswith("$Nstep"):
            P.Nstep = get_value()
        elif line.startswith("$temperature"):
            P.temperature = get_value()
        elif line.startswith("$dt"):
            P.dt = get_value()
        elif line.startswith("$Isimulation"):
            P.Isimulation = get_value()
        # elif line.startswith("$Nref"):
        #     P.Nref = get_value()
        elif line.startswith("$out_step"):
            P.out_step = get_value()
        # elif line.startswith("$Nys"):
        #     P.Nys = get_value()
        elif line.startswith("$Nnhc"):
            P.Nnhc = get_value()
        elif line.startswith("$Ncent"):
            P.Ncent = get_value()
        elif line.startswith("$Lperiodic"):
            P.Lperiodic = get_bool()
        elif line.startswith("$Lrestart"):
            P.Lrestart = get_bool()
        elif line.startswith("$Coords"):
            coords_start = i + 1
        elif line.startswith("$end Coords"):
            coords_end = i
        # elif line.startswith("$end parameter"):
        elif line.startswith("$end Coords"):
            break
        i += 1
    # print("Lperiodic",P.Lperiodic)
    # print("Langstrom",P.Langstrom)
    # exit(1)

    if coords_start >= 0 and coords_end > coords_start:
        P.Natom = coords_end - coords_start
    elif coords_start >= 0 and coords_end <= coords_start:
        raise ValueError("ERROR!!: '$Coords' exists but '$end Coords' not found properly.")

    if "$end parameter" not in ''.join(lines):
        raise ValueError('ERROR!!: There is no "$end parameter"')


def check_input(output_file=None):
    """
    Check and print input parameters to an output file.
    """
    if output_file is None:
        output_file = P.Fout if hasattr(P, "Fout") else "output.log"

    with open(output_file, "a") as f:
        f.write("+++++ Input Check +++++\n")
        f.write(f"{'+++++ Isimulation':26s}{P.Isimulation:12d}\n")
        # f.write(f"{'+++++ Simulation type':26s}{P.name_simulation}\n")
        f.write(f"{'+++++ Number of Atoms':26s}{P.Natom:12d}\n")
        f.write(f"{'+++++ Number of Beads':26s}{P.Nbead:12d}\n")
        f.write(f"{'+++++ Number of Steps':26s}{P.Nstep:12d}\n")
        f.write(f"{'+++++ Given Temperature':26s}{P.temperature:12.4f}\n")
        f.write(f"{'+++++ Given time step':26s}{P.dt / P.fs2AU:12.4f}\n")
        f.write(f"{'+++++ Output step (fs)':26s}{P.dt / P.fs2AU * P.out_step:12.4f}\n")
        f.write(f"{'+++++ Method of Centr NHC':26s}{P.Ncent:12d}\n")
        f.write(f"{'+++++ Length of Centr NHC':26s}{P.Nnhc:12d}\n")
        f.write(f"{'+++++ Flag for Restart':26s}{str(P.Lrestart):>12}\n")
        # f.write(f"{'+++++ Seed for Random No.1':26s}{P.Iseeds:12d}\n")
        f.write(f"{'+++++ Address of Result':26s}{P.dir_result}\n")
        f.write(f"{'+++++ Address of Scratch':26s}{P.dir_scr}\n")
        f.write("\n+++++ Atomic Label, Mass, and Coords +++++\n")

        label = P.alabel # arrays["alabel"]
        mass = P.physmass/P.factmass
        temp = P.alabel

        for i in range(P.Natom):
            coord = P.ur[:,i,0]*P.AUtoAng  # ur(3, Natom, Nbead) → 1st bead
        #     if P.Langstrom:
        #         coord = coord * P.AUtoAng
            f.write(f"{label[i]:4s} {mass[i]:12.6f} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n")

        f.write("  " + "+" * 42 + "\n\n")
