#!/usr/bin/env python3
import numpy as np
import parameters as P
import os
from pathlib import Path
import shutil

def restart_read():
    filepath = Path(P.dir_result) / "restart.dat"

    with open(filepath, "r") as f:
        # Irestep を読み込み
        P.Irestep = int(f.readline())

        for j in range(P.Nbead):
            for i in range(P.Natom):
                P.ur[:, i, j] = np.fromstring(f.readline(), sep=" ")

        for j in range(P.Nbead):
            for i in range(P.Natom):
                pur = np.fromstring(f.readline(), sep=" ")
                P.vur[:, i, j] = pur / np.sqrt(P.fictmass[i, j])

        for j in range(P.Nbead):
            for i in range(P.Natom):
                P.fur[:, i, j] = np.fromstring(f.readline(), sep=" ")

        for j in range(P.Nbead):
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    P.rbath[:, i, inhc, j] = np.fromstring(f.readline(), sep=" ")

        for j in range(P.Nbead):
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    P.vrbath[:, i, inhc, j] = np.fromstring(f.readline(), sep=" ")

        for j in range(P.Nbead):
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    P.frbath[:, i, inhc, j] = np.fromstring(f.readline(), sep=" ")

        if P.Ncent == 3:
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    P.rbc31[:, i, inhc] = np.fromstring(f.readline(), sep=" ")
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    P.vrbc31[:, i, inhc] = np.fromstring(f.readline(), sep=" ")
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    P.frbc31[:, i, inhc] = np.fromstring(f.readline(), sep=" ")



def restart_write(istep):
    restart_path = os.path.join(P.dir_result, 'restart.dat')
    restart1_path = os.path.join(P.dir_result, 'restart1.dat')

    if istep > 1 and os.path.exists(restart_path):
        shutil.copyfile(restart_path, restart1_path)

    with open(restart_path, 'w') as f:
        f.write(f"{istep:10d}\n")

        # ur
        for j in range(P.Nbead):
            for i in range(P.Natom):
                f.write(" ".join(f"{x:.15e}" for x in P.ur[:, i, j]) + "\n")

        # vur * sqrt(mass)
        for j in range(P.Nbead):
            for i in range(P.Natom):
                adjusted_vur = P.vur[:, i, j] * np.sqrt(P.fictmass[i, j])
                f.write(" ".join(f"{x:.15e}" for x in adjusted_vur) + "\n")

        # fur
        for j in range(P.Nbead):
            for i in range(P.Natom):
                f.write(" ".join(f"{x:.15e}" for x in P.fur[:, i, j]) + "\n")

        # rbath
        for j in range(P.Nbead):
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    f.write(" ".join(f"{x:.15e}" for x in P.rbath[:, i, inhc, j]) + "\n")

        # vrbath
        for j in range(P.Nbead):
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    f.write(" ".join(f"{x:.15e}" for x in P.vrbath[:, i, inhc, j]) + "\n")

        # frbath
        for j in range(P.Nbead):
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    f.write(" ".join(f"{x:.15e}" for x in P.frbath[:, i, inhc, j]) + "\n")

        # centroid thermostat if Ncent == 3
        if P.Ncent == 3:
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    f.write(" ".join(f"{x:.15e}" for x in P.rbc31[:, i, inhc]) + "\n")
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    f.write(" ".join(f"{x:.15e}" for x in P.vrbc31[:, i, inhc]) + "\n")
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    f.write(" ".join(f"{x:.15e}" for x in P.frbc31[:, i, inhc]) + "\n")



def restart_read_cl():
    path = f"{P.dir_result}/restart.dat"
    with open(path, 'r') as f:
        P.Irestep = int(f.readline())

        for i in range(P.Natom):
            P.ur[:, i, 0] = np.fromstring(f.readline(), sep=' ')

        for i in range(P.Natom):
            P.vur[:, i, 0] = np.fromstring(f.readline(), sep=' ')

        for i in range(P.Natom):
            P.fur[:, i, 0] = np.fromstring(f.readline(), sep=' ')

        if P.Ncent == 3:
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    P.rbc31[:, i, inhc] = np.fromstring(f.readline(), sep=' ')
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    P.vrbc31[:, i, inhc] = np.fromstring(f.readline(), sep=' ')
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    P.frbc31[:, i, inhc] = np.fromstring(f.readline(), sep=' ')



def restart_write_cl(istep):
    if istep > 1:
        path_old = f"{P.dir_result}/restart.dat"
        path_new = f"{P.dir_result}/restart1.dat"
        shutil.copyfile(path_old, path_new)

    with open(f"{P.dir_result}/restart.dat", 'w') as f:
        f.write(f"{istep:10d}\n")

        for i in range(P.Natom):
            f.write(" ".join(f"{x:.16e}" for x in P.ur[:, i, 0]) + "\n")

        for i in range(P.Natom):
            f.write(" ".join(f"{x:.16e}" for x in P.vur[:, i, 0]) + "\n")

        for i in range(P.Natom):
            f.write(" ".join(f"{x:.16e}" for x in P.fur[:, i, 0]) + "\n")

        if P.Ncent == 3:
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    f.write(" ".join(f"{x:.16e}" for x in P.rbc31[:, i, inhc]) + "\n")
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    f.write(" ".join(f"{x:.16e}" for x in P.vrbc31[:, i, inhc]) + "\n")
            for inhc in range(P.Nnhc):
                for i in range(P.Natom):
                    f.write(" ".join(f"{x:.16e}" for x in P.frbc31[:, i, inhc]) + "\n")
