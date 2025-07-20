#!/usr/bin/env python3
import numpy as np
import parameters as P
from utility import atom2mass

def read_structure(filename='input.inp'):
    """
    input.inpファイルから原子座標と元素ラベルを読み込み、urとphysmassを更新
    """

    with open(filename, 'r') as f:
        lines = f.readlines()

    istart = next(i for i, line in enumerate(lines) if line.strip() == '$Coords') + 1
    iend = next(i for i, line in enumerate(lines) if line.strip() == '$end Coords')

    coord_lines = lines[istart:iend]

    if len(coord_lines) != P.Natom:
        raise ValueError(f"Natom={P.Natom} but {len(coord_lines)} coordinate lines found.")

    # 一時配列
    ur = np.zeros((3, P.Natom, P.Nbead))
    labels = []
    physmass = np.zeros(P.Natom)

    for i, line in enumerate(coord_lines):
        tokens = line.split()
        if len(tokens) != 4:
            raise ValueError(f"Invalid coordinate line: {line}")
        label, x, y, z = tokens
        labels.append(label)
        for ib in range(P.Nbead):
            ur[:, i, ib] = [float(x), float(y), float(z)]

        physmass[i] = atom2mass(label)

    # グローバル変数に反映
    P.ur[:] = ur * P.AngtoAU
    P.alabel = labels
    P.physmass[:] = physmass * P.factmass


