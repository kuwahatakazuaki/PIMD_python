#!/usr/bin/env python3
import os
import parameters as P
from utility import program_abort


def err_lattice():
    lattice_template = [
        "H2O molecule",
        "1.0",
        "        7.9376997948         0.0000000000         0.0000000000",
        "        0.0000000000         7.9376997948         0.0000000000",
        "        0.0000000000         0.0000000000         7.9376997948",
        "    O    H",
        "    1    2",
        "Cartesian"
    ]

    with open("LATTICE", "w") as f:
        for line in lattice_template:
            f.write(line.strip() + "\n")

    print('ERROR!! "LATTICE" is NOT exist!')
    print('Please use "LATTICE" template')
    program_abort('')

def set_nnp_matlantis():
    # name_file = 'run_matlantis.py'
    # if not os.path.exists(name_file):
    #     program_abort(f'ERROR!!! There is no "{name_file}"')

    # 周期系であれば LATTICE の存在チェック
    if P.Lperiodic:
        if not os.path.exists("LATTICE"):
            err_lattice()

