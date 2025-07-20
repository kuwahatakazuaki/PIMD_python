#!/usr/bin/env python3
# from ase import Atoms
from ase.io import read
import numpy as np
import sys
import parameters as P

# # === For Matlantis ===
# import pfp_api_client
# from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
# from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

# # === EstimatorCalcMode : CRYSTAL, CRYSTAL_U0, CRYSTAL_PLUS_D3, MOLECULE ===
# estimator = Estimator(calc_mode=EstimatorCalcMode.MOLECULE)
# calculator = ASECalculator(estimator)
# === For Matlantis ===

# # === For Effective Medium Theory ===
from ase.calculators.emt import EMT
calculator = EMT()
# # === For Effective Medium Theory ===

def run_cal():
    Fforce = 'forces.out'
    Fenergy = 'energy.out'

    # (path_dir, Ista, Iend, Lperi) = sys.argv[1:]
    # Ista = int(Ista)
    # Iend = int(Iend)

    if P.Lperiodic == True:
        extension = "vasp"
    elif P.Lperiodic == False:
        extension = "xyz"
    else:
        print("ERROR!!! Bad statement for periodic condition in run_matlantis.py")
        sys.exit()

    all_forces = np.array([])
    all_energy = np.array([])

    for i in range(P.Nbead):
        path_inp=''
        Fname = 'str{0:05}.'.format(i+1) + extension
        path_inp = P.addresstmp + Fname

        atoms = read(path_inp)
        atoms.calc = calculator

        energy = atoms.get_total_energy()
        forces = atoms.get_forces()
        # charges = atoms.get_charges()

        if i == 0:
            all_forces = forces
        else:
            all_forces = np.concatenate([all_forces,forces])
        all_energy = np.append(all_energy,energy)


    np.savetxt(Fforce,all_forces)
    np.savetxt(Fenergy,all_energy)


