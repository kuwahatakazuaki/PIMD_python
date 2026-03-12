#!/usr/bin/env python3
import numpy as np
from ase import Atoms
import parameters as P
import run_matlantis


def prepare_ase_atoms():
    """Build one ASE Atoms object for each bead from the in-memory coordinates."""
    atoms_list = np.empty(P.Nbead, dtype=object)
    # LATTICE is read directly from input.inp in Angstrom, so do not convert it again.
    cell = P.lattice.T if P.Lperiodic else None
    pbc = [True, True, True] if P.Lperiodic else [False, False, False]

    for imode in range(P.Nbead):
        positions = P.r[:, :, imode].T * P.AUtoAng
        atoms_list[imode] = Atoms(
            symbols=P.alabel,
            positions=positions,
            cell=cell,
            pbc=pbc,
        )

    P.ase_atoms = atoms_list



def force_nnp_matlantis():
    prepare_ase_atoms()
    run_matlantis.run_cal()
