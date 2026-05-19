#!/usr/bin/env python3
import numpy as np
import parameters as P

_CALCULATOR = None


def _get_matlantis_calc_mode(mode_name, estimator_modes):
    mode_key = (mode_name or "MOLECULE").strip().upper()
    mode = getattr(estimator_modes, mode_key, None)
    if mode is None:
        valid_modes = (
            "CRYSTAL",
            "CRYSTAL_U0",
            "CRYSTAL_PLUS_D3",
            "MOLECULE",
        )
        raise ValueError(
            f"Unsupported EstimatorCalcMode: {mode_name!r}. "
            f"Choose from {', '.join(valid_modes)}."
        )
    return mode

def _build_calculator():
    force_module = P.force_module.lower()
    model_path = P.model_path or None
    device = P.device

    if force_module == "emt":
        from ase.calculators.emt import EMT
        return EMT()

    if force_module == "matlantis":
        import pfp_api_client
        from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
        from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

        calc_mode = _get_matlantis_calc_mode(
            getattr(P, "estimator_calc_mode", "MOLECULE"),
            EstimatorCalcMode,
        )
        estimator = Estimator(model_version="v7.0.0", calc_mode=calc_mode)
        return ASECalculator(estimator)

    if force_module == "mattersim":
        # Prefer macer's factory path to match macer PIMD behavior.
        try:
            from macer.calculator.factory import get_calculator
            return get_calculator(ff_name="mattersim", model_path=model_path, device=device)
        except Exception:
            from mattersim.forcefield import MatterSimCalculator
            if model_path:
                return MatterSimCalculator(device=device, load_path=model_path)
            return MatterSimCalculator(device=device)

    raise ValueError(f"Unsupported force field: {force_module}")


def _get_calculator():
    # Build the ASE calculator once and reuse it for all force evaluations.
    global _CALCULATOR
    if _CALCULATOR is None:
        _CALCULATOR = _build_calculator()
    return _CALCULATOR


def run_cal():
    calculator = _get_calculator()
    atoms_bead = P.ase_atoms
    if atoms_bead is None:
        raise ValueError("P.ase_atoms is not prepared. Call prepare_ase_atoms() before run_cal().")

    for i, atoms in enumerate(atoms_bead):
        atoms.calc = calculator

        energy = atoms.get_total_energy()
        forces = atoms.get_forces()

        P.Eenergy[i] = energy
        for j in range(P.Natom):
            P.fr[:, j, i] = forces[j, :]

    P.fr *= P.eVAng2AU * P.dp_inv
    P.Eenergy *= P.eVtoAU
    P.potential = np.sum(P.Eenergy) * P.dp_inv

