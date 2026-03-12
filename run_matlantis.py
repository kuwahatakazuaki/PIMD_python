#!/usr/bin/env python3
import numpy as np
import parameters as P

_CALCULATOR = None
_CALCULATOR_KEY = None

# # === For Matlantis ===
# import pfp_api_client
# from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
# from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

# # === EstimatorCalcMode : CRYSTAL, CRYSTAL_U0, CRYSTAL_PLUS_D3, MOLECULE ===
# estimator = Estimator(calc_mode=EstimatorCalcMode.MOLECULE)
# calculator = ASECalculator(estimator)
# === For Matlantis ===

# # # === For Effective Medium Theory ===
# from ase.calculators.emt import EMT
# calculator = EMT()
# # # === For Effective Medium Theory ===


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
    ff_name = getattr(P, "force_module", "emt").lower()
    model_path = getattr(P, "model_path", "") or None
    device = getattr(P, "device", "cpu")

    if ff_name == "emt":
        from ase.calculators.emt import EMT
        return EMT()

    if ff_name == "matlantis":
        import pfp_api_client
        from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
        from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode

        calc_mode = _get_matlantis_calc_mode(
            getattr(P, "estimator_calc_mode", "MOLECULE"),
            EstimatorCalcMode,
        )
        estimator = Estimator(model_version="v7.0.0", calc_mode=calc_mode)
        return ASECalculator(estimator)

    if ff_name == "mattersim":
        # Prefer macer's factory path to match macer PIMD behavior.
        try:
            from macer.calculator.factory import get_calculator
            return get_calculator(ff_name="mattersim", model_path=model_path, device=device)
        except Exception:
            from mattersim.forcefield import MatterSimCalculator
            if model_path:
                return MatterSimCalculator(device=device, load_path=model_path)
            return MatterSimCalculator(device=device)

    raise ValueError(f"Unsupported force field: {ff_name}")


def _get_calculator():
    global _CALCULATOR, _CALCULATOR_KEY
    key = (
        getattr(P, "force_module", "emt"),
        getattr(P, "model_path", "") or "",
        getattr(P, "device", "cpu"),
    )
    if _CALCULATOR is None or _CALCULATOR_KEY != key:
        _CALCULATOR = _build_calculator()
        _CALCULATOR_KEY = key
    return _CALCULATOR


def run_cal():
    calculator = _get_calculator()
    atoms_list = getattr(P, "ase_atoms", None)
    if atoms_list is None:
        raise ValueError("P.ase_atoms is not prepared. Call prepare_ase_atoms() before run_cal().")

    for i, atoms in enumerate(atoms_list):
        atoms.calc = calculator

        energy = atoms.get_total_energy()
        forces = atoms.get_forces()

        P.Eenergy[i] = energy
        for j in range(P.Natom):
            P.fr[:, j, i] = forces[j, :]

    P.fr *= P.eVAng2AU * P.dp_inv
    P.Eenergy *= P.eVtoAU
    P.potential = np.sum(P.Eenergy) * P.dp_inv
