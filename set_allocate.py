#!/usr/bin/env python3
import numpy as np
import parameters as P

def allocate_arrays():
    """
    各変数を NumPy 配列としてグローバルに初期化する（Fortranの allocate に相当）
    """

    P.atom_num   = np.zeros(P.Natom, dtype=int)
    P.r          = np.zeros((3, P.Natom, P.Nbead))
    P.ur         = np.zeros((3, P.Natom, P.Nbead))
    P.vur        = np.zeros((3, P.Natom, P.Nbead))
    P.fr         = np.zeros((3, P.Natom, P.Nbead))
    P.fur        = np.zeros((3, P.Natom, P.Nbead))
    P.physmass   = np.zeros(P.Natom)
    P.dnmmass    = np.zeros((P.Natom, P.Nbead))
    P.fictmass   = np.zeros((P.Natom, P.Nbead))

    # ラベル・エネルギー
    P.alabel     = np.empty(P.Natom, dtype=object)
    P.Eenergy    = np.zeros(P.Nbead)
    P.ysweight   = np.zeros(P.Nys)

    # パス積分用の追加配列（古典でないとき）
    if P.Isimulation != 10:
        P.qmass     = np.zeros(P.Nbead)
        P.fur_ref   = np.zeros((3, P.Natom, P.Nbead))
        P.tnm       = np.zeros((P.Nbead, P.Nbead))
        P.tnminv    = np.zeros((P.Nbead, P.Nbead))
        P.u         = np.zeros((P.Nbead, P.Nbead))
        P.uinv      = np.zeros((P.Nbead, P.Nbead))
        P.rbath     = np.zeros((3, P.Natom, P.Nnhc, P.Nbead))
        P.vrbath    = np.zeros((3, P.Natom, P.Nnhc, P.Nbead))
        P.frbath    = np.zeros((3, P.Natom, P.Nnhc, P.Nbead))

    # 中心点の選択肢による配列確保
    if P.Ncent == 3:
        P.rbc31     = np.zeros((3, P.Natom, P.Nnhc))
        P.vrbc31    = np.zeros((3, P.Natom, P.Nnhc))
        P.frbc31    = np.zeros((3, P.Natom, P.Nnhc))
        P.qmcent31  = np.zeros(P.Nnhc)


# import numpy as np
# from typing import Dict, Any

# def allocate_arrays(Natom: int, Nbead: int, Nnhc: int, Nys: int, Ncent: int, Isimulation: int) -> Dict[str, Any]:
#     """
#     Allocate and return all required simulation arrays as a dictionary.
#     """

#     arrays = {
#         "atom_num": np.zeros(Natom, dtype=int),
#         "r":        np.zeros((3, Natom, Nbead)),
#         "ur":       np.zeros((3, Natom, Nbead)),
#         "vur":      np.zeros((3, Natom, Nbead)),
#         "fr":       np.zeros((3, Natom, Nbead)),
#         "fur":      np.zeros((3, Natom, Nbead)),
#         "physmass": np.zeros(Natom),
#         "dnmmass":  np.zeros((Natom, Nbead)),
#         "fictmass": np.zeros((Natom, Nbead)),
#         "alabel":   np.empty(Natom, dtype=object),
#         "Eenergy":  np.zeros(Nbead),
#         "ysweight": np.zeros(Nys),
#     }

#     if Isimulation != 10:
#         arrays.update({
#             "qmass":    np.zeros(Nbead),
#             "fur_ref":  np.zeros((3, Natom, Nbead)),
#             "tnm":      np.zeros((Nbead, Nbead)),
#             "tnminv":   np.zeros((Nbead, Nbead)),
#             "u":        np.zeros((Nbead, Nbead)),
#             "uinv":     np.zeros((Nbead, Nbead)),
#             "rbath":    np.zeros((3, Natom, Nnhc, Nbead)),
#             "vrbath":   np.zeros((3, Natom, Nnhc, Nbead)),
#             "frbath":   np.zeros((3, Natom, Nnhc, Nbead)),
#         })

#     if Ncent == 3:
#         arrays.update({
#             "rbc31":    np.zeros((3, Natom, Nnhc)),
#             "vrbc31":   np.zeros((3, Natom, Nnhc)),
#             "frbc31":   np.zeros((3, Natom, Nnhc)),
#             "qmcent31": np.zeros(Nnhc),
#         })

#     return arrays

