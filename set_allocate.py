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
    # P.lattice    = np.zeros((3,3))

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


