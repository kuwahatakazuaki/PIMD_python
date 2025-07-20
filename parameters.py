#!/usr/bin/env python3
import numpy as np
# グローバルに共有される物理定数・計算条件・配列の定義モジュール

# ----------------------------
# 物理定数・単位変換係数など
# ----------------------------
# pi       = 3.141592653589793
# fs2AU    = 41.3413745758      # fs → 原子単位
# factmass = 1822.888486209     # 質量のamu → 電子質量
# eVtoAU   = 0.0367493221757    # eV → 原子単位
# KtoAU    = 3.1668114e-6       # K → 原子単位エネルギー

pi = np.pi
fs2AU = 1.0 / 0.024188843
factmass = 1.6605402e-27 / 9.1093897e-31
eVtoAU = 1.0 / 27.21162
AngtoAU = 1.0 / 0.529177249
AUtoAng = 0.529177249
AUtoJ = 4.35974434e-18
KtoAU = 8.617333262145e-5 / 27.211396132
eVAng2AU = (1.0 / 27.21162) * 0.529177249


# ----------------------------
# パラメータ（実行前に read_parameter 等で代入）
# ----------------------------
Natom       = -1   # 原子数
Nbead       = -1   # ビーズ数（パス積分の数）
Nnhc        = -1   # NHC熱浴の自由度
Nys         = 5   # YS法の重み数
Ncent       = -1   # 中心点設定（例: 3で中心重心など）
Nstep       = -1   # MDステップ数
temperature = -1.0 # 温度（K）
dt          = -1.0 # タイムステップ（fs）
Isimulation = -1   # シミュレーションタイプ（古典 or PIMD）
Nref        = 5
gamma1      = 1.0
gamma2      = 0.0
istepsv     = 0
Irestep     = 0


# ----------------------------
# 配列（allocate_arrays で代入される）
# ----------------------------
atom_num   = None
r          = None
ur         = None
vur        = None
fr         = None
fur        = None
physmass   = None
dnmmass    = None
fictmass   = None
alabel     = None
Eenergy    = None
ysweight   = None

# 経路積分用
qmass      = None
fur_ref    = None
tnm        = None
tnminv     = None
u          = None
uinv       = None
rbath      = None
vrbath     = None
frbath     = None

# 中心点モード用
rbc31      = None
vrbc31     = None
frbc31     = None
qmcent31   = None

# ----------------------------
# 定数計算（calc_constant で代入）
# ----------------------------
beta       = None
omega_p2   = None
omega_sys  = None
gkt        = None
gnkt       = None

# --- File names and paths ---
address0   = ""
addresstmp = ""
Finp = "input.inp"
Fout = "std.out"
Ferr = "std.err"
# name_simulation: str = ""
dir_result = "./Result"
dir_scr = "./Scr"


Lsave_force = False
Langstrom = True
# Lperiodic = False
Lperiodic = True
Lrestart = False
Lrandom_coor = False



# import numpy as np
# from dataclasses import dataclass, field
# from typing import ClassVar, Optional

# @dataclass
# class Parameters:
#     # --- Physical constants ---
#     pi: ClassVar[float] = np.pi
#     fs2AU: ClassVar[float] = 1.0 / 0.024188843
#     factmass: ClassVar[float] = 1.6605402e-27 / 9.1093897e-31
#     eVtoAU: ClassVar[float] = 1.0 / 27.21162
#     AngtoAU: ClassVar[float] = 1.0 / 0.529177249
#     AUtoAng: ClassVar[float] = 0.529177249
#     AUtoJ: ClassVar[float] = 4.35974434e-18
#     KtoAU: ClassVar[float] = 8.617333262145e-5 / 27.211396132
#     eVAng2AU: ClassVar[float] = (1.0 / 27.21162) * 0.529177249

#     # --- Integer parameters ---
#     Natom: int = 0
#     Nbead: int = 0
#     Nstep: int = 0
#     Isimulation: int = 0
#     Nref: int = 0
#     Nys: int = 0
#     Nnhc: int = 0
#     out_step: int = 1
#     Ncent: int = 0
#     Irestep: int = 0
#     Iseeds: int = 0
#     laddress: int = 0
#     istepsv: int = 0

#     # --- Real scalar variables ---
#     gamma1: float = 1.0
#     gamma2: float = 0.0
#     omega_system: float = 0.0
#     omega_p2: float = 0.0
#     omega2: float = 0.0
#     gkt: float = 0.0
#     gnkt: float = 0.0
#     dp_inv: float = 0.0
#     E_Virial: float = 0.0
#     ebath: float = 0.0
#     ebath_cent: float = 0.0
#     dkinetic: float = 0.0
#     qkinetic: float = 0.0
#     beta: float = 0.0
#     temperature: float = 0.0
#     dt: float = 0.0
#     dt_ref: float = 0.0
#     potential: float = 0.0
#     hamiltonian: float = 0.0
#     temp: float = 0.0

#     # --- File names and paths ---
#     Finp: str = "input.inp"
#     Fout: str = "std.out"
#     Ferr: str = "std.err"
#     name_simulation: str = ""
#     dir_result: str = ""
#     dir_scr: str = ""
#     address0: str = ""
#     addresstmp: str = ""

#     # --- Logical flags ---
#     Lsave_force: bool = False
#     Langstrom: bool = True
#     Lperiodic: bool = False
#     Lrestart: bool = False
#     Lrandom_coor: bool = False

#     # --- Arrays (to be initialized later) ---
#     r: Optional[np.ndarray] = field(default=None)
#     fr: Optional[np.ndarray] = field(default=None)
#     ur: Optional[np.ndarray] = field(default=None)
#     vur: Optional[np.ndarray] = field(default=None)
#     fur: Optional[np.ndarray] = field(default=None)
#     fur_ref: Optional[np.ndarray] = field(default=None)
#     rbath: Optional[np.ndarray] = field(default=None)
#     vrbath: Optional[np.ndarray] = field(default=None)
#     frbath: Optional[np.ndarray] = field(default=None)
#     rbc31: Optional[np.ndarray] = field(default=None)
#     vrbc31: Optional[np.ndarray] = field(default=None)
#     frbc31: Optional[np.ndarray] = field(default=None)
#     tnm: Optional[np.ndarray] = field(default=None)
#     tnminv: Optional[np.ndarray] = field(default=None)
#     u: Optional[np.ndarray] = field(default=None)
#     uinv: Optional[np.ndarray] = field(default=None)
#     pot: Optional[np.ndarray] = field(default=None)
#     physmass: Optional[np.ndarray] = field(default=None)
#     dnmmass: Optional[np.ndarray] = field(default=None)
#     fictmass: Optional[np.ndarray] = field(default=None)
#     qmass: Optional[np.ndarray] = field(default=None)
#     ysweight: Optional[np.ndarray] = field(default=None)
#     qmcent31: Optional[np.ndarray] = field(default=None)
#     atom_num: Optional[np.ndarray] = field(default=None)
#     Eenergy: Optional[np.ndarray] = field(default=None)
#     alabel: Optional[np.ndarray] = field(default=None)

