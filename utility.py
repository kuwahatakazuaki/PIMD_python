#!/usr/bin/env python3

def atom2mass(symbol: str) -> float:
    """
    Convert atomic symbol to atomic mass (in amu).
    
    Parameters
    ----------
    symbol : str
        Atomic symbol (e.g. "H", "He", "O", "Fe", ...)

    Returns
    -------
    mass : float
        Atomic mass in atomic mass units (amu)
    """
    symbol = symbol.strip().lower()

    mass_table = {
        "h": 1.00782503223, "d": 2.01410177812, "mu": 0.1134289257,
        "he": 4.00260325413, "li": 7.0160034366, "be": 9.012183065,
        "b": 6.9675, "c": 12.0, "n": 14.00307400443, "o": 15.99491461957,
        "f": 18.99840316273, "ne": 19.9924401762, "na": 22.9897692820,
        "mg": 24.3055, "al": 26.9815386, "si": 28.0855, "p": 30.97376199842,
        "s": 32.065, "cl": 35.453, "ar": 39.948, "k": 39.0983,
        "ca": 40.078, "sc": 44.955912, "ti": 47.867, "v": 50.9415,
        "cr": 51.9961, "mn": 54.938045, "fe": 55.845, "co": 58.933195,
        "ni": 58.6934, "cu": 63.546, "zn": 65.38, "ga": 69.723,
        "ge": 72.64, "as": 74.92160, "se": 78.96, "br": 79.904,
        "kr": 83.798, "rb": 85.4678, "sr": 87.62, "y": 88.90585,
        "zr": 91.224, "nb": 92.90638, "mo": 95.96, "tc": 98.0,
        "ru": 101.07, "rh": 102.90550, "pd": 106.42, "ag": 107.8682,
        "cd": 112.411, "in": 114.818, "sn": 118.710, "sb": 121.760,
        "te": 127.60, "i": 126.90447, "xe": 131.294, "cs": 132.9054519,
        "ba": 137.327, "la": 138.90547, "ce": 140.116, "pt": 195.084,
        "au": 196.966569, "hg": 200.59
    }

    if symbol in mass_table:
        return mass_table[symbol]
    else:
        raise ValueError(f'ERROR!! "{symbol}" is not defined in atom2mass')

# utility.py
def program_abort(message: str):
    """
    標準出力にメッセージを表示して異常終了する。
    """
    print(message)
    raise SystemExit(1)  # または raise RuntimeError(message)

