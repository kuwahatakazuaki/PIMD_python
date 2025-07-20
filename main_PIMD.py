#!/usr/bin/env python3
import parameters as P
from print_start_end import print_start, print_end
from read_parameter import read_parameter, check_input
from read_structure import read_structure
from calc_constant import calc_constant
import set_allocate
from simulation import simulation_qm, simulation_cl

if __name__ == "__main__":
    read_parameter("input.inp")
    print_start("std.out")
    set_allocate.allocate_arrays()
    read_structure("input.inp")
    
    calc_constant()
    check_input(output_file="std.out")

    # シミュレーション分岐
    if 0 <= P.Isimulation <= 2:
        simulation_qm()
    elif P.Isimulation == 10:
        simulation_cl()
    else:
        print(f'"Simulation" is {P.Isimulation}')
        raise ValueError('ERROR!!! Wrong "Simulation" option')

    # set_deallocate()
    print_end("std.out")

