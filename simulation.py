#!/usr/bin/env python3
import parameters as P
from neural_network import set_nnp_matlantis
from utility import program_abort
from restart import restart_read, restart_write, restart_read_cl, restart_write_cl
from print_start_end import print_ini, print_result_qm, print_ham
from force_nnp_matlantis import force_nnp_matlantis
from mod_md_subroutine import (
    temp_ctr, nmtrans_fr2fur, nmtrans_ur2r, get_force_ref,
    Uupdate, Vupdate, Vupdate_Ref)
from calc_energy import ham_temp
from nhc_integrate import nhc_integrate_cent3, nhc_integrate
from set_ini_md import (
    setup_time_mass, normal_mode, init_mass, nm_position,
    init_velocity, init_bath, init_bath_cl, ham_temp_cl
    )

def simulation_qm():
    setup_time_mass()
    normal_mode()
    init_mass()
    set_nnp_matlantis()

    if P.Lrestart:
        restart_read()
    else:
        print_ini()
        nm_position()
        init_velocity()
        init_bath()
        temp_ctr()
        nmtrans_ur2r()
        force_nnp_matlantis()
        if P.istepsv % P.out_step == 0:
            print_result_qm()
        nmtrans_fr2fur()
        ham_temp()
        print_ham(P.Irestep)

    get_force_ref()

    for istepsv in range(P.Irestep + 1, P.Nstep + 1):
        P.istepsv = istepsv

        if P.Ncent == 3:
            nhc_integrate_cent3()
        Vupdate()

        for iref in range(1, P.Nref + 1):
            if P.Ncent == 0:
                Vupdate_Ref()
                Uupdate()
                get_force_ref()
                Vupdate_Ref()
            else:
                nhc_integrate()
                Vupdate_Ref()
                Uupdate()
                get_force_ref()
                Vupdate_Ref()
                nhc_integrate()

        nmtrans_ur2r()
        force_nnp_matlantis()
        if istepsv % P.out_step == 0:
            print_result_qm()
        nmtrans_fr2fur()
        Vupdate()

        if P.Ncent == 3:
            nhc_integrate_cent3()

        ham_temp()
        print_ham(istepsv)

        if istepsv % P.out_step == 0:
            restart_write(istepsv)

    #     if istepsv % 10 == 0:
    #         exit_program()

    with open(P.Fout, "a") as fout:
        fout.write(" " + "*" * 121 + "\n")



def simulation_cl():
    setup_time_mass()
    init_mass()
    set_nnp_matlantis()

    # for i in range(P.Natom):
    #     print(P.r[:,i,0])
    # exit(1)

    if P.Lrestart:
        restart_read_cl()
    else:
        print_ini()
        init_velocity()
        init_bath_cl()
        P.r[:,:,0] = P.ur[:,:,0]
        force_nnp_matlantis()
        P.fur[:,:,0] = P.fr[:,:,0]
        if P.istepsv % P.out_step == 0:
            print_result_qm()
        ham_temp_cl()
        print_ham(P.Irestep)

    for istepsv in range(P.Irestep + 1, P.Nstep + 1):
        P.istepsv = istepsv

        if P.Ncent == 3:
            nhc_integrate_cent3()

        Vupdate()
        Uupdate()
        P.r[:,:,0] = P.ur[:,:,0]
        force_nnp_matlantis()
        P.fur[:,:,0] = P.fr[:,:,0]

        if istepsv % P.out_step == 0:
            print_result_qm()

        Vupdate()

        if P.Ncent == 3:
            nhc_integrate_cent3()

        ham_temp_cl()
        print_ham(istepsv)

        if istepsv % P.out_step == 0:
            restart_write_cl(istepsv)

    with open(P.Fout, "a") as fout:
        fout.write(" " + "*" * 95 + "\n")

