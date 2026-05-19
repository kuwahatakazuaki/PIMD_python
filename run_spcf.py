#!/usr/bin/env python3
import numpy as np
import parameters as P


def _validate_spcf_system():
    if P.Lperiodic:
        raise ValueError('SPC/F is only supported for non-periodic systems.')

    if P.Natom % 3 != 0:
        raise ValueError('SPC/F requires Natom to be a multiple of 3 (O, H, H ordering).')

    for iatom in range(0, P.Natom, 3):
        labels = [str(P.alabel[iatom + offset]).strip().upper() for offset in range(3)]
        if labels != ["O", "H", "H"]:
            raise ValueError(
                'SPC/F requires atoms to be ordered as O, H, H, O, H, H, ...'
            )


def run_cal():
    _validate_spcf_system()

    au_length = 0.529177249e-10
    au_energy = 4.3597482e-18
    boltz = 0.316682968e-5

    rho_w = 2.361e10 * au_length
    d_w = 0.708e-18 / au_energy

    b_oh = 1.000e-10 / au_length
    b_hh = b_oh * np.sin(P.pi * 108.0 / 2.0 / 180.0) * 2.0

    b_const = 1.803e02 / au_energy * au_length**2
    c_const = -1.469e02 / au_energy * au_length**2
    d_const = 0.776e02 / au_energy * au_length**2

    sigma = 3.165e-10 / au_length
    epsilon = 78.22 * boltz

    q_o = -0.82
    q_h = 0.41
    charges = np.array([q_o, q_h, q_h], dtype=float)

    es6 = epsilon * sigma**6
    es12 = epsilon * sigma**12

    P.Eenergy[:] = 0.0
    P.fr[:, :, :] = 0.0

    for ibead in range(P.Nbead):
        for iatom in range(0, P.Natom, 3):
            io = iatom
            ih1 = iatom + 1
            ih2 = iatom + 2

            rr_oh1 = P.r[:, io, ibead] - P.r[:, ih1, ibead]
            rr_oh2 = P.r[:, io, ibead] - P.r[:, ih2, ibead]
            rr_hh = P.r[:, ih1, ibead] - P.r[:, ih2, ibead]

            r_oh1 = np.linalg.norm(rr_oh1)
            r_oh2 = np.linalg.norm(rr_oh2)
            r_hh = np.linalg.norm(rr_hh)

            P.Eenergy[ibead] += (
                rho_w * rho_w * d_w * (r_oh1 - b_oh) * (r_oh1 - b_oh)
                + rho_w * rho_w * d_w * (r_oh2 - b_oh) * (r_oh2 - b_oh)
                + 0.5 * b_const * (r_hh - b_hh) * (r_hh - b_hh)
                + c_const * (r_oh1 + r_oh2 - 2.0 * b_oh) * (r_hh - b_hh)
                + d_const * (r_oh1 - b_oh) * (r_oh2 - b_oh)
            )

            P.fr[:, io, ibead] += (
                -2.0 * rho_w * rho_w * d_w * rr_oh1 / r_oh1 * (r_oh1 - b_oh)
                -2.0 * rho_w * rho_w * d_w * rr_oh2 / r_oh2 * (r_oh2 - b_oh)
                -c_const * rr_oh1 / r_oh1 * (r_hh - b_hh)
                -c_const * rr_oh2 / r_oh2 * (r_hh - b_hh)
                -d_const * rr_oh1 / r_oh1 * (r_oh2 - b_oh)
                -d_const * rr_oh2 / r_oh2 * (r_oh1 - b_oh)
            )

            P.fr[:, ih1, ibead] += (
                2.0 * rho_w * rho_w * d_w * rr_oh1 / r_oh1 * (r_oh1 - b_oh)
                -b_const * rr_hh / r_hh * (r_hh - b_hh)
                + c_const * rr_oh1 / r_oh1 * (r_hh - b_hh)
                - c_const * (r_oh1 + r_oh2 - 2.0 * b_oh) * rr_hh / r_hh
                + d_const * rr_oh1 / r_oh1 * (r_oh2 - b_oh)
            )

            P.fr[:, ih2, ibead] += (
                2.0 * rho_w * rho_w * d_w * rr_oh2 / r_oh2 * (r_oh2 - b_oh)
                + b_const * rr_hh / r_hh * (r_hh - b_hh)
                + c_const * rr_oh2 / r_oh2 * (r_hh - b_hh)
                + c_const * (r_oh1 + r_oh2 - 2.0 * b_oh) * rr_hh / r_hh
                + d_const * rr_oh2 / r_oh2 * (r_oh1 - b_oh)
            )

        for iatom in range(0, P.Natom, 3):
            for jatom in range(iatom + 3, P.Natom, 3):
                rr = P.r[:, iatom, ibead] - P.r[:, jatom, ibead]
                r1 = np.linalg.norm(rr)
                rinv = 1.0 / r1
                rinv2 = rinv * rinv
                rinv6 = rinv2 * rinv2 * rinv2
                rinv12 = rinv6 * rinv6

                P.Eenergy[ibead] += -4.0 * es6 * rinv6 + 4.0 * es12 * rinv12
                dvdr = 24.0 * es6 * rinv6 * rinv - 48.0 * es12 * rinv12 * rinv

                P.fr[:, iatom, ibead] -= dvdr * rr * rinv
                P.fr[:, jatom, ibead] += dvdr * rr * rinv

        for iatom in range(P.Natom):
            qi = charges[iatom % 3]
            mol_i = iatom // 3

            for jatom in range((mol_i + 1) * 3, P.Natom):
                qj = charges[jatom % 3]
                rr = P.r[:, iatom, ibead] - P.r[:, jatom, ibead]
                r1 = np.linalg.norm(rr)
                rinv = 1.0 / r1
                rinv2 = rinv * rinv

                P.Eenergy[ibead] += qi * qj * rinv
                dvdr = -qi * qj * rinv2

                P.fr[:, iatom, ibead] -= dvdr * rr * rinv
                P.fr[:, jatom, ibead] += dvdr * rr * rinv

    P.fr *= P.dp_inv
    P.potential = np.sum(P.Eenergy) * P.dp_inv
