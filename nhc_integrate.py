#!/usr/bin/env python3
import numpy as np
import parameters as P

def nhc_integrate_cent3():
    """
    Python implementation of the Fortran subroutine nhc_integrate_cent3.
    Integrates Nosé-Hoover chain thermostats attached to the centroid (mode 1).
    """

    for iys in range(P.Nys):
        dt_ys = P.dt * P.ysweight[iys]

        for iatom in range(P.Natom):
            # Kinetic energy (element-wise squared)
            dkinr = P.fictmass[iatom, 0] * P.vur[:, iatom, 0]**2

            # Update the force for the first thermostat
            P.frbc31[:, iatom, 0] = (dkinr - P.gkt) / P.qmcent31[0]

            # Higher chain levels
            for inhc in range(1, P.Nnhc):
                prev_v = P.vrbc31[:, iatom, inhc - 1]
                P.frbc31[:, iatom, inhc] = (
                    P.qmcent31[inhc - 1] * prev_v**2 - P.gkt
                ) / P.qmcent31[inhc]

            # Update thermostat velocities (first half)
            P.vrbc31[:, iatom, P.Nnhc - 1] += 0.25 * P.frbc31[:, iatom, P.Nnhc - 1] * dt_ys

            for inhc in range(P.Nnhc - 2, -1, -1):
                v = P.vrbc31[:, iatom, inhc + 1]
                vrfact = np.exp(-0.125 * v * dt_ys)
                P.vrbc31[:, iatom, inhc] = (
                    P.vrbc31[:, iatom, inhc] * vrfact**2
                    + 0.25 * P.frbc31[:, iatom, inhc] * vrfact * dt_ys
                )

            # Update particle velocity scaling
            pvrfact = np.exp(-0.5 * P.vrbc31[:, iatom, 0] * dt_ys)

            # Update force for first thermostat again
            P.frbc31[:, iatom, 0] = (pvrfact**2 * dkinr - P.gkt) / P.qmcent31[0]

            # Update thermostat positions
            for inhc in range(P.Nnhc):
                P.rbc31[:, iatom, inhc] += 0.5 * P.vrbc31[:, iatom, inhc] * dt_ys

            # Update thermostat velocities (second half)
            for inhc in range(P.Nnhc - 1):
                v = P.vrbc31[:, iatom, inhc + 1]
                vrfact = np.exp(-0.125 * v * dt_ys)
                P.vrbc31[:, iatom, inhc] = (
                    P.vrbc31[:, iatom, inhc] * vrfact**2
                    + 0.25 * P.frbc31[:, iatom, inhc] * vrfact * dt_ys
                )

                P.frbc31[:, iatom, inhc + 1] = (
                    P.qmcent31[inhc] * P.vrbc31[:, iatom, inhc]**2 - P.gkt
                ) / P.qmcent31[inhc + 1]

            # Final velocity update for last chain level
            P.vrbc31[:, iatom, P.Nnhc - 1] += 0.25 * P.frbc31[:, iatom, P.Nnhc - 1] * dt_ys

            # Final particle velocity update
            P.vur[:, iatom, 0] *= pvrfact


def nhc_integrate():
    """
    Velocity-Verlet形式のNose-Hooverチェーン（非中心モード）を時間発展。
    """
    for iys in range(P.Nys):
        dt_ys = P.dt_ref * P.ysweight[iys]

        for imode in range(1, P.Nbead):  # Fortranの2:Nbead → Pythonの1:Nbead-1
            for iatom in range(P.Natom):

                dkinr = P.fictmass[iatom, imode] * P.vur[:, iatom, imode]**2

                # update force on first NHC
                P.frbath[:, iatom, 0, imode] = (dkinr - P.gkt) / P.qmass[imode]

                for inhc in range(1, P.Nnhc):
                    v_sq = P.vrbath[:, iatom, inhc-1, imode]**2
                    P.frbath[:, iatom, inhc, imode] = (P.qmass[imode] * v_sq - P.gkt) / P.qmass[imode]

                # update thermostat velocities (1st half)
                P.vrbath[:, iatom, P.Nnhc-1, imode] += 0.25 * P.frbath[:, iatom, P.Nnhc-1, imode] * dt_ys

                for inhc in range(P.Nnhc-2, -1, -1):
                    vrfact = np.exp(-0.125 * P.vrbath[:, iatom, inhc+1, imode] * dt_ys)
                    P.vrbath[:, iatom, inhc, imode] = (
                        P.vrbath[:, iatom, inhc, imode] * vrfact**2 +
                        0.25 * P.frbath[:, iatom, inhc, imode] * vrfact * dt_ys
                    )

                # update particle velocities
                pvrfact = np.exp(-0.5 * P.vrbath[:, iatom, 0, imode] * dt_ys)

                # re-update force
                P.frbath[:, iatom, 0, imode] = (pvrfact**2 * dkinr - P.gkt) / P.qmass[imode]

                # update thermostat positions
                for inhc in range(P.Nnhc):
                    P.rbath[:, iatom, inhc, imode] += 0.5 * P.vrbath[:, iatom, inhc, imode] * dt_ys

                # update thermostat velocities (2nd half)
                for inhc in range(P.Nnhc-1):
                    vrfact = np.exp(-0.125 * P.vrbath[:, iatom, inhc+1, imode] * dt_ys)
                    P.vrbath[:, iatom, inhc, imode] = (
                        P.vrbath[:, iatom, inhc, imode] * vrfact**2 +
                        0.25 * P.frbath[:, iatom, inhc, imode] * vrfact * dt_ys
                    )
                    P.frbath[:, iatom, inhc+1, imode] = (
                        P.qmass[imode] * P.vrbath[:, iatom, inhc, imode]**2 - P.gkt
                    ) / P.qmass[imode]

                P.vrbath[:, iatom, P.Nnhc-1, imode] += 0.25 * P.frbath[:, iatom, P.Nnhc-1, imode] * dt_ys

                # final velocity update
                P.vur[:, iatom, imode] *= pvrfact
