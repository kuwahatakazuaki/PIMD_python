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
        vur = P.vur[:, :, 0]
        vr = P.vrbc31
        fr = P.frbc31
        rb = P.rbc31

        # Kinetic energy for all atoms at the centroid mode.
        dkinr = P.fictmass[:, 0][None, :] * vur**2

        # Update the force for the first thermostat.
        fr[:, :, 0] = (dkinr - P.gkt) / P.qmcent31[0]

        # Higher chain levels remain sequential in chain index, but vectorized over atoms.
        for inhc in range(1, P.Nnhc):
            prev_v = vr[:, :, inhc - 1]
            fr[:, :, inhc] = (
                P.qmcent31[inhc - 1] * prev_v**2 - P.gkt
            ) / P.qmcent31[inhc]

        # Update thermostat velocities (first half).
        vr[:, :, P.Nnhc - 1] += 0.25 * fr[:, :, P.Nnhc - 1] * dt_ys

        for inhc in range(P.Nnhc - 2, -1, -1):
            v = vr[:, :, inhc + 1]
            vrfact = np.exp(-0.125 * v * dt_ys)
            vr[:, :, inhc] = (
                vr[:, :, inhc] * vrfact**2
                + 0.25 * fr[:, :, inhc] * vrfact * dt_ys
            )

        # Update particle velocity scaling.
        pvrfact = np.exp(-0.5 * vr[:, :, 0] * dt_ys)

        # Update force for the first thermostat again.
        fr[:, :, 0] = (pvrfact**2 * dkinr - P.gkt) / P.qmcent31[0]

        # Update thermostat positions.
        rb += 0.5 * vr * dt_ys

        # Update thermostat velocities (second half).
        for inhc in range(P.Nnhc - 1):
            v = vr[:, :, inhc + 1]
            vrfact = np.exp(-0.125 * v * dt_ys)
            vr[:, :, inhc] = (
                vr[:, :, inhc] * vrfact**2
                + 0.25 * fr[:, :, inhc] * vrfact * dt_ys
            )

            fr[:, :, inhc + 1] = (
                P.qmcent31[inhc] * vr[:, :, inhc]**2 - P.gkt
            ) / P.qmcent31[inhc + 1]

        # Final velocity update for last chain level.
        vr[:, :, P.Nnhc - 1] += 0.25 * fr[:, :, P.Nnhc - 1] * dt_ys

        # Final particle velocity update.
        vur *= pvrfact


# def nhc_integrate():
#     """
#     Newer version of Velocity-Verlet
#     """
#     for iys in range(P.Nys):
#         dt_ys = P.dt_ref * P.ysweight[iys]

#         for imode in range(1, P.Nbead):
#             qmass = P.qmass[imode]

#             for iatom in range(P.Natom):
#                 # --- Step 1: dkinr and initial force on first NHC ---
#                 vur = P.vur[:, iatom, imode]
#                 dkinr = P.fictmass[iatom, imode] * vur**2
#                 P.frbath[:, iatom, 0, imode] = (dkinr - P.gkt) / qmass

#                 # --- Step 2: update force on rest of chain ---
#                 vr = P.vrbath[:, iatom, :, imode]
#                 fr = P.frbath[:, iatom, :, imode]

#                 for inhc in range(1, P.Nnhc):
#                     prev_v_sq = vr[:, inhc - 1]**2
#                     fr[:, inhc] = (qmass * prev_v_sq - P.gkt) / qmass

#                 # --- Step 3: update vr (1st half) ---
#                 vr[:, -1] += 0.25 * fr[:, -1] * dt_ys
#                 for inhc in reversed(range(P.Nnhc - 1)):
#                     v = vr[:, inhc + 1]
#                     vrfact = np.exp(-0.125 * v * dt_ys)
#                     vr[:, inhc] = vr[:, inhc] * vrfact**2 + 0.25 * fr[:, inhc] * vrfact * dt_ys

#                 # --- Step 4: particle velocity update ---
#                 pvrfact = np.exp(-0.5 * vr[:, 0] * dt_ys)
#                 vur *= pvrfact

#                 # --- Step 5: re-update fr (first NHC) ---
#                 dkinr = P.fictmass[iatom, imode] * vur**2
#                 fr[:, 0] = (pvrfact**2 * dkinr - P.gkt) / qmass

#                 # --- Step 6: update thermostat positions ---
#                 P.rbath[:, iatom, :, imode] += 0.5 * vr * dt_ys

#                 # --- Step 7: vr (2nd half) ---
#                 for inhc in range(P.Nnhc - 1):
#                     v = vr[:, inhc + 1]
#                     vrfact = np.exp(-0.125 * v * dt_ys)
#                     vr[:, inhc] = vr[:, inhc] * vrfact**2 + 0.25 * fr[:, inhc] * vrfact * dt_ys
#                     fr[:, inhc + 1] = (qmass * vr[:, inhc]**2 - P.gkt) / qmass

#                 # --- Step 8: final velocity update ---
#                 vr[:, -1] += 0.25 * fr[:, -1] * dt_ys

#                 # --- Step 9: apply updated velocity ---
#                 P.vur[:, iatom, imode] = vur


def nhc_integrate():
    """
    Velocity-Verlet形式のNose-Hooverチェーン（非中心モード）を時間発展。
    """
    for iys in range(P.Nys):
        dt_ys = P.dt_ref * P.ysweight[iys]

        for imode in range(1, P.Nbead):  # Fortranの2:Nbead → Pythonの1:Nbead-1
            qmass = P.qmass[imode]
            vur = P.vur[:, :, imode]
            vr = P.vrbath[:, :, :, imode]
            fr = P.frbath[:, :, :, imode]
            rb = P.rbath[:, :, :, imode]
            dkinr = P.fictmass[:, imode][None, :] * vur**2

            # Update force on first NHC.
            fr[:, :, 0] = (dkinr - P.gkt) / qmass

            for inhc in range(1, P.Nnhc):
                v_sq = vr[:, :, inhc - 1]**2
                fr[:, :, inhc] = (qmass * v_sq - P.gkt) / qmass

            # Update thermostat velocities (1st half).
            vr[:, :, P.Nnhc - 1] += 0.25 * fr[:, :, P.Nnhc - 1] * dt_ys

            for inhc in range(P.Nnhc - 2, -1, -1):
                vrfact = np.exp(-0.125 * vr[:, :, inhc + 1] * dt_ys)
                vr[:, :, inhc] = (
                    vr[:, :, inhc] * vrfact**2
                    + 0.25 * fr[:, :, inhc] * vrfact * dt_ys
                )

            # Update particle velocities.
            pvrfact = np.exp(-0.5 * vr[:, :, 0] * dt_ys)

            # Re-update force while preserving the original algebra.
            fr[:, :, 0] = (pvrfact**2 * dkinr - P.gkt) / qmass

            # Update thermostat positions.
            rb += 0.5 * vr * dt_ys

            # Update thermostat velocities (2nd half).
            for inhc in range(P.Nnhc - 1):
                vrfact = np.exp(-0.125 * vr[:, :, inhc + 1] * dt_ys)
                vr[:, :, inhc] = (
                    vr[:, :, inhc] * vrfact**2
                    + 0.25 * fr[:, :, inhc] * vrfact * dt_ys
                )
                fr[:, :, inhc + 1] = (
                    qmass * vr[:, :, inhc]**2 - P.gkt
                ) / qmass

            vr[:, :, P.Nnhc - 1] += 0.25 * fr[:, :, P.Nnhc - 1] * dt_ys

            # Final velocity update.
            vur *= pvrfact
