r"""Di Russo-inspired controller adapted to the converted gait2392 model.
This is NOT validated against the original SCONE implementation. The
author's source is https://github.com/DiRussoAndrea/Spinal_controller and
https://gitlab.com/simgait/SCONE. Remaining differences include individual
edge weights, length offsets, delays, normalization, and balance-state gates.
Earlier adaptation specification:
replication/dirusso_connectome_draft.json (arbitration-verified 2026-10-24:
two-oscillator pair + FIVE raised-cosine primitives).

Implemented structure: phase-oscillator pair with hard heel-strike reset; 5 raised
cosines (sigma 0.2, mu 0.1..0.9) with signed per-muscle W; receptors
rIa=65/200*sqrt(max(0,v_norm)), rII=l_norm, rIb=f_norm; leaky-integrator
neurons tau=10 ms; MN f=clip(x,0,1); Ia mono + reciprocal (IaIN) + IaIN
mutual; II excitatory(same)/inhibitory(antagonist) INs; Ib autogenic
inhibition + IbIN mutual + extensor positive-force reversal; Renshaw
drive->own MN and own IaIN + RC<->antagonist-RC mutual (Methods p.6
reading; the Table-1-caption conflict stays flagged); NO phase gating;
Ong-style PD trunk balance on ILPSO/GMAX/HAMS in stance.

Adaptations (labeled per the B3 rule): sagittal-locked converted gait2392;
muscle groups
mapped (MAP); normalizations from MuJoCo actuator quantities over the
converted lengthrange; receptor delays = one 2 ms step (hers by proximity,
not numerically extracted).

Per Ben's 2026-10-25 conductance ruling: EVERY weight is an individual
named quantity (arrays are per-muscle; no family collapse anywhere).
"""
from __future__ import annotations

import numpy as np

TAU = 0.01          # s, interneuron membrane (paper)
TAU_A = TAU         # s, MN integration, paper equations 1 and 3
MUSCLES = ["ILPSO", "GMAX", "HAMS", "BFSH", "RF", "VAS", "GAS", "SOL", "TA"]
NM = len(MUSCLES)
MIDX = {m: i for i, m in enumerate(MUSCLES)}

MAP = {
    "ILPSO": ["iliacus_{s}", "psoas_{s}"],
    "GMAX": ["glut_max1_{s}", "glut_max2_{s}", "glut_max3_{s}"],
    "HAMS": ["semimem_{s}", "semiten_{s}", "bifemlh_{s}"],
    "BFSH": ["bifemsh_{s}"],
    "RF": ["rect_fem_{s}"],
    "VAS": ["vas_med_{s}", "vas_int_{s}", "vas_lat_{s}"],
    "GAS": ["med_gas_{s}", "lat_gas_{s}"],
    "SOL": ["soleus_{s}"],
    "TA": ["tib_ant_{s}"],
}
ANTAG = {
    "ILPSO": ["GMAX", "HAMS"], "GMAX": ["ILPSO"], "HAMS": ["VAS", "RF"],
    "RF": ["HAMS", "BFSH"], "BFSH": ["VAS", "RF"], "VAS": ["HAMS", "BFSH"],
    "GAS": ["TA"], "SOL": ["TA"], "TA": ["GAS", "SOL"],
}
ANT_I = [[MIDX[a] for a in ANTAG[m]] for m in MUSCLES]
EXTENSOR = {"GMAX", "HAMS", "RF", "VAS", "GAS", "SOL"}
BALANCE_MUSCLES = ["ILPSO", "GMAX", "HAMS"]
MU = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
SIGMA = 0.2


def primitives(phi: float) -> np.ndarray:
    d = ((phi - MU + 0.5) % 1.0) - 0.5
    return np.where(np.abs(d) <= SIGMA,
                    0.5 * (1.0 + np.cos(np.pi * d / SIGMA)), 0.0)


def default_gains() -> dict:
    return dict(
        ff=1.0,
        ia_aff=np.ones(NM), ia_mono=np.ones(NM), ia_iain=np.ones(NM),
        ia_recip=np.ones(NM), ia_iain_mut=np.ones(NM),
        ii_aff=np.ones(NM), ii_exc=np.ones(NM), ii_inh=np.ones(NM),
        ib_aff=np.ones(NM), ib_auto=np.ones(NM), ib_rev=np.ones(NM),
        ibin_mut=np.ones(NM),
        rc_drive=np.ones(NM), rc_mn=np.ones(NM), rc_iain=np.ones(NM),
        rc_mut=np.ones(NM),
        cut_aff=1.0, cut_exc=np.ones(NM),
    )


class DiRussoController:
    def __init__(self, W_ff: np.ndarray, dt: float,
                 omega: float = 2 * np.pi / 1.1, gamma: float = 1.0,
                 gains: dict | None = None,
                 balance_gains: dict | None = None):
        self.dt, self.omega, self.gamma = dt, omega, gamma
        self.W = np.asarray(W_ff, dtype=float).copy()   # [2*NM, 5] signed
        self.g = default_gains()
        if gains:
            for k, v in gains.items():
                self.g[k] = float(v) if np.isscalar(v) \
                    else np.asarray(v, dtype=float)
        self.bal = balance_gains or dict(kp=8.0, kd=1.5, target=0.0)
        self.phi = np.array([0.0, 0.5])
        self.mn = np.zeros((2, NM))
        self.mn_state = np.zeros((2, NM))
        self.trunk_prev = None
        for s in "LR":
            for nm_ in ("ia", "ii", "ib", "iain", "rc", "ibin", "ibe",
                        "iie", "iii"):
                setattr(self, f"st_{nm_}_{s}", np.zeros(NM))

    @staticmethod
    def receptors(l_norm, v_norm, f_norm):
        return ((65.0 / 200.0) * np.sqrt(np.maximum(0.0, v_norm)),
                l_norm.copy(), f_norm.copy())

    def step(self, sens, heel_strike, trunk_pitch) -> np.ndarray:
        """sens: {'L': (ria,rii,rib,cut), 'R': (...)}; heel_strike {'L':bool,
        'R':bool}; trunk_pitch rad. Returns MN excitations [2, NM]."""
        k = self.dt / TAU
        # phi is stored in cycles, while omega and equation 6 use radians.
        # Update both oscillators simultaneously from the previous state.
        previous_phi = self.phi.copy()
        for i, s in enumerate(("L", "R")):
            self.phi[i] = (previous_phi[i] + self.dt * (
                self.omega - self.gamma * np.sin(2 * np.pi * (
                    previous_phi[i] - previous_phi[1 - i]) - np.pi)) /
                (2 * np.pi)) % 1.0
            if heel_strike[s]:
                self.phi[i] = 0.0
        out = np.zeros((2, NM))
        e = trunk_pitch - self.bal["target"]
        de = 0.0 if self.trunk_prev is None else (e - self.trunk_prev) / self.dt
        self.trunk_prev = e
        u_bal = -self.bal["kp"] * e - self.bal["kd"] * de
        for i, s in enumerate(("L", "R")):
            ria, rii, rib, cut = sens[s]
            g = self.g
            ia = getattr(self, f"st_ia_{s}")
            ii = getattr(self, f"st_ii_{s}")
            ib = getattr(self, f"st_ib_{s}")
            iain = getattr(self, f"st_iain_{s}")
            rc = getattr(self, f"st_rc_{s}")
            ibin = getattr(self, f"st_ibin_{s}")
            ibe = getattr(self, f"st_ibe_{s}")
            iie = getattr(self, f"st_iie_{s}")
            iii = getattr(self, f"st_iii_{s}")
            # Equation 1: SN outputs rectify; IN outputs saturate. Keep
            # integrated states separate from these firing-rate outputs.
            ia_o, ii_o, ib_o = (np.maximum(z, 0.0) for z in (ia, ii, ib))
            iain_o, rc_o, ibin_o = (np.clip(z, 0.0, 1.0)
                                   for z in (iain, rc, ibin))
            ia_n = ia + k * (g["ia_aff"] * ria - ia)
            ii_n = ii + k * (g["ii_aff"] * rii - ii)
            ib_n = ib + k * (g["ib_aff"] * rib - ib)
            # IaIN: own Ia drive, RC inhibition, antagonist-IaIN mutual inh
            iain_drive = g["ia_iain"] * ia_o - g["rc_iain"] * rc_o
            for m in range(NM):
                for a in ANT_I[m]:
                    iain_drive[m] -= g["ia_iain_mut"][a] * iain_o[a]
            iain_n = iain + k * (iain_drive - iain)
            # RC: MN collateral drive, antagonist-RC mutual inhibition
            rc_drive = g["rc_drive"] * self.mn[i]
            for m in range(NM):
                for a in ANT_I[m]:
                    rc_drive[m] -= g["rc_mut"][a] * rc_o[a]
            rc_n = rc + k * (rc_drive - rc)
            ibin_drive = ib_o.copy()
            for m in range(NM):
                for a in ANT_I[m]:
                    ibin_drive[m] -= g["ibin_mut"][a] * ibin_o[a]
            ibin_n = ibin + k * (ibin_drive - ibin)
            ibe_n = ibe + k * (ib_o - ibe)
            iie_n = iie + k * (ii_o - iie)
            iii_n = iii + k * (ii_o - iii)
            # MN input: signed weighted sums (per-muscle individual weights)
            p = primitives(self.phi[i])
            x = g["ff"] * (self.W[i * NM:(i + 1) * NM] @ p)
            for m in range(NM):
                x[m] += g["ia_mono"][m] * max(ia_n[m], 0.0) + g["ii_exc"][m] * np.clip(iie_n[m], 0.0, 1.0)
                x[m] -= g["ib_auto"][m] * np.clip(ibin_n[m], 0.0, 1.0)
                x[m] -= g["rc_mn"][m] * np.clip(rc_n[m], 0.0, 1.0)
                x[m] += g["cut_exc"][m] * g["cut_aff"] * cut
                if MUSCLES[m] in EXTENSOR:
                    x[m] += g["ib_rev"][m] * np.clip(ibe_n[m], 0.0, 1.0)
                for a in ANT_I[m]:
                    x[m] -= g["ia_recip"][a] * np.clip(iain_n[a], 0.0, 1.0)
                    x[m] -= g["ii_inh"][a] * np.clip(iii_n[a], 0.0, 1.0)
            setattr(self, f"st_ia_{s}", np.clip(ia_n, -50, 50))
            setattr(self, f"st_ii_{s}", np.clip(ii_n, -50, 50))
            setattr(self, f"st_ib_{s}", np.clip(ib_n, -50, 50))
            setattr(self, f"st_iain_{s}", np.clip(iain_n, -50, 50))
            setattr(self, f"st_rc_{s}", np.clip(rc_n, 0, 50))
            setattr(self, f"st_ibin_{s}", np.clip(ibin_n, -50, 50))
            setattr(self, f"st_ibe_{s}", np.clip(ibe_n, 0, 50))
            setattr(self, f"st_iie_{s}", np.clip(iie_n, 0, 50))
            setattr(self, f"st_iii_{s}", np.clip(iii_n, -50, 50))
            out[i] = x
        # Equation 3 clips the OUTPUT, not the integrated state. Equation 4
        # adds balance after MN integration, so inhibition cannot suppress it.
        self.mn_state += (self.dt / TAU_A) * (out - self.mn_state)
        self.mn[:] = np.clip(self.mn_state, 0.0, 1.0)
        for i, s in enumerate(("L", "R")):
            for bi, mm in enumerate(BALANCE_MUSCLES):
                signal = u_bal if np.ndim(u_bal) == 0 else u_bal[bi]
                self.mn[i, MIDX[mm]] += 0.10 * signal * (
                    1.0 if sens[s][3] > 0.2 else 0.0)
        np.clip(self.mn, 0.0, 1.0, out=self.mn)
        return self.mn
