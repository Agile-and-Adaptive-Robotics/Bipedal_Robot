r"""Di Russo et al. 2023 JNE controller — faithful implementation for the
gait2392 plant (Ben's directive 2, 2026-10-25). Spec of record:
replication/dirusso_connectome_draft.json (arbitration-verified 2026-10-24:
two-oscillator pair + FIVE raised-cosine primitives).

Paper-faithful: phase-oscillator pair with hard heel-strike reset; 5 raised
cosines (sigma 0.2, mu 0.1..0.9) with signed per-muscle W; receptors
rIa=65/200*sqrt(max(0,v_norm)), rII=l_norm, rIb=f_norm; leaky-integrator
neurons tau=10 ms; MN f=clip(x,0,1); Ia mono + reciprocal (IaIN) + IaIN
mutual; II excitatory(same)/inhibitory(antagonist) INs; Ib autogenic
inhibition + IbIN mutual + extensor positive-force reversal; Renshaw
drive->own MN and own IaIN + RC<->antagonist-RC mutual (Methods p.6
reading; the Table-1-caption conflict stays flagged); NO phase gating;
Ong-style PD trunk balance on ILPSO/GMAX/HAMS in stance.

Adaptations (labeled per the B3 rule): plant = sagittal-locked converted
gait2392 (shared Delp lineage, not her custom 2D walker); muscle groups
mapped (MAP); normalizations from MuJoCo actuator quantities over the
converted lengthrange; receptor delays = one 2 ms step (hers by proximity,
not numerically extracted).

Per Ben's 2026-10-25 conductance ruling: EVERY weight is an individual
named quantity (arrays are per-muscle; no family collapse anywhere).
"""
from __future__ import annotations

import numpy as np

TAU = 0.01          # s, interneuron membrane (paper)
TAU_A = 0.02        # s, MN activation smoothing
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
EXTENSOR = {"GMAX", "VAS", "SOL"}
BALANCE_MUSCLES = ["ILPSO", "GMAX", "HAMS"]
MU = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
SIGMA = 0.2


def primitives(phi: float) -> np.ndarray:
    d = ((phi - MU + 0.5) % 1.0) - 0.5
    return np.exp(-(d ** 2) / (2 * SIGMA ** 2))


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
        for i, s in enumerate(("L", "R")):
            self.phi[i] = (self.phi[i] + self.dt * (
                self.omega - self.gamma * np.sin(
                    self.phi[i] - self.phi[1 - i] - np.pi))) % 1.0
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
            ia_n = ia + k * (g["ia_aff"] * ria - ia)
            ii_n = ii + k * (g["ii_aff"] * rii - ii)
            ib_n = ib + k * (g["ib_aff"] * rib - ib)
            # IaIN: own Ia drive, RC inhibition, antagonist-IaIN mutual inh
            iain_drive = g["ia_iain"] * ia - g["rc_iain"] * rc
            for m in range(NM):
                for a in ANT_I[m]:
                    iain_drive[m] -= g["ia_iain_mut"][a] * iain[a]
            iain_n = iain + k * (iain_drive - iain)
            # RC: MN collateral drive, antagonist-RC mutual inhibition
            rc_drive = g["rc_drive"] * self.mn[i]
            for m in range(NM):
                for a in ANT_I[m]:
                    rc_drive[m] -= g["rc_mut"][a] * rc[a]
            rc_n = rc + k * (rc_drive - rc)
            ibin_drive = ib.copy()
            for m in range(NM):
                for a in ANT_I[m]:
                    ibin_drive[m] -= g["ibin_mut"][a] * ibin[a]
            ibin_n = ibin + k * (ibin_drive - ibin)
            ibe_n = ibe + k * (ib - ibe)
            iie_n = iie + k * (ii - iie)
            iii_n = iii + k * (ii - iii)
            # MN input: signed weighted sums (per-muscle individual weights)
            p = primitives(self.phi[i])
            x = g["ff"] * (self.W[i * NM:(i + 1) * NM] @ p)
            for m in range(NM):
                x[m] += g["ia_mono"][m] * ia_n[m] + g["ii_exc"][m] * iie_n[m]
                x[m] -= g["ib_auto"][m] * np.maximum(ibin_n[m], 0.0)
                x[m] -= g["rc_mn"][m] * rc_n[m]
                x[m] += g["cut_exc"][m] * g["cut_aff"] * cut
                if MUSCLES[m] in EXTENSOR:
                    x[m] += g["ib_rev"][m] * np.maximum(ibe_n[m], 0.0)
                for a in ANT_I[m]:
                    x[m] -= g["ia_recip"][a] * np.maximum(iain_n[a], 0.0)
                    x[m] -= g["ii_inh"][a] * np.maximum(iii_n[a], 0.0)
            for mm in BALANCE_MUSCLES:
                x[MIDX[mm]] += 0.10 * u_bal * (1.0 if cut > 0.2 else 0.0)
            setattr(self, f"st_ia_{s}", ia_n)
            setattr(self, f"st_ii_{s}", ii_n)
            setattr(self, f"st_ib_{s}", ib_n)
            setattr(self, f"st_iain_{s}", iain_n)
            setattr(self, f"st_rc_{s}", rc_n)
            setattr(self, f"st_ibin_{s}", ibin_n)
            setattr(self, f"st_ibe_{s}", ibe_n)
            setattr(self, f"st_iie_{s}", iie_n)
            setattr(self, f"st_iii_{s}", iii_n)
            out[i] = x
        self.mn += (self.dt / TAU_A) * (out - self.mn)
        return np.clip(self.mn, 0.0, 1.0)
