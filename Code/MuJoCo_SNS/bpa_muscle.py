"""
Custom BPA (braided pneumatic actuator / pneumatic artificial muscle) force
model for MuJoCo + SNS-Toolbox.

Python port of Ben Bolen's MATLAB model:

  festo4.m                          normalized force surface (sfit fit)
      F_norm(rel, P) = a0*(exp(-a1*rel) - 1) + P*exp(-a3*rel^2)
      rel = relative strain = contraction / KMAX, P = pressure / 620 kPa
      clamps: F=0 for rel>1; F=0 for F<0.  Result is a fraction of Fmax.

  maxBPAforce.m                     max force at zero strain
      10/20 mm: Fmax = P * a1 * atan(a2 * (L_rest - 0.0075 m) * P)
      40 mm:    Fmax = 6398.4 N (Festo tool, length-independent)

  MonoPamDataExplicit_balanceX3.m   stiffness-aware contraction & equilibrium
      KMAX   = (rest - kmax) / rest           max contraction fraction
      Xi0    constant length offset (m)
      Xi3    bend-loss: delta_L = Xi3 * bend_measure * comp^2,
             comp = max(0, 1 - relstrain_first_pass)
      series equilibrium (fortz): the muscle force balances the series
      stiffness (tendon spring rate + bracket compliance) through the
      stretch r:
          festo4(dia, (rest-(Lm-r))/rest/KMAX, P) * mif  =  k_eff * r
      solved on r in [0, Lm - kmax]; force is re-evaluated at the solved r.

The MuJoCo wiring treats each BPA as a site-based <tendon> with a zero-gain
<general> actuator so MuJoCo exposes actuator_length / actuator_velocity /
actuator_moment; the actual force is computed HERE each step and applied via
qfrc_applied (see bpa_mujoco.py). This keeps the force law entirely in Ben's
model instead of MuJoCo's Hill-type muscle.

Velocity: Ben's model is quasi-static (no damping term). An optional small
series damping can be enabled for simulation stability; it is an engineering
addition, not part of the identified model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
from scipy.optimize import brentq

__all__ = ["festo4", "maxBPAforce", "tendon_spring_rate", "BPAMuscle"]

# (a0, a1, a3) of the sfit surface per diameter, from FestoLookup.mat via
# export_festo_lookup.m (festo_lookup_coeffs.json); x=rel strain, y=P/620
FESTO_COEFFS = {
    10: (0.568207874671, 4.25442545542, 0.55972777762),
    20: (0.257852586017, 6.4766142989, 1.32087718059),
    40: (0.122366088343, 10.4714342293, 2.02328814095),
}

# (a1, a2) of maxBPAforce.m per diameter; 40 mm is a constant 6398.4 N
MAXBPA_COEFFS = {10: (0.4895, 0.03068), 20: (1.4877, 0.0248)}
MAXBPA_40MM_N = 6398.4

# tendon spring rate constants from Spr() in MonoPamDataExplicit_balanceX3.m
CABLE_AREA_M2 = 1.51e-6   # 19-strand steel cable effective area
CABLE_E_PA = 193e9        # Young's modulus
_CABLE_MULT = {10: 2, 20: 6}   # default wrap multipliers (40 mm -> 2)


def festo4(dia: float, rel, pres_kpa) -> np.ndarray:
    """Normalized BPA force (fraction of Fmax), port of festo4.m."""
    rel = np.asarray(rel, dtype=float)
    pres_kpa = np.asarray(pres_kpa, dtype=float)
    a0, a1, a3 = FESTO_COEFFS[int(dia)]
    p = pres_kpa / 620.0
    f = a0 * (np.exp(-a1 * rel) - 1.0) + p * np.exp(-a3 * rel * rel)
    f = np.where(rel > 1.0, 0.0, f)
    return np.maximum(f, 0.0)


def maxBPAforce(rest_m: float, dia: float, pres_kpa: float = 620.0) -> float:
    """Maximum force (N) at zero strain, port of maxBPAforce.m."""
    if int(dia) == 40:
        return MAXBPA_40MM_N
    a1, a2 = MAXBPA_COEFFS[int(dia)]
    return pres_kpa * a1 * np.arctan(a2 * (rest_m - 0.0075) * pres_kpa)


def tendon_spring_rate(tendon_len_m: float, dia: float, wraps: int = 1,
                       bpa_count: int = 1) -> float:
    """Effective tendon/cable spring rate (N/m), port of Spr()."""
    if tendon_len_m <= 0:
        return np.inf
    mult = _CABLE_MULT.get(int(dia), 2) * (wraps if wraps else 1)
    return bpa_count * mult * CABLE_AREA_M2 * CABLE_E_PA / tendon_len_m


@dataclass
class BPAMuscle:
    """One BPA muscle-tendon unit (optionally several BPAs in parallel).

    Set pressure via set_activation(a) with a in [0,1] -> a*pressure_max kPa,
    then call force(musculotendon_length) each simulation step.
    """

    name: str
    diameter: float                  # 10, 20 or 40 (mm)
    resting_length: float            # BPA sleeve resting length (m)
    kmax_length: float               # fully contracted BPA length (m)
    fitting_length: float = 0.0      # one end fitting length (m)
    tendon_length: float = 0.0       # cable length (m)
    bpa_count: int = 1               # parallel BPAs
    pressure_max_kpa: float = 620.0  # activation 1.0 -> this pressure
    wraps: int = 1                   # cable wraps (scales tendon stiffness)

    # stiffness-aware terms from the Xi pipeline (all optional)
    xi0: float = 0.0                 # constant length offset (m)
    xi3: float = 0.0                 # bend-loss scale factor
    bend_measure: Callable[[], float] = None   # R*alpha (m) at current pose
    series_stiffness: Optional[float] = None   # N/m; default = tendon rate
    extra_damping: float = 0.0       # N per (m/s), optional, NOT in Ben's model

    pressure_kpa: float = field(default=0.0, init=False, repr=False)

    def __post_init__(self):
        self.fmax = maxBPAforce(self.resting_length, self.diameter)
        if self.series_stiffness is None:
            self.series_stiffness = tendon_spring_rate(
                self.tendon_length, self.diameter, self.wraps, self.bpa_count)
        self.set_activation(0.0)

    # ------------------------------------------------------------------ state
    @property
    def kmax(self) -> float:
        """Max contraction fraction, KMAX in the MATLAB code."""
        return (self.resting_length - self.kmax_length) / self.resting_length

    def set_activation(self, a: float) -> float:
        """Map SNS/neural activation in [0,1] to supply pressure (kPa)."""
        a = float(np.clip(a, 0.0, 1.0))
        self.pressure_kpa = a * self.pressure_max_kpa
        return self.pressure_kpa

    # ------------------------------------------------------------------ model
    def _delta_l(self, lm, relstrain):
        """Xi3 bend-loss (m); comp from the first-pass relative strain."""
        if self.xi3 == 0.0 or self.bend_measure is None:
            return 0.0
        comp = np.clip(1.0 - relstrain, 0.0, None)
        return self.xi3 * self.bend_measure() * comp ** 2

    def force(self, lmt: float, velocity: float = 0.0) -> float:
        """Total BPA force (N) at musculotendon length lmt (m).

        lmt is the full route length MuJoCo reports for the tendon
        (actuator_length / ten_length). Replicates the balanceX3 pipeline:
        first-pass Xi3 loss -> series-equilibrium stretch -> final strain.
        """
        rest = self.resting_length
        l0 = lmt - self.tendon_length - 2.0 * self.fitting_length - self.xi0

        # first pass (no deformation): strain for the Xi3 comp factor
        relstrain0 = (rest - l0) / rest / self.kmax
        delta_l = self._delta_l(l0, relstrain0)

        # series equilibrium (fortz): muscle force vs series stiffness
        lm = l0 - delta_l
        mif = self.bpa_count * self.fmax          # max total force
        k_eff = self.series_stiffness

        contraction0 = (rest - lm) / rest
        relstrain = contraction0 / self.kmax
        if relstrain >= 1.0:
            # would need to be shorter than fully contracted: no force
            return 0.0

        def balance(r):
            rel = (rest - (lm - r)) / rest / self.kmax
            return float(festo4(self.diameter, rel, self.pressure_kpa)) * mif - k_eff * r

        r = 0.0
        f0 = balance(0.0)
        if f0 > 0.0 and k_eff > 0.0 and lm > self.kmax_length:
            try:
                r = brentq(balance, 0.0, lm - self.kmax_length)
            except ValueError:
                r = 0.0

        rel_f = (rest - (lm - r)) / rest / self.kmax
        f_mag = float(festo4(self.diameter, rel_f, self.pressure_kpa)) * mif
        return max(f_mag, 0.0) - self.extra_damping * velocity
