"""SKELETON: sensor -> SNS afferent encoding. These are the SAME formulas the
spinal model uses (spinal/DESIGN.md conversion map, Deng-style figure):
    Ia rate ~ (dL/dt) / Lmax_dot      sign = stretch of own muscle
    II rate ~ (L - Lmid) / Lhalf      tonic length signal
    Ib rate ~ F / Fmax                autogenic force (festo4(P, L))
Contact and vestibular are events/tonic signals feeding the existing
HEEL_c/TOE_c and VEST input ports of the built network.

Constants (Lmid, Lhalf, Fmax, Lmax_dot) come from bpa_actuators_27.json +
maxBPAforce via config_loader below. Nothing in this file may hardcode a
muscle-specific number.
"""
from __future__ import annotations
import json


def load_config(path_to_json: str) -> dict:
    with open(path_to_json) as f:
        return json.load(f)


class AfferentEncoder:
    def __init__(self, muscle_cfg: dict):
        """muscle_cfg: per-muscle dict with Fmax [N], Lmid [m], Lhalf [m],
        Lmax_dot [m/s]. Built from bpa_actuators_27.json + design rest lengths."""
        self.m = muscle_cfg
        self.out: dict[str, float] = {}

    def update(self, name: str, L_m: float, dLdt: float, F_N: float) -> None:
        c = self.m[name]
        ia = max(-1.0, min(1.0, dLdt / c["Lmax_dot"]))          # signed: stretch +
        ii = (L_m - c["Lmid"]) / c["Lhalf"]                      # tonic, unsigned
        ib = max(0.0, min(1.0, F_N / c["Fmax"]))                 # normalized force
        self.out[f"Ia_{name}"] = ia
        self.out[f"II_{name}"] = ii
        self.out[f"Ib_{name}"] = ib
        return None

    # ---- contact / vestibular pass-throughs (already-port ports) -----------
    def contact_event(self, zone: str, pressed: bool) -> dict:
        """Heel/toe edges drive S2W phase reset: heel ON -> stance trigger of the
        ipsilateral RG; toe -> dorsiflexion inhibition chain (per Ben's rules)."""
        return {"contact": zone, "on": int(pressed)}

    def vestibular(self, roll_deg: float, pitch_deg: float, gyro: tuple) -> dict:
        """Trunk/pelvis attitude -> VEST cells (stage-4 balance inputs)."""
        return {"roll": roll_deg, "pitch": pitch_deg, "gyro": gyro}

    # ---- SNS port scaling ---------------------------------------------------
    # The built network's afferent input ports expect nA-level currents. The
    # scale below is the LAST tuning knob before the net: keep it in config.
    def to_port_current(self, rate: float, gain_nA: float) -> float:
        return gain_nA * rate
