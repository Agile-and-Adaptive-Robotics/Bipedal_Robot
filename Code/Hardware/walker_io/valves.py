"""SKELETON: Festo valve manifold driver (fill/vent per BPA line).

Two modes, both proven on the bench:
  bang_bang           closed-loop hysteresis on pressure error (default, safest)
  proportional_taped  short PWM packets; duty ~= activation (requires the
                      bench duty->dP/dt calibration table filled in below)

Hardware access points are marked TODO(hw): swap the stub Driver for the real
GPIO/PWM layer (PCA9685 or 74HC595 + MOSFET boards on the Orin).
"""
from __future__ import annotations
import time
from protocol import now


class ValveChannel:
    def __init__(self, name: str, fill_idx: int, vent_idx: int,
                 p_max_kpa: float = 620.0):
        self.name = name
        self.fill_idx, self.vent_idx = fill_idx, vent_idx
        self.p_max = p_max_kpa
        self.state = "neutral"          # fill | vent | neutral
        self.duty = 0.0
        self.last_change = 0.0


class ValveDriver:
    """One instance owns the whole manifold (channels = BPAs)."""

    def __init__(self, channels: list[tuple[str, int, int]], mode: str = "bang_bang",
                 pwm_hz: float = 200.0):
        self.mode = mode
        self.pwm_hz = pwm_hz
        self.ch = {c[0]: ValveChannel(*c) for c in channels}
        # TODO(hw): open the PWM/GPIO device here (board config per G3 schematic)
        self.outputs = {i: 0.0 for c in channels for i in c[1:3]}

    # ---- calibration (fill from the bench bring-up, valves.py companion csv) --
    # duty -> dP/dt (kPa/s) for the taped-proportional mode, per channel family
    DUTY_TABLE = {
        # "20mm": {0.2: 40.0, 0.4: 90.0, 0.6: 150.0, 0.8: 210.0, 1.0: 260.0},
    }

    def apply_activation(self, name: str, a: float, p_meas: float, p_cmd: float) -> str:
        """One control tick for one BPA. a in [0,1] from the SNS; p_cmd from the
        pressure mapping a->pressure (or manual override). Returns new state."""
        v = self.ch[name]
        if self.mode == "bang_bang":
            err = p_cmd - p_meas
            if p_meas > v.p_max:
                new = "vent"
            elif err > 15.0:
                new = "fill"
            elif err < -15.0:
                new = "vent"
            else:
                new = "neutral"
        else:  # proportional_taped: duty drives PWM packet density
            v.duty = max(0.0, min(1.0, a))
            new = "fill" if v.duty > 0.02 else "neutral"
            # TODO(hw): schedule PWM packets at self.pwm_hz with duty v.duty
        if new != v.state and now() - v.last_change > 0.005:   # 5 ms debounce
            v.state, v.last_change = new, now()
            self._write(v)
        return v.state

    def vent_all(self) -> None:
        """Safety neutral: everything to vent (exhaust) immediately."""
        for v in self.ch.values():
            v.state = "vent"
            self._write(v)

    def _write(self, v: ValveChannel) -> None:
        # TODO(hw): set the physical fill/vent outputs for this channel here.
        self.outputs[v.fill_idx] = 1.0 if v.state == "fill" else 0.0
        self.outputs[v.vent_idx] = 1.0 if v.state == "vent" else 0.0
