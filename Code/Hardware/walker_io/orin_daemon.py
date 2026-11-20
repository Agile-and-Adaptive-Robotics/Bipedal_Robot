"""SKELETON: Orin Nano main loop — sensors -> afferents -> SNS -> valves,
telemetry to the lab PC, command overrides, safety watchdogs.

Loop budget (1 kHz): sensor reads ~100 us, afferents ~10 us, SNS step (numpy,
pools only) ~200-400 us, valve tick ~10 us. If the SNS step won't fit, drop the
SNS to 500 Hz and hold activation between ticks (valves stay 1 kHz).
"""
from __future__ import annotations
import time
from protocol import (UdpLink, Telemetry, Command, TELEM_PORT, CMD_PORT,
                      LAB_IP, TELEM_HZ, SENSOR_STALE_S, now)
from valves import ValveDriver
from afferents import AfferentEncoder, load_config
# from sns_rt import SnsNetwork            # TODO(hw): plug the exported SNS step
from sensors import (PressureSensor, Encoder, LiquidWire, Imu, InsoleArray)  # noqa


class Watchdogs:
    """Deliberately dumb and always-on."""

    def __init__(self, valve_driver: ValveDriver):
        self.vd = valve_driver
        self.tripped = None

    def check(self, last_sensor_t: float, pressures: dict) -> str | None:
        if now() - last_sensor_t > SENSOR_STALE_S:
            self.tripped = "sensor stale"
        elif any(p > 650.0 for p in pressures.values()):
            self.tripped = "overpressure"
        if self.tripped:
            self.vd.vent_all()
        return self.tripped


def main(config_dir: str = "."):
    telem = UdpLink(TELEM_PORT, remote=(LAB_IP, TELEM_PORT))
    cmd = UdpLink(CMD_PORT)

    # --- build IO from the frozen config (G3) -------------------------------
    cfg = load_config(f"{config_dir}/walker_io_config.json")     # channels, maps
    vd = ValveDriver(cfg["valve_channels"], mode=cfg.get("valve_mode", "bang_bang"))
    enc = {j: Encoder(**k) for j, k in cfg["encoders"].items()}          # TODO(hw)
    liq = {m: LiquidWire(**k) for m, k in cfg["liquid_wire"].items()}    # TODO(hw)
    pres = {p: PressureSensor(**k) for p, k in cfg["pressure"].items()}  # TODO(hw)
    enc_aff = AfferentEncoder(cfg["muscles"])

    overrides: dict = {"G": {}, "stim": []}          # ephemeral, from the console
    last_sensor_t = now()
    last_telem, last_hb = 0.0, 0.0
    wd = Watchdogs(vd)

    # sns = SnsNetwork(cfg["sns"], overrides["G"])  # TODO(hw): SNS step function
    t_next = now()
    while True:
        t_next += 0.001
        # 1) sensors ----------------------------------------------------------
        # joints  = {j: e.read() for j, e in enc.items()}            # TODO(hw)
        # lengths = {m: lw.read_m() for m, lw in liq.items()}        # TODO(hw)
        # press   = {p: ps.read_kpa() for p, ps in pres.items()}     # TODO(hw)
        # imu     = imu_pelvis.read()                                # TODO(hw)
        # contact = insole.read_contact()                            # TODO(hw)
        last_sensor_t = now()

        # 2) afferents ---------------------------------------------------------
        # for name, L in lengths.items():
        #     F = festo4(press[name], L) * maxBPAforce(...)
        #     enc_aff.update(name, L, dLdt(name), F)                 # TODO(hw)
        # sns.set_afferents(enc_aff.out, contact, vestib)            # TODO(hw)

        # 3) SNS step ----------------------------------------------------------
        # a = sns.step(dt=0.001)   # pool activations in [0,1]        # TODO(hw)
        a = {}

        # 4) valves -------------------------------------------------------------
        # for name, ai in a.items():
        #     vd.apply_activation(name, ai, press[name], p_cmd(name, ai))
        if wd.check(last_sensor_t, {}) or "stop" in overrides.get("flags", {}):
            break

        # 5) commands (non-blocking) --------------------------------------------
        while (m := cmd.recv(0.0)) is not None and m.get("k") == "cmd":
            if m["c"] == "override":
                overrides["G"].update(m["d"]["G"])
                # sns.apply_overrides(overrides["G"])                # TODO(hw)
            elif m["c"] == "stim":
                overrides["stim"].append(m["d"])
                # sns.inject(**m["d"])                               # TODO(hw)
            elif m["c"] == "valve_override":
                vd.apply_activation(m["d"]["name"], m["d"]["duty"],
                                    0.0, 999.0)
            elif m["c"] == "stop":
                overrides.setdefault("flags", {})["stop"] = True

        # 6) telemetry -----------------------------------------------------------
        if now() - last_telem > 1.0 / TELEM_HZ:
            last_telem = now()
            T = Telemetry(t=now(), activations=a,
                          afferents=dict(enc_aff.out), pressures={})
            telem.send(T.to_json())

    vd.vent_all()
    print("daemon exit:", wd.tripped or "clean")


if __name__ == "__main__":
    main()
