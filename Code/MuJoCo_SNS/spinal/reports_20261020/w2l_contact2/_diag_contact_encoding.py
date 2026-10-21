r"""CONTACT-ENCODER DIAGNOSTIC for the tuned winner (Ben's ask 2026-10-20:
"is your heel and crossed edge weak BECAUSE THE MODEL IS IN THE AIR AND
NOT MAKING CONTACT? Test that.").

Runs the round_best_r030 tuned winner (obj 91.36, survive30 fresh) for
20 s via test_w2l_ground.run_phase IN-PROCESS with the campaign's exact
knob recipe (EV.BASE protocol knobs + the winner's params), then answers,
with numbers over the walk window (t >= hold+lower+0.5):

  - heel/toe PORT CURRENTS (heelin/toein logs, nA): mean / max per side;
  - fraction of neural steps ACTIVE (0 < I < ctn_amp), SATURATED
    (F_plate >= cref, the encoder's clip ceiling), SILENT (I == 0);
  - plate-force percentiles vs cref (how deep into the encoder's range
    the physics actually drives it);
  - contact-NEURON voltages (heelV/toeV logs, mV): mean/max per side
    (are the encoders' postsynaptic membranes even moving?).

Verbatim baseline behavior: no code edits ride this script - it reads
the existing logs only.
"""
from __future__ import annotations

import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
W2L = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
for p in (W2L, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import numpy as np                     # noqa: E402

_argv = sys.argv[:]
sys.argv = ["test_w2l_ground.py"]      # keep the test's argv parser quiet
import test_w2l_ground as T            # noqa: E402
sys.argv = _argv

WINNER = json.load(open(os.path.join(HERE, "round_best_r030.json"),
                        encoding="utf-8"))
BASE = dict(legcoll=0.0, sink=0.06, rig=0.25, walker_yield=1.5)


def stats(x: np.ndarray) -> dict:
    return dict(mean=float(x.mean()), max=float(x.max()),
                p50=float(np.percentile(x, 50)),
                p90=float(np.percentile(x, 90)),
                p99=float(np.percentile(x, 99)))


def main() -> int:
    kv = dict(BASE, net="bilateralrg2", dur_walk=20.0)
    kv.update({k: float(v) for k, v in WINNER["params"].items()})
    T.KV.update(kv)
    print(f"knobs: {kv}", flush=True)
    wlog = T.run_phase("walk", T.KV["dur_walk"], T.KV["rig"], None)

    t0 = T.KV["hold"] + T.KV["lower"] + 0.5      # the walk metric window
    w = wlog["t"] >= t0
    amp, cref = T.KV["ctn_amp"], T.KV["cref"]
    n_steps = int(w.sum())
    print(f"\nfinite={wlog['finite']} "
          f"fall_t={wlog['fall_t']} window t>={t0} s "
          f"({n_steps} physics rows, neural-step-sampled encoders)")
    print(f"encoder: I = ctn_amp*clip(F/cref, 0, 1), "
          f"ctn_amp={amp:.3f} nA, cref={cref:.2f} N")
    out = dict(knobs=kv, finite=bool(wlog["finite"]),
               fall_t=wlog["fall_t"], ctn_amp=amp, cref=cref)
    for kind, LOGI, LOGV, LOGF in (("heel", "heelin", "heelV", "heelF"),
                                   ("toe", "toein", "toeV", "toeF")):
        cur = wlog[LOGI][w]              # nA at neural steps (held fwd)
        vol = wlog[LOGV][w]              # mV
        force = wlog[LOGF][w]            # N (per physics step)
        for k, s in enumerate(("L", "R")):
            c, v, f = cur[:, k], vol[:, k], force[:, k]
            sat = float((f >= cref).mean())        # encoder clip ceiling hit
            act = float(((c > 1e-9) & (f < cref)).mean())
            sil = float((c <= 1e-9).mean())
            d = dict(port_current=stats(c), neuron_V=stats(v),
                     plate_force_N=stats(f),
                     frac_steps_saturated=sat, frac_steps_partial=act,
                     frac_steps_silent=sil,
                     frac_force_gt5N=float((f > 5.0).mean()))
            key = f"{kind}_{s}"
            out[key] = d
            print(f"\n{kind} {s}: port current nA mean {c.mean():.3f} "
                  f"max {c.max():.3f} (p50 {np.percentile(c,50):.3f}, "
                  f"p90 {np.percentile(c,90):.3f}, "
                  f"p99 {np.percentile(c,99):.3f}) | ceiling {amp:.3f}")
            print(f"  steps: SATURATED (F>=cref) {sat:.3f} | PARTIAL "
                  f"(0<I<ceiling) {act:.3f} | SILENT {sil:.3f}")
            print(f"  plate force N: mean {f.mean():.1f} max {f.max():.1f} "
                  f"(p50 {np.percentile(f,50):.1f}, p90 "
                  f"{np.percentile(f,90):.1f}, p99 {np.percentile(f,99):.1f})"
                  f" | F>5N {(f>5.0).mean():.3f} | F/cref at p99 "
                  f"{np.percentile(f,99)/cref:.3f}")
            print(f"  contact-neuron V mV: mean {v.mean():.3f} "
                  f"max {v.max():.3f} p99 {np.percentile(v,99):.3f} "
                  f"(0..5 mV scale, drive ceiling at ~{amp}*tau/C ratio)")
    path = os.path.join(HERE, "diag_contact_encoding.json")
    json.dump(out, open(path, "w", encoding="utf-8"), indent=1)
    print(f"\nwritten: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
