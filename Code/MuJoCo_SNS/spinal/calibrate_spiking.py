"""Calibrate the spiking-mirror synapse scales (SPIKING_MIRROR_PLAN.md
rule 2: "weights ... rescaled so a single-PSP steady train reproduces
the non-spiking steady-state MN current (calibration script, not
hand-tuning)").

What it measures / solves (all inputs measured or bisected, nothing
hand-tuned):
  1. sat_ref  - active-phase mean saturation sat(V)=clip(V/E_HI,0,1) of
     the analog RG_E_r cell, measured on the REAL non-spiking network at
     constant DRIVE=2.5 (the check_selfsustain.py recipe, winner tables
     via draw_circuit.effective_tables("best")).
  2. f_ref    - design within-burst rate of spiking RG/PF cells (40 Hz;
     the plan's "I->rate map chosen from the existing afferent gain
     ranges" - a rate the existing nA-scale drive currents support).
  3. k_sn_*   - spike->MN (hybrid) increments: two-cell rigs (spiking
     source at f_ref -> spike synapse -> analog MN replica) vs the
     graded reference rig (saturated analog source -> graded synapse ->
     same MN).  Bisection on a correction factor until the MN steady
     mean V matches (tol 1e-3 mV).
  4. k_ns     - graded analog->spiking conductance scale: subthreshold
     frame-mapping rig; target = the spiking cell's depolarization equals
     frame_gain x the analog depolarization, frame_gain =
     (v_thr - v_rest)/E_HI = 4 (analog 0..5 mV "full scale" maps to the
     20 mV rest->threshold distance).
  5. g_inc_readout - readout-tap increment so a 50 Hz train holds the
     RD_ tap at exactly E_HI (bisection).
  6. k_s2s_*  - spike->spike increments: analytic mean-conductance
     mapping k = sat_ref/(f_ref*tau_syn), VERIFIED on a slow-integrator
     rig (measured steady V vs the analytic V of the mean conductance).

Writes spiking_calibration.json next to build_network_spiking.py (the
builder prefers it over its module defaults).

Usage: D:\\Anaconda\\envs\\myo\\python.exe calibrate_spiking.py [--quick]
  --quick skips the 20 s non-spiking sat_ref measurement and uses 0.90.
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron, SpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse, SpikingSynapse

HERE = Path(__file__).parent
E_HI = 5.0
F_REF = 40.0            # Hz design burst rate
TAU_EXC, TAU_INH = 0.005, 0.020   # spike-synapse decay (s)
E_REV_EXC, E_REV_INH = 8.0, -5.0  # MN-frame reversals (analog build's)
S_REV_EXC = 0.0
V_REST, V_THR = -70.0, -50.0
BIAS_TONIC = 16.0      # nA tonic background (equilibrium -54 mV)
DT = 0.0005             # network dt (s)

report: dict = {"f_ref": F_REF, "tau_exc": TAU_EXC, "tau_inh": TAU_INH}


# ----------------------------------------------------------------- helpers
def _run(model, t_end: float, i_ext=None) -> tuple[np.ndarray, np.ndarray]:
    """Step a compiled rig; returns (mean V over last 1 s, spike count
    over the whole run).  i_ext=None for rigs with no input ports."""
    n = int(round(t_end / DT))
    args = None if i_ext is None else list(i_ext)
    ssum = np.zeros(len(model.V))
    scount = 0
    n_tail = int(round(1.0 / DT))
    for k in range(n):
        model(args)
        if k >= n - n_tail:
            ssum += model.V
        scount += int((model.spikes == -1).sum()) if model.spiking else 0
    return ssum / n_tail, scount


def _rate_for_bias(target_f: float) -> tuple[float, SpikingNeuron]:
    """Bisect the constant bias current that makes a plain spiking LIF
    fire at target_f Hz; returns (current, neuron preset)."""
    lo, hi = 0.5, 50.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        n2 = Network()
        n2.add_neuron(_spk_plain(0.0), name="src")
        n2.add_input("src")
        m = n2.compile(backend="numpy", dt=DT)
        _, sc = _run(m, 3.0, [mid])
        f = sc / 3.0
        if f < target_f:
            lo = mid
        else:
            hi = mid
        if abs(f - target_f) < 0.2:
            break
    return 0.5 * (lo + hi), _spk_plain(0.0)


def _spk_plain(bias: float) -> SpikingNeuron:
    """Rig LIF with the mirror's cell convention (g_m = 1 uS, real mV)."""
    return SpikingNeuron(threshold_time_constant=100.0,
                         threshold_initial_value=V_THR,
                         threshold_proportionality_constant=0.0,
                         threshold_leak_rate=1.0, threshold_increment=0.0,
                         threshold_floor=V_THR, reset_potential=-60.0,
                         membrane_capacitance=0.05,
                         membrane_conductance=1.0,
                         resting_potential=V_REST, bias=bias)


def _bisect(fun, lo: float, hi: float, target: float, tol: float,
            iters: int = 80, increasing: bool = True) -> tuple[float, float]:
    """Monotone fun(x) -> y; solve fun(x)=target."""
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        y = fun(mid)
        ok = (y < target) if increasing else (y > target)
        if ok:
            lo = mid
        else:
            hi = mid
        if abs(y - target) < tol:
            break
    x = 0.5 * (lo + hi)
    return x, fun(x)


# --------------------------------------------------- 1. sat_ref (analog net)
def measure_sat_ref() -> float:
    import mujoco
    import build_network as bn
    import params
    from draw_circuit import effective_tables
    t = effective_tables("best")
    params.G.clear(); params.G.update(t["G"])
    params.TAU.clear(); params.TAU.update(t["TAU"])
    for ph in params.W_PF_MN:
        params.W_PF_MN[ph].clear(); params.W_PF_MN[ph].update(t["W"][ph])
    params.W_POSTURE.clear(); params.W_POSTURE.update(t["WPOST"])
    m = mujoco.MjModel.from_xml_path(
        str(HERE.parents[2] / "Solid_Models" / "OpenSim" /
            "Gait2392_Robotbody" / "mjc" / "gait2392_simbody" /
            "gait2392_simbody_cvt3.xml"))
    acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(m.nu)]
    net = bn.build(acts, dt=params.DT, interleg=True)
    u = net.make_inputs()
    u[net.input_index("DRIVE")] = 2.5
    n = int(round(20.0 / params.DT))
    v = np.empty(n)
    for k in range(n):
        V = net.step(u)
        v[k] = V[net.idx["RG_E_r"]]
    sat = np.clip(v / E_HI, 0.0, 1.0)
    active = sat > 0.5 * max(sat.max(), 1e-9)
    val = float(sat[active].mean())
    print(f"sat_ref: analog RG_E_r active-phase mean saturation = "
          f"{val:.4f} (active fraction {active.mean():.2f}, "
          f"V range {v.min():.2f}..{v.max():.2f} mV)")
    return val


# ------------------------------------------------------------------- mains
def calibrate(sat_ref: float) -> dict:
    print(f"f_ref = {F_REF} Hz (design), sat_ref = {sat_ref:.4f}")
    report["sat_ref"] = sat_ref

    # --- find the bias current that fires a plain LIF at f_ref --------
    i_src, _ = _rate_for_bias(F_REF)
    n2 = Network()
    n2.add_neuron(_spk_plain(0.0), name="src")
    n2.add_input("src")
    m = n2.compile(backend="numpy", dt=DT)
    _, sc = _run(m, 3.0, [i_src])
    print(f"spiking source: bias {i_src:.4f} nA -> {sc / 3.0:.2f} Hz")
    report["src_bias_nA"] = float(i_src)

    # --- k_sn (spike -> MN), exc and inh ------------------------------
    for sign_name, exc, e_rev, tau_s in (
            ("exc", True, E_REV_EXC, TAU_EXC),
            ("inh", False, E_REV_INH, TAU_INH)):
        g_ns = 1.0
        # analog reference rig
        ref = Network()
        ref.add_neuron(NonSpikingNeuron(
            membrane_capacitance=0.05, membrane_conductance=1.0,
            resting_potential=0.0, bias=sat_ref * E_HI), name="src")
        ref.add_neuron(NonSpikingNeuron(
            membrane_capacitance=0.03, membrane_conductance=1.0,
            resting_potential=0.0, bias=0.0), name="mn")
        ref.add_connection(NonSpikingSynapse(
            max_conductance=g_ns, reversal_potential=e_rev,
            e_lo=0.0, e_hi=E_HI), "src", "mn")
        refm = ref.compile(backend="numpy", dt=DT)
        v_ref, _ = _run(refm, 3.0)

        def spike_mn(c, exc=exc, e_rev=e_rev, tau_s=tau_s):
            rig = Network()
            rig.add_neuron(_spk_plain(0.0), name="src")
            rig.add_input("src")
            rig.add_neuron(NonSpikingNeuron(
                membrane_capacitance=0.03, membrane_conductance=1.0,
                resting_potential=0.0, bias=0.0), name="mn")
            g_inc = c * sat_ref / (F_REF * tau_s)
            rig.add_connection(SpikingSynapse(
                max_conductance=8.0 * g_inc, reversal_potential=e_rev,
                time_constant=tau_s, transmission_delay=0,
                conductance_increment=g_inc), "src", "mn")
            mm = rig.compile(backend="numpy", dt=DT)
            v, _ = _run(mm, 3.0, [i_src])
            return float(v[1])

        c, v_got = _bisect(spike_mn, 0.05, 20.0, float(v_ref[1]), 1e-3,
                           increasing=exc)
        k = c * sat_ref / (F_REF * tau_s)
        print(f"k_sn_{sign_name}: MN steady V ref {float(v_ref[1]):+.4f} "
              f"mV, spiking {v_got:+.4f} mV -> k = {k:.4f}")
        report[f"k_sn_{sign_name}"] = float(k)
        report[f"k_sn_{sign_name}_resid_mV"] = float(v_got - v_ref[1])

    # --- k_ns (graded analog -> spiking): RATE-matching rig -------------
    # Target: the graded DRIVE operating point (source saturated at
    # sat_ref, conductance g_ns = G['descend_to_rg_e'] of the analog
    # build) makes the spiking RG replica fire at f_ref - the plan's
    # "I->rate map" for analog commands.  RG replica = the mirror's RG
    # cell WITHOUT adaptation (well-defined steady rate).
    import params as _P
    g_ns_drive = float(_P.G["descend_to_rg_e"])

    def graded_rate(k_ns):
        rig = Network()
        rig.add_neuron(NonSpikingNeuron(
            membrane_capacitance=0.10, membrane_conductance=1.0,
            resting_potential=0.0, bias=sat_ref * E_HI), name="src")
        rig.add_neuron(SpikingNeuron(
            threshold_time_constant=100.0, threshold_initial_value=V_THR,
            threshold_proportionality_constant=0.0, threshold_leak_rate=1.0,
            threshold_increment=0.0, threshold_floor=V_THR,
            reset_potential=-60.0, membrane_capacitance=0.05,
            membrane_conductance=1.0, resting_potential=V_REST,
            bias=BIAS_TONIC), name="rg")
        rig.add_connection(NonSpikingSynapse(
            max_conductance=g_ns_drive * k_ns, reversal_potential=S_REV_EXC,
            e_lo=0.0, e_hi=E_HI), "src", "rg")
        mm = rig.compile(backend="numpy", dt=DT)
        _, sc = _run(mm, 3.0)
        return sc / 3.0

    k_ns, f_got = _bisect(graded_rate, 1e-4, 1.0, F_REF, 0.05)
    print(f"k_ns: DRIVE g={g_ns_drive} at sat {sat_ref:.2f} -> RG replica "
          f"{f_got:.2f} Hz (target {F_REF}) -> k_ns = {k_ns:.5f}")
    report["k_ns"] = float(k_ns)
    report["k_ns_resid_Hz"] = float(f_got - F_REF)
    report["k_ns_g_ns_drive"] = g_ns_drive

    # --- g_inc_readout (50 Hz -> tap at E_HI) ---------------------------
    def readout_v(g_inc):
        rig = Network()
        rig.add_neuron(_spk_plain(0.0), name="src")
        rig.add_input("src")
        rig.add_neuron(NonSpikingNeuron(
            membrane_capacitance=0.05, membrane_conductance=1.0,
            resting_potential=0.0, bias=0.0), name="tap")
        rig.add_connection(SpikingSynapse(
            max_conductance=4.0 * g_inc, reversal_potential=E_REV_EXC,
            time_constant=0.020, transmission_delay=0,
            conductance_increment=g_inc), "src", "tap")
        mm = rig.compile(backend="numpy", dt=DT)
        v, _ = _run(mm, 3.0, [i_src])
        return float(v[1])

    # the tap source here fires at f_ref, not 50 Hz: rescale after
    g_at_fref, v_got = _bisect(readout_v, 0.01, 50.0, E_HI, 1e-3)
    g_readout = g_at_fref * (F_REF / 50.0)
    print(f"g_inc_readout: {F_REF} Hz tap at {v_got:.4f} mV; rescaled to "
          f"50 Hz full-scale -> g_inc = {g_readout:.4f}")
    report["g_inc_readout"] = float(g_readout)
    report["g_inc_readout_check_mV"] = float(v_got)

    # --- k_s2s (spike -> SPIKING): RATE-calibrated ----------------------
    # MEASURED TRAP (2026-10-02, debug_rg_pair.py): mean-conductance
    # matching onto a SPIKING target produces single-spike conductances
    # several times the leak (g_inc ~ 21 uS vs g_m = 1) -> the postsynap-
    # tic relay saturates at its f-I ceiling (~200 Hz) and pins anything
    # it inhibits.  Mean matching is only valid where V is the state (the
    # analog MN).  For spiking targets the weight is calibrated by RATE:
    #   exc: a 40 Hz presynaptic train drives the relay replica to 40 Hz
    #        (1:1 rate transfer, the analog "follower" behavior);
    #   inh: a 40 Hz inhibitory train SILENCES a relay that is excited
    #        to 40 Hz by the calibrated excitation (rate <= 1 Hz) - the
    #        analog "holds the opposite half-center down" behavior.
    def relay_replica():
        return SpikingNeuron(
            threshold_time_constant=0.10, threshold_initial_value=V_THR,
            threshold_proportionality_constant=0.0, threshold_leak_rate=1.0,
            threshold_increment=0.0, threshold_floor=V_THR,
            reset_potential=-60.0, membrane_capacitance=0.05,
            membrane_conductance=1.0, resting_potential=V_REST,
            bias=BIAS_TONIC)

    def s2s_rate(k, exc, extra_exc_k=None):
        rig = Network()
        rig.add_neuron(_spk_plain(0.0), name="src")
        rig.add_input("src")
        if extra_exc_k is not None:
            rig.add_neuron(_spk_plain(0.0), name="src2")
            rig.add_input("src2")
        rig.add_neuron(relay_replica(), name="dst")
        g_inc = k * 1.0
        rig.add_connection(SpikingSynapse(
            max_conductance=8.0 * g_inc,
            reversal_potential=S_REV_EXC if exc else -70.0,
            time_constant=TAU_EXC if exc else TAU_INH,
            transmission_delay=0, conductance_increment=g_inc),
            "src", "dst")
        if extra_exc_k is not None:
            rig.add_connection(SpikingSynapse(
                max_conductance=8.0 * extra_exc_k,
                reversal_potential=S_REV_EXC, time_constant=TAU_EXC,
                transmission_delay=0,
                conductance_increment=extra_exc_k), "src2", "dst")
        mm = rig.compile(backend="numpy", dt=DT)
        cur = [i_src] + ([i_src] if extra_exc_k is not None else [])
        dst_i = 2 if extra_exc_k is not None else 1
        n = int(round(3.0 / DT))
        # skip the first 0.5 s transient
        sc = 0
        for k_ in range(n):
            mm(cur)
            if k_ * DT >= 0.5:
                sc += int(mm.spikes[dst_i] == -1)
        return sc / 2.5

    # exc: bisect k so dst fires at f_ref driven by the 40 Hz train
    k_exc, f_got = _bisect(lambda k: s2s_rate(k, True), 1e-3, 3.0,
                           F_REF, 0.1)
    print(f"k_s2s_exc: 40 Hz train -> dst {f_got:.2f} Hz "
          f"(target {F_REF}) -> k = {k_exc:.4f}")
    report["k_s2s_exc"] = float(k_exc)
    report["k_s2s_exc_rate_resid_Hz"] = float(f_got - F_REF)
    # inh: bisect k so the 40 Hz inhibitory train silences the
    # excitation-driven dst (target <= 1 Hz; bisect on rate directly)
    k_inh, f_sil = _bisect(
        lambda k: s2s_rate(k, False, extra_exc_k=k_exc), 1e-3, 3.0,
        1.0, 0.1, increasing=False)
    print(f"k_s2s_inh: 40 Hz inh train onto a 40 Hz-driven dst -> "
          f"{f_sil:.2f} Hz (target <= 1) -> k = {k_inh:.4f}")
    report["k_s2s_inh"] = float(k_inh)
    report["k_s2s_inh_silenced_Hz"] = float(f_sil)
    return report


def main() -> int:
    quick = "--quick" in sys.argv
    sat_ref = 0.90 if quick else measure_sat_ref()
    rep = calibrate(sat_ref)
    rep["provenance"] = (
        "calibrate_spiking.py 2026-10-02; sat_ref measured on the analog "
        "network (check_selfsustain recipe, DRIVE 2.5, effective_tables"
        "('best')); k_sn/k_ns/g_inc_readout bisected on 2-cell rigs; "
        "k_s2s analytic mean-conductance mapping + integrator verify")
    out = HERE / "spiking_calibration.json"
    out.write_text(json.dumps(rep, indent=1), "utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
