"""TEMPLATE-MODEL VERIFICATION: Rybak et al. 2006a reduced two-level CPG
(activity-based reinterpretation; see rybak2006_standalone.py for circuit,
parameter provenance and flags F-A..F-I).

LITERATURE CLAIMS TESTED (quotes from the local full text,
spinal/lit_rybak2006_fulltext.txt):

  EXP A - DRIVE ASYMMETRY (Fig 2B/C): "In B, the MLR drive to RG-E
      (d rg-e) population is larger than to RG-F population
      (d rg-f = 0.43; d rg-e = 0.5), and the model generates a rhythm
      with a longer duration extensor phase (T_E > T_F). In C, the RG-F
      population receives a larger drive (d rg-f = 0.51; d rg-e = 0.45),
      and the model generates a flexion-dominated rhythm (T_F > T_E)."
      GATES: regime (0.5, 0.43) -> T_E > T_F; regime (0.45, 0.51) ->
      T_F > T_E; rhythm present in both (>= 5 RG episodes/side,
      period cv <= 0.05).

  EXP B - DELETION DIAGNOSTIC (the paper's central two-level claim):
      B1 NON-RESETTING at the PF level (Figs 8-11; Fig 10Ba drives
      d_rg-e = 0.55, d_rg-f = 0.44): "a deletion of flexor activity was
      produced by increased excitatory drive to the PF-E population ...
      As in the experimental recordings, the phase of the locomotor
      rhythm is maintained after the deletion." (Fig 10Aa: "two bursts
      in the iTA trace were omitted".)
      GATES: (i) flexor MN bursts deleted inside the pulse window
      while extensor MN activity stays sustained; (ii) the RG-E onset
      grid is unchanged through and after the pulse (max residual
      <= 0.15 period; the RG receives no PF input - Table 2); (iii) the
      first post-deletion Mn-F burst reappears ON the extrapolated
      pre-deletion grid (residual <= 0.15 period).
      B2 RESETTING at the RG level (Fig 6A; caption drives d_rg-e = 0.5,
      d_rg-f = 0.46): "the deletion of flexor activity was produced by a
      temporary 400% increase in excitatory drive to the RG-E
      population ... The increased drive changed RG-E population
      activity from phasic to sustained ... when the applied
      'perturbation' ends, the post-deletion rhythm re-appears with a
      phase shift in respect to the pre-deletion rhythm."
      GATES: (i) RG-E phasic -> sustained inside the window (no falling
      edge) and Mn-F deleted; (ii) the first post-pulse Mn-F burst is
      phase-SHIFTED (residual > 0.15 period vs the extrapolated grid);
      (iii) the rhythm re-establishes afterwards (post-pulse period
      within 15% of pre-pulse).

  Pulses (flag F-H): extra excitatory drive on the target population,
  started 0.1 period into an extensor burst (mid-extensor, the paper's
  phasic->sustained regime), lasting 2.0 periods (the paper's 1-3-cycle
  deletions; two flexor bursts omitted as in Fig 10Aa). Amplitudes are
  shell-scaled: PF-E +1.0 d-unit, RG-E +5.0 d-units (paper percentages
  quoted in the standalone docstring; the activity-based conductance
  regime differs from the spiking original).

Protocol per run: deterministic (no noise), dt = 0.05 ms, 10 s settling
skip (paper: 20 s settling, 0.1 ms exponential Euler). Grid residuals
are measured against each channel's own pre-pulse onset grid, as the
paper's Fig 6 dotted lines do ("plotted at intervals representing the
step cycle period preceding the deletion").

Verdict: RYBAK2006_TEMPLATE PASS iff EXP A (both signs) and EXP B
(B1 i-iii, B2 i-iii) all pass.
"""
from __future__ import annotations

import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np                     # noqa: E402

from rybak2006_standalone import (      # noqa: E402
    DRIVE_REGIMES, build_circuit, onsets, regime_stats, simulate)

TOL = 0.15          # phase residual tolerance, fraction of the period
SKIP = 10.0         # settling skip (s)
PULSE_DUR_T = 2.0   # deletion length in periods (paper range 1-3 cycles)


def grid_residual(t_onset, last_pre_onset, period):
    """Phase distance of t_onset from the grid last_pre_onset + k*period,
    in [0, 0.5]."""
    ph = ((t_onset - last_pre_onset) / period) % 1.0
    return min(ph, 1.0 - ph)


def exp_a(c):
    print("EXP A - drive asymmetry (Fig 2B/C)")
    checks = []
    stats = {}
    for name, te_claim in (("fig2b_E_dominant", "T_E > T_F"),
                           ("fig2c_F_dominant", "T_F > T_E")):
        dr = DRIVE_REGIMES[name]
        t, out, fout, _ = simulate(**dr, t_end_s=24.0, circuit=c)
        st = regime_stats(t, fout, c, skip_s=SKIP)
        stats[name] = st
        te, tf = st["RG-E"]["dur_mean"], st["RG-F"]["dur_mean"]
        ok_rhythm = (np.isfinite(st["period"])
                     and st["RG-E"]["n"] >= 5 and st["RG-F"]["n"] >= 5
                     and st["period_cv"] <= 0.05)
        if te_claim == "T_E > T_F":
            ok_sign = ok_rhythm and te > tf
        else:
            ok_sign = ok_rhythm and tf > te
        checks.append((f"A[{name}] rhythm present (period "
                       f"{st['period']:.3f} s, cv {st['period_cv']:.3f}, "
                       f"{st['RG-E']['n']}/{st['RG-F']['n']} RG eps)",
                       ok_rhythm))
        checks.append((f"A[{name}] {te_claim}: T_E {te:.3f} vs "
                       f"T_F {tf:.3f} s (MN-E {st['Mn-E']['dur_mean']:.3f} / "
                       f"MN-F {st['Mn-F']['dur_mean']:.3f})", ok_sign))
    # determinism / repeatability (standing directive)
    dr = DRIVE_REGIMES["fig2b_E_dominant"]
    _, out1, _, _ = simulate(**dr, t_end_s=24.0, circuit=c)
    _, out2, _, _ = simulate(**dr, t_end_s=24.0, circuit=c)
    checks.append(("A[repeat] identical traces on re-run (results repeat)",
                   bool(np.array_equal(out1, out2))))
    return checks


def deletion_test(c, regime_name, target, delta_d, label):
    """Baseline + pulsed run; deletion measurements. Returns (checks, info)."""
    dr = DRIVE_REGIMES[regime_name]
    t_end = 36.0
    t, out, fout, _ = simulate(**dr, t_end_s=t_end, circuit=c)
    st = regime_stats(t, fout, c, skip_s=SKIP)
    iRGE, iMnF, iMnE = c["idx"]["RG-E"], c["idx"]["Mn-F"], c["idx"]["Mn-E"]
    ons_rge = onsets(t, fout[:, iRGE])
    ons_mnf = onsets(t, fout[:, iMnF])
    T = st["period"]
    pre_rge = ons_rge[(ons_rge >= SKIP) & (ons_rge <= 20.0)]
    anchor = pre_rge[np.searchsorted(pre_rge, 18.0)]
    p_start = anchor + 0.1 * T
    p_end = p_start + PULSE_DUR_T * T
    pulse = dict(target=target, t_start_s=p_start,
                 dur_s=PULSE_DUR_T * T, delta_d=delta_d)
    t2, out2, fout2, _ = simulate(**dr, pulse=pulse, t_end_s=t_end,
                                  circuit=c)
    ons_rge2 = onsets(t2, fout2[:, iRGE])
    ons_mnf2 = onsets(t2, fout2[:, iMnF])

    # pre-pulse determinism: the pulse cannot act before its window
    npre = np.searchsorted(t, p_start - 1e-9)
    ident = bool(np.array_equal(out[:npre], out2[:npre]))

    last_pre_mnf = ons_mnf2[ons_mnf2 < p_start][-1]
    last_pre_rge = ons_rge2[ons_rge2 < p_start][-1]

    # deletion evidence inside the window
    base_in = ons_mnf[(ons_mnf >= p_start) & (ons_mnf <= p_end)]
    del_in = ons_mnf2[(ons_mnf2 >= p_start) & (ons_mnf2 <= p_end)]
    f_mne_win = fout2[(t2 >= p_start) & (t2 <= p_end), iMnE]
    f_rge_win = fout2[(t2 >= p_start) & (t2 <= p_end), iRGE]
    mne_sustained = f_mne_win.min() >= 0.05
    deleted = len(del_in) == 0 and len(base_in) >= 1

    # first post-deletion bursts vs the extrapolated pre-deletion grid
    post_mnf = ons_mnf2[ons_mnf2 >= p_end]
    res_mnf = (grid_residual(post_mnf[0], last_pre_mnf, T)
               if len(post_mnf) else float("nan"))
    # RG grid integrity: every post-skip RG-E onset (through + after pulse)
    post_rge2 = ons_rge2[ons_rge2 >= p_start]
    res_rge_all = ([grid_residual(x, last_pre_rge, T) for x in post_rge2]
                   if len(post_rge2) else [float("nan")])
    res_rge_max = max(res_rge_all) if res_rge_all else float("nan")
    # post-pulse rhythm re-establishment (skip the first transient cycle)
    late = ons_rge2[ons_rge2 >= p_end + 1.5 * T]
    per_post = (float(np.mean(np.diff(late))) if len(late) >= 3
                else float("nan"))

    info = dict(T=T, p_start=p_start, p_end=p_end,
                base_bursts_in_window=len(base_in),
                deleted_bursts_in_window=len(del_in),
                deleted=deleted, mne_sustained=mne_sustained,
                mne_min=float(f_mne_win.min()),
                rge_min_in_window=float(f_rge_win.min()),
                res_mnf=res_mnf, res_rge_max=res_rge_max,
                per_post=per_post, n_post_mnf=len(post_mnf),
                pre_pulse_identical=ident)
    print(f"  {label}: period {T:.3f} s; pulse [{p_start:.2f}, {p_end:.2f}] s "
          f"({PULSE_DUR_T:.1f} T); baseline flexor bursts in window "
          f"{len(base_in)} -> pulsed {len(del_in)}; Mn-E min f in window "
          f"{info['mne_min']:.2f}; RG-E min f in window "
          f"{info['rge_min_in_window']:.2f}; first post-deletion Mn-F "
          f"residual {res_mnf:.3f} T; RG grid max residual "
          f"{res_rge_max:.3f} T; post-pulse period "
          f"{per_post if np.isfinite(per_post) else float('nan'):.3f} s; "
          f"pre-pulse trace identical: {ident}")
    return info


def main() -> int:
    c = build_circuit()
    all_checks = []

    print("=" * 72)
    a_checks = exp_a(c)
    for nm, ok in a_checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {nm}")
    all_checks += a_checks

    print("=" * 72)
    print("EXP B1 - NON-RESETTING deletion, pulse at the PF level "
          "(Fig 10Ba regime, extra drive to PF-E)")
    i1 = deletion_test(c, "fig10ba_nonresetting", "PF-E", 1.0, "B1 PF-pulse")
    b1 = [
        (f"B1(i) flexor MN bursts deleted in window "
         f"({i1['base_bursts_in_window']} -> {i1['deleted_bursts_in_window']})"
         f" with extensor MN sustained (min f {i1['mne_min']:.2f})",
         i1["deleted"] and i1["mne_sustained"]),
        (f"B1(ii) RG-E rhythm grid maintained through/after the PF pulse "
         f"(max residual {i1['res_rge_max']:.3f} T <= {TOL})",
         i1["res_rge_max"] <= TOL),
        (f"B1(iii) first post-deletion flexor MN burst ON schedule "
         f"(residual {i1['res_mnf']:.3f} T <= {TOL})",
         np.isfinite(i1["res_mnf"]) and i1["res_mnf"] <= TOL),
        ("B1(iv) pre-pulse trace identical to baseline (determinism)",
         i1["pre_pulse_identical"]),
    ]
    for nm, ok in b1:
        print(f"  [{'PASS' if ok else 'FAIL'}] {nm}")
    all_checks += b1

    print("=" * 72)
    print("EXP B2 - RESETTING deletion, pulse at the RG level "
          "(Fig 6A regime, extra drive to RG-E)")
    i2 = deletion_test(c, "fig6_resetting", "RG-E", 5.0, "B2 RG-pulse")
    b2 = [
        (f"B2(i) RG-E phasic->sustained in window (min f "
         f"{i2['rge_min_in_window']:.2f} >= 0.25) and flexor MN deleted "
         f"({i2['base_bursts_in_window']} -> {i2['deleted_bursts_in_window']})",
         i2["rge_min_in_window"] >= 0.25 and i2["deleted"]),
        (f"B2(ii) first post-pulse flexor MN burst PHASE-SHIFTED "
         f"(residual {i2['res_mnf']:.3f} T > {TOL})",
         np.isfinite(i2["res_mnf"]) and i2["res_mnf"] > TOL),
        (f"B2(iii) rhythm re-establishes after the perturbation "
         f"(post period {i2['per_post']:.3f} s vs {i2['T']:.3f} s)",
         np.isfinite(i2["per_post"])
         and abs(i2["per_post"] - i2["T"]) <= TOL * i2["T"]),
        ("B2(iv) pre-pulse trace identical to baseline (determinism)",
         i2["pre_pulse_identical"]),
    ]
    for nm, ok in b2:
        print(f"  [{'PASS' if ok else 'FAIL'}] {nm}")
    all_checks += b2

    print("=" * 72)
    verdict = all(ok for _, ok in all_checks)
    print("RYBAK2006_TEMPLATE " + ("PASS" if verdict else "FAIL"))
    if not verdict:
        for nm, ok in all_checks:
            if not ok:
                print(f"  FAILED: {nm}")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
