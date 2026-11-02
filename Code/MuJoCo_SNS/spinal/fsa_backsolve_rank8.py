"""Rank-8 what-if variant of the FSA activation-to-circuit backsolve.

Companion to :mod:`fsa_backsolve`, which selected the shared NMF synergy
count 6 by the smallest centered-VAF>=0.90 criterion (right minimum 5,
left 6).  This script re-runs the identical pipeline with the synergy
count fixed at 8 for both legs and writes every output to
``fsa_results/rank8/``.  The recorded rank-6 outputs in ``fsa_results/``
are not modified.  Fixing the count at 8 is an override of the selection
criterion, not a new selection of record.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import fsa_backsolve as F

HERE = Path(__file__).resolve().parent
OUT = HERE / "fsa_results" / "rank8"
RECORDED = HERE / "fsa_results" / "fsa_backsolve.json"
RANK = 8

# Redirect every fsa_backsolve figure/json/npz write into rank8/.
F.OUT = OUT


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = np.load(HERE / "bsolve_out.npz", allow_pickle=True)
    time = np.asarray(data["t"], dtype=float)
    side_raw = {side: F.side_data(data, side) for side in ("r", "l")}

    curves = {}
    for side, (names, activation, _) in side_raw.items():
        target = np.clip(F.smooth_matrix(time, activation), 0.0, 1.0)
        curves[side], _ = F.rank_curve(target)

    results = {side: F.fit_side(time, names, activation, RANK)
               for side, (names, activation, _) in side_raw.items()}
    metrics = {side: F.metric_block(result)
               for side, result in results.items()}

    F.plot_rank(curves, RANK)
    phase_results = {side: F.plot_side(time, result, side)
                     for side, result in results.items()}

    recorded = json.loads(RECORDED.read_text(encoding="utf-8"))
    rec_metrics = recorded["metrics"]
    rec_curves = recorded["rank_selection"]["curves"]

    def curve_vaf(side, rank):
        row = next(c for c in rec_curves[side] if c["rank"] == rank)
        return row["vaf_centered"]

    payload = {
        "source": "bsolve_out.npz['acts'] (converted-MuJoCo ridge/NNLS)",
        "variant": "synergy count fixed at 8; override of the centered-VAF "
                   "selection criterion, not a selection of record",
        "time_s": [float(time[0]), float(time[-1])],
        "frames": int(len(time)),
        "rank_selection": {
            "criterion": "overridden; shared synergy count fixed at 8",
            "recorded_shared_synergy_count":
                recorded["rank_selection"]["shared_synergy_count"],
            "curves": curves,
        },
        "neuron_model": recorded["neuron_model"],
        "metrics": metrics,
        "sides": {},
    }
    arrays = {"time": time}
    for side, result in results.items():
        payload["sides"][side] = {
            "muscle_names": result["names"].tolist(),
            "pf_to_mn_dynamic_g_uS": result["g_fit"].tolist(),
            "mn_bias_nA": result["bias"].tolist(),
            "pf_upstream_exc_g_uS": result["pf"][1].tolist(),
            "pf_upstream_inh_g_uS": result["pf"][2].tolist(),
            "phase_normalization": {
                "convention": "heel strike=0, toe off=50, next heel strike=100",
                "n_complete_cycles": int(len(phase_results[side]["cycles"])),
                "events_s": phase_results[side]["events"],
                "measured_duty": phase_results[side]["duty"].tolist(),
            },
        }
        for key in ("target", "source", "gain", "nmf_pred", "g_analytic",
                    "analytic_pred", "g_fit", "bias", "fitted_pred"):
            arrays[f"{side}_{key}"] = result[key]
        arrays[f"{side}_names"] = result["names"]
        arrays[f"{side}_pf_src_exc"] = result["pf"][3]
        arrays[f"{side}_pf_src_inh"] = result["pf"][4]
        arrays[f"{side}_phase_grid"] = phase_results[side]["grid"]
        arrays[f"{side}_pf_phase_cycles"] = phase_results[side]["cycles"]
        arrays[f"{side}_pf_phase_mean"] = phase_results[side]["mean"]
        arrays[f"{side}_pf_phase_std"] = phase_results[side]["std"]

    (OUT / "fsa_backsolve.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")
    np.savez_compressed(OUT / "fsa_backsolve.npz", **arrays)

    lines = [
        "# FSA backsolve with 8 synergies (rank override)",
        "",
        f"What-if companion to `fsa_backsolve_report.md` (selected 6 by the "
        f"smallest centered-VAF>=90% criterion: right 5, left 6). This run "
        f"fixes the shared synergy count at **8 for both legs**. Source: "
        f"`{payload['source']}`; {len(time)} frames, "
        f"{time[0]:.3f}-{time[-1]:.3f} s, identical inputs to the rank-6 run.",
        "",
        "The recorded rank-6 outputs in `fsa_results/` were not modified; "
        "every file of this variant lives in `fsa_results/rank8/`.",
        "",
        "## Why 8 is not a selection of record",
        "",
        f"The centered-VAF criterion was already satisfied at rank 6 "
        f"(right {curve_vaf('r', 6):.4f}, left {curve_vaf('l', 6):.4f}). "
        f"At rank 8 the curve reaches "
        f"{curve_vaf('r', 8):.4f} (right) and {curve_vaf('l', 8):.4f} "
        f"(left), so components S7-S8 buy "
        f"+{1e2 * (curve_vaf('r', 8) - curve_vaf('r', 6)):.2f} / "
        f"+{1e2 * (curve_vaf('l', 8) - curve_vaf('l', 6)):.2f} percentage "
        f"points of centered VAF. Only one complete measured stride per "
        f"side exists, so the extra components are unconstrained between "
        f"cycles and carry no functional labels. NMF components also carry "
        f"scaling and permutation ambiguities; the S1-S8 ordering is not "
        f"identified beyond this run.",
        "",
    ]
    for side in ("r", "l"):
        m = metrics[side]
        lines += [
            f"## Side {side} (8 synergies)", "",
            f"- NMF centered VAF: {m['nmf_vaf_centered']:.4f}",
            f"- Dynamic FSA centered VAF: {m['dynamic_fsa_vaf_centered']:.4f}",
            f"- Dynamic FSA activation RMSE: "
            f"{m['dynamic_fsa_rmse_activation']:.4f}",
            f"- Analytical one-synapse gains outside Eexc/R limit: "
            f"{m['invalid_analytic_synapses']}",
            f"- Mean frames needing net negative MN current: "
            f"{m['mean_negative_required_fraction']:.3f}", "",
        ]
    lines += [
        "## Rank 6 (recorded) vs rank 8 (this run)", "",
        "| metric | r rank 6 | r rank 8 | l rank 6 | l rank 8 |",
        "| --- | --- | --- | --- | --- |",
    ]
    keys = [("nmf_vaf_centered", "NMF centered VAF", "{:.4f}"),
            ("dynamic_fsa_vaf_centered", "Dynamic FSA centered VAF", "{:.4f}"),
            ("dynamic_fsa_rmse_activation", "Dynamic FSA RMSE", "{:.4f}"),
            ("invalid_analytic_synapses", "Invalid analytic synapses", "{}"),
            ("mean_negative_required_fraction", "Negative-current frames", "{:.3f}")]
    for key, label, fmt in keys:
        cells = [fmt.format(rec_metrics["r"][key]), fmt.format(metrics["r"][key]),
                 fmt.format(rec_metrics["l"][key]), fmt.format(metrics["l"][key])]
        lines.append("| " + label + " | " + " | ".join(cells) + " |")
    lines += [
        "",
        "## Interpretation", "",
        "The upstream PF excitation/inhibition traces are exact inverse "
        "requirements for the PF leaky membranes. They are not evidence "
        "that the current two-half-center RG generates those waveforms. "
        "That forward closure must be tested separately. S7-S8 have no "
        "provisional functional labels, unlike S1-S6. The direct "
        "S-to-MN FSA fit remains an exploratory reduced mapping.",
        "",
        "References: Szczecinski, Hunt & Quinn (2017), DOI "
        "10.3389/fnbot.2017.00037; Szczecinski, Quinn & Hunt (2020), DOI "
        "10.3389/fnbot.2020.577804.",
    ]
    (OUT / "fsa_rank8_report.md").write_text(
        "\n".join(lines), encoding="utf-8")

    print(f"source frames={len(time)}, t={time[0]:.3f}-{time[-1]:.3f}s")
    for side in ("r", "l"):
        print(side, json.dumps(metrics[side], sort_keys=True))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
