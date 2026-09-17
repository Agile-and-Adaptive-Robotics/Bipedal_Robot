"""Backsolve MuJoCo gait activations through a nonspiking FSA circuit.

Input
-----
``bsolve_out.npz['acts']``: the converted-MuJoCo ridge/NNLS muscle
activation solution produced by :mod:`bsolve_ik` from subject01 IK + GRF.

Method
------
1. Select a per-leg muscle-synergy count by nonnegative matrix
   factorization (smallest rank with centered VAF >= 0.90; one PF per
   retained synergy).
2. Map activation to the implemented MN potential, V_MN = E_HI * a.
3. Compute the analytical FSA signal-transmission conductance

       g = k * R * Gm / (DeltaE - k * R)

   (Szczecinski et al. 2017, Eq. 18) for each PF->MN gain k.
4. Fit constant nonnegative PF->MN conductances directly to the complete
   LIF dynamics C*dV/dt = -Gm*V + I_bias + sum g*u*(Eexc-V).  This second
   result quantifies the correction required when multiple inputs overlap
   and the project's Eexc/R ratio is modest (8/5 rather than 194/20).
5. Invert each PF membrane to the excitatory/inhibitory upstream source
   fractions required to produce its target waveform.  These are RG/sensory
   requirements, not a claim that the current RG already produces them.

Outputs are written to ``fsa_results/``.  This script does not mutate the
running network parameters.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import lsq_linear
from scipy.signal import butter, filtfilt
from sklearn.decomposition import NMF

from muscle_map import classify
from params import E_HI, TAU
from build_network import E_REV_EXC, E_REV_INH
from gait_phase import phase_normalize


HERE = Path(__file__).parent
OUT = HERE / "fsa_results"
G_M = 1.0
VAF_TARGET = 0.90
MAX_RANK = 8
LOWPASS_HZ = 6.0
PF_LABELS = (
    "hip ext. / biartic. knee flx.",
    "knee ext. / ankle DF",
    "posterior chain / ankle PF",
    "hip flex. + frontal stabilization",
    "trunk-pelvis / ankle DF / hip abd.",
    "hip ext. / knee flx. residual",
)


def centered_vaf(target, pred):
    sse = float(np.sum((target - pred) ** 2))
    sst = float(np.sum((target - target.mean()) ** 2))
    return 1.0 - sse / max(sst, 1e-12)


def uncentered_vaf(target, pred):
    return 1.0 - float(np.sum((target - pred) ** 2)) / \
        max(float(np.sum(target ** 2)), 1e-12)


def smooth_matrix(time, x):
    dt = float(np.median(np.diff(time)))
    nyq = 0.5 / dt
    if LOWPASS_HZ >= 0.95 * nyq:
        return x.copy()
    b, a = butter(4, LOWPASS_HZ / nyq)
    return filtfilt(b, a, x, axis=0)


def nmf_model(x, rank):
    model = NMF(n_components=rank, init="nndsvda", max_iter=5000,
                tol=1e-7, random_state=42)
    coeff = model.fit_transform(np.maximum(x, 0.0))
    weight = model.components_
    return coeff, weight, coeff @ weight


def rank_curve(x):
    curve = []
    for rank in range(1, MAX_RANK + 1):
        _, _, pred = nmf_model(x, rank)
        curve.append({
            "rank": rank,
            "vaf_centered": centered_vaf(x, pred),
            "vaf_uncentered": uncentered_vaf(x, pred),
            "rmse": float(np.sqrt(np.mean((x - pred) ** 2))),
        })
    selected = next((row["rank"] for row in curve
                     if row["vaf_centered"] >= VAF_TARGET), MAX_RANK)
    return curve, selected


def normalize_components(coeff, weight):
    """Scale PF coefficients to [0,1], moving scale into spatial gains."""
    peak = np.maximum(coeff.max(axis=0), 1e-12)
    return coeff / peak, weight * peak[:, None]


def neuron_tau(name):
    info = classify(name)
    return TAU["mn"] * (1.0 + 0.5 * float(info.biarticular if info else 0))


def simulate_mn(time, source, conductance, bias, tau):
    """Forward Euler at the recorded sample times (stable for dt/tau<1)."""
    v = np.zeros((len(time), conductance.shape[1]))
    v[0] = np.clip(bias, 0.0, E_HI)
    for i in range(len(time) - 1):
        dt = time[i + 1] - time[i]
        syn = np.sum(source[i, :, None] * conductance *
                     (E_REV_EXC - v[i])[None, :], axis=0)
        dv = (-G_M * v[i] + bias + syn) / tau
        v[i + 1] = v[i] + dt * dv
    return v


def analytical_conductance(gain):
    denom = E_REV_EXC - gain * E_HI
    valid = denom > 1e-9
    g = np.full_like(gain, np.nan)
    g[valid] = gain[valid] * E_HI * G_M / denom[valid]
    return g, valid


def dynamic_conductance_fit(time, target_a, source, names):
    target_v = E_HI * target_a
    d_v = np.gradient(target_v, time, axis=0)
    conductance = np.zeros((source.shape[1], target_a.shape[1]))
    bias = np.zeros(target_a.shape[1])
    current_rmse = np.zeros(target_a.shape[1])
    negative_required = np.zeros(target_a.shape[1])
    for j, name in enumerate(names):
        c_m = neuron_tau(name) * G_M
        rhs = c_m * d_v[:, j] + G_M * target_v[:, j]
        negative_required[j] = np.mean(rhs < 0.0)
        design = np.column_stack([
            np.ones(len(time)),
            source * (E_REV_EXC - target_v[:, j, None]),
        ])
        # Light Tikhonov regularization resolves collinearity without
        # changing the nonnegative FSA sign constraint.
        scale = max(float(np.linalg.norm(design, ord="fro")), 1e-9)
        ridge = 1e-5 * scale
        aug_a = np.vstack([design, ridge * np.eye(design.shape[1])])
        aug_b = np.concatenate([rhs, np.zeros(design.shape[1])])
        sol = lsq_linear(aug_a, aug_b, bounds=(0.0, np.inf),
                         tol=1e-12, max_iter=1000)
        bias[j] = sol.x[0]
        conductance[:, j] = sol.x[1:]
        current_rmse[j] = np.sqrt(np.mean((design @ sol.x - rhs) ** 2))
    tau = np.array([neuron_tau(name) for name in names])
    pred_v = simulate_mn(time, source, conductance, bias, tau)
    return conductance, bias, pred_v / E_HI, current_rmse, negative_required


def pf_inverse(time, source):
    """Exact FSA source fractions for each PF target waveform."""
    v = E_HI * source
    dv = np.gradient(v, time, axis=0)
    required = TAU["pf"] * G_M * dv + G_M * v
    exc_req_g = np.maximum(required, 0.0) / \
        np.maximum(E_REV_EXC - v, 1e-9)
    inh_req_g = np.maximum(-required, 0.0) / \
        np.maximum(v - E_REV_INH, 1e-9)
    g_exc = np.maximum(exc_req_g.max(axis=0), 1e-12)
    g_inh = np.maximum(inh_req_g.max(axis=0), 1e-12)
    src_exc = exc_req_g / g_exc
    src_inh = inh_req_g / g_inh
    current_recon = (src_exc * g_exc * (E_REV_EXC - v) +
                     src_inh * g_inh * (E_REV_INH - v))
    current_rmse = np.sqrt(np.mean((current_recon - required) ** 2, axis=0))
    return required, g_exc, g_inh, src_exc, src_inh, current_rmse


def side_data(data, side):
    names = np.asarray([str(x) for x in data["act_names"]])
    fmax = np.asarray(data["Fmax"], dtype=float)
    activation = np.asarray(data["acts"], dtype=float)
    mask = np.array([name.endswith("_" + side) for name in names]) & \
        (fmax > 5.0) & (activation.max(axis=0) > 0.05)
    return names[mask], activation[:, mask], mask


def fit_side(time, names, activation, rank):
    target = np.clip(smooth_matrix(time, activation), 0.0, 1.0)
    coeff, weight, nmf_pred = nmf_model(target, rank)
    source, gain = normalize_components(coeff, weight)
    g_analytic, valid = analytical_conductance(gain)
    # Invalid gains cannot be realized by one excitatory FSA synapse with
    # Eexc=8 mV. Keep them NaN for the audit and zero for simulation.
    g_for_sim = np.nan_to_num(g_analytic, nan=0.0, posinf=0.0)
    analytic_pred = simulate_mn(
        time, source, g_for_sim, np.zeros(len(names)),
        np.array([neuron_tau(name) for name in names])) / E_HI
    g_fit, bias, fitted_pred, i_rmse, neg_frac = dynamic_conductance_fit(
        time, target, source, names)
    pf = pf_inverse(time, source)
    return {
        "names": names, "target": target, "source": source, "gain": gain,
        "nmf_pred": nmf_pred, "g_analytic": g_analytic,
        "analytic_valid": valid, "analytic_pred": analytic_pred,
        "g_fit": g_fit, "bias": bias, "fitted_pred": fitted_pred,
        "current_rmse": i_rmse, "negative_required": neg_frac,
        "pf": pf,
    }


def metric_block(result):
    target = result["target"]
    return {
        "nmf_vaf_centered": centered_vaf(target, result["nmf_pred"]),
        "nmf_vaf_uncentered": uncentered_vaf(target, result["nmf_pred"]),
        "analytic_fsa_vaf_centered": centered_vaf(target, result["analytic_pred"]),
        "dynamic_fsa_vaf_centered": centered_vaf(target, result["fitted_pred"]),
        "dynamic_fsa_rmse_activation": float(np.sqrt(np.mean(
            (target - result["fitted_pred"]) ** 2))),
        "invalid_analytic_synapses": int(np.size(result["analytic_valid"]) -
                                         np.count_nonzero(result["analytic_valid"])),
        "mean_negative_required_fraction": float(np.mean(
            result["negative_required"])),
        "pf_inverse_current_rmse_max_nA": float(np.max(result["pf"][-1])),
    }


def save_figure(fig, stem):
    """Write review PNG plus vector dissertation PDF."""
    fig.savefig(OUT / f"{stem}.png", dpi=220)
    fig.savefig(OUT / f"{stem}.pdf")


def plot_rank(curves, chosen):
    fig, ax = plt.subplots(figsize=(5.8, 3.5))
    for side, curve in curves.items():
        ax.plot([x["rank"] for x in curve],
                [x["vaf_centered"] for x in curve], "o-", label=side)
    ax.axhline(VAF_TARGET, color="0.4", ls="--", lw=1, label="90% criterion")
    ax.axvline(chosen, color="tab:red", ls=":", lw=1.2,
               label=f"chosen synergy count = {chosen}")
    ax.set(xlabel="NMF synergy count", ylabel="centered VAF",
           ylim=(0, 1.02), xticks=range(1, MAX_RANK + 1))
    ax.grid(alpha=0.25); ax.legend(fontsize=8)
    fig.tight_layout(); save_figure(fig, "fsa_rank_selection")
    plt.close(fig)


def plot_side(time, result, side):
    rank = result["source"].shape[1]
    phase = phase_normalize(time, result["source"], side)
    phase_grid = phase["grid"]
    phase_mean = phase["mean"]
    phase_std = phase["std"]
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.2),
                             gridspec_kw={"height_ratios": [1.1, 1.7]})
    for k in range(rank):
        line, = axes[0].plot(phase_grid, phase_mean[:, k],
                             label=f"S{k+1}")
        if len(phase["cycles"]) > 1:
            axes[0].fill_between(
                phase_grid, phase_mean[:, k] - phase_std[:, k],
                phase_mean[:, k] + phase_std[:, k],
                color=line.get_color(), alpha=0.14, linewidth=0)
    axes[0].axvspan(0, 50, color="0.92", zorder=-5)
    axes[0].axvline(50, color="0.35", ls="--", lw=0.9)
    axes[0].text(25, 1.02, "stance", ha="center", va="bottom",
                 transform=axes[0].get_xaxis_transform(), fontsize=8)
    axes[0].text(75, 1.02, "swing", ha="center", va="bottom",
                 transform=axes[0].get_xaxis_transform(), fontsize=8)
    axes[0].set(
        xlabel="stance-rescaled gait phase (%)",
        xlim=(0, 100), ylabel="normalized PF voltage",
        title=f"side {side}: cycle mean NMF synergy coefficients "
              f"(n={len(phase['cycles'])} complete stride)",
    )
    axes[0].grid(alpha=0.25); axes[0].legend(ncol=min(rank, 4), fontsize=7)
    im = axes[1].imshow(result["gain"], aspect="auto", cmap="viridis")
    labels = [f"S{k+1}  {PF_LABELS[k] if k < len(PF_LABELS) else ''}"
              for k in range(rank)]
    axes[1].set(yticks=range(rank), yticklabels=labels,
                xticks=range(len(result["names"])),
                xticklabels=result["names"],
                title="NMF spatial gains (extracted synergies; "
                      "functional labels provisional)")
    axes[1].tick_params(axis="y", labelsize=6.2)
    axes[1].tick_params(axis="x", rotation=90, labelsize=5.5)
    fig.colorbar(im, ax=axes[1], label="activation gain")
    fig.tight_layout(); save_figure(fig, f"fsa_pf_synergies_{side}")
    plt.close(fig)

    # Six muscles spanning proximal/distal flexor and extensor functions.
    desired = [f"psoas_{side}", f"glut_max1_{side}", f"vas_med_{side}",
               f"semimem_{side}", f"soleus_{side}", f"tib_ant_{side}"]
    idx = [list(result["names"]).index(name) for name in desired
           if name in result["names"]]
    fig, axes = plt.subplots(len(idx), 1, figsize=(9, 1.55 * len(idx)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, j in zip(axes, idx):
        ax.plot(time, result["target"][:, j], color="black", lw=1.5,
                label="MN target")
        ax.plot(time, result["nmf_pred"][:, j], color="tab:blue", ls="--",
                label="NMF")
        ax.plot(time, result["fitted_pred"][:, j], color="tab:orange",
                label="dynamic FSA")
        ax.set_ylabel(str(result["names"][j]), fontsize=7)
        ax.grid(alpha=0.2)
    axes[0].legend(ncol=3, fontsize=7); axes[-1].set_xlabel("time (s)")
    fig.tight_layout(); save_figure(fig, f"fsa_mn_reconstruction_{side}")
    plt.close(fig)

    req, g_exc, g_inh, src_exc, src_inh, _ = result["pf"]
    fig, axes = plt.subplots(rank, 1, figsize=(9, 1.5 * rank), sharex=True)
    axes = np.atleast_1d(axes)
    for k, ax in enumerate(axes):
        ax.plot(time, src_exc[:, k], color="tab:green", label="upstream exc")
        ax.plot(time, src_inh[:, k], color="tab:red", label="upstream inh")
        ax.set_ylabel(f"S{k+1}", fontsize=7); ax.grid(alpha=0.2)
        ax.text(1.005, 0.5, f"g+={g_exc[k]:.2f}\ng-={g_inh[k]:.2f} uS",
                transform=ax.transAxes, va="center", fontsize=6)
    axes[0].legend(ncol=2, fontsize=7); axes[-1].set_xlabel("time (s)")
    fig.tight_layout(); save_figure(fig, f"fsa_pf_upstream_requirements_{side}")
    plt.close(fig)
    return phase


def main():
    OUT.mkdir(exist_ok=True)
    data = np.load(HERE / "bsolve_out.npz", allow_pickle=True)
    time = np.asarray(data["t"], dtype=float)
    side_raw = {side: side_data(data, side) for side in ("r", "l")}
    curves, initial_rank = {}, {}
    for side, (names, activation, _) in side_raw.items():
        target = np.clip(smooth_matrix(time, activation), 0.0, 1.0)
        curves[side], initial_rank[side] = rank_curve(target)
    rank = max(initial_rank.values())
    results = {side: fit_side(time, names, activation, rank)
               for side, (names, activation, _) in side_raw.items()}
    metrics = {side: metric_block(result) for side, result in results.items()}

    plot_rank(curves, rank)
    phase_results = {}
    for side, result in results.items():
        phase_results[side] = plot_side(time, result, side)

    payload = {
        "source": "bsolve_out.npz['acts'] (converted-MuJoCo ridge/NNLS)",
        "time_s": [float(time[0]), float(time[-1])],
        "frames": int(len(time)),
        "rank_selection": {"criterion": "smallest centered VAF >= 0.90",
                           "per_side_minimum": initial_rank,
                           "shared_synergy_count": rank,
                           "curves": curves},
        "neuron_model": {"R_mV": E_HI, "Gm_uS": G_M,
                         "Eexc_mV": E_REV_EXC, "Einh_mV": E_REV_INH,
                         "tau_pf_s": TAU["pf"], "tau_mn_s": TAU["mn"]},
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
        "# FSA activation-to-circuit backsolve",
        "",
        f"Source: `{payload['source']}`; {len(time)} frames, "
        f"{time[0]:.3f}-{time[-1]:.3f} s.",
        "",
        f"Selected **{rank} NMF synergies per leg**. "
        f"Per-side minimum ranks at centered VAF >= {VAF_TARGET:.0%}: "
        f"right {initial_rank['r']}, left {initial_rank['l']}.",
        "",
        "These S1--S6 labels are NMF components, not PF neurons and not the "
        "existing ZCode E1/E2/F1/F2 early/late stance/swing channels. The "
        "direct S-to-MN FSA fit is an exploratory reduced mapping, not a "
        "claim of one PF neuron or one PF layer per synergy. The cycle plot maps "
        "heel strike to 0%, measured toe-off to 50%, and the next heel "
        "strike to 100%. Only one complete measured stride per side is "
        "available in this recording, so no between-cycle variance can yet "
        "be estimated.",
        "",
        "The target mapping is `V_MN = 5 mV * activation`. Analytical "
        "conductances use Szczecinski et al. (2017) Eq. 18; the dynamic "
        "fit uses the implemented conductance-based LIF equation with "
        "nonnegative excitatory PF synapses and a nonnegative tonic bias.",
        "",
    ]
    for side in ("r", "l"):
        m = metrics[side]
        lines += [
            f"## Side {side}", "",
            f"- NMF centered VAF: {m['nmf_vaf_centered']:.4f}",
            f"- Dynamic FSA centered VAF: {m['dynamic_fsa_vaf_centered']:.4f}",
            f"- Dynamic FSA activation RMSE: {m['dynamic_fsa_rmse_activation']:.4f}",
            f"- Analytical one-synapse gains outside Eexc/R limit: "
            f"{m['invalid_analytic_synapses']}",
            f"- Mean frames needing net negative MN current: "
            f"{m['mean_negative_required_fraction']:.3f}", "",
        ]
    lines += [
        "## Interpretation", "",
        "The upstream PF excitation/inhibition traces are exact inverse "
        "requirements for the PF leaky membranes. They are not evidence "
        "that the current two-half-center RG generates those waveforms. "
        "That forward closure must be tested separately.", "",
        "References: Szczecinski, Hunt & Quinn (2017), DOI "
        "10.3389/fnbot.2017.00037; Szczecinski, Quinn & Hunt (2020), DOI "
        "10.3389/fnbot.2020.577804.",
    ]
    (OUT / "fsa_backsolve_report.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"source frames={len(time)}, t={time[0]:.3f}-{time[-1]:.3f}s")
    print(f"rank minima: {initial_rank}; shared synergy count={rank}")
    for side in ("r", "l"):
        print(side, json.dumps(metrics[side], sort_keys=True))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
