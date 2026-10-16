"""Generate per-model joint-angle overlay figures with the OpenSim reference
repeated to match the simulation duration. Replaces the current single-cycle
overlays in the slideshow folder."""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
SLIDESHOW = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\slideshow")
JOINTS = [("hip_flexion_r", "Hip"), ("knee_angle_r", "Knee"), ("ankle_angle_r", "Ankle")]
COLORS = {"sim": "#0072B2", "ref": "#E69F00"}  # Tol: blue, orange


def load_ref():
    sys.path.insert(0, str(SPINAL))
    from kine_ref import load_reference
    return load_reference()


def make_figure(label, npz_path, ref, model_suffix=""):
    d = np.load(npz_path, allow_pickle=True)
    if "q" not in d or "key_joints" not in d:
        print(f"  SKIP {label}: npz missing q/key_joints")
        return
    q = d["q"]
    t = d["t"] if "t" in d else np.arange(len(q)) * 0.002
    kj = [str(x) for x in d["key_joints"]]
    T_ref = ref.get("T_r", 1.23)  # reference cycle period

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    fig.suptitle(f"{label}{model_suffix}: joint angles vs OpenSim reference", fontsize=12)

    for ax, (col_name, jname) in zip(axes, JOINTS):
        if col_name not in kj:
            ax.set_ylabel(f"{jname}\n(missing)")
            continue
        col = kj.index(col_name)
        sim = q[:, col]
        sim_t = t - t[0]  # start at 0

        # Plot simulation
        ax.plot(sim_t, sim, color=COLORS["sim"], linewidth=1.2, label="SNS model")

        # Plot reference, repeated to cover sim duration
        ref_cycle = np.asarray(ref["r"][jname.lower()])
        n_ref = len(ref_cycle)
        ref_dt = T_ref / n_ref
        sim_duration = sim_t[-1]
        n_cycles = int(np.ceil(sim_duration / T_ref)) + 1
        ref_t_full = []
        ref_v_full = []
        for c in range(n_cycles):
            ref_t_full.append(np.arange(n_ref) * ref_dt + c * T_ref)
            ref_v_full.append(ref_cycle)
        ref_t = np.concatenate(ref_t_full)
        ref_v = np.concatenate(ref_v_full)
        mask = ref_t <= sim_duration
        ax.plot(ref_t[mask], ref_v[mask], color=COLORS["ref"], linewidth=1.0,
                alpha=0.7, linestyle="--", label="OpenSim reference")

        ax.set_ylabel(f"{jname} (deg)")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time (s)")
    plt.tight_layout()

    out = SLIDESHOW / f"{label}_joints.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out.name} ({n_cycles} ref cycles over {sim_duration:.1f}s)")


def main():
    ref = load_ref()
    models = [
        ("s3k_mutfix", "spinal_run_s3k_mutfix_best_20261003.npz"),
        ("syn6", "spinal_run_syn6_best_20261003.npz"),
        ("s3kpruned", "spinal_run_pruned_best_20261003.npz"),
        ("w2lvar", "spinal_run_w2lvar_s1_best_20261003.npz"),
        ("spiking_mirror_s1", "spinal_run_spk_s1_best_20261003.npz"),
    ]
    for label, npz_name in models:
        npz = SPINAL / npz_name
        if npz.exists():
            print(f"Generating {label}...")
            make_figure(label, npz, ref)
        else:
            print(f"SKIP {label}: {npz_name} not found")


if __name__ == "__main__":
    main()
