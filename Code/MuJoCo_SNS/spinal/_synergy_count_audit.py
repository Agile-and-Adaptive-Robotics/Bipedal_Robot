"""Per-leg NMF rank audit for the PF-population count.

The original exploratory script factorized both legs together, allowing
left/right phase to consume components.  If one PF population is intended per
muscle synergy, the relevant first check is the rank curve for each leg.
"""
from pathlib import Path

import numpy as np
from sklearn.decomposition import NMF


OSIM_DIR = Path(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim"
                r"\Gait2392_Robotbody\ResultsBSolve")


def read_sto(path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    end = next(i for i, line in enumerate(lines) if line.strip() == "endheader")
    names = lines[end + 1].split()[1:]
    data = np.array([[float(x) for x in line.split()]
                     for line in lines[end + 2:] if line.split()])
    return data[:, 0], names, data[:, 1:]


def audit(label, x, names, time):
    print(f"\n{label}: {x.shape[0]} frames x {x.shape[1]} muscles")
    print(f"{'n':>3s} {'VAF':>8s} {'delta':>8s} {'RMSE':>9s} {'BIC-ish':>10s}")
    prev = 0.0
    t, m = x.shape
    four = None
    for n in range(1, 9):
        model = NMF(n_components=n, init="nndsvda", max_iter=2000,
                    random_state=42)
        w = model.fit_transform(x)
        xr = w @ model.components_
        mse = float(np.mean((x - xr) ** 2))
        vaf = 1.0 - mse / max(float(np.var(x)), 1e-12)
        k = n * (t + m)
        bic = t * m * np.log(max(mse, 1e-12)) + k * np.log(t * m)
        print(f"{n:3d} {vaf:8.3f} {vaf - prev:8.3f} "
              f"{np.sqrt(mse):9.4f} {bic:10.0f}")
        if n == 4:
            four = (w, model.components_)
        prev = vaf

    w, h = four
    print("  four-component composition/timing:")
    for component in range(4):
        order = np.argsort(h[component])[::-1]
        top = ", ".join(names[i] for i in order[:6])
        peak_index = int(np.argmax(w[:, component]))
        peak_phase = float(np.mod(time[peak_index], 1.23) / 1.23)
        print(f"    S{component + 1}: peak phase={peak_phase:.2f}; {top}")


def main():
    time, names, activation = read_sto(
        OSIM_DIR / "zz_bsolve_StaticOptimization_activation.sto")
    print(f"source: {time[0]:.2f}-{time[-1]:.2f} s ({len(time)} frames)")
    for side in ("r", "l"):
        mask = np.array([name.endswith(f"_{side}") for name in names])
        side_names = np.asarray(names)[mask]
        active = activation[:, mask]
        keep = active.max(axis=0) > 0.05
        active = active[:, keep]
        audit(f"side {side}", active, side_names[keep], time)


if __name__ == "__main__":
    main()
