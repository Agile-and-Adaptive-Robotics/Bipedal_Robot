import numpy as np
z = np.load("probe_w01.npz")
t = z["t"]; w = t >= 4.5
heel, toe = z["heelF"][w], z["toeF"][w]
for k, s in enumerate(("L", "R")):
    h5 = heel[:, k] > 5; t5 = toe[:, k] > 5
    any5 = h5 | t5
    both = heel[:, k] + toe[:, k]
    print(f"{s}: heel>5N {h5.mean():.3f} toe>5N {t5.mean():.3f} "
          f"combined>5N {any5.mean():.3f} totalF>20N {(both>20).mean():.3f} "
          f"mean total F {both.mean():.1f} N")
