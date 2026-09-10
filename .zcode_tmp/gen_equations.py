import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"D:\GitHub\Bipedal_Robot\.zcode_tmp\deck_assets\eqs"
os.makedirs(OUT, exist_ok=True)
C = "#2B2B2B"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.family"] = "STIXGeneral"

EQS = {
    "strain":  r"$\varepsilon \,=\, \dfrac{l_0 - l_m}{l_0} \qquad \varepsilon_{620} \,=\, \dfrac{l_0 - l_{620}}{l_0} \qquad \varepsilon^* \,=\, \dfrac{\varepsilon}{\varepsilon_{620}}$",
    "fstar":   r"$F^*(\varepsilon^*, P^*) \,=\, c_0\left(e^{-c_1\varepsilon^*} - 1\right) \,+\, P^*\, e^{-c_2(\varepsilon^*)^2} \qquad (F^* \leq 0 \Rightarrow 0)$",
    "fmax":    r"$F_{620_{10}} \,=\, 303.5\,\mathrm{N} \cdot \arctan\!\left(19.03\,\mathrm{m}^{-1}(l_0 - 7.5\,\mathrm{mm})\right)$",
    "fmax20":  r"$F_{620_{20}} \,=\, 922.4\,\mathrm{N} \cdot \arctan\!\left(15.37\,\mathrm{m}^{-1}(l_0 - 13\,\mathrm{mm})\right)$",
    "fmax3d":  r"$F_{max}(l_0, P) \,=\, a_1 \, P \, \arctan\!\left(a_2 \, P \, (l_0 - 7.5\,\mathrm{mm})\right)$",
    "force":   r"$F(\varepsilon^*, P^*, l_0) \,=\, F^*(\varepsilon^*, P^*) \cdot F_{620}(l_0)$",
    "lm":      r"$l_m \,=\, l_{LMT} \,-\, 2\,l_{fitting} \,-\, \chi_0 \;[-\, \Delta l]$",
    "kbr":     r"$\mathbf{K}_{br} \,=\, \mathrm{diag}(\chi_1, \chi_2, \chi_1) \qquad k_{\hat{u}} \,=\, \left(\hat{u}^{\,\mathsf{T}} \mathbf{K}_{br}^{-1} \hat{u}\right)^{-1}$",
    "series":  r"$k_{eq} \,=\, \left(k_{\hat{u}}^{-1} + k_{ten}^{-1}\right)^{-1} \qquad F_{620}(l_0)\, F^*\!\left(\varepsilon^*(\delta), P^*\right) - k_{eq}\,\delta \,=\, 0$",
    "wrap":    r"$\Delta l \,=\, \chi_3\, R\, \left|\theta_{wrap} - \theta_k\right| \zeta^2$",
    "torque":  r"$M \,=\, r_k \times F \qquad \vec{r} \,=\, \vec{p} - \hat{u}\left(\hat{u} \cdot \vec{p}\right)$",
    "xi_table": r"$\chi_0\ \mathrm{[mm]} \quad \chi_1, \chi_2\ \mathrm{[N/m]} \quad \chi_3\ \mathrm{[-]}$",
}

for name, tex in EQS.items():
    fig = plt.figure(figsize=(0.1, 0.1))
    t = fig.text(0, 0, tex, fontsize=30, color=C)
    fig.canvas.draw()
    bb = t.get_window_extent()
    fig.set_size_inches(bb.width / fig.dpi + 0.08, bb.height / fig.dpi + 0.08)
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("ok", name)
