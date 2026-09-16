"""Patch the figure block in _muscle_force_compare.py to matched activation."""
import io

p = "_muscle_force_compare.py"
t = io.open(p, encoding="utf-8").read()

old_fig = """    # --- figure: 6 representative muscles, F/Fmax_osim vs % cycle ---
    # Fmax_osim from the audit parse: reuse raw tendon-force normalization
    # (plot normalized to each curve's own peak for shape comparison)
    fig, axes = plt.subplots(3, 2, figsize=(9.5, 8.0), sharex=True)"""
new_fig = """    # --- figure: matched-activation comparison, normalized to own peak ---
    fig, axes = plt.subplots(3, 2, figsize=(9.5, 8.0), sharex=True)"""
assert old_fig in t, "fig header not found"
t = t.replace(old_fig, new_fig)

old_plot = """        ax.plot(pc, F_mj[:, j] / max(F_mj[:, j].max(), 1e-9), lw=1.8,
                color="#d55e00", label="MuJoCo (rigid tendon)")"""
new_plot = """        ax.plot(pc, F_mj_so[:, j] / max(F_mj_so[:, j].max(), 1e-9),
                lw=1.8, color="#d55e00", label="MuJoCo (rigid tendon)")"""
assert old_plot in t, "plot line not found"
t = t.replace(old_plot, new_plot)

old_title = 'ax.set_title(m, fontsize=9)'
new_title = ('ax.set_title(f"{m} (R2 {rows[j][1]:.2f}, '
             'lag {rows[j][3]}f)", fontsize=8)')
assert old_title in t, "title not found"
t = t.replace(old_title, new_title)

old_sup = ('    fig.suptitle("Isometric muscle force along subject01 walk "\n'
           '                 "(normalized to own peak; act = 1)", fontsize=10)')
new_sup = ('    fig.suptitle("Muscle force along subject01 walk at matched '
           'OpenSim SO\\nactivations (normalized to own peak)", fontsize=10)')
assert old_sup in t, "suptitle not found"
t = t.replace(old_sup, new_sup)

io.open(p, "w", encoding="utf-8").write(t)
print("figure block patched")
