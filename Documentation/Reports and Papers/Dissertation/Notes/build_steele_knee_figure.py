"""Build the Steele biomimetic-knee background figure for the dissertation.

Panel A: crossed four-bar schematic adapted from A.G. Steele, "Biomimetic
Design and Construction of a Bipedal Walking Robot", PSU thesis 2018
(Fig. 52 linkage): link lengths [L1,L2,L3,L4] = [1.85, 0.92, 2.03, 1.57] in;
link 1 = PCL analog, link 2 = femoral end, link 3 = ACL analog,
link 4 = tibial head (ground). ICR = intersection of the crossed link lines.
Drawn at three flexion angles with the ICR path traced.

Panel B: ICR migration path in the sagittal plane (femur frame, relative to
full extension). Robot (Steele knee, robot-body model) from
buildKneeFlexorContext20mm.m (knee_x_Pam/knee_y_Pam). Human: true ICR
computed from the OpenSim Gait2392 tibiofemoral translation splines
(Yamaguchi & Zajac 1989 planar knee) - identical knots to KneeFunction.m.

Panel C: photo placeholder (test-stand photo inserted by B. Bolen).

Run: myoconv python. Outputs steeleknee.pdf/.png next to this script's
--outdir.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.interpolate import CubicSpline
import os, argparse

IN2MM = 25.4
L1, L2, L3, L4 = 1.85 * IN2MM, 0.92 * IN2MM, 2.03 * IN2MM, 1.57 * IN2MM  # mm

# ---------------- four-bar kinematics ----------------
def seg_cross(p1, p2, p3, p4):
    """True if segment p1p2 crosses segment p3p4 (proper interior crossing)."""
    def side(a, b, p):
        return np.sign((b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0]))
    s1, s2 = side(p1, p2, p3), side(p1, p2, p4)
    s3, s4 = side(p3, p4, p1), side(p3, p4, p2)
    return s1 * s2 < 0 and s3 * s4 < 0

def line_intersect(p1, d1, p2, d2):
    """Intersection of lines p1+t*d1 and p2+s*d2."""
    M = np.array([[d1[0], -d2[0]], [d1[1], -d2[1]]])
    if abs(np.linalg.det(M)) < 1e-12:
        return None
    t, s = np.linalg.solve(M, p2 - p1)
    return p1 + t * d1

# Ground pivots D=(0,0) (posterior, PCL attachment) and A=(L4,0) (anterior,
# ACL attachment). Link 1 (PCL analog) sweeps from D; C solves the closure
# |BC| = L2, |AC| = L3 with the crossed branch selected. The crossed family
# exists for th1 in [37, 105.5] deg and spans ~181 deg of femoral-end
# rotation (verified numerically; see fb_diag probes).
Ap = np.array([L4, 0.0]); Dp = np.zeros(2)

sel = []
for th1 in np.arange(37.0, 105.51, 0.1):
    t1 = np.radians(th1)
    B = Dp + L1 * np.array([np.cos(t1), np.sin(t1)])
    d = np.linalg.norm(B - Ap)
    a = (L3**2 - L2**2 + d**2) / (2 * d)
    h2 = L3**2 - a**2
    if h2 < 0:
        continue
    h = np.sqrt(h2)
    u = (B - Ap) / d
    P = Ap + a * u
    perp = np.array([-u[1], u[0]])
    for C in (P + h * perp, P - h * perp):
        if not seg_cross(Dp, B, Ap, C):
            continue
        icr = line_intersect(Dp, B - Dp, Ap, C - Ap)
        if icr is None:
            continue
        psi = np.degrees(np.arctan2(C[1] - B[1], C[0] - B[0]))
        sel.append((th1, psi, B.copy(), C.copy(), icr.copy()))
        break  # one crossed branch per th1

psi_u = np.degrees(np.unwrap(np.radians([s[1] for s in sel])))
# flexion measured from the pose where the femoral-end bar is horizontal
flex = psi_u - 180.0
keep = (flex >= -8.0) & (flex <= 106.0)
sel = [s for s, k in zip(sel, keep) if k]
flex = flex[keep]
print(f"four-bar: {len(sel)} crossed configs, flexion {flex.min():.1f}"
      f"..{flex.max():.1f} deg of femoral-end rotation")
assert len(sel) > 100, "crossed configuration family too small"

# ---------------- ICR data (dissertation model) ----------------
ka_r = np.array([0.17,0.09,0.03,0.00,-0.09,-0.17,-0.26,-0.52,-0.79,-1.05,-1.31,-1.57,-1.83,-2.09,-2.36,-2.62])
xr = np.array([23.30,22.22,21.55,21.09,19.91,18.70,17.48,13.82,10.44,7.60,5.52,4.35,4.16,5.01,7.04,10.47])/1000
yr = np.array([-416.65,-417.03,-417.19,-417.28,-417.41,-417.41,-417.30,-416.28,-414.36,-411.72,-408.62,-405.32,-402.08,-399.16,-396.85,-395.66])/1000
ka_hx = np.array([-2.0944,-1.74533,-1.39626,-1.0472,-0.698132,-0.349066,-0.174533,0.197344,0.337395,0.490178,1.52146,2.0944])
hx    = np.array([-0.0032,0.00179,0.00411,0.0041,0.00212,-0.001,-0.0031,-0.005227,-0.005435,-0.005574,-0.005435,-0.00525])
ka_hy = np.array([-2.0944,-1.22173,-0.523599,-0.349066,-0.174533,0.159149,2.0944])
hy    = np.array([-0.4226,-0.4082,-0.399,-0.3976,-0.3966,-0.395264,-0.396])

def cs(x, y):
    o = np.argsort(x); return CubicSpline(np.asarray(x)[o], np.asarray(y)[o])

fx_r, fy_r = cs(ka_r, xr), cs(ka_r, yr)
fx_h, fy_h = cs(ka_hx, hx), cs(ka_hy, hy)

thR = np.linspace(-2.62, 0.0, 200)                 # robot 0..150 deg flexion
rx, ry = fx_r(thR), fy_r(thR)
drx = (rx - fx_r(0.0)) * 1000
dry = (ry - fy_r(0.0)) * 1000
flexR = -np.degrees(thR)

thH = np.linspace(-2.0944, 0.0, 200)               # human 0..120 deg
px, py = fx_h(thH), fy_h(thH)
vx, vy = fx_h(thH, 1), fy_h(thH, 1)
cxh, cyh = px - vy, py + vx                        # true ICR, femur frame
dcx = (cxh - (fx_h(0.0) - fy_h(0.0, 1))) * 1000
dcy = (cyh - (fy_h(0.0) + fx_h(0.0, 1))) * 1000
flexH = -np.degrees(thH)
print(f"robot ICR migration @150deg: ({drx[0]:.1f}, {dry[0]:.1f}) mm; "
      f"human ICR @120deg: ({dcx[0]:.1f}, {dcy[0]:.1f}) mm")

# ---------------- figure ----------------
plt.rcParams.update({
    "font.family": "serif", "font.size": 8.0,
    "axes.linewidth": 0.7, "xtick.direction": "out", "ytick.direction": "out",
    "pdf.fonttype": 42,
})
fig = plt.figure(figsize=(7.0, 2.75))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 0.75],
                      wspace=0.34, left=0.055, right=0.985, top=0.93, bottom=0.15)

# ---- Panel A: four-bar schematic ----
axA = fig.add_subplot(gs[0])
CB = {"l1": "#0072B2", "l2": "#000000", "l3": "#D55E00", "l4": "#555555",
      "icr": "#CC79A7"}

def cfg_at(target_deg):
    i = int(np.argmin(np.abs(flex - target_deg)))
    return sel[i]

icr_path = np.array([c[4] for c in sel])
icr_lo = int(np.argmin(np.abs(flex - 2.0)))
icr_hi = int(np.argmin(np.abs(flex - 102.0)))
icr_path = icr_path[icr_lo:icr_hi + 1]
draw_cfgs = [cfg_at(f) for f in (2.0, 102.0)]
alphas = [0.95, 0.50]

def femur_stub(B, C):
    """shaft direction: normal to the femoral-end bar, pointing away from
    the ground pivots at the extension pose"""
    m = (B + C) / 2
    n = np.array([-(C - B)[1], (C - B)[0]])
    n = n / np.linalg.norm(n)
    if np.dot(n, m - np.array([L4 / 2, 0.0])) < 0:
        n = -n
    return m, m + 20.0 * n

for k, ((th1, psi, B, C, icr), al) in enumerate(zip(draw_cfgs, alphas)):
    # femur shaft stub (drawn first, behind the links)
    s0, s1 = femur_stub(B, C)
    axA.plot([s0[0], s1[0]], [s0[1], s1[1]], color=CB["l2"], lw=3.2,
             alpha=al * 0.85, zorder=3, solid_capstyle="round")
    # link 1: D->B (PCL analog)
    axA.plot([Dp[0], B[0]], [Dp[1], B[1]], color=CB["l1"], lw=2.6, alpha=al,
             zorder=4, solid_capstyle="round")
    # link 3: A->C (ACL analog)
    axA.plot([Ap[0], C[0]], [Ap[1], C[1]], color=CB["l3"], lw=2.6, alpha=al,
             zorder=4, solid_capstyle="round")
    # link 2: femoral end B->C
    axA.plot([B[0], C[0]], [B[1], C[1]], color=CB["l2"], lw=6.5, alpha=al,
             zorder=5, solid_capstyle="round")
    # pivots
    for P in (Dp, Ap, B, C):
        axA.plot(*P, marker="o", ms=3.4, mfc="white", mec="black", mew=0.8,
                 alpha=min(1.0, al + 0.15), zorder=7)
    if k == 0:
        axA.plot(*icr, marker="+", ms=9, mew=1.6, color=CB["icr"], zorder=8)
# ground link (tibial head) on top of everything, single draw
axA.plot([0, L4], [0, 0], color=CB["l4"], lw=5.0, solid_capstyle="butt",
         zorder=6)
for xb in np.linspace(-2, L4 + 2, 8):
    axA.plot([xb, xb - 3.0], [0, -3.0], color=CB["l4"], lw=0.8, zorder=6)
# ICR path
axA.plot(icr_path[:, 0], icr_path[:, 1], "--", color=CB["icr"], lw=1.3,
         zorder=6)
# labels (data-coordinate placement against measured geometry)
th1, psi, B, C, icr = draw_cfgs[0]
axA.annotate("link 1\n(PCL analog)", xy=(4.5, 6.2), xytext=(-46, -3),
             fontsize=7, color=CB["l1"], ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=CB["l1"], lw=0.7,
                             shrinkA=2, shrinkB=2))
axA.annotate("link 3\n(ACL analog)", xy=(27.0, 13.4), xytext=(43, 1),
             fontsize=7, color=CB["l3"], ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=CB["l3"], lw=0.7,
                             shrinkA=2, shrinkB=2))
axA.annotate("femoral end\n(link 2)", xy=(26.0, 38.5), xytext=(38, 46),
             fontsize=7, ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color="black", lw=0.7,
                             shrinkA=2, shrinkB=2))
axA.annotate("tibial head\n(link 4, ground)", xy=(12.0, -0.8), xytext=(-46, -17),
             fontsize=7, color=CB["l4"], ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=CB["l4"], lw=0.7,
                             shrinkA=2, shrinkB=2))
s0, s1 = femur_stub(*draw_cfgs[0][2:4])
axA.annotate("femur", xy=(s1[0] + 2, s1[1] + 1), fontsize=8, style="italic",
             ha="left", va="bottom")
axA.annotate("ICR", xy=(25.0, 19.5), fontsize=7.5,
             color=CB["icr"], fontweight="bold", ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=CB["icr"], lw=0.7,
                             shrinkA=2, shrinkB=3))
axA.annotate("ICR path", xy=(-3.0, 21.9), xytext=(-45, 11),
             fontsize=7.5, color=CB["icr"], ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=CB["icr"], lw=0.7,
                             shrinkA=2, shrinkB=2))
_ext, _deep = draw_cfgs[0], draw_cfgs[1]
axA.add_patch(FancyArrowPatch(femur_stub(*_ext[2:4])[1],
                              femur_stub(*_deep[2:4])[1],
                              connectionstyle="arc3,rad=-0.3",
                              arrowstyle="-|>", mutation_scale=9,
                              color="0.35", lw=0.9, zorder=2))
axA.annotate("flexion", xy=(-30.0, 26.0), fontsize=7.5, color="0.25",
             ha="right", va="top")
axA.set_xlim(-50, 72)
axA.set_ylim(-21, 68)
axA.set_title("(A)  Crossed four-bar knee linkage", fontsize=8.5, pad=4)
axA.set_aspect("equal")
axA.axis("off")

# ---- Panel B: ICR migration paths ----
axB = fig.add_subplot(gs[1])
axB.plot(drx, dry, "-", color=CB["l1"], lw=1.6, zorder=3,
         label="Steele knee (robot model)")
axB.plot(dcx, dcy, "--", color=CB["l3"], lw=1.6, zorder=3,
         label="Human knee (Gait2392)")
# angle markers every 30 deg
for f in range(30, 160, 30):
    if f <= 150:
        i = np.argmin(np.abs(flexR - f))
        axB.plot(drx[i], dry[i], "o", ms=3.4, mfc="white", mec=CB["l1"],
                 mew=1.0, zorder=4)
    if f <= 120:
        j = np.argmin(np.abs(flexH - f))
        axB.plot(dcx[j], dcy[j], "^", ms=3.8, mfc="white", mec=CB["l3"],
                 mew=1.0, zorder=4)
# direction arrows
i = np.argmin(np.abs(flexR - 60))
axB.annotate("", xy=(drx[i+3], dry[i+3]), xytext=(drx[i], dry[i]),
             arrowprops=dict(arrowstyle="-|>", color=CB["l1"], lw=1.0))
j = np.argmin(np.abs(flexH - 60))
axB.annotate("", xy=(dcx[j+3], dcy[j+3]), xytext=(dcx[j], dcy[j]),
             arrowprops=dict(arrowstyle="-|>", color=CB["l3"], lw=1.0))
axB.plot(0, 0, "k*", ms=7, zorder=5)
axB.annotate("full\nextension", xy=(0, 0), xytext=(6, -14),
             textcoords="offset points", fontsize=7)
axB.annotate("flexion\nincreases", xy=(-11, 14), fontsize=7, color="0.25",
             ha="left")
axB.set_xlabel("posterior migration  (mm)")
axB.set_ylabel("proximal migration  (mm)")
axB.set_title("(B)  ICR migration, sagittal plane", fontsize=8.5, pad=4)
axB.legend(fontsize=6.5, loc="lower right", frameon=False, handlelength=1.8)
axB.axhline(0, color="0.85", lw=0.6, zorder=1)
axB.axvline(0, color="0.85", lw=0.6, zorder=1)
axB.tick_params(labelsize=7)

# ---- Panel C: photo placeholder ----
axC = fig.add_subplot(gs[2])
axC.add_patch(plt.Rectangle((0.03, 0.03), 0.94, 0.94, fill=True,
              facecolor="0.96", edgecolor="0.45", lw=1.0, ls=(0, (4, 3))))
axC.text(0.5, 0.56, "test-stand photo", ha="center", va="center",
         fontsize=8, color="0.35")
axC.text(0.5, 0.42, "(insert)", ha="center", va="center", fontsize=7.5,
         color="0.45")
axC.set_title("(C)  Physical knee", fontsize=8.5, pad=4)
axC.set_xlim(0, 1); axC.set_ylim(0, 1)
axC.set_xticks([]); axC.set_yticks([])
for s in axC.spines.values():
    s.set_visible(False)

out = os.path.abspath(os.path.join(os.path.dirname(__file__),
                    "..", "..", "..", "..",
                    "Documentation", "Reports and Papers", "Dissertation"))
if not os.path.isdir(out):
    out = os.path.dirname(os.path.abspath(__file__))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(out, f"steeleknee.{ext}"), dpi=600 if ext == "png" else None)
print("wrote", out)
