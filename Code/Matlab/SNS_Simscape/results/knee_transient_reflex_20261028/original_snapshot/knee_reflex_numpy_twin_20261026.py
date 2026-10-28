# NumPy twin of KneeReflexDemo (cross-software comparison for the figure).
# Equations/constants transcribed EXACTLY from demos/sns_build_demo.m and
# sns_build_library.m (2026-09-22 architecture). Euler dt = 1 ms, 5 s.
import numpy as np
from scipy.io import savemat

dt = 1e-3
T = 5.0
n = int(T / dt)

# mechanics + drives (sns_build_demo.m PreLoadFcn)
I_knee, b_knee, K_knee, th0 = 0.06, 0.5, 0.5, 0.26
ROM, r_arm, Tload = 1.5708, 0.035, 0.5
desc_ext, desc_flex = 4.2, 2.5
Fmax_ext, Fmax_flex, epsMax, epsScale = 500.0, 450.0, 0.25, 0.15

# neurons: Vrest, Gm(uS), Cm(nF), Thr, Slope ; currents nA, V mV
SN = dict(Vrest=-52.0, Gm=0.4, Cm=2.0, Thr=-55.0, Slope=1.0)
MN = dict(Vrest=-52.0, Gm=0.5, Cm=2.5, Thr=-45.0, Slope=1.0)
# synapses onto MNs: (gmax uS, Esyn mV); Sat window ThrPre -45, SlopePre 0.5
SYN = dict(ThrPre=-45.0, SlopePre=0.5)
ia_ext_exc = (0.18, 0.0)
ia_ext_inh_on_flex = (0.15, -72.0)
ib_ext_inh = (0.20, -72.0)
ia_flex_exc = (0.18, 0.0)
ia_flex_inh_on_ext = (0.15, -72.0)
ib_flex_inh = (0.20, -72.0)

# afferents
Imax, Wl, Wv = 10.0, 6.0, 8.0
Kf_ext, Kf_flex = 0.02, 0.0222
tauAct = 50.0  # ms

V = {k: SN["Vrest"] for k in ("Ia_ext", "Ia_flex", "Ib_ext", "Ib_flex")}
Vm = {"ext": MN["Vrest"], "flex": MN["Vrest"]}
A = {"ext": 0.0, "flex": 0.0}
th, thd = th0, 0.0

t_arr = np.arange(n) * dt
th_arr = np.zeros(n)


def sat(x, lo, hi):
    return np.clip(x, lo, hi)


def syn_g(Vpre, gmax):
    return gmax * sat((Vpre - SYN["ThrPre"]) / SYN["SlopePre"], 0, 1)


def step(Vp, par, Iapp, gsyn):
    Cm, Gm = par["Cm"], par["Gm"]
    Isyn = sum(g * (E - Vp) for g, E in gsyn)
    dV = (Gm * (par["Vrest"] - Vp) + Iapp + Isyn) / Cm  # nA/nF = mV/ms
    return Vp + dV * dt * 1000.0  # dt in s -> ms


for i in range(n):
    thn = th / ROM
    thdn = thd / ROM
    eps_ext = epsScale * thn
    eps_flex = epsScale * (1 - thn)
    F_ext = sat(Fmax_ext * A["ext"] * max(0.0, epsMax - eps_ext), 0, Fmax_ext)
    F_flex = sat(Fmax_flex * A["flex"] * max(0.0, epsMax - eps_flex), 0, Fmax_flex)
    # afferent currents
    I_ia_ext = sat(Wl * thn + Wv * thdn, 0, Imax)         # ext stretch + vel
    I_ia_flex = sat(Wl * (1 - thn) + Wv * (-thdn), 0, Imax)
    I_ib_ext = sat(Kf_ext * F_ext, 0, Imax)
    I_ib_flex = sat(Kf_flex * F_flex, 0, Imax)
    # sensory neurons
    for k, Iin in (("Ia_ext", I_ia_ext), ("Ia_flex", I_ia_flex),
                   ("Ib_ext", I_ib_ext), ("Ib_flex", I_ib_flex)):
        V[k] = step(V[k], SN, Iin, [])
    # motoneurons: Iapp = descending drive; synapses from SN voltages
    gE = [syn_g(V["Ia_ext"], ia_ext_exc[0] * 0 + ia_ext_exc[0] * 0 + 0) ]
    gE = []
    gE.append((syn_g(V["Ia_ext"], ia_ext_exc[0]), ia_ext_exc[1]))
    gE.append((syn_g(V["Ib_ext"], ib_ext_inh[0]), ib_ext_inh[1]))
    gE.append((syn_g(V["Ia_flex"], ia_flex_inh_on_ext[0]), ia_flex_inh_on_ext[1]))
    Vm["ext"] = step(Vm["ext"], MN, desc_ext, gE)
    gF = []
    gF.append((syn_g(V["Ia_flex"], ia_flex_exc[0]), ia_flex_exc[1]))
    gF.append((syn_g(V["Ib_flex"], ib_flex_inh[0]), ib_flex_inh[1]))
    gF.append((syn_g(V["Ia_ext"], ia_ext_inh_on_flex[0]), ia_ext_inh_on_flex[1]))
    Vm["flex"] = step(Vm["flex"], MN, desc_flex, gF)
    # activation + mechanics
    for m in ("ext", "flex"):
        S = sat((Vm[m] - MN["Thr"]) / MN["Slope"], 0, 1)
        A[m] = sat(A[m] + (S - A[m]) * 1000.0 * dt / tauAct, 0, 1)
    T_net = r_arm * F_flex - r_arm * F_ext + Tload \
        - b_knee * thd - K_knee * (th - th0)
    thd += T_net / I_knee * dt
    th += thd * dt
    th_arr[i] = th

print("numpy twin: theta end %.3f rad = %.2f deg (Simulink settle ~43.5 deg)"
      % (th, np.degrees(th)))
savemat(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\knee_twin_ref.mat",
        {"t": t_arr, "theta": th_arr}, do_compression=True)
print("saved results/knee_twin_ref.mat")
