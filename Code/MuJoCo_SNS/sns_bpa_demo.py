"""
End-to-end demo: SNS-Toolbox neural network -> custom BPA muscles -> MuJoCo.

Builds a minimal 1-DOF knee (shank on a hinge) with two antagonist BPA
tendons (20 mm Festo BPAs, Ben's force law), drives them with a tiny
Synthetic Nervous System (command neuron -> antagonist motor neurons with
reciprocal inhibition), and simulates the coupled loop.

The network commands a flexion burst, then an extension burst; forces,
pressure (activation) and the joint angle are logged and plotted to
sns_bpa_demo.png.

Run inside the `myoconv` env:
    conda run --live-stream -n myoconv python sns_bpa_demo.py
"""

from pathlib import Path

import mujoco
import numpy as np

from bpa_muscle import BPAMuscle
from bpa_mujoco import BPAMuscleSystem

HERE = Path(__file__).parent

# --------------------------------------------------------------------- MJCF
MJCF = """
<mujoco model="bpa_knee_demo">
  <compiler angle="radian"/>
  <option timestep="0.001" gravity="0 0 -9.81"/>
  <worldbody>
    <!-- thigh welded to world; shank on a hinge (positive angle = flexion) -->
    <body name="thigh" pos="0 0 0.45">
      <geom type="box" size="0.02 0.02 0.20" pos="0 0 -0.10" rgba="0.6 0.6 0.6 1"/>
      <site name="flex_org" pos=" 0.05 -0.005 -0.10"/>
      <site name="ext_org"  pos="-0.05 -0.005 -0.10"/>
      <body name="shank" pos="0 0 -0.20">
        <joint name="knee" type="hinge" pos="0 0 0" axis="0 -1 0" limited="true"
               range="-2.2 0" damping="0.5" armature="0.01"/>
        <geom type="capsule" size="0.02" fromto="0 0 0 0 0 -0.20" rgba="0.8 0.4 0.4 1"/>
        <geom type="sphere" size="0.03" pos="0 0 -0.22" density="1200"/>
        <site name="flex_ins"  pos=" 0.02 -0.005 -0.10"/>
        <site name="ext_ins"   pos="-0.02 -0.005 -0.10"/>
      </body>
    </body>
  </worldbody>
  <tendon>
    <spatial name="knee_flex_tendon">
      <site site="flex_org"/>
      <site site="flex_ins"/>
    </spatial>
    <spatial name="knee_ext_tendon">
      <site site="ext_org"/>
      <site site="ext_ins"/>
    </spatial>
  </tendon>
  <actuator>
    <general name="knee_flex_bpa" tendon="knee_flex_tendon"
             gaintype="fixed" gainprm="0 0 0" biastype="none" biasprm="0 0 0"
             ctrlrange="0 1" ctrllimited="true"/>
    <general name="knee_ext_bpa" tendon="knee_ext_tendon"
             gaintype="fixed" gainprm="0 0 0" biastype="none" biasprm="0 0 0"
             ctrlrange="0 1" ctrllimited="true"/>
  </actuator>
</mujoco>
"""

# Typical 20 mm BPA parameters (Ben's MonoPam contexts use these orders)
# 20 mm BPA sized to the ~0.20 m straight-line routes (rest ~= route - 2*fit)
BPA = dict(diameter=20, resting_length=0.155, kmax_length=0.116,  # 25% contraction
           fitting_length=0.025, tendon_length=0.0, bpa_count=2,
           pressure_max_kpa=620.0, wraps=1)


def build_network(dt):
    """Command neuron -> two antagonist motor neurons with recip. inhibition."""
    from sns_toolbox.networks import Network
    from sns_toolbox.connections import NonSpikingSynapse
    from sns_toolbox.neurons import NonSpikingNeuron
    net = Network(name="BPA knee")
    # F = flexor MN, E = extensor MN (resting 0 mV, e_hi ~ 5 mV);
    # separate command inputs per MN + reciprocal inhibition for switching
    net.add_neuron(NonSpikingNeuron(membrane_capacitance=0.02,
                                    membrane_conductance=1.0), name="F")
    net.add_neuron(NonSpikingNeuron(membrane_capacitance=0.02,
                                    membrane_conductance=1.0), name="E")
    net.add_input("F")
    net.add_input("E")
    # whichever MN is driven harder suppresses the other
    net.add_connection(NonSpikingSynapse(max_conductance=3.0, reversal_potential=-5.0), "F", "E")
    net.add_connection(NonSpikingSynapse(max_conductance=3.0, reversal_potential=-5.0), "E", "F")
    return net.compile(dt=dt, backend="numpy")


def main(sim_time=3.0, dt=0.001):
    model = mujoco.MjModel.from_xml_string(MJCF)
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0) if model.nkey > 0 else None

    muscles = {
        "knee_flex_bpa": BPAMuscle(name="knee_flex_bpa", **BPA),
        "knee_ext_bpa": BPAMuscle(name="knee_ext_bpa", **BPA),
    }
    system = BPAMuscleSystem(muscles)
    system.attach(model)
    mujoco.set_mjcb_control(system.cb)

    net = build_network(dt)
    qadoi = model.jnt_qposadr[model.joint("knee").id]

    steps = int(sim_time / dt)
    log = {k: np.zeros(steps) for k in
           ("t", "q", "F_flex", "F_ext", "a_flex", "a_ext")}

    for k in range(steps):
        t = k * dt
        # command: rest (0-1 s) -> flex (1-2 s) -> extend (2-3 s)
        # command: rest (0-1 s) -> flex (1-2 s) -> extend (2-3 s)
        i_flex = 8.0 if 1.0 <= t < 2.0 else 0.0
        i_ext = 8.0 if t >= 2.0 else 0.0
        out = net.forward([i_flex, i_ext])    # SNS network step
        v = net.V                             # membrane potentials (mV)
        a_flex = float(np.clip(v[0] / 5.0, 0.0, 1.0))   # normalize e_hi=5 mV
        a_ext = float(np.clip(v[1] / 5.0, 0.0, 1.0))
        system.set_activation("knee_flex_bpa", a_flex)
        system.set_activation("knee_ext_bpa", a_ext)

        mujoco.mj_step(model, data)

        log["t"][k] = t
        log["q"][k] = data.qpos[qadoi]
        log["F_flex"][k], log["F_ext"][k] = (system.forces(data)["knee_flex_bpa"],
                                             system.forces(data)["knee_ext_bpa"])
        log["a_flex"][k], log["a_ext"][k] = a_flex, a_ext

    mujoco.set_mjcb_control(None)   # release the callback

    print(f"final knee angle: {np.degrees(log['q'][-1]):7.2f} deg "
          f"(0 = straight, negative = flexed)")
    print(f"peak flex force: {log['F_flex'].max():7.1f} N, "
          f"peak ext force: {log['F_ext'].max():7.1f} N")

    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
        ax[0].plot(log["t"], log["a_flex"], label="flex activation (-> kPa)")
        ax[0].plot(log["t"], log["a_ext"], label="ext activation")
        ax[0].set_ylabel("activation [0-1]")
        ax[0].legend(); ax[0].grid(True)
        ax[1].plot(log["t"], log["F_flex"], label="flex BPA force")
        ax[1].plot(log["t"], log["F_ext"], label="ext BPA force")
        ax[1].set_ylabel("force [N]"); ax[1].legend(); ax[1].grid(True)
        ax[2].plot(log["t"], np.degrees(log["q"]), "k")
        ax[2].set_ylabel("knee angle [deg]"); ax[2].set_xlabel("t [s]")
        ax[2].grid(True)
        fig.suptitle("SNS -> BPA muscles -> MuJoCo knee")
        fig.tight_layout()
        out = HERE / "sns_bpa_demo.png"
        fig.savefig(out, dpi=140)
        print(f"plot written to {out}")
    except Exception as e:      # plotting is optional
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
