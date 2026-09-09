"""
MuJoCo glue for the custom BPA muscle (bpa_muscle.BPAMuscle).

Pattern: each BPA is a site-based <tendon> with a zero-gain <general>
actuator. MuJoCo therefore applies no force itself but exposes the tendon
length/velocity and the exact moment arm (actuator_moment) every step.
The force is computed HERE from Ben's BPA model and applied to the joints
via data.qfrc_applied, so the actuated dynamics use Ben's force law
instead of MuJoCo's Hill-type muscle.

Usage
-----
    system = BPAMuscleSystem({"knee_flex_bpa": BPAMuscle(...), ...})
    system.attach(model)                     # resolve actuator ids once
    mujoco.set_mjcb_control(system.cb)       # forces applied inside mj_step
    ...
    system.set_activation("knee_flex_bpa", 0.7)   # from SNS neuron output
    mujoco.mj_step(model, data)

Note: the control callback owns data.qfrc_applied and zeroes it every step.
If you also need external applied generalized forces, add them inside the
same callback (or extend BPAMuscleSystem with an external-forces hook).
"""

from __future__ import annotations

import numpy as np


class BPAMuscleSystem:
    """Attaches BPAMuscle instances to MuJoCo actuators and applies forces."""

    def __init__(self, muscles):
        self.muscles = dict(muscles)      # actuator name -> BPAMuscle
        self.actuator_ids = {}
        self._transposed_moment = None    # mujoco version dependent

    def attach(self, model):
        """Resolve actuator ids and moment-matrix orientation once."""
        for name in self.muscles:
            self.actuator_ids[name] = model.actuator(name).id
        # data.actuator_moment is (nv, nu) in mujoco 2.x, (nu, nv) in 3.x;
        # resolve after first data exists, so defer to first cb call.
        self._transposed_moment = None

    # ----------------------------------------------------------------- state
    def set_activation(self, name, a):
        """Set BPA pressure from an SNS/neural activation in [0, 1]."""
        self.muscles[name].set_activation(a)

    def forces(self, data):
        """Current BPA forces (N), keyed like self.muscles."""
        return {name: self.muscles[name].force(data.actuator_length[aid],
                                               data.actuator_velocity[aid])
                for name, aid in self.actuator_ids.items()}

    # ------------------------------------------------------------------ force
    def cb(self, model, data):
        """mjcb_control callback: apply each BPA force along its tendon."""
        moment = data.actuator_moment
        if self._transposed_moment is None:
            self._transposed_moment = moment.shape[0] == model.nv

        # this callback owns qfrc_applied; clear last step's contribution
        data.qfrc_applied[:] = 0.0

        for name, aid in self.actuator_ids.items():
            mus = self.muscles[name]
            f = mus.force(data.actuator_length[aid],
                          data.actuator_velocity[aid])
            col = moment[:, aid] if self._transposed_moment else moment[aid, :]
            data.qfrc_applied[:] += col * f
