"""Evaluate the force reported by OpenSim static optimization on MuJoCo routes.

StaticOptimizationTarget uses the unit-activation inextensible-tendon active
capacity at the current fiber length and speed, multiplied by recruitment.
This diagnostic does not replace the actuator's total force in forward dynamics.
Upstream sources: OpenSim 4.3 StaticOptimizationTarget::prepareToOptimize
and Thelen2003Muscle::calcInextensibleTendonActiveFiberForce.
https://github.com/opensim-org/opensim-core/blob/4.3/OpenSim/Analyses/StaticOptimizationTarget.cpp
https://github.com/opensim-org/opensim-core/blob/4.3/OpenSim/Actuators/Thelen2003Muscle.cpp
"""
from dataclasses import dataclass

import numpy as np
import opensim

from thelen_dynamic import force_velocity


@dataclass(frozen=True)
class StaticOptimizationForce:
    actuator_ids: np.ndarray
    parameters: np.ndarray

    @classmethod
    def from_models(cls, mujoco_model, opensim_model, names):
        """Read parameters from the model that was actually converted."""
        muscles = opensim_model.getMuscles()
        rows = []
        for name in names:
            muscle = opensim.Thelen2003Muscle.safeDownCast(muscles.get(name))
            if muscle is None or not np.isclose(muscle.get_fv_linear_extrap_threshold(), .95):
                raise ValueError(f"{name}: expected the Thelen force–velocity law with 0.95 extrapolation")
            optimal = muscle.getOptimalFiberLength()
            rows.append((
                muscle.getMaxIsometricForce(), optimal, muscle.getTendonSlackLength(),
                optimal * np.sin(muscle.getPennationAngleAtOptimalFiberLength()),
                muscle.getPennationModel().getMinimumFiberLength(),
                muscle.get_KshapeActive(), muscle.get_Af(), muscle.get_Flen(),
                muscle.getMaxContractionVelocity(),
            ))
        return cls(np.array([mujoco_model.actuator(n).id for n in names]), np.array(rows))

    def from_lengths(self, lengths, velocities, activation):
        """Return active tension (N) for muscle–tendon lengths (m) and speeds (m/s)."""
        lengths, velocities, activation = map(np.asarray, (lengths, velocities, activation))
        expected = (len(self.parameters),)
        if any(arr.shape != expected for arr in (lengths, velocities, activation)):
            raise ValueError("One length, velocity, and recruitment are required per muscle")
        if not np.all(np.isfinite([lengths, velocities, activation])) or np.any((activation < 0) | (activation > 1)):
            raise ValueError("Lengths and speeds must be finite, and recruitment must lie in [0, 1]")
        force = np.zeros_like(lengths, dtype=float)
        for i, (fmax, optimal, slack, thickness, minimum, shape, af, flen, vmax) in enumerate(self.parameters):
            along_tendon = lengths[i] - slack
            fiber = np.hypot(along_tendon, thickness)
            if along_tendon <= 0 or fiber <= minimum:
                continue
            cosine = along_tendon / fiber
            normalized_length = fiber / optimal
            fiber_speed = velocities[i] * cosine / optimal
            force_length = np.exp(-((normalized_length - 1.)**2) / shape)
            force_speed, _ = force_velocity(fiber_speed, 1., af, flen, vmax)
            force[i] = activation[i] * fmax * force_length * force_speed * cosine
        return force

    def evaluate(self, data, activation):
        """Read an already realized MuJoCo state, without mutating any force or state."""
        return self.from_lengths(data.actuator_length[self.actuator_ids],
                                 data.actuator_velocity[self.actuator_ids], activation)
