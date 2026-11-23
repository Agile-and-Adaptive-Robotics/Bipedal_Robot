"""Record each physical foot–ground contact once, separately from the legacy heel/toe controller channels."""
from pathlib import Path

import mujoco
import numpy as np


class FootGroundLoads:
    def __init__(self, model, nsteps):
        self.ground = model.geom_type == mujoco.mjtGeom.mjGEOM_PLANE
        self.side = np.full(model.ngeom, -1, dtype=int)
        self.region = {}
        for index, body in enumerate(model.geom_bodyid):
            name = model.body(int(body)).name or ''
            if name in ('calcn_r', 'toes_r', 'calcn_l', 'toes_l'):
                self.side[index] = 0 if name.endswith('_r') else 1
                self.region[index] = ('heel' if name.startswith('calcn_') else 'toe', name[-1])
        gravity = np.asarray(model.opt.gravity)
        self.up = -gravity / np.linalg.norm(gravity)
        self.normal = np.zeros((nsteps, 2))
        self.upward = np.zeros_like(self.normal)
        self.count = np.zeros((nsteps, 2), dtype=np.int32)

    def add(self, step, geom1, geom2, frame, wrench):
        if self.ground[geom1] and self.side[geom2] >= 0:
            foot, direction = geom2, 1.0
        elif self.ground[geom2] and self.side[geom1] >= 0:
            foot, direction = geom1, -1.0
        else:
            return
        side = self.side[foot]
        self.normal[step, side] += max(float(wrench[0]), 0.0)
        force_on_geom2 = np.asarray(frame).reshape(3, 3).T @ wrench[:3]
        self.upward[step, side] += direction * float(self.up @ force_on_geom2)
        self.count[step, side] += 1
        return self.region[foot]

    def save(self, destination, t):
        with Path(destination).open('wb') as out:
            np.savez_compressed(
                out,
                t=t,
                normal_N=self.normal,
                upward_N=self.upward,
                contact_count=self.count,
                sides=np.array(['right', 'left']),
            )
