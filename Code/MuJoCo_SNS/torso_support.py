"""Torso-mounted support for early curriculum tests of the trunk load path."""

import numpy as np
import mujoco


TRUNK_ROUTES = tuple(f"{name}_{side}" for side in ("r", "l")
                     for name in ("ercspn", "intobl", "extobl"))
LUMBAR_JOINTS = ("lumbar_extension", "lumbar_bending", "lumbar_rotation", "lumbar_ball")
PELVIS_JOINTS = ("pelvis_tx", "pelvis_ty", "pelvis_tz", "pelvis_tilt", "pelvis_list",
                 "pelvis_rotation", "pelvis_free", "root")


def _id(model, kind, name):
    return mujoco.mj_name2id(model, kind, name)


def torso_roster(model):
    """Report anatomy and route support before allowing a trunk curriculum gate."""
    body_kind = mujoco.mjtObj.mjOBJ_BODY
    joint_kind = mujoco.mjtObj.mjOBJ_JOINT
    actuator_kind = mujoco.mjtObj.mjOBJ_ACTUATOR
    torso = next((name for name in ("torso", "lumbar_torso")
                  if _id(model, body_kind, name) >= 0), None)
    lumbar = [name for name in LUMBAR_JOINTS if _id(model, joint_kind, name) >= 0]
    available = [mujoco.mj_id2name(model, actuator_kind, i) for i in range(model.nu)]
    routes = {name: next((a for a in (name, name + "_ref", name + "_bpa")
                          if a in available), None) for name in TRUNK_ROUTES}
    missing = [name for name, channel in routes.items() if channel is None]
    return dict(torso=torso, lumbar_joints=lumbar, routes=routes,
                missing_routes=missing, eligible=bool(torso and lumbar and not missing))


def restore_joint_mechanics(model, source, names):
    """Remove test-rig springs/dampers while retaining the source model mechanics."""
    kind = mujoco.mjtObj.mjOBJ_JOINT
    for name in names:
        joint = _id(model, kind, name)
        if joint < 0:
            continue
        original = _id(source, kind, name)
        if original < 0 or source.jnt_type[original] != model.jnt_type[joint]:
            raise ValueError(f"Source joint does not match {name}")
        if any(active and typ == mujoco.mjtEq.mjEQ_JOINT and joint in (a, b)
               for active, typ, a, b in zip(model.eq_active, model.eq_type,
                                            model.eq_obj1id, model.eq_obj2id)):
            raise ValueError(f"Release the equality restraint on {name} before torso training")
        model.jnt_stiffness[joint] = source.jnt_stiffness[original]
        dof_count = {mujoco.mjtJoint.mjJNT_FREE: 6, mujoco.mjtJoint.mjJNT_BALL: 3}.get(
            model.jnt_type[joint], 1)
        dst, src = model.jnt_dofadr[joint], source.jnt_dofadr[original]
        model.dof_damping[dst:dst+dof_count] = source.dof_damping[src:src+dof_count]
        q_count = {mujoco.mjtJoint.mjJNT_FREE: 7, mujoco.mjtJoint.mjJNT_BALL: 4}.get(
            model.jnt_type[joint], 1)
        dst, src = model.jnt_qposadr[joint], source.jnt_qposadr[original]
        model.qpos_spring[dst:dst+q_count] = source.qpos_spring[src:src+q_count]


def trunk_gate(summary, excursion_limit=15.):
    if not 0 < excursion_limit <= 30:
        raise ValueError("Early-stage lumbar excursion gate must be within (0, 30] degrees")
    return bool(summary["trunk_exercised"]
                and summary["max_lumbar_excursion_deg"] <= excursion_limit)


class TorsoSupport:
    """Apply and measure a Cartesian harness wrench at the torso COM."""

    def __init__(self, model, data, *, stage, scale=1., fraction=1.,
                 lateral_stiffness=1600., vertical_stiffness=2500.,
                 angular_stiffness=120., linear_damping=150., angular_damping=25.):
        if stage not in (1, 2) or not 0 <= scale <= 1 or not 0 <= fraction <= 1:
            raise ValueError("Use curriculum stage 1 or 2 and support settings in [0, 1]")
        if min(lateral_stiffness, vertical_stiffness, angular_stiffness,
               linear_damping, angular_damping) < 0:
            raise ValueError("Harness gains must be nonnegative")
        self.roster = torso_roster(model)
        if not self.roster["eligible"]:
            raise ValueError(f"Cannot test all trunk routes: {self.roster}")
        self.stage, self.scale, self.fraction = stage, float(scale), float(fraction)
        self.body = model.body(self.roster["torso"]).id
        root = self.body
        while model.body_parentid[root]:
            root = model.body_parentid[root]
        g = np.linalg.norm(model.opt.gravity)
        if g <= 0:
            raise ValueError("Torso weight support requires a gravity direction")
        self.weight = float(model.body_subtreemass[root] * g)
        self.vertical_cap_bw = 1.1
        self.up = -model.opt.gravity / g
        self.position = data.xipos[self.body].copy()
        self.rotation = data.xmat[self.body].reshape(3, 3).copy()
        self.pelvis = next((model.body(name).id for name in ("pelvis", "pelvis_link")
                            if _id(model, mujoco.mjtObj.mjOBJ_BODY, name) >= 0), root)
        self.lumbar_ref = (data.xmat[self.pelvis].reshape(3, 3).T @ self.rotation)
        self.ids = [model.actuator(self.roster["routes"][route]).id for route in TRUNK_ROUTES]
        self.act_addresses = model.actuator_actadr[self.ids]
        self.k = self.scale * (lateral_stiffness * np.eye(3)
                               + (vertical_stiffness-lateral_stiffness)*np.outer(self.up,self.up))
        self.kr = self.scale * angular_stiffness
        self.c = np.sqrt(self.scale) * linear_damping
        self.cr = np.sqrt(self.scale) * angular_damping
        self.velocity = np.zeros(6)
        self.wrench = np.zeros(6)
        self.peak_tension = np.zeros(len(self.ids))
        self.peak_driven_tension = np.zeros(len(self.ids))
        self.max_lumbar_excursion = 0.
        self.sum_upward = 0.
        self.sum_horizontal_squared = 0.
        self.sum_moment_squared = 0.
        self.samples = 0

    def apply(self, model, data):
        mujoco.mj_objectVelocity(model, data, mujoco.mjtObj.mjOBJ_BODY, self.body,
                                 self.velocity, 0)
        rot = data.xmat[self.body].reshape(3, 3)
        rotation_error = .5*np.cross(rot.T, self.rotation.T).sum(axis=0)
        force = (self.k @ (self.position - data.xipos[self.body])
                 - self.c*self.velocity[3:] + self.up*self.fraction*self.weight)
        upward = force @ self.up
        force += self.up * (np.clip(upward, 0., self.vertical_cap_bw*self.weight) - upward)
        self.wrench[:3] = force
        self.wrench[3:] = self.kr*rotation_error - self.cr*self.velocity[:3]
        data.xfrc_applied[self.body] = self.wrench

    def record(self, data, tension=None):
        """For an external BPA bank, supply its six measured model tensions."""
        if tension is None:
            tension = np.abs(data.actuator_force[self.ids])
            driven = data.ctrl[self.ids] > .01
            stateful = self.act_addresses >= 0
            driven[stateful] = data.act[self.act_addresses[stateful]] > .01
        else:
            tension = np.asarray(tension)
            driven = tension > 1.
        if tension.shape != (len(self.ids),) or not np.all(np.isfinite(tension) & (tension >= 0)):
            raise ValueError("Expected finite tensile force for all six trunk routes")
        self.peak_tension = np.maximum(self.peak_tension, tension)
        self.peak_driven_tension = np.maximum(self.peak_driven_tension, np.where(driven, tension, 0.))
        current = (data.xmat[self.pelvis].reshape(3, 3).T
                   @ data.xmat[self.body].reshape(3, 3))
        cosine = .5*(np.trace(self.lumbar_ref.T @ current)-1.)
        angle = np.degrees(np.arccos(np.clip(cosine, -1., 1.)))
        self.max_lumbar_excursion = max(self.max_lumbar_excursion, float(angle))
        upward = float(self.wrench[:3] @ self.up)
        horizontal = self.wrench[:3] - upward*self.up
        self.sum_upward += upward
        self.sum_horizontal_squared += float(horizontal @ horizontal)
        self.sum_moment_squared += float(self.wrench[3:] @ self.wrench[3:])
        self.samples += 1

    def summary(self):
        return dict(attachment=self.roster["torso"], stage=self.stage, scale=self.scale,
                    body_weight_N=self.weight, nominal_upward_support_bw=self.fraction,
                    mean_upward_support_bw=self.sum_upward/(max(1,self.samples)*self.weight),
                    rms_horizontal_support_N=float(np.sqrt(self.sum_horizontal_squared/max(1,self.samples))),
                    rms_orientation_support_Nm=float(np.sqrt(self.sum_moment_squared/max(1,self.samples))),
                    peak_tension_N=dict(zip(TRUNK_ROUTES,self.peak_tension.tolist())),
                    peak_tension_while_driven_N=dict(zip(TRUNK_ROUTES,self.peak_driven_tension.tolist())),
                    max_lumbar_excursion_deg=self.max_lumbar_excursion, samples=self.samples,
                    trunk_exercised=bool(self.samples and np.min(self.peak_driven_tension) > 1.),
                    roster=self.roster)
