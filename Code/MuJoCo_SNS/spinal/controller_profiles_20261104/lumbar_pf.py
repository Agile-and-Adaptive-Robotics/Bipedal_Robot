"""Fit and run three antagonistic lumbar PF pairs on native SNS motoneurons."""

import numpy as np
import mujoco
from scipy.optimize import lsq_linear

AXES = ("lumbar_extension", "lumbar_bending", "lumbar_rotation")
LABELS = ("LUM_Z", "LUM_X", "LUM_Y")
ROUTES = tuple(f"{group}_{side}" for side in ("r", "l")
               for group in ("ercspn", "intobl", "extobl"))
MIRROR = np.r_[3:6, 0:3]
MIRROR_TASK = np.array([1., -1., -1.])
EXCITATORY = 8.
INHIBITORY = -5.
FULL_SCALE = 5.
PHASE_FILTERS_S = np.array([.035,.12,.30,.65])


def fit(source, coactivation=.20, torque_step=22.):
    """Identify static inverse torque primitives and fit their MN conductances."""
    matrix = np.asarray(source["torque_at_full_activation_Nm"], dtype=float)
    bias = np.asarray(source["required_bias_Nm"], dtype=float)
    if tuple(source["routes"]) != ROUTES or tuple(source["joint_names"]) != AXES:
        raise ValueError("Lumbar actuator or axis order differs")
    if np.linalg.matrix_rank(matrix) != 3:
        raise ValueError("Trunk routes lack three-dimensional torque authority")
    tied = .5 * (matrix + MIRROR_TASK[:, None] * matrix[:, MIRROR])
    target_bias = np.array([bias[0], 0., 0.])
    pairing = np.vstack((np.eye(3), np.eye(3)))
    ridge = 6.
    pair_bias = lsq_linear(np.vstack((tied @ pairing, ridge*np.eye(3))),
                           np.r_[target_bias, np.full(3, ridge*coactivation)],
                           bounds=(.025, .8), method="bvls").x
    base = pairing @ pair_bias
    directions = (np.eye(3)[:, None, :] * np.array([1., -1.])[None, :, None]).reshape(6, 3)
    targets = []
    for i, direction in enumerate(directions):
        if i in (3, 5):
            targets.append(targets[-1][MIRROR])
            continue
        activation = lsq_linear(np.vstack((tied, ridge*np.eye(6))),
                                np.r_[target_bias + torque_step*direction, ridge*base],
                                bounds=(.005, .85), method="bvls").x
        targets.append(activation)
    targets = np.array(targets)
    change = targets-base
    reversal = np.where(change >= 0, EXCITATORY, INHIBITORY)
    voltage = targets*FULL_SCALE
    conductance = (voltage-base*FULL_SCALE)/(reversal-voltage)
    training = np.random.default_rng(219).uniform(-.5,.5,(160,3))
    training = np.vstack((directions,training))
    source = np.stack((training.clip(0),(-training).clip(0)),axis=2).reshape(-1,6)
    target_activation = np.array([lsq_linear(
        np.vstack((tied,ridge*np.eye(6))),
        np.r_[target_bias+torque_step*command,ridge*base],
        bounds=(.005,.85),method="bvls").x for command in training])
    target_v = target_activation*FULL_SCALE
    weights = np.ones(len(training))
    weights[:6] = 4.
    direction_mirror = [0,1,3,2,5,4]
    for j in range(3):
        other = MIRROR[j]
        a = source*(reversal[:,j]-target_v[:,j,None])
        b = source*(reversal[:,other]-target_v[:,other,None])
        design = np.vstack((weights[:,None]*a,weights[:,None]*b[:,direction_mirror],
                            .1*np.eye(6)))
        rhs = np.r_[weights*(target_v[:,j]-FULL_SCALE*base[j]),
                    weights*(target_v[:,other]-FULL_SCALE*base[other]),.1*conductance[:,j]]
        fitted = lsq_linear(design,rhs,bounds=(0.,2.),method="bvls").x
        conductance[:,j] = fitted
        conductance[:,other] = fitted[direction_mirror]
    # Holdout samples combine partial commands about different axes. Native
    # conductances sum, so these samples need not follow a linear torque map.
    rng = np.random.default_rng(961)
    holdout = rng.uniform(-.45, .45, (96,3))
    input_fraction = np.stack((holdout.clip(0), (-holdout).clip(0)),axis=2).reshape(-1,6)
    pred_v = (base*FULL_SCALE + input_fraction @ (conductance*reversal))/(1 + input_fraction @ conductance)
    observed = (pred_v/FULL_SCALE-base) @ matrix.T
    endpoint_v = (base*FULL_SCALE+conductance*reversal)/(1+conductance)
    return dict(schema="lumbar_pf_fsa_v1", routes=ROUTES, axes=AXES,
                torque_step_Nm=torque_step, baseline_activation=base.tolist(),
                conductance_uS=conductance.tolist(), reversal_mV=reversal.tolist(),
                trained_torque_increment_Nm=((endpoint_v/FULL_SCALE-base) @ matrix.T).tolist(),
                desired_torque_increment_Nm=(torque_step*directions).tolist(),
                current_fit_training_samples=len(training),
                baseline_residual_Nm=(matrix@base-bias).tolist(),
                mixed_command_rmse_Nm=np.sqrt(np.mean((observed-torque_step*holdout)**2,axis=0)).tolist())


def attach(net, spec):
    """Reuse Hunt's lumbar pairs, or add the same signed topology to another family."""
    import build_network as BN

    if spec["schema"] != "lumbar_pf_fsa_v1" or tuple(spec["routes"]) != ROUTES:
        raise ValueError("Unexpected lumbar PF fit")
    aliases = dict(getattr(net,"_aliases",{}))
    for name in aliases:
        net.idx.pop(name,None)
    cells, relays, ports = [], [], []
    for label in LABELS:
        existing = getattr(net,"extra_pairs",{}).get(label+"_mid")
        pair = list(existing or (f"PF_{label}_{sign}_mid" for sign in ("P","N")))
        inhibition = [f"PF_IN_{label}_{sign}_mid" for sign in ("P","N")]
        if existing:
            inhibition = [name.replace("PF_FUNC_","PF_IN_FUNC_",1) for name in pair]
        else:
            for name in pair+inhibition:
                net._add(name,.045,net.net)
            for i in range(2):
                net.net.add_connection(BN._syn(5./(EXCITATORY-5.),True),pair[i],inhibition[i])
                net.net.add_connection(BN._syn(.35,False),inhibition[i],pair[1-i])
        cells.extend(pair)
        relays.extend(inhibition)
        for cell in pair:
            port = f"LUMBAR_COMMAND_{len(ports)}"
            net.net.add_input(cell)
            net.inputs.append(port)
            ports.append(port)
    names = [p["name"] for p in net.net.populations]
    existing_edges = {(names[c["source"]],names[c["destination"]]) for c in net.net.connections}
    weights = {}
    for i, (cell, relay) in enumerate(zip(cells,relays)):
        for j, route in enumerate(ROUTES):
            gain = spec["conductance_uS"][i][j]
            if gain <= 1e-12:
                continue
            excitatory = spec["reversal_mV"][i][j] > 0
            pre, post = (cell if excitatory else relay), "MN_"+route
            if (pre,post) not in existing_edges:
                net.net.add_connection(BN._syn(gain,excitatory),pre,post)
            weights[pre,post] = gain
    net._n_true = len(net.idx)
    net.idx.update(aliases)
    net.lumbar_spec = spec
    net.lumbar_cells = cells
    net.lumbar_relays = relays
    net.lumbar_ports = ports
    net.lumbar_weights = weights


class Runtime:
    def __init__(self, net, stage, stiffness, damping, load_compensation=False,
                 precharge=False, phase_teacher=False, dt=None, pose_recruitment=False,
                 feedback_scale=1.):
        self.net, self.stage = net, stage
        self.k, self.c = float(stiffness), float(damping)
        self.load_compensation = load_compensation
        self.precharge = precharge
        self.feedback_enabled = stage >= 2 or phase_teacher
        self.feedback_scale = feedback_scale
        self.pose_recruitment = pose_recruitment
        self.ports = [net.inputs.index(name) for name in net.lumbar_ports]
        self.limit = float(net.lumbar_spec["torque_step_Nm"])
        names = [p["name"] for p in net.net.populations]
        c = net.compiled
        self.motors = [net.idx["MN_"+route] for route in ROUTES]
        self.trunk_slots = np.flatnonzero(np.isin(c.post,self.motors))
        # Hunt's compiled layer retains its phase-bank slots. The experiment
        # isolates the new signed commands by masking inherited trunk drive.
        self.bank_slots = np.array([i for i,(pre,post) in enumerate(zip(c.pre,c.post))
                                    if names[post] in net.lumbar_cells
                                    and names[pre].startswith("HUNT_BANK")],dtype=int)
        self.weights = np.array([net.lumbar_weights.get((names[c.pre[i]],names[c.post[i]]),0.)
                                 for i in self.trunk_slots])
        self.bias = FULL_SCALE*np.asarray(net.lumbar_spec["baseline_activation"])
        if pose_recruitment:
            self.moment_gain = np.empty((3, len(ROUTES)))
            self.conductance = np.asarray(net.lumbar_spec["conductance_uS"])
            self.driving = self.conductance*np.asarray(net.lumbar_spec["reversal_mV"])
        self.positive = np.array([net.idx[name] for name in net.lumbar_cells])
        self.relay = np.array([net.idx[name] for name in net.lumbar_relays])
        self.relay_slots = np.array([np.flatnonzero((c.pre == pf) & (c.post == relay)).item()
                                     for pf, relay in zip(self.positive,self.relay)])
        self.goal = None
        self.command_peak = np.zeros(6)
        self.history = []
        self.rg_ids = [net.idx[f"RG_{phase}_{side}"] for side in ("r","l") for phase in ("E","F")]
        phase = net.lumbar_spec.get("phase_readout")
        self.phase_weights = None if phase is None else np.asarray(phase["weights"])
        self.phase_state = np.zeros((len(PHASE_FILTERS_S),4))
        dt = net.command_runtime.dt if dt is None else dt
        self.phase_blend = 1.-np.exp(-dt/PHASE_FILTERS_S)
        template = net.lumbar_spec.get("phase_template")
        self.phase_tables = None
        if template is not None:
            self.phase_tables = np.asarray(template["commands"],dtype=float)
            self.phase_durations = np.asarray(template["duration_s"],dtype=float)
            if (self.phase_weights is not None or self.phase_tables.ndim != 3
                    or self.phase_tables.shape[1] < 2 or self.phase_tables.shape[2] != 3
                    or len(self.phase_tables) != len(self.phase_durations)
                    or not np.isfinite(self.phase_tables).all()
                    or not np.isfinite(self.phase_durations).all()
                    or np.any(self.phase_durations <= 0.)):
                raise ValueError("Invalid lumbar phase template")
            self.phase_threshold = float(template["right_flexor_threshold_mV"])
            if not np.isfinite(self.phase_threshold):
                raise ValueError("Invalid lumbar phase threshold")
            self.phase_segment = 0
            self.phase_age = 0.
            self.phase_above = False
            self.phase_dt = dt
        tau = net.lumbar_spec.get("tau_pf_s")
        if tau is not None:
            c.c.time_factor_membrane[np.r_[self.positive,self.relay]] = dt/float(tau)

    def apply_weights(self):
        c = self.net.compiled
        c.gmax[self.trunk_slots] = self.weights
        c.gmax[self.bank_slots] = 0.
        c.c.i_b[self.motors] = self.bias
        # The inhibitory relay reaches the PF full-scale voltage at 5 mV.
        # Some inherited Hunt relays have different conductances.
        c.gmax[self.relay_slots] = 5./(EXCITATORY-5.)
        self.net.sync()

    def update(self, u, model, data):
        if self.goal is None:
            self.ball = mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,"lumbar_ball")
            if self.ball >= 0:
                start = model.jnt_qposadr[self.ball]
                self.position = np.arange(start,start+4)
                self.ball_order = np.array([2,0,1])
                self.velocity = model.jnt_dofadr[self.ball]+self.ball_order
                self.rotation_error = np.zeros(3)
            else:
                self.position = np.array([model.joint(name).qposadr[0] for name in AXES])
                self.velocity = np.array([model.joint(name).dofadr[0] for name in AXES])
            self.goal = data.qpos[self.position].copy()
            self.actuators = [model.actuator(route).id for route in ROUTES]
            self.torso = model.body("lumbar_torso" if self.ball >= 0 else "torso").id
            self.external = np.zeros(model.nv)
            if self.precharge:
                data.act[model.actuator_actadr[self.actuators]] = self.bias/FULL_SCALE
                self.net.compiled.V[self.motors] = self.bias
        command = np.zeros(3)
        if self.ball >= 0:
            mujoco.mju_subQuat(self.rotation_error,data.qpos[self.position],self.goal)
            error = self.rotation_error[self.ball_order]
        else:
            error = data.qpos[self.position]-self.goal
        if self.feedback_enabled:
            torque = -self.k*error-self.c*data.qvel[self.velocity]
            if self.load_compensation:
                torque += self.net_load(model,data)
            torque *= self.feedback_scale
            command = self.recruit(torque) if self.pose_recruitment else torque/self.limit
        rg = self.net.compiled.V[self.rg_ids]
        if self.phase_weights is not None:
            self.phase_state += self.phase_blend[:,None]*(np.clip(rg/FULL_SCALE,0.,1.)-self.phase_state)
            command += np.r_[1.,self.phase_state.ravel()] @ self.phase_weights
        if self.phase_tables is not None:
            command += self.phase_command(rg[1])
        command = np.clip(command,-1.,1.)
        signed = np.stack((command.clip(0),(-command).clip(0)),axis=1).reshape(-1)
        u[self.ports] = FULL_SCALE*signed
        self.command_peak = np.maximum(self.command_peak,signed)
        self.history.append(np.r_[data.time,command,rg,error])

    def phase_command(self, right_flexor):
        above = right_flexor >= self.phase_threshold
        if above and not self.phase_above:
            self.phase_segment = min(self.phase_segment+1,len(self.phase_tables)-1)
            self.phase_age = 0.
        self.phase_above = above
        table = self.phase_tables[self.phase_segment]
        index = min(len(table)-1,self.phase_age/self.phase_durations[self.phase_segment]*(len(table)-1))
        lower = min(int(index),len(table)-2)
        fraction = index-lower
        self.phase_age += self.phase_dt
        return (1.-fraction)*table[lower]+fraction*table[lower+1]

    def recruit(self, torque):
        """Invert the fitted PF-to-torque response at the current pose."""
        command = np.clip(torque/self.limit,-1.,1.)
        best, best_error = command, np.inf
        regularization = (0.08*self.limit)**2*np.eye(3)
        for _ in range(5):
            signs = np.where(command < 0.,-1.,1.)
            rows = 2*np.arange(3)+(signs < 0.)
            g, e = self.conductance[rows], self.driving[rows]
            fraction = np.abs(command)
            denominator = 1.+fraction@g
            voltage = (self.bias+fraction@e)/denominator
            residual = torque-self.moment_gain@(voltage-self.bias)/FULL_SCALE
            error = residual@residual
            if error < best_error:
                best, best_error = command, error
            jacobian = (self.moment_gain@((e-g*voltage)/denominator).T)*signs/FULL_SCALE
            step = np.linalg.solve(jacobian.T@jacobian+regularization,jacobian.T@residual)
            command = np.clip(command+np.clip(step,-.35,.35),-1.,1.)
        return best

    def net_load(self, model, data):
        self.external[:] = 0.
        force, moment = np.split(data.xfrc_applied[self.torso],2)
        mujoco.mj_applyFT(model,data,force,moment,data.xipos[self.torso],self.torso,self.external)
        torque = (data.qfrc_bias-data.qfrc_passive-data.qfrc_applied-self.external)[self.velocity]
        for j,(activation,i) in enumerate(zip(self.bias/FULL_SCALE,self.actuators)):
            if model.actuator_gaintype[i] == mujoco.mjtGain.mjGAIN_FIXED:
                gain, passive = model.actuator_gainprm[i,0], 0.
            else:
                length, speed = data.actuator_length[i], data.actuator_velocity[i]
                lengthrange, acc = model.actuator_lengthrange[i], model.actuator_acc0[i]
                passive = mujoco.mju_muscleBias(length,lengthrange,acc,model.actuator_biasprm[i,:9])
                gain = mujoco.mju_muscleGain(length,speed,lengthrange,acc,model.actuator_gainprm[i,:9])
            moment_arm = data.actuator_moment[i,self.velocity]
            torque -= moment_arm*(passive+activation*gain)
            if self.pose_recruitment:
                self.moment_gain[:,j] = moment_arm*gain
        return torque
