    def run(self, te=3.0, tf=4.0, wsup=0.7, dur=10.0,
            stimN=0.0, stimT=8.0, stimDur=1.0,
            hold=3.0, lower=1.0, sink=0.01, drop=0.026,
            cap=0.5, acap=-1, rig=0.25, settle=0.4) -> dict:
        """One trial with the full drop protocol, compiled neural step."""
        m, d = self.model, self.data
        mujoco = self.mujoco
        dt = self.dt
        n_steps = int(dur / dt)

        # reset
        mujoco.mj_resetDataKeyframe(m, d, 0)
        self.V[:] = 0.0
        self.b_gate[:] = 0.0
        self.c_gate[:] = 0.0

        # ---- drop-protocol vertical harness (from test_w2l_ground) ----
        jadr = m.jnt_qposadr[mujoco.mj_name2id(
            m, mujoco.mjtObj.mjOBJ_JOINT, "root")]
        z_air = float(d.qpos[jadr + 2])
        z_contact = z_air - drop - sink
        W = float(np.sum(m.body_mass)) * (-m.opt.gravity[2])
        K_V = 4.0 * W / max(drop, 0.01)
        M_EFF = 42.0
        C_V = 2.0 * float(np.sqrt(K_V * M_EFF))

        # ---- ankle range fix (from test_w2l_ground line ~464) ----
        for anm in ("ankle_L", "ankle_R"):
            jid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, anm)
            m.jnt_range[jid] = np.array([np.radians(-20.0), np.radians(5.0)])

        # ---- drive caps ----
        cap = 0.5 if cap <= 0 else cap
        acap_v = cap if acap < 0 else acap

        # ---- port indices ----
        P_HEEL = [self.net.input_index(p) for p in
                  ("L heel contact", "R heel contact")]
        P_TOE = [self.net.input_index(p) for p in
                 ("L toe contact", "R toe contact")]
        P_E = [self.net.input_index(p) for p in
               ("TONIC L RG ext", "TONIC R RG ext")]
        P_F = [self.net.input_index(p) for p in
               ("TONIC L RG flx", "TONIC R RG flx")]
        P_S1 = self.net.input_index("Stimulus_1")
        P_S2 = self.net.input_index("Stimulus_2 (10nA)")
        P_IB = [[self.net.input_index(p) for p in self.net.inputs
                 if p.startswith("PORT Ib " + s)] for s in ("L", "R")]

        EXT_ACTS = {s: [self.act_idx[k] for k in self.act_idx
                        if "ext" in k.lower() and s.lower() in k.lower()]
                    for s in ("L", "R")}
        FMAX = {k: max(float(m.actuator_gainprm[a][2]), 1.0)
                for k, a in self.act_idx.items()}

        toe_geoms = {s: mujoco.mj_name2id(
            m, mujoco.mjtObj.mjOBJ_GEOM, "toe_" + s.lower() + "_contact")
            for s in ("L", "R")}

        ctn_amp, cref, ibnA = 4.0, 50.0, 1.0
        blend = 0.6

        heelF = np.zeros((n_steps, 2))
        toeF = np.zeros((n_steps, 2))
        com_x = np.zeros(n_steps)
        pel_z = np.zeros(n_steps)
        qbuf = np.zeros((n_steps, 6))
        f6 = np.zeros(6)

        cf = {g: 0.0 for g in list(self.heel_geoms.values()) +
              list(toe_geoms.values())}

        for i in range(n_steps):
            t = i * dt

            # --- sensors (per-plate, like the runner)
            mujoco.mj_forward(m, d)
            for g in cf:
                cf[g] = 0.0
            for k in range(d.ncon):
                c_con = d.contact[k]
                mujoco.mj_contactForce(m, d, k, f6)
                fz = max(f6[2], 0.0)
                if c_con.geom1 in cf:
                    cf[c_con.geom1] += fz
                elif c_con.geom2 in cf:
                    cf[c_con.geom2] += fz

            heelF[i] = [cf[self.heel_geoms["L"]], cf[self.heel_geoms["R"]]]
            toeF[i] = [cf[toe_geoms["L"]], cf[toe_geoms["R"]]]
            com_x[i] = d.subtree_com[0][0]
            pel_z[i] = d.xpos[self.root_body, 2]
            for j in range(6):
                qbuf[i, j] = d.qpos[self.qadr[j]]

            # --- rig scale
            if t < settle:
                S_now = 1.0
            elif t < settle + blend:
                S_now = 1.0 + (rig - 1.0) * (t - settle) / blend
            else:
                S_now = rig

            # --- input currents
            self.u[:] = 0.0
            if t < 0.01:
                self.u[P_S1] = 10.0
                self.u[P_S2] = 10.0
            self.u[P_E[0]] = te
            self.u[P_E[1]] = te
            self.u[P_F[0]] = tf
            self.u[P_F[1]] = tf
            if stimN > 0.0 and stimT <= t < stimT + stimDur:
                self.u[P_E[0]] = max(0.0, te - stimN)
                self.u[P_E[1]] = max(0.0, te - stimN)
            for k in range(2):
                self.u[P_HEEL[k]] = ctn_amp * min(heelF[i, k] / cref, 1.0)
                self.u[P_TOE[k]] = ctn_amp * min(toeF[i, k] / cref, 1.0)
            for si, s in enumerate(("L", "R")):
                if EXT_ACTS[s]:
                    fnorm = np.mean([
                        max(d.actuator_force[a] / FMAX[k2], 0.0)
                        for k2, a in self.act_idx.items()
                        if a in EXT_ACTS[s]])
                    for p in P_IB[si]:
                        self.u[p] = ibnA * fnorm

            # --- neural step (NUMBA-COMPILED: 5.4 us)
            self._step_compiled(
                self.V, self.V_last, self.i_syn, self.u,
                self.syn_pre, self.syn_post, self.syn_gmax,
                self.syn_del_e, self.syn_e_lo, self.syn_range,
                self.in_pre, self.in_post, self.in_gmax,
                self.g_m, self.v_rest, self.i_b, self.tf_mem,
                self.gated_flag, self.k_a, self.k_b, self.k_c,
                self.slope_a, self.slope_b, self.slope_c,
                self.e_a, self.e_b, self.e_c, self.e_ion,
                self.g_ion, self.pow_a, self.pow_b, self.pow_c,
                self.tau_max_b, self.tau_max_c,
                self.b_gate, self.c_gate,
                dt, self.n, len(self.syn_pre))

            # --- muscle activation (with caps)
            ctrl = self.net.muscle_ctrl(self.V)
            for template_label, act_i in self.act_idx.items():
                c_val = ctrl.get(template_label, 0.0)
                c_cap = acap_v if "ankle" in template_label.lower() else cap
                d.ctrl[act_i] = min(c_val, c_cap)

            # --- harness: vertical PD + horizontal leash + tilt
            d.xfrc_applied[:] = 0.0
            if t < hold:
                z_t, wff = z_air, 1.0
            elif t < hold + lower:
                z_t = z_air + (z_contact - z_air) * (t - hold) / max(lower, 1e-6)
                wff = 1.0
            else:
                z_t, wff = z_contact, wsup
            Fz = max(wff * W + K_V * (z_t - float(d.qpos[jadr + 2]))
                     - C_V * d.qvel[jadr + 2], 0.0)
            d.xfrc_applied[self.root_body, 2] = Fz

            # --- physics
            mujoco.mj_step(m, d)

        return dict(heelF=heelF, toeF=toeF, com_x=com_x,
                    pel_z=pel_z, q=qbuf)
