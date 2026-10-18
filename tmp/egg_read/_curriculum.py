"""Staged tuning curriculum (Ben 2026-09-15): deafferented-air ->
afferented-air -> STANDING -> ground walk, each stage seeded by the
previous, only the NEW pathway gains searched per stage.
(2026-09-24 night, Ben: "standing should be stage 3 and walking
stage 4" -- the ladder is stand-before-walk.)

Stage 1 (air_deaff): deafferented air-stepping. Tunes: drive/rg_adapt/
    desc_e/rg_to_pf (rhythm core).
Stage 2 (air_aff): afferented air (heel/toe + load signals active, feet
    free). Adds: phase_reset_e/f, heel_rge, toe_rge.
Stage 3 (balance): STANDING-BALANCE stage modeled on SCONE Tutorial 3a
    "Balance" (proprioceptive autogenic length feedback + vestibular
    torso-point PD; Ben's request 2026-09-23, MOVED BEFORE WALKING by
    Ben 2026-09-24 night). Adds: vest_ext / vest_flex_inh (vestibular-
    analog VEST cells -> extensor tone / flexor inhibition) and
    vest_prop (stance-gated II length-loop boost), plus rig_scale
    (wean the support rig so the stage measures real balance, not the
    springs). Objective = COM sway radius + tilt envelope + contact
    symmetry + no falls over an 8 s standing eval (runner --stand-eval).
    All new gains default 0 = previous behavior.
Stage 4 (ground walk): supported ground walking (the s3b..s3k lineage;
    renumbered from 3 on 2026-09-24 night). Adds: ib_rge (stance-Ib
    prolonger), ia_in (IaIN pathway), ankle_post_walk_trim, etc.,
    with full_rules + no_cross pinned on. Seeds from the s3k winner.
Stage 5 (pfvariant, 2026-09-24 night): the per-PF-layer contact
    variant on the ground eval (joint-layer PF fixed ON; only
    heel_pf_layer / toe_df_inh / heel_in_f_exc / ia_pf_f / ii_pf_f
    searched; everything else pinned at the s3k winner).

Usage: python _curriculum.py <stage> [n_trials]
Each stage writes curriculum_stage<N>.json (its best) so the next stage
seeds from it. Resumable per-stage via distinct optuna study names.
"""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

import optuna_walk_v10 as OW
import params as P
import runner as R

DB = "sqlite:///optuna_walk.db"
BASE_MUL = None  # v10 winner multipliers; loaded in main, merged per trial

# searched keys per stage (single source of truth for the objective, the
# seed builder, and the prev-winner filter)
KEYS1 = ("drive", "rg_nap_h", "desc_e", "desc_f", "rg_to_pf")
KEYS2 = KEYS1 + ("heel_rge", "toe_rge", "ib_e_central", "ia_f_central",
                 "ii_f_central", "ii_e_central", "c1_gain", "v3_gain")
# 2026-09-20 v3-objective additions (Ben's critique): pelvis height
# (AARL_PELVIS_TY env; lower the walker so BOTH feet work) and the
# swing-phase ankle-PF inhibition (dorsiflexion in swing - the toe drag).
# 2026-09-21 (t54 diagnosis): pf_gain (v10 left it at 0.41 - PF drive
# alone sits below MN threshold, so an UNLOADED leg's muscles go silent:
# the t54 right leg was a passive flail) and contra_swing (crossed swing
# trigger: opposite foot's heel strike inhibits this side's RG-E /
# disinhibits RG-F; breaks the frozen-stance latch).
KEYS3 = KEYS2 + ("ib_rge", "ia_in", "ia_f_contra_f", "v3_to_ibexc",
                 "ankle_post_walk_trim", "contact_onset",
                 "pelvis_ty", "f1_anklepf_inh",
                 "pf_gain", "contra_swing", "contra_kinh",
                 "pm_gain", "pm_T", "pm_ws", "ky_scale", "pm_add",
                 "pm_aff", "full_rules", "no_cross")
# 2026-09-23 goal2 stage 4 (balance): vestibular-analog tone + stance-gated
# II length loop + support-rig wean. vest_* are new params.G knobs (JSON
# RULE: params default, set_stage loader, runner --best loader branch,
# suggest surface here); rig_scale is a runner arg, not a G knob.
KEYS4 = ("vest_ext", "vest_flex_inh", "vest_prop", "rig_scale")
# 2026-09-24 night: BEN'S REORDER -- standing balance is stage 3 and
# ground walking is stage 4 ("standing should be stage 3 and walking
# stage 4"). The old curr_s4_balance study stays in the DB under its
# name; the renumbered stages use fresh study names (curr_s3_balance /
# curr_s4_nocross) so no study name changes meaning.
KEYS_BAL = KEYS4                      # stage 3: standing balance
KEYS_WALK = KEYS3                     # stage 4: ground walk (s3 lineage)
# 2026-09-24 night stage 5 (Ben: "update the working models with the
# correct rules and run simulations to tune"): PER-PF-LAYER CONTACT
# VARIANT on the ground eval. joint_pf is a FIXED topology switch
# (the PF micro-layers the new edges target exist only in that build);
# the five new gains enter SMALL ([0, 0.5], toe up to the drawn 5).
# Everything else is pinned at the s3k production winner.
KEYS5 = ("heel_pf_layer", "toe_df_inh", "heel_in_f_exc",
         "ia_pf_f", "ii_pf_f", "joint_pf")
STAGE_KEYS = {1: KEYS1, 2: KEYS2, 3: KEYS_BAL, 4: KEYS_WALK, 5: KEYS5}


def set_stage(stage, p):
    OW.load_fitted_baseline()
    OW.set_params(p)
    P.G["renshaw"] = 0.5
    P.G["ankle_post_walk_trim"] = float(p.get("ankle_post_walk_trim", 1.0))
    P.TAU["rg_nap_h"] = float(p.get("rg_nap_h", 0.35))
    if stage in (2, 4, 5):
        P.G["heel_rge"] = float(p.get("heel_rge", 0.0))
        P.G["toe_rge"] = float(p.get("toe_rge", 0.0))
        P.G["ib_e_central"] = float(p.get("ib_e_central", 0.0))
        P.G["ia_f_central"] = float(p.get("ia_f_central", 0.0))
        P.G["ii_f_central"] = float(p.get("ii_f_central", 0.0))
        P.G["ii_e_central"] = float(p.get("ii_e_central", 0.0))
        P.G["c1_gain"] = float(p.get("c1_gain", 1.0))
        P.G["v3_gain"] = float(p.get("v3_gain", 0.0))
    if stage in (4, 5):
        P.G["ib_rge"] = float(p.get("ib_rge", 0.0))
        P.G["ia_in"] = float(p.get("ia_in", 0.0))
        P.G["ia_f_contra_f"] = float(p.get("ia_f_contra_f", 0.0))
        P.G["v3_to_ibexc"] = float(p.get("v3_to_ibexc", 0.0))
        # contact-EVENT transients on the heel/toe ports (runner-side;
        # 2026-09-20 stage-3 latch hypothesis)
        P.G["contact_onset"] = float(p.get("contact_onset", 0.0))
        # 2026-09-20 v3 knobs: swing-phase ankle-PF inhibition (toe-drag)
        P.G["f1_anklepf_inh"] = float(p.get("f1_anklepf_inh", 0.0))
        # 2026-09-21 knobs: crossed swing trigger (runner-side)
        P.G["contra_swing"] = float(p.get("contra_swing", 0.0))
        # crossed KINH drive (build_network conditional edge)
        P.G["contra_kinh"] = float(p.get("contra_kinh", 0.0))
        # per-side contact-reset phase machine (runner-side)
        P.G["pm_gain"] = float(p.get("pm_gain", 0.0))
        P.G["pm_T"] = float(p.get("pm_T", 1.2))
        P.G["pm_ws"] = float(p.get("pm_ws", 0.0))
        P.G["pm_add"] = float(p.get("pm_add", 0.0))
        P.G["pm_aff"] = float(p.get("pm_aff", 0.0))
        P.G["full_rules"] = float(p.get("full_rules", 0.0))
        # s3k: SEVERED-L/R switch (the only config that walked
        # bilaterally) - pin ALL crossed communication off
        if float(p.get("no_cross", 0.0)) > 0.0:
            P.G["contra_swing"] = 0.0
            P.G["contra_kinh"] = 0.0
            P.G["ia_f_contra_f"] = 0.0
            P.G["c1_gain"] = 0.1
            P.G["v3_gain"] = 0.0
        # lateral rig spring scale (runner env AARL_KY)
        import os
        if "ky_scale" in p:
            os.environ["AARL_KY"] = repr(float(p["ky_scale"]))
        else:
            os.environ.pop("AARL_KY", None)
        # pelvis height: runner env (lower the walker onto the ground;
        # the rig anchors at this height)
        import os
        if "pelvis_ty" in p:
            os.environ["AARL_PELVIS_TY"] = repr(float(p["pelvis_ty"]))
        else:
            os.environ.pop("AARL_PELVIS_TY", None)
    if stage == 5:
        # 2026-09-24 night per-PF-layer contact variant (Ben's drawing,
        # Circuit_rules_CONNECTOME_md__connectome.json): joint-layer PF
        # FIXED ON + the five default-0 keys. JSON-RULE loaders.
        P.G["joint_pf"] = float(p.get("joint_pf", 0.0))
        P.G["heel_pf_layer"] = float(p.get("heel_pf_layer", 0.0))
        P.G["toe_df_inh"] = float(p.get("toe_df_inh", 0.0))
        P.G["heel_in_f_exc"] = float(p.get("heel_in_f_exc", 0.0))
        P.G["ia_pf_f"] = float(p.get("ia_pf_f", 0.0))
        P.G["ii_pf_f"] = float(p.get("ii_pf_f", 0.0))
    if stage == 3:
        # goal2 balance stage (2026-09-23; MOVED to stage 3 on
        # 2026-09-24 night, Ben: stand before walk): JSON-RULE loaders
        # for the vest_* G knobs (defaults 0 = previous behavior).
        P.G["vest_ext"] = float(p.get("vest_ext", 0.0))
        P.G["vest_flex_inh"] = float(p.get("vest_flex_inh", 0.0))
        P.G["vest_prop"] = float(p.get("vest_prop", 0.0))


def objective(stage):
    def obj(trial):
        sug = dict(
            # drive ub raised 3.2 -> 4.0 (2026-09-20): the stage-1 winner
            # sits AT the old 3.2 bound (many plateau trials at 3.1-3.2)
            drive=trial.suggest_float("drive", 1.2, 4.0),
            rg_nap_h=trial.suggest_float("rg_nap_h", 0.15, 0.90),
            desc_e=trial.suggest_float("desc_e", 0.8, 1.8),
            desc_f=trial.suggest_float("desc_f", 0.7, 2.2),
            rg_to_pf=trial.suggest_float("rg_to_pf", 1.8, 3.0),
        )
        if stage == 2:
            # NEW GAINS ENTER SMALL (2026-09-20): the 09-18 s2 study
            # sampled these in [0, 1] and EVERY trial lost the stage-1
            # rhythm (trial 0 = stage-1 winner + gains 0.2-0.95 -> static
            # -4.2). Small gains first; widen only if the front pins at
            # 0.5.
            sug["heel_rge"] = trial.suggest_float("heel_rge", 0.0, 0.5)
            sug["toe_rge"] = trial.suggest_float("toe_rge", 0.0, 0.5)
            sug["ib_e_central"] = trial.suggest_float("ib_e_central",
                                                      0.0, 0.5)
            sug["ia_f_central"] = trial.suggest_float("ia_f_central",
                                                      0.0, 0.5)
            sug["ii_f_central"] = trial.suggest_float("ii_f_central",
                                                      0.0, 0.5)
            sug["ii_e_central"] = trial.suggest_float("ii_e_central",
                                                      0.0, 0.5)
            sug["c1_gain"] = trial.suggest_float("c1_gain", 0.1, 0.8)
            sug["v3_gain"] = trial.suggest_float("v3_gain", 0.0, 0.3)
        if stage == 4:
            # GROUND WALK (renumbered from 3 on 2026-09-24 night):
            # the s3-lineage key set, seeded from the s3k winner.
            sug["heel_rge"] = trial.suggest_float("heel_rge", 0.0, 0.5)
            sug["toe_rge"] = trial.suggest_float("toe_rge", 0.0, 0.5)
            sug["ib_e_central"] = trial.suggest_float("ib_e_central",
                                                      0.0, 0.5)
            sug["ia_f_central"] = trial.suggest_float("ia_f_central",
                                                      0.0, 0.5)
            sug["ii_f_central"] = trial.suggest_float("ii_f_central",
                                                      0.0, 0.5)
            sug["ii_e_central"] = trial.suggest_float("ii_e_central",
                                                      0.0, 0.5)
            sug["c1_gain"] = trial.suggest_float("c1_gain", 0.1, 0.8)
            sug["v3_gain"] = trial.suggest_float("v3_gain", 0.0, 0.3)
            sug["ib_rge"] = trial.suggest_float("ib_rge", 0.0, 1.0)
            sug["ia_in"] = trial.suggest_float("ia_in", 0.0, 1.2)
            sug["ia_f_contra_f"] = trial.suggest_float("ia_f_contra_f",
                                                       0.0, 1.0)
            sug["v3_to_ibexc"] = trial.suggest_float("v3_to_ibexc",
                                                     0.0, 1.0)
            sug["ankle_post_walk_trim"] = trial.suggest_float(
                "ankle_post_walk_trim", 0.05, 1.0)
            sug["contact_onset"] = trial.suggest_float("contact_onset",
                                                       0.0, 1.0)
            sug["pelvis_ty"] = trial.suggest_float("pelvis_ty",
                                                   0.88, 0.93)
            sug["f1_anklepf_inh"] = trial.suggest_float(
                "f1_anklepf_inh", 0.0, 1.2)
            # 2026-09-21: PF->MN gain (log-sampled; v10 pinned it at
            # 0.41 which leaves unloaded-leg MNs silent) and the
            # crossed swing trigger
            sug["pf_gain"] = trial.suggest_float("pf_gain", 0.3, 3.0,
                                                 log=True)
            sug["contra_swing"] = trial.suggest_float("contra_swing",
                                                      0.0, 1.5)
            sug["contra_kinh"] = trial.suggest_float("contra_kinh",
                                                     0.0, 2.0)
            # 2026-09-21: per-side contact-reset phase machine
            sug["pm_gain"] = trial.suggest_float("pm_gain", 0.0, 1.0)
            sug["pm_T"] = trial.suggest_float("pm_T", 0.8, 2.0)
            sug["pm_ws"] = trial.suggest_float("pm_ws", 0.0, 1.0)
            # v3: ADDITIVE flexor burst in the swing window
            sug["pm_add"] = trial.suggest_float("pm_add", 0.0, 1.0)
            # v4: afferent disfacilitation in the swing window
            sug["pm_aff"] = trial.suggest_float("pm_aff", 0.0, 1.0)
            # s3j: FULL Deng-style connectome FIXED on (topology, not a
            # free parameter)
            sug["full_rules"] = trial.suggest_categorical(
                "full_rules", [1.0])
            # s3k: SEVERED-L/R fixed on - tune the bilateral config
            sug["no_cross"] = trial.suggest_categorical("no_cross",
                                                        [1.0])
            # lateral rig compliance: ky=5e5 anchors the pelvis and
            # blocks weight transfer (s3g diagnosis) - scale it
            sug["ky_scale"] = trial.suggest_float("ky_scale", 0.02,
                                                  1.0, log=True)
        if stage == 5:
            # 2026-09-24 night per-PF-layer variant: ONLY the five new
            # gains searched (curriculum philosophy: new pathways only;
            # everything else pinned at the s3k winner). Small-enter
            # ranges per the 09-18/09-20 lesson; toe_df_inh tops at the
            # drawn g 5. joint_pf FIXED ON (topology, not a parameter).
            sug["heel_pf_layer"] = trial.suggest_float("heel_pf_layer",
                                                       0.0, 0.5)
            sug["toe_df_inh"] = trial.suggest_float("toe_df_inh",
                                                    0.0, 5.0)
            sug["heel_in_f_exc"] = trial.suggest_float("heel_in_f_exc",
                                                       0.0, 0.5)
            sug["ia_pf_f"] = trial.suggest_float("ia_pf_f", 0.0, 0.5)
            sug["ii_pf_f"] = trial.suggest_float("ii_pf_f", 0.0, 0.5)
            sug["joint_pf"] = trial.suggest_categorical("joint_pf",
                                                        [1.0])
        if stage == 3:
            # STANDING BALANCE (moved before walking, Ben 2026-09-24
            # night): vestibular-analog tone + reciprocal flexor
            # inhibition + stance-gated II loop; gains enter SMALL
            # ([0, 0.5/0.3/1.0] - the 09-18 full-range lesson). rig_scale
            # weans the support rig so the objective measures the
            # controllers, not the springs (AGENTS: support boundary
            # S~0.8-1.0 "until the pelvis-balance piece exists" - this
            # stage is that piece).
            sug["vest_ext"] = trial.suggest_float("vest_ext", 0.0, 0.5)
            sug["vest_flex_inh"] = trial.suggest_float("vest_flex_inh",
                                                       0.0, 0.3)
            sug["vest_prop"] = trial.suggest_float("vest_prop", 0.0, 1.0)
            sug["rig_scale"] = trial.suggest_float("rig_scale", 0.05,
                                                   1.0)
        # searched keys override; everything else pinned at v10 winner
        p = {**BASE_MUL, **sug}
        set_stage(stage, p)
        if stage <= 2:
            # stages 1-2: AIR stepping (stage 1 deafferented + interleg
            # off; stage 2 AFFERENTED + interleg ON, exercising the c1/V3
            # commissurals in the Ivanenko air-stepping prep)
            args = ["--no-ground", "--no-afferents", "--no-interleg",
                    "--time", "14", "--drive", repr(p["drive"])]
            if stage == 2:
                args = ["--no-ground", "--time", "14",
                        "--drive", repr(p["drive"])]
            m = R.main(args)
            import numpy as np
            z = __import__("numpy").load("spinal_run.npz", allow_pickle=True)
            t, q, neuro = z["t"], z["q"], z["neuro"]
            m = (t >= 5.0) & (t <= 17.0)
            if not np.all(np.isfinite(q[m])) or \
                    not np.all(np.isfinite(neuro[m])):
                return -200.0
            knee = q[m, 4]
            if not (-360.0 < float(knee.min()) < 360.0) or \
                    not (-360.0 < float(knee.max()) < 360.0):
                return -200.0  # unphysical RoM (finite but exploded)
            rge = neuro[m, 2]
            on = rge > 0.5 * max(rge.max(), 1e-9)
            rises = int(np.sum(np.diff(on.astype(int)) == 1))
            knee = q[m, 4]
            if rises > 30:
                return -200.0  # runaway flutter, not stepping
            # RHYTHM GATE (2026-09-20): the 09-18 s2 exploit - a static
            # deep-flexion pose scores 0.5*(-knee_min) with rises=0 (the
            # recorded 36.906 was EXACTLY 0.5*73.81, zero bursts) and it
            # won the study while never stepping. Static poses must
            # order below every genuine rhythm (>=3 bursts = >=~0.25 Hz
            # in the 12 s window, the Ivanenko-preferred slow end).
            if rises < 3 or (float(rge.max()) - float(rge.min())) < 1.0:
                return -10.0 + 0.05 * (-float(knee.min()))
            # air objective: rhythmic + deep knee swing flexion
            score = 3.0 * rises + 0.5 * (-float(knee.min()))
            if not np.isfinite(score):
                return -200.0
            return float(score)
        if stage == 3:
            # STANDING-BALANCE eval (SCONE Tutorial-3a analog):
            # 8 s standing (--stand-eval 8, DRIVE 0 throughout) at this
            # trial's rig wean. Sentinels: NaN -200 < fall -150 < any
            # stander. Terms: COM sway radius (m, 400/m), tilt envelope
            # (deg, 1/deg), foot-load symmetry deviation from even
            # (bal_contact_sym = 0.5 = even split; ABS - either side
            # overloading costs the same; 20 max cost at sym 0 or 1).
            m = R.main(["--stand-eval", "8",
                        "--rig-scale", repr(float(p["rig_scale"]))])
            if m["nan"]:
                return -200.0
            if m["bal_fell"]:
                return -150.0
            return (100.0
                    - 400.0 * float(m["bal_sway"])
                    - 1.0 * float(m["bal_tilt_max"])
                    - 40.0 * abs(0.5 - float(m["bal_contact_sym"])))
        m = R.main(["--eval", "--drive", repr(p["drive"])])
        if m["nan"]:
            return -400.0
        if m.get("kine") is None:
            # no cycles on either leg. MUST sit below the worst genuine
            # walker: under kine_ref v2 (both legs + guards) the s3b
            # one-legged "winner" scores -279, and a frozen walker is
            # worse than a bad stepper but must still beat this sentinel.
            return -320.0
        # clip so pathological walkers stay above the sentinels
        score = max(float(m["kine_score"]), -315.0)
        if m["kz"] < 0.62:
            score -= 20.0
        if m["tilt_max"] > 40.0:
            score -= 10.0
        return score
    return obj


def main():
    global BASE_MUL
    stage = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    # 2026-09-20: fresh study names. The curr_s2_air_aff study was won
    # by the static-pose exploit (best 36.906 = 0.5*73.81, rises 0) and
    # curr_s3_ground (50 trials, all -100) seeded from it; both archived
    # to curriculum_exploit_archive_20260920.json and deleted.
    # curr_s3b_ground was then won under the v1 (right-leg, DC-removed,
    # neural-duty) objective by a ONE-LEGGED gait (left foot planted
    # 100% of frames - Ben's critique); it scores -279 under kine_ref
    # v2. curr_s3c_ground = the both-leg objective retune (80 trials;
    # corrected best t54 -181.6, ALL trials frozen-left). curr_s3d adds
    # pf_gain + contra_swing after the t54 passive-flail diagnosis
    # (best -192.0, still frozen-left: RG-level crossed kicks do not
    # release a loaded jammed leg). curr_s3e adds contra_kinh (crossed
    # KINH; best -191.8 = its seed, still frozen-left). curr_s3f adds
    # the per-side CONTACT-RESET PHASE MACHINE (heel-strike reset,
    # antiphase coupling, swing-window MN gating) - the Di Russo
    # the per-side CONTACT-RESET PHASE MACHINE (heel-strike reset,
    # antiphase coupling, swing-window MN gating) - the Di Russo
    # eq-7 analog (s3f best -164.4 t32, right leg cycles, left frozen).
    # curr_s3g added WEIGHT-SHIFT (pm_ws + load-gated window +
    # ky_scale lateral rig compliance): seed still best (-182.4), left
    # frozen. curr_s3h adds pm_add: ADDITIVE flexor burst in the swing
    # window - the v2 multiplicative boost multiplied a ~0 ctrl (MNs
    # below threshold) and never actually drove anything. s3h: search
    # (break the re-latch loop at the source). curr_s3j: FULL
    # Deng-style connectome FIXED ON (full_rules=1) - retune the
    # operating point for the corrected topology: best -164.9 (t7),
    # left frozen in EVERY trial. curr_s3k: SEVERED-L/R (no_cross
    # fixed on, full_rules kept on) - retune the ONLY configuration
    # that produced bilateral stepping.
    name = {1: "curr_s1_air_deaff", 2: "curr_s2b_air_aff",
            3: "curr_s3_balance", 4: "curr_s4_nocross",
            5: "curr_s5_pfvariant"}[stage]
    prev = json.loads(open("best_walk_params_v10.json",
                           encoding="utf-8").read())
    BASE_MUL = dict(prev["multipliers"])
    BASE_MUL["renshaw"] = 0.5
    if stage in (4, 5):
        # pin the non-searched keys at the s3k PRODUCTION winner (the
        # full_rules/no_cross operating point stages 4-5 retune)
        s3k = json.loads(open(
            "reports_20260923/s3k_trial34_full_params.json",
            encoding="utf-8").read())["params"]
        BASE_MUL.update({k: v for k, v in s3k.items()})
        if s3k.get("pf_gain") is None:
            BASE_MUL.pop("pf_gain", None)
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction="maximize", storage=DB,
                                study_name=name, load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(
                                    seed=21 + stage, n_startup_trials=8))
    if len(study.trials) == 0:
        # FULL-dict seed (JSON-rule lesson: missing keys get SAMPLED by
        # optuna, so the seed must pin every searched key explicitly).
        # Stage keys absent from the previous winner start at 0 = the
        # previous stage reproduced exactly (trial 0 = stage-1 rhythm
        # with all new pathways OFF).
        sk = STAGE_KEYS[stage]
        seed = {k: 0.0 for k in sk}
        if stage == 5:
            # joint-layer PF is a FIXED-ON topology switch in stage 5;
            # the seed carries ONLY the stage keys -- the s3k operating
            # point comes from the BASE_MUL merge above
            seed["joint_pf"] = 1.0
        elif stage == 3:
            # balance stage: vest gains seed at 0 (defaults-off =
            # today's standing reproduced exactly); rig_scale seeds FULLY
            # SUPPORTED (1.0) so trial 0 is the incumbent standing config.
            seed["rig_scale"] = 1.0
        elif stage in (1, 2):
            seed["c1_gain"] = 0.1  # floor of its range
            # sensible mid-range defaults for keys the previous winner
            # cannot carry (a 0.0 pelvis_ty would put the walker at ground
            # level!)
            seed["pelvis_ty"] = 0.905
            seed["f1_anklepf_inh"] = 0.3
            seed["pf_gain"] = 1.0      # 2026-09-21: raise from v10's 0.41
            seed["contra_swing"] = 0.5
            seed["contra_kinh"] = 0.5
            seed["pm_gain"] = 0.5
            seed["pm_T"] = 1.23
            seed["pm_ws"] = 0.4
            seed["pm_add"] = 0.4
            seed["pm_aff"] = 0.8
            seed["ky_scale"] = 0.05
        try:
            # stages 4-5 (walk + pfvariant) chain from the s3k
            # PRODUCTION winner directly (their numeric predecessors
            # are the air/balance stages, which carry no walk keys)
            if stage in (4, 5):
                prevw = json.loads(open(
                    "reports_20260923/s3k_trial34_full_params.json",
                    encoding="utf-8").read())["params"]
            else:
                prevw = json.loads(open(
                    f"curriculum_stage{stage}.json",
                    encoding="utf-8").read())["params"]
        except FileNotFoundError:
            try:
                prevw = json.loads(open(
                    f"curriculum_stage{stage-1}.json",
                    encoding="utf-8").read())["params"]
            except FileNotFoundError:
                prevw = {}
        for k in sk:
            if k in prevw:
                seed[k] = float(prevw[k])
        if "c1_gain" in sk:
            seed["c1_gain"] = max(seed["c1_gain"], 0.1)
        # full_rules is a FIXED topology switch in s3j: force the seed
        # to match the categorical even when the previous winner
        # predates the connectome change
        if "full_rules" in sk:
            seed["full_rules"] = 1.0
        if "no_cross" in sk:
            seed["no_cross"] = 1.0
        study.enqueue_trial(seed)
        print("seeded", json.dumps(seed), flush=True)
    study.optimize(objective(stage), n_trials=n, gc_after_trial=True)
    best = study.best_trial
    print(f"== stage {stage} best {best.value:.3f} (trial {best.number})")
    print(json.dumps(best.params, indent=1))
    with open(f"curriculum_stage{stage}.json", "w", encoding="utf-8") as f:
        json.dump({"stage": stage, "score": best.value,
                   "params": best.params, "trial": best.number,
                   "study": name}, f, indent=2)
    print(f"saved curriculum_stage{stage}.json")


if __name__ == "__main__":
    main()
