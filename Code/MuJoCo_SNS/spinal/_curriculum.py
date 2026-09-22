"""Staged tuning curriculum (Ben 2026-09-15): deafferented-air ->
afferented-air -> supported ground, each stage seeded by the previous,
only the NEW pathway gains searched per stage.

Stage 1 (air_deaff): deafferented air-stepping. Tunes: drive/rg_adapt/
    desc_e/rg_to_pf (rhythm core).
Stage 2 (air_aff): afferented air (heel/toe + load signals active, feet
    free). Adds: phase_reset_e/f, heel_rge, toe_rge.
Stage 3 (ground): supported ground walk. Adds: ib_rge (stance-Ib
    prolonger), ia_in (IaIN pathway), ankle_post_walk_trim.

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
                 "pm_aff")
STAGE_KEYS = {1: KEYS1, 2: KEYS2, 3: KEYS3}


def set_stage(stage, p):
    OW.load_fitted_baseline()
    OW.set_params(p)
    P.G["renshaw"] = 0.5
    P.G["ankle_post_walk_trim"] = float(p.get("ankle_post_walk_trim", 1.0))
    P.TAU["rg_nap_h"] = float(p.get("rg_nap_h", 0.35))
    if stage >= 2:
        P.G["heel_rge"] = float(p.get("heel_rge", 0.0))
        P.G["toe_rge"] = float(p.get("toe_rge", 0.0))
        P.G["ib_e_central"] = float(p.get("ib_e_central", 0.0))
        P.G["ia_f_central"] = float(p.get("ia_f_central", 0.0))
        P.G["ii_f_central"] = float(p.get("ii_f_central", 0.0))
        P.G["ii_e_central"] = float(p.get("ii_e_central", 0.0))
        P.G["c1_gain"] = float(p.get("c1_gain", 1.0))
        P.G["v3_gain"] = float(p.get("v3_gain", 0.0))
    if stage >= 3:
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
        if stage >= 2:
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
        if stage >= 3:
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
            # lateral rig compliance: ky=5e5 anchors the pelvis and
            # blocks weight transfer (s3g diagnosis) - scale it
            sug["ky_scale"] = trial.suggest_float("ky_scale", 0.02,
                                                  1.0, log=True)
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
    # picked pm_add 0.12 only, right leg 12 cycles, left frozen.
    # curr_s3i adds pm_aff: AFFERENT DISFACILITATION - silence the
    # swing leg's load-afferent inputs during its own swing window
    # (break the re-latch loop at the source).
    name = {1: "curr_s1_air_deaff", 2: "curr_s2b_air_aff",
            3: "curr_s3i_ground"}[stage]
    prev = json.loads(open("best_walk_params_v10.json",
                           encoding="utf-8").read())
    BASE_MUL = dict(prev["multipliers"])
    BASE_MUL["renshaw"] = 0.5
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
            # stage-3 retunes chain from the previous GROUND winner when
            # one exists (s3b); fresh chains fall back to the prior stage
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
        seed["c1_gain"] = max(seed["c1_gain"], 0.1)
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
