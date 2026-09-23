# Goal 3 — MuJoCo muscle fidelity: series-elastic tendons, non-linear contact damping, adaptive integration

**Prepared:** 2026-09-23 (EB475WS4, myo env, MuJoCo 2.3.7 — the pinned version).
**Prompt:** https://hyfydy.com/hyfydy-vs-mujoco/ and its bottom-of-page citations.
**Scope:** ANALYSIS ONLY. No model XML or code was edited. Every empirical number below was
produced this session in throwaway scripts under `%TEMP%` (kept: `goal3_see_demo.py` family).

---

## TL;DR — what Ben must decide

1. **All three Hyfydy differentiators are real gaps in MuJoCo, and MuJoCo's own documentation
   admits two of them in words.** None of the three can be added natively in 2.3.7.
2. **Nothing here should touch the SNS/spinal pipeline before the dissertation deadline.**
   The bit-exact reproduction contract (regression gate `--fitted --best --eval` =
   −61.1741850610298; bit-exact reproductions of v6/v7/v9/v10 winners; Simulink E2 parity at
   fixed 2 ms) is a property of the *exact* model text + code + stepping semantics. ANY change
   to muscle parameters, contact parameters, or the integrator **invalidates every recorded
   winner and study in the optuna DB** (ground_walk_v3…v11, curriculum s1–s3b, basin-gate
   state dumps, and the Simulink E1/E2 results). Retuning is measured in days of campaign
   time, not an evening.
3. **Within 2.3.7 the honest menu is:** (a) SEE — engine-level NO; parameter-level Zajac
   approximation is possible today without new states; a structural extra-DOF workaround is
   *expressible* (demonstrated below) but demands ≥1 kg virtual tendon masses or a 4× smaller
   timestep — both invasive. (b) non-linear contact damping — *partially* native via
   `solimp`'s depth-dependent impedance (effective k **and** b grow with penetration —
   Hunt-Crossley-like in effect, not in law); exact HC and viscous friction are impossible
   because contact forces come out of the constraint solver. (c) adaptive integration —
   impossible natively; hand-rolled outer-loop stepping only, and it fights both the 2 ms
   SNS lockstep and cross-run determinism.
4. **Recommended path:** keep 2.3.7 + `implicitfast` + 2 ms for everything dissertation-bound;
   buy fidelity confidence with an *offline convergence ladder* (re-run winner configs at
   1 ms / 0.5 ms and compare ensemble statistics — pointwise comparison is meaningless past
   ~0.4 s per the bifurcation findings already banked in DESIGN.md). Revisit SCONE+Hyfydy
   **after** the deadline as the standing-balance validation engine: it is already installed
   and verified on this machine; the gate is a free license key, and the port cost is a
   co-sim rewrite + `.hfd` conversion + full retune (details §7).

---

## 1. What Hyfydy claims, and which citations support each claim

Source: `WebFetch` of https://hyfydy.com/hyfydy-vs-mujoco/ (2026-09-23). The page is a
marketing comparison; **it contains no quantitative benchmarks** — it says Hyfydy and MuJoCo
are "similar in speed" and both are "orders of magnitude faster than OpenSim", nothing more.

| # | Hyfydy claim (page's words) | Citation at bottom of page | My read |
|---|---|---|---|
| 1 | Muscles "accurately simulate tendon elasticity"; MuJoCo "does not model this phenomenon" | Millard, Uchida, Seth & Delp 2013 (*J Biome Eng* 135(2):021005, flexing computational muscle); Blazevich & Fletcher 2023 (*Biol Rev*, long-Achilles benefits) | **True.** MuJoCo docs: "We assume inelastic tendons while OpenSim can model tendon elasticity." Claim↔citation mapping is the page's bottom-list ordering; I treat Millard 2013 → MT-dynamics model and Blazevich & Fletcher → *why elasticity matters* as the page's intent (inference — the fetched text did not anchor each citation inline). |
| 2 | Contacts include "realistic non-linear damping" plus "dynamic and viscous friction coefficients"; MuJoCo contacts "do not model these properties" | Hunt & Crossley 1975 (*J Appl Mech* 42(2):440, restitution as damping in vibroimpact) | **Largely true.** MuJoCo has no HC force law and no viscous/velocity-dependent friction (verified against docs, §4.2). Caveat: MuJoCo's `solimp` impedance ramp does give *depth-dependent* effective stiffness **and** damping — a parameter-level cousin of HC damping — so "does not model these properties" is fair as stated (no damping *law*, no viscous friction) but not the whole story. |
| 3 | "Error-controlled integration, which adapts the integration step size"; MuJoCo uses fixed-step with no error control | (no dedicated citation on the page) | **True.** Verified in the 2.3.7 headers and docs (§5). |
| 4 | (context) Both engines fast; Hyfydy proprietary, uneven terrain; OpenSim slow | Schumacher et al. 2025 (*iScience* 112203, learned bipedal walking — a Hyfydy user); Seth et al. 2018 (*PLoS CB*, OpenSim) | Schumacher 2025 = evidence the engine is usable for learned bipedal gait; Seth 2018 = the OpenSim baseline. |

Note the irony worth having in the dissertation's methods section: MuJoCo's own muscle
activation dynamics (`timeconst` defaults 0.01/0.04 s) are documented "per Millard et al.
(2013)" — the same Millard paper Hyfydy cites for its *elastic-tendon* MT model. The
difference is the tendon state, not the activation curve family.

---

## 2. Our current stack — measured facts (all run this session)

Model of record: `Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml`
(referred to as cvt3; loaded and inspected with `mujoco` 2.3.7 from the myo env).

| Item | Value | Citation |
|---|---|---|
| `<option>` | `timestep="0.005" collision="predefined"` — no integrator attribute → **Euler** | cvt3.xml:4 |
| Muscle default class | `dyntype="muscle" gaintype="muscle" biastype="muscle"`, `dynprm="0.01 0.04 …"` | cvt3.xml:16 (default class block) |
| Muscle state count | **`actnum = 1` per muscle (activation only)** — no tendon-stretch state exists | measured on compiled cvt3 + simbridge2 (`vas_med_r`: dyntype 3, actnum 1) |
| Per-muscle gainprm (example `vas_med_r`) | `[0.49, 1.4095, 1215.67, 1, 0, 2, …]` → range 0.49–1.41, **Fmax = gainprm[2] = 1215.7 N** | measured; matches the skill/AGENTS rule `gainprm[:,2]` = Fmax |
| Tendons | **92 spatial tendons, `tendon_stiffness` unique = {0.0}, `tendon_damping` unique = {0.0}** — pure kinematic transmissions | measured on compiled model |
| Contacts | 19 predefined ground-vs-body `<pair>` lines, **no solref/solimp/friction/condim overrides** → compiled values: `solref=[0.02 1]`, `solimp=[0.9 0.95 0.001 0.5 2]`, `friction=[1 1 0.005 1e-4 1e-4]`, `condim=3` | cvt3.xml:629-647; measured `pair_solref/solimp/friction/dim` |
| The one custom `solimp` | `0.9999 0.9999 0.001 0.5 2` — on the 82 **equality-coupled pathpoint followers**, not on any contact | cvt3.xml:650-731 |
| Production option line | **`timestep 0.002 + integrator="implicitfast"`**, chosen because "stiff joint springs + explicit Euler at dt=5 ms explode; implicitfast handles the joint-damping/spring terms stably" | runner.py:176-180 |
| Contact surgery in the runner | 19 pairs → keep the 6 foot-ground pairs (talus/calcn/toes); air mode deletes the section; predefined pairs ignore contype/conaffinity | runner.py:352-368 |
| Rig / ligaments | spring+damper specs on pelvis/lumbar/hip-rotation (e.g. 2e5 N/m, d 2000) and 10–30 N·m/rad ligament surrogates; `--rig-scale S` scales k×S, c×√S | runner.py:392-464 |
| Simulink bridge model | `simbridge2.xml` = same rewrite `0.005 → 0.002 + implicitfast` so the bridge plant steps at the SAME 2 ms as the SNS network (single-rate co-sim, exact parity) | make_simbridge_xml.py:43-54 (Code/Matlab/SNS_Simscape/mujoco_bridge/matlab/) |
| MuJoCo pin | 2.3.7 — the Simulink blockset is mexed against it; 3.3.6 cannot parse `collision="predefined"` | AGENTS.md bridge section; skill |
| Conversion kept tendon data but discarded elasticity | Step3 pkls carry `tendon_sla` (soleus 0.25 m, med_gas 0.39 m, tib_ant 0.223 m) + `fiber_opt`, `penna_opt`, `fmax` — then bake into the rigid-tendon MuJoCo muscle | `Step3_muscleKinetics/soleus_l.pkl` etc., read this session |

**Reading:** the muscles where SEE matters most physiologically (soleus, gastrocnemius — the
Achilles group; Blazevich & Fletcher 2023's whole subject) are exactly the ones with the
longest tendon slack lengths sitting unused as elasticity in our model.

---

## 3. Feature (a): series-elastic tendon (SEE)

### 3.1 Native in 2.3.7 — NO (and still NO in 3.x)

- MuJoCo docs (stable, fetched 2026-09-23): *"We assume inelastic tendons while OpenSim can
  model tendon elasticity."* / *"We assume that the biological tendon is inelastic, with
  constant length L_T, while the biological muscle length L_M varies over time."* The
  suggested in-model approximation is Zajac-1989-style: *"their effect can be captured
  approximately by stretching the FL curve"* (shortening the muscle operating range). For a
  genuine elastic tendon the docs point at **user callbacks**: custom models *"could then
  simulate elastic tendons."*
- Compiled model evidence: every muscle has `actnum = 1` (activation only) and all 92 tendon
  elements have `stiffness = damping = 0`. There is no per-tendon force state anywhere.
- Upgrading MuJoCo does **not** buy (a): the stable docs (3.x-era, incl. the new Sept-2026
  `discrete` integrator) still say inelastic tendon. **Migrating off 2.3.7 for SEE alone
  would break the Simulink bridge for nothing.**

### 3.2 Workarounds inside 2.3.7, with exact knobs

**(i) Parameter-level — the Zajac approximation (cheap, contract-breaking but structurally safe).**
Widen/shift the muscle operating range and FLV shape so the lumped compliance appears as a
softer force-length slope: `gainprm[0..1]` (range, e.g. `0.75 1.05`) on the `<general
class="muscle">` elements, plus `lmin/lmax/fpmax` (gainprm[4,5,7]). This is what MuJoCo's
docs recommend, it adds **no states and no stiffness**, so 2 ms stepping is untouched.
What it does NOT give: energy storage-and-release dynamics (no rebound, no SEE work loop),
no force filtering — it is a static lumping. This is also exactly the knob MyoConverter's
Step3 fit already exercises, so re-fitting ranges is a pipeline we already own
(`Step3_muscleKinetics/*.pkl` → `opt_results.res_opt`).

**(ii) Structural — extra DOF + spring tendon (demonstrated this session).**
A real SEE needs an independent stretch coordinate. In MJCF that is: an intermediate body
with a slide joint along the muscle line, the actuator pulling *that* body, and a
`<spatial>` tendon with `stiffness`/`damping`/`springlength` from the intermediate body to
the skeleton site. Knobs (toy model kept at `%TEMP%\goal3_see_demo.py`, MuJoCo 2.3.7):

```xml
<body name="ce" pos="0 0 0.5">           <!-- virtual CE body -->
  <joint name="ce_z" type="slide" axis="0 0 1" damping="0.1"/>
  <geom type="box" size="0.03 0.03 0.03" mass="1.0"/>   <!-- the knob that matters -->
  <site name="ce_bot" pos="0 0 -0.05"/>
</body>
<tendon>
  <spatial name="muscle_t"><site site="anchor"/><site site="ce_top"/></spatial>
  <spatial name="see" stiffness="3.2e5" damping="100" springlength="0.4 0.4">
    <site site="ce_bot"/><site site="load_site"/>       <!-- Achilles-like k -->
  </spatial>
</tendon>
```

Measured behavior (1000 N muscle step at t=0.1 s, 80 kg load, `implicitfast`):

| config | dt | stable? | max SEE force on load | verdict |
|---|---|---|---|---|
| rigid baseline (motor on load directly) | 2 ms | yes | direct 1000 N | current-world reference |
| SEE, CE mass **1.0 kg** | 2 ms | yes | **16,300 N** (16× step overshoot, ringing) | expressible, ugly transient |
| SEE, CE mass **0.1 kg** (pathpoint-like) | 2 ms | **no** — MuJoCo itself warns "simulation is unstable"; forces reach 4.6×10⁸ N | garbage | **unusable at the contract timestep** |
| SEE, CE mass 0.1 kg | **0.5 ms** | yes | 10,760 N | needs 4× smaller dt |

Two structural traps compound this for *our* model: (1) MuJoCo tendon springs are functions
of the **kinematic** tendon length — with stiffness 0 they are transmissions (what we have
now), and the actuator force and any spring force act **in parallel** on the same
coordinate; the series property only appears once you donate a real DOF to the stretch, as
above. (2) The converted model's pathpoint bodies are already the smallest-mass danger zone
(`boundmass 0.1` was itself a hard-won fight, runner.py:181-193); 92 extra ~1 kg
high-frequency DOFs would be a new numerical regime for the equality-coupled mess, and
0.5 ms stepping quadruples every campaign.

**(iii) Callback-level — custom actuator dynamics without an engine change.**
The 2.3.7 Python bindings expose `set_mjcb_act_dyn / set_mjcb_act_gain / set_mjcb_act_bias /
set_mjcb_passive / set_mjcb_control / set_mjcb_sensor` (verified: `dir(mujoco)` this
session). A dyntype="user" actuator could therefore carry a tendon-stretch state computed
in a Python callback — no C plugin strictly required for a prototype. But `mj_step`
integrates that state with the same explicit/semi-implicit scheme at 2 ms, so a
3.2×10⁵ N/m tendon state hits exactly the instability the toy demonstrated; a stable
version needs an internal quasi-static/equilibrium solve (which degenerates back toward
rigid tendon) or sub-stepping. Inference (mine, not measured): performance is acceptable
(one callback per step over actuator arrays, not per muscle), but this is research code,
not a knob.

**(iv) Full Millard/Thelen with elastic tendon = engine/plugin change** (C plugin with
actuator states, or a different engine). Nothing in the 2.3.7→3.x line provides it natively.

### 3.3 Effort/risk summary (a)

| option | effort | timestep survives 2 ms? | bit-exact winners survive? |
|---|---|---|---|
| (i) range/FLV re-fit | hours (re-run Step3 fits) | yes | **no** (model text changes) |
| (ii) extra-DOF SEE ×92 | days–weeks, new numerical regime | **no** (0.5 ms for light masses) | no |
| (iii) act_dyn callback prototype | days, research-grade | marginal (stiffness trap) | no |
| (iv) plugin/engine | weeks | n/a | no |

---

## 4. Feature (b): non-linear contact damping / realistic friction

### 4.1 What we run today

The 6 foot-ground pairs (runner.py:363-368) inherit MuJoCo defaults: `solref=[0.02 1]`,
`solimp=[0.9 0.95 0.001 0.5 2]`, `friction=[1 1 0.005 1e-4 1e-4]`, `condim=3`, Newton
solver (`opt.solver=2`), `impratio=1`, `noslip_iterations=0` (all measured this session).
Nothing anywhere tunes contact — the only contact-adjacent damping in the stack is joint
damping (leg hinges, rig dampers) and the `damping="5"`-style DOF damping.

### 4.2 Native in 2.3.7 — partially, at the parameter level

- **Depth-dependent (non-linear) stiffness AND damping: YES, via `solimp`.** Docs (fetched):
  impedance d(r) interpolates sigmoidally from `d0` at r=0 to `dw` at r=width, shaped by
  `midpoint` and `power` (power ≥ 1; power=1 linear). Effective stiffness is d·k and
  effective damping d·b, so **both grow non-linearly with penetration** — structurally the
  same idea as Hunt-Crossley damping (damping growing with indentation) though not the same
  law (HC: F = k·δⁿ(1 + (3/2)α·δ̇)).
- **Direct k,b control: YES, via negative `solref`.** `solref="-stiffness -damping"` sets
  them directly with fixed impedance scaling — docs recommend this convention for system
  identification; it also permits perfectly elastic (b=0) contacts.
- **Separate friction-constraint softness: YES, `solreffriction`** — exists at **pair**
  level in 2.3.7 (`model.pair_solreffriction`, verified; NOT on geoms in this version).
- **Friction-cone and slip knobs:** elliptic vs pyramidal cones, `impratio` (verified
  present), `noslip_iterations` (docs recipe against slip: elliptic cones + large impratio
  + Newton + tight tolerance, then NoSlip).
- **`timeconst` floor:** docs require `timeconst ≥ 2× timestep` unless `refsafe` is
  disabled — at 2 ms that means `timeconst ≥ 0.004`; we run 0.02 (10× dt), i.e. plenty of
  headroom to tighten to 0.01 or 0.005 if Ben ever wants crisper foot contacts.

### 4.3 Impossible without an engine change

- **An exact Hunt-Crossley (or any custom) contact force law.** MuJoCo contact forces are
  the solution of a constrained optimization (Newton/CG/PGS over soft constraints), not a
  per-contact force function; there is **no custom contact-force callback** (docs + callback
  inventory above — the closest hooks are `mjcb_contactfilter`, which only filters pairs).
  You can shape the effective law with `solref`/`solimp`; you cannot write one.
- **Viscous friction (tangential force ∝ tangential velocity) and dynamic-vs-static
  friction split.** The friction constraints have r ≡ 0 and k = 0 (docs) — first-order
  exponential velocity decay inside the Coulomb limit, single coefficient set per contact.
  Nearest proxies: DOF damping on the joints (models joint viscous friction, not surface),
  global `opt.density`/`opt.viscosity` fluid drag (verified present; global, not
  per-contact), or the NoSlip post-pass.

### 4.4 What I would change if Ben wants "more realistic foot contacts" (not this round)

Exact knobs, all on the 6 kept pairs (or the geoms' class): tighten `solref` timeconst
toward `0.01 1`; consider negative-solref with k,b set from shoe-floor stiffness measurements;
raise `power` in `solimp` (e.g. 2 → 1.5/3 changes the depth-hardening curve shape); set
`solreffriction` on the pairs; switch to elliptic cone + `impratio 10` if slip artifacts
appear in standing balance. Cost: every one of these changes the model → **full retune,
contract void** (§6). Effort to *evaluate*: an hour of swept re-runs; effort to *adopt*:
a retune campaign.

---

## 5. Feature (c): error-controlled / adaptive integration

### 5.1 Native in 2.3.7 — NO

- `mjmodel.h:127-132` (shipped with the pinned wheel): `mjINT_EULER` (semi-implicit Euler),
  `mjINT_RK4`, `mjINT_IMPLICIT`, `mjINT_IMPLICITFAST`. Nothing else. Verified again via
  `dir(mujoco.mjtIntegrator)` this session.
- Docs (stable, fetched): RK4 is *"the fixed-step 4th-order Runge-Kutta method"*; the
  timestep *"is perhaps the single most important parameter that the user can adjust"*;
  `implicitfast` is *"the recommended integrator… best tradeoff of stability and
  performance"*; RK4 is recommended for near-energy-conserving systems, while *"in the
  presence of large velocity-dependent forces, if the chosen single-step method integrates
  those forces implicitly, single-step methods can be significantly more stable than RK4."*
  Our model (springs + dampers everywhere) is exactly the velocity-dependent-force case —
  which is why runner.py already chose `implicitfast` (runner.py:176-180).
- The new `discrete` integrator documented on the stable docs page is a Sept-2026 / MuJoCo
  3.x feature and is itself fixed-step. **Upgrading does not buy (c) either.**

### 5.2 Workaround — hand-rolled outer-loop step control (and why I don't recommend it)

`m.opt.timestep` is mutable at runtime, so an outer loop can implement step-doubling /
embedded error estimates (e.g., compare an RK4 step against two half-steps, adjust dt to
hold an error tolerance). Costs, all of which land on things we care about:

1. **It dissolves the SNS lockstep.** The network solves at 2 ms and the Simulink bridge is
   single-rate 2 ms *by design* (make_simbridge_xml.py:43-54, E2). Variable MuJoCo stepping
   needs interpolation or a fixed communication step anyway — at which point the physics is
   adaptive but the coupling is not, and E2 parity is gone.
2. **It dissolves determinism.** We already live with chaos (pointwise agreement between any
   two integrators dies at ~0.4 s; only seeded identical binaries reproduce bit-exactly —
   DESIGN.md/AGENTS.md). An adaptive controller multiplies the branch-on-floating-point
   surface: runs become irreproducible across BLAS orders and machines, which is worse than
   the current controlled chaos.
3. **It is the wrong tool for our failure modes.** Our blowups have been structural
   (mass-matrix singularities, constraint fights, undamped rig springs), not truncation
   error — each was fixed by model repairs, and `implicitfast` handles the stiff
   velocity-dependent terms, which is precisely where error control would otherwise help.

### 5.3 The right fidelity tool: an offline convergence ladder

Re-run frozen winner configs at 1 ms and 0.5 ms (same integrator) and compare **ensemble
statistics** (duty, cadence, amplitudes, kine_score), not pointwise traces. This quantifies
how much the 2 ms discretization biases the gait metrics without touching the production
contract. Cheap (one `--eval` sweep per dt), zero risk, and it is the honest answer to a
reviewer asking "is 2 ms enough?". Not run this round — it is a campaign-scale run, listed
as the recommended next step.

---

## 6. Risk to the bit-exact reproduction contract — explicit

**ANY change to the model XML (muscle params, contact params, extra DOFs), the option line
(timestep, integrator), or stepping semantics invalidates, at minimum:**

- the regression gate `runner --fitted --best --eval --drive 2.5868` = −61.1741850610298;
- every recorded winner reproduction (`--best6/--best7/--best10`, v8/v9/v10/v11 numbers,
  curriculum s1–s3b winner trial 52 = −85.22);
- the optuna studies' comparability (ground_walk_v3…v11, curr_s2b/s3b — old trials were
  scored on the old physics);
- `basin_gate_results.json` (state dumps reference the current param set);
- the Simulink E1/E2 results and the bridge's bit-exact sensor parity (bridge mex + model
  are pinned to the same physics);
- every figure in `Dissertation\CPG_airstepping_figs` regenerated from spinal_run.npz.

This is not a "re-check" cost, it is a **re-tune-everything** cost. If Ben adopts any model
change, the clean process is: freeze the old world (tag the current cvt3+patch as a named
configuration), make the change a patch_xml branch behind a default-off flag (the
conditional-topology pattern the codebase already uses for zero-gain synapses), re-run
sentinels, then re-open studies seeded from the old winners.

---

## 7. SCONE + Hyfydy for standing-balance work

### 7.1 What it buys

- **All three features natively**: elastic-tendon muscle models (SCONE exposes
  `tendon_slack_length`, `stiffness_multiplier` as `MuscleModifier` parameters — i.e. SEE
  is a first-class, *optimizable* property), Hunt-Crossley-style non-linear contact damping
  + dynamic/viscous friction, and error-controlled integration. Uneven terrain is
  Hyfydy-only.
- **Speed**: skill-verified claim ~10-50× faster than the OpenSim backends; vs our MuJoCo
  pipeline it is "similar" per the Hyfydy page (no numbers published).
- **Purpose-built objective machinery for balance**: `BalanceMeasure`, perturbation training
  (`PerturbationController` pushes), NoiseController robustness — the standing-balance
  curriculum we are hand-building as optuna objectives already exists there as maintained
  measure types.
- **The Geyer-Herr 2010 controller lineage** lives there (`GaitStateController`,
  `ControllerGH2010v12.scone` with paper-value comments) — a direct literature-faithful
  baseline to compare our SNS network against (supervisor-review friendly).

### 7.2 What porting costs (grounded in the installed 2.4.4 on THIS machine)

1. **License gate (the only hard blocker)**: Hyfydy scenarios currently fail with "no
   active license key found". Free non-commercial key: `sconecmd --hyfydy id` → request at
   hyfydy.com → `sconecmd --hyfydy <KEY>`. Until then only the OpenSim 3/4 backends run.
2. **Model**: SCONE loads `.osim` directly (our Gait2392 family *should* load — skill marks
   it UNTESTED; stock Gait2392_Simbody first). Hyfydy needs `.hfd` via `hfdmodeltool.exe`.
   Our robotbody repairs (hip-axis flips, rect_fem patella reroute, prune list, Fmax fixes)
   would have to be re-applied in the .osim/.hfd world — they currently live as MuJoCo-side
   `patch_xml` surgery, not in any .osim.
3. **Controller**: our 410-neuron SNS network has no SCONE controller type. The realistic
   path is the SconePy RL-loop (`load_model`, `set_actuator_inputs`,
   `advance_simulation_to(t)` — verified working on this machine in the py3.9 `scone` env):
   port runner.py's numpy SNS↔backend loop to drive SCONE actuators at a 2 ms
   `fixed_control_step_size`. That is a moderate rewrite of the co-sim harness; the network
   itself (numpy) ports as-is. The alternative — re-expressing the network as
   ReflexController/Lua — is a rewrite, not a port, and loses the connectome tooling.
4. **Physics ≠ physics**: different contact law, different MT model, error-controlled vs
   fixed stepping ⇒ **every tuned gain is void**; a full re-run of the curriculum. Also
   irreproducibility across machines is already documented for SCONE (floating-point
   accumulation, skill gotcha) — the bit-exact discipline we have in MuJoCo does not
   transfer.
5. **No Simulink story**: our MuJoCo↔Simulink bridge investment (E1/E2, bit-exact sensor
   parity) is MuJoCo-specific; a Hyfydy port orphans it.

### 7.3 Verdict

For **standing-balance science** (does realistic contact damping + tendon elasticity change
balance strategy?) SCONE+Hyfydy is the better engine, and it is one license key away on this
machine. For **the dissertation-deadline SNS pipeline** it is strictly a cost. Recommended
sequencing: (1) now — convergence ladder in MuJoCo (§5.3); (2) post-deadline — get the key,
load stock Gait2392 in Hyfydy, run their Tutorial-3 balance scenarios to sanity-check the
backend, then port the SNS loop via SconePy and A/B the same balance task across engines;
(3) cite the cross-engine difference as a robustness check of the SNS conclusions, not as a
failure of MuJoCo.

---

## 8. Verified / not verified

**Verified this session (commands run, myo env, MuJoCo 2.3.7):**

- `WebFetch https://hyfydy.com/hyfydy-vs-mujoco/` — claims + 5 bottom citations extracted;
  no benchmarks on page.
- `WebFetch` MuJoCo stable docs (modeling.html; computation/index.html; solver-parameters
  section) — inelastic-tendon quotes, solref/solimp semantics, integrator statements.
- XML inspection of cvt3 + simbridge2 (script `goal3_inspect.py`): option lines, muscle
  default class, 92 named muscles, 19 pairs, tendon/equality inventory.
- Compiled-model checks (`goal3_env_check.py`): mujoco 2.3.7; integrator enums; cvt3 =
  Euler@5 ms, simbridge2 = implicitfast@2 ms; `vas_med_r` dyntype 3 / actnum 1 / Fmax
  1215.7 N; all tendon stiffness/damping = 0; pair solref/solimp/friction/condim = MuJoCo
  defaults; solver Newton, noslip 0.
- `findstr` on `mjmodel.h` (2.3.7 wheel): the four fixed-step integrators.
- Step3 pkl reads (`goal3_pkl_check*.py`): soleus/med_gas/tib_ant `tendon_sla` = 0.25/0.39/
  0.223 m, `fiber_opt`, `penna_opt`, `fmax` retained in conversion data.
- SEE toy experiment (`goal3_see_demo.py`): rigid baseline stable@2 ms; extra-DOF SEE
  stable@2 ms with 1 kg CE mass (16× force overshoot), unstable@2 ms with 0.1 kg
  (4.6×10⁸ N, MuJoCo instability warning), stable@0.5 ms with 0.1 kg; Python callback
  setters present (`set_mjcb_act_dyn` et al.).
- Contact-knob presence (`goal3_contact_knobs.py`): `pair_solreffriction` True,
  `geom_solreffriction` False, `opt.impratio` 1.0, `opt.density/viscosity` 0.0.
- Local-code reads: runner.py:166-465 (patch_xml), make_simbridge_xml.py:1-54,
  cvt3.xml:4/16/629-731; skills `mujoco` and `scone` (installed/verified 2026-09-23).

**Not verified / not run (and why):**

- **No quantitative Hyfydy-vs-MuJoCo benchmark** exists on the fetched page; I did not run
  Hyfydy (license key absent), so all Hyfydy capability statements rest on the vendor page
  + the SCONE skill's verified install notes, not on my own runs.
- **Convergence ladder (1 ms / 0.5 ms re-runs)** — recommended, not run (campaign-scale,
  out of analysis-only scope).
- **SCONE loading our Gait2392-family .osim** — skill says untested; I did not test it.
- **`mjcb_act_dyn` custom-muscle prototype** — callback availability verified, the actual
  elastic-tendon actuator not implemented (out of scope).
- Exact claim↔citation anchoring on the Hyfydy page (which ref supports which sentence) —
  the fetched text lists citations at the bottom only; the mapping in §1 beyond
  Hunt-Crossley→contact and Millard→musculotendon is my inference, labeled as such.
