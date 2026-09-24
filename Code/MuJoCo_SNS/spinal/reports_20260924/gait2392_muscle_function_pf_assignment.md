# Gait2392 muscle-function literature → PF-layer assignment design

**Date:** 2026-09-24 (EB475WS4) · **Task:** gait2392 muscle-function study
**Purpose:** turn the gait2392 muscle-function literature Ben pointed at into (1) a
per-muscle-group reference table, (2) a concrete design for Ben's requested
frontal-plane PF pair (hip adductors = agonist, hip abductors = antagonist),
and (3) a 3-DoF foot force-sensor design grounded in John et al. 2012.

**Claim provenance tags used throughout:**
- **[FETCHED]** — from a source fetched and read in this session (see §6 Sources for URL + locator).
- **[REPO]** — from a repo file read in this session (cited `path:line`).
- **[REPO-MEMORY]** — standing repo memory (AGENTS.md / DESIGN.md summaries) not re-verified line-by-line here.
- **[LIT-STD]** — standard gait-analysis knowledge (Perry 1992 / Liu 2008 school) NOT contained in the fetched sources; used only where the fetched sources are silent, and marked as such.

---

## 1. What the fetched sources establish

### 1.1 The model page (OpenSim Confluence, "Gait 2392 and 2354 Models") [FETCHED]

- Gait2392 = 3-D, 23-DoF model, **92 musculotendon actuators representing 76 muscles**
  (lower extremities + torso); Gait2354 is the reduced-muscle variant (patella removed by
  Seth; quadriceps insertions as moving points in the tibial frame). Default unscaled
  subject ≈ 1.8 m, 75.16 kg. Joints: hip ball-and-socket; knee = Yamaguchi–Zajac planar
  model; ankle / subtalar / MTP = frictionless revolutes (Delp 1990 axes; Inman 1976; MTP
  axis rotated −8°). Seven segments per leg: pelvis, femur, patella, tibia/fibula, talus,
  foot (calcaneus+…), toes.
- The page's inline PDFs are Confluence attachments (resolved via the public attachment
  API, pageId 53086215): the **"what muscles are included in each of the model"** PDF is
  the attachment **"Gait 2392 vs. Gait 2354.pdf"** (att53093029), and the **"maximum
  isometric muscle forces … along with the scale factors"** PDF is **"MuscleIsometricForces.pdf"**
  (att53092274). Both downloaded and text-extracted this session (local copies in `D:\temp\`).
- Full attachment inventory of the page (titles from the API response):
  `Anderson2001.pdf`, `Anderson&Pandy1999.pdf`, `Yamaguchi1989.pdf`, `Delp1990.pdf`,
  `Gait 2392 vs. Gait 2354.pdf`, `Gait2392ComparisonResultsCMC.pdf`,
  `MuscleIsometricForces.pdf`, `MuscleIsometricForces 2.pdf` (duplicate), 3 figure PNGs
  (body frames / knee geometry / ankle-subtalar-MTP axes). **Fetched: "Gait 2392 vs. Gait
  2354.pdf" and "MuscleIsometricForces.pdf". NOT fetched (inventoried only): the four
  journal PDFs, the CMC-comparison PDF, the duplicate, the PNGs.**
- The page also links the experimental-data paper — the same **John et al. 2012** (DOI
  10.1080/10255842.2011.627560) studied in §1.3 — and simtk.org model-distribution links
  (listed, not fetched).

### 1.2 Muscle inventory (92 actuators; "Gait 2392 vs. Gait 2354.pdf" + "MuscleIsometricForces.pdf") [FETCHED]

The appendix lists **46 actuators per side (r/l) = 92 total**: 43 Delp leg muscles plus
erector spinae, internal oblique, external oblique (the "Back+abs" group John 2012 names;
verified by re-parsing the extracted text: 92 rows with a Gait2392 X, 46 per side,
`ercspn_r/l`, `intobl_r/l`, `extobl_r/l` present). Cross-check:
`findstr /r "general.*name.*class=.muscle" gait2392_simbody_cvt3.xml | find /c /v ""` → **92**
on the converted MJCF the spinal runner loads (`runner.py:55`).

Per-muscle primary actions + Carhart-2000 (=Gait2392 column) max isometric forces, from
`MuscleIsometricForces.pdf` (selected; full text in `D:\temp\muscle_isometric_forces_text.txt`):

| muscle (base) | action (as printed in the PDF) | Fmax [N] |
|---|---|---|
| glut_med_1 / 2 / 3 | flex, abd, inrot / abd / ext, abd, exrot | 819 / 573 / 653 |
| glut_min_1 / 2 / 3 | flex, abd, inrot / abd / ext, abd, exrot | 270 / 285 / 323 |
| glut_max_1 / 2 / 3 | ext, abd / ext / ext | 573 / 819 / 552 |
| add_long / add_brev / pect | flex, ext(?), add / flex, add / flex, add | 627 / 429 / 266 |
| add_mag_1 / 2 / 3 | ext, add / ext, add / ext, add | 381 / 343 / 488 |
| quad_fem / gem / piriformis (peri) | exrot / exrot / abd, exrot | 381 / 164 / 444 |
| iliacus / psoas | flex, inrot / flex, inrot | 1073 / 1113 |
| rect_fem | h_flex, k_ext | 1169 |
| gracilis / sartorius / tfl | h_flex,h_add,k_flex / h_flex,h_abd,k_flex / h_flex,h_abd,h_inrot | 162 / 156 / 233 |
| semimem / semiten / bi_fem_lh | h_ext, h_add, k_flex (×3) | 1288 / 410 / 896 |
| bi_fem_sh | flex (knee) | 804 |
| vast_med / vast_int / vast_lat | ext (×3) | 1294 / 1365 / 1871 |
| med_gas / lat_gas | k_flex, a_pf | 1558 / 683 |
| soleus | pf | 3549 |
| tib_post | pf, inv | 1588 |
| flex_dig / flex_hal | pf, inv | 310 / 322 |
| per_brev / per_long | pf, ev | 435 / 943 |
| per_tert | **df, ev** | 180 |
| tib_ant | df, inv | 905 |
| ext_dig | df, ev | 512 |
| ext_hal | df, inv | 162 |

(The trunk muscles are not in this PDF — it covers the 43 Delp leg muscles per side.)

### 1.3 John et al. 2012 — muscles and mediolateral GRF [FETCHED]

John, C.T., Seth, A., Schwartz, M.H., Delp, S.L., "Contributions of muscles to
mediolateral ground reaction force over a range of walking speeds," *J Biomech*
45:2438–2443. Muscle-driven simulations (OpenSim CMC, 19-DoF / 92-actuator generic model
scaled per subject; 8 subjects × 4 speeds 0.54–1.56 m/s; data at
simtk.org/home/mspeedwalksims). Findings, with page locators:

1. **Three mediolateral (ML) windows** (p. 2439–2440, Fig. 1–2):
   - **Early stance 0–6 % GC:** peak GRF on the **leading** foot is **LATERAL**, grows with speed (p<0.05).
   - **Early single support 14–30 % GC:** peak GRF on the stance foot is **MEDIAL**, grows with speed (p<0.01).
   - **Late stance 40–56 % GC:** peak GRF on the **trailing** foot is **MEDIAL** (no speed effect on the 8-subject sample; in an 80-subject check the late-stance medial peak *decreases* with speed, Fig. 6).
   - Magnitudes are single-digit %BW: Fig. 2's axis spans −5…+10 %BW; Fig. 4's group-contribution axes span −10…+20 %BW.
2. **Muscles contribute ≥92 % of the ML GRF** on average over all speeds; gravity and
   velocity-related forces are small (p. 2440, Fig. 3).
3. **Abductors produce the largest MEDIAL contribution at all speeds in all three periods**
   (p. 2440, Fig. 4 A–C). Adductors/vasti/gastroc/soleus produce the opposing **LATERAL** contributions.
4. Window-by-window largest contributors (p. 2440, Fig. 4):
   - **0–6 %:** lateral = hamstrings + **contralateral** abductors (all speeds), plus vasti,
     adductors, gluteus maximus (free/fast, increasing with speed); medial = abductors
     (largest), dorsiflexors second (free/fast).
   - **14–30 %:** lateral = gastroc+soleus (slow), **soleus+vasti** (free/fast); medial =
     abductors, then dorsiflexors and **contralateral back+abs** (free/fast).
   - **40–56 %:** lateral = **gastrocnemius, soleus, adductors** at all speeds (adductor lateral contribution grows with speed).
5. **Double-support weight transfer** (p. 2441–2442, Fig. 5; trailing right → leading left):
   the leading leg's GRF flips **lateral → medial** and the trailing leg's flips **medial →
   lateral** before toe-off, driven mainly by the **rise of the leading abductors' medial
   contribution** and the **fall of the trailing abductors' medial contribution**, against
   the leading hamstrings/glut-max/vasti (early-stance support) and the trailing
   plantarflexors+adductors (late-stance lateral push).
6. Activation-timing anchors the paper itself cites: vasti active early-to-mid stance
   (Perry 1992); plantarflexors most active in the latter half of stance (Cappellini 2006);
   abductor activity increases with speed (Cappellini 2006) (p. 2441).
7. **Subtalar limitation, directly relevant to us:** the simulations ran **without subtalar
   motion**; tibialis-posterior/peroneal EMG was not recorded, so inverter/everter
   contributions are unmeasured; level-ground subtalar motion is typically <5° and the
   inverter/everter ML contributions are likely smaller than abductors'/gastroc/soleus' (p. 2442, citing Jenkyn 2010).
8. Group definitions used (Table 1, p. 2439): ABD = glut med/min (all compartments)+TFL;
   DF = tib_ant + ext_hal; HAMS = semimem+semiten+bflh+bfsh; VAS = 3 vasti;
   ADD = add_long+add_brev+add_mag(3)+pect+grac; GMAX = 3 compartments; GAS = med+lat
   gastroc; SOL; Back+abs = erector spinae + internal + external oblique.

---

## 2. Per-muscle-group table (BOTH legs; gait2392 × our SNS)

Pool naming is side-suffixed (`*_r` / `*_l`); everything below applies identically to both
legs — "our SNS target pool" column gives the exact actuator names per side from
`muscle_map.py` (all read this session). Network pools = 92 actuators − 6 pruned =
**86 MN pools** (40 leg + 3 trunk per side); the prune set is `PRUNE_MUSCLES` at
`runner.py:70-71`: **quad_fem_r/l, gem_r/l, peri_r/l** (matches the ask's "ignores 3
muscles per side").

Peak-use phases: Perry windows in brackets where used — IC (0–2 %), LR (2–12 %), MSt
(12–31 %), TSt (31–50 %), PSw (50–62 %), ISw/MSw/TSw (62–100 %).

| Group | Muscles per side (SNS pool, from `muscle_map.py:24-91`) | Primary action [FETCHED MuscleIsometricForces.pdf] | Gait phase of peak use | Our current SNS routing [REPO] |
|---|---|---|---|---|
| **Hip extensors** | glut_max1, glut_max2, glut_max3 (primary); **add_mag3** (primary hip_ext, secondary hip_add); ~~quad_fem, gem, peri~~ pruned | glut_max1 "ext, abd", 2/3 "ext"; add_mag1/2/3 "ext, add" | **Early stance**: GMAX among the large **lateral**-GRF contributors at 0–6 % (free/fast, ↑ with speed) [FETCHED John p.2440]; body-weight support in early stance during weight transfer [FETCHED p.2441]. Classic hip-ext torque peaks at IC–LR and TSt [LIT-STD] | group `hip_ext`; W_PF_MN E1 0.45 / E2 0.50 (`params.py:328-331`); POSTURE 0.22; BAL_DF route (`build_network.py:886-888`); IBEXC stance group; VEST sagittal target (`build_network.py:915-917`); in EXTENSOR_STANCE_GROUPS (`muscle_map.py:125`) |
| **Hip flexors** | iliacus, psoas, sar, tfl (primary); rect_fem (secondary, half-weight) | iliacus/psoas "flex, inrot"; sar "h_flex,h_abd,k_flex"; tfl "h_flex,h_abd,h_inrot"; rect_fem "h_flex,k_ext" | **Pre-swing → initial swing** (~50–75 % GC) for swing initiation [LIT-STD — the fetched sources do not time hip flexors; they are absent from John 2012's contributor table [FETCHED — Table 1, p.2439]] | group `hip_flex`; W_PF_MN F1 0.75 / F2 0.15 (`params.py:332-333`); POSTURE 0.10; BAL_PF route (`build_network.py:883-885`). NOTE: sar/tfl/grac knee-flex & hip-abd ride-alongs were deliberately REMOVED (2026-09-11 swing-leg splay; comments at `muscle_map.py:38-49,57-61`) — tension with the fetched action table (tfl IS an abductor there) |
| **Hip abductors** | glut_med1, glut_med2, glut_med3, glut_min1, glut_min2, glut_min3 | med/min1 "flex, abd, inrot"; 2 "abd"; 3 "ext, abd, exrot" | **All of stance, all three windows**: largest **MEDIAL** GRF contributor at 0–6, 14–30, 40–56 % at every speed [FETCHED John p.2440]; the decisive weight-transfer muscle (medial contribution rises on the leading leg, falls on the trailing one) [FETCHED p.2441-2442]; ↑ with speed [FETCHED p.2441] | group `hip_abd`; W_PF_MN E1 0.05 / E2 0.05 — near-absent from the PF drive (`params.py:328-331`); POSTURE 0.05; **BAL_LAT_R/L → MN_glut_*_** (frontal balance, `build_network.py:889-891`); in EXTENSOR_STANCE_GROUPS (Ib stance reversal); central heel/toe path (`build_network.py:534-541`). **No PF agonist pair — this is the gap Ben's assignment targets** |
| **Hip adductors** | add_long, add_brev, add_mag1, add_mag2, pect, grac (primary); add_mag3 (secondary) | add_long "flex, ext(?), add"; add_brev/pect "flex, add"; add_mag1/2/3 "ext, add"; grac "h_flex,h_add,k_flex" | **Late stance 40–56 %**: among the largest **LATERAL** contributors with gastroc+soleus at all speeds (↑ with speed) [FETCHED John p.2440, Fig.4C]; also lateral contributors at 0–6 % (free/fast) [FETCHED]; (classic adductor EMG: stance double peak IC–LR + PSw [LIT-STD]) | group `hip_add`; **W_PF_MN F1 hip_add = 0.0 — the ONLY explicit zero in the PF table** (`params.py:332`); no other PF phase lists hip_add; POSTURE 0.05 (`params.py:341-343`); NOT in EXTENSOR_STANCE_GROUPS; antagonist loops exist: ANTAGONIST hip_add↔hip_abd (`build_network.py:192-193`). **Essentially un-driven by the CPG today** |
| **Knee extensors** | vas_med, vas_int, vas_lat, rect_fem (primary) | vasti "ext" ×3; rect_fem "h_flex,k_ext" | **Early–mid stance**: lateral-GRF contributor at 0–6 % (↑ with speed) and, with soleus, the largest at 14–30 % (free/fast) [FETCHED John p.2440]; active early-mid stance per Perry as cited in the paper [FETCHED p.2441] | group `knee_ext`; W_PF_MN E1 0.10 / E2 0.15 / F2 0.0 (`params.py:328-333`); POSTURE 0.62 (strongest); **KINH swing suppression** `f1_kneext_inh` (`build_network.py:836-857`, v6 lever); IBEXC; VEST |
| **Knee flexors (hamstrings)** | semimem, semiten, bifemlh, bifemsh (primary); med_gas, lat_gas (secondary) | semimem/semiten/bflh "h_ext,h_add,k_flex"; bfsh "flex"; (gas "k_flex,a_pf") | **Early stance 0–6 %**: among the largest **LATERAL** contributors at all speeds [FETCHED John p.2440]; classic double burst TSw + IC–LR [LIT-STD] | group `knee_flex`; W_PF_MN **F1 1.80** (largest PF weight) / E1 0.05; hip_ext secondary at half-weight; POSTURE 0.08; IaIN antagonist machinery; VEST flexor list (`build_network.py:917-922`) |
| **Ankle plantarflexors** (triceps surae + deep PF + peronei, invertors/everters folded in) | soleus, med_gas, lat_gas, tib_post, flex_dig, flex_hal, per_brev, per_long, per_tert | soleus "pf" (3549 N — biggest in the model); gas "k_flex,a_pf"; tib_post/fl_dig/fl_hal "pf, inv"; per_brev/long "pf, ev"; **per_tert "df, ev" — NOTE: our map routes it ankle_pf, the fetched table calls it a dorsiflexor+everter (flagged §5)** | **Late stance**: GAS+SOL largest lateral at 14–30 % (slow) / SOL+VAS (free/fast); **GAS, SOL largest lateral at 40–56 %**; PF most active in latter half of stance (Cappellini, cited) [FETCHED John p.2440-2441] | group `ankle_pf`; W_PF_MN E2 0.25; POSTURE 0.30 with per-muscle POSTURE_OVERRIDE (soleus 0.55, tib_post 0.35, gas 0.08 …; `params.py:348-353`) + `ankle_post_walk_trim`; KINH `f1_anklepf_inh` (`build_network.py:861-868`); IBEXC |
| **Ankle dorsiflexors** | tib_ant, ext_dig, ext_hal | tib_ant "df, inv" (905 N); ext_dig "df, ev"; ext_hal "df, inv" (162 N; was a 1-N converter bug, fixed to stock — [REPO-MEMORY AGENTS 2026-09-13]) | **Early stance 0–6 %**: second-largest **MEDIAL** contributor at free/fast [FETCHED John p.2440]; same at 14–30 % (free/fast); classic swing activity + IC–LR lowering [LIT-STD] | group `ankle_df`; W_PF_MN E1 0.10 / F1 0.55 / F2 0.55 (`params.py:328-333`); POSTURE 0.15; BAL_DF route |
| **Subtalar (inverters/everters)** | **no dedicated group** — tib_post/fl_dig/fl_hal/tib_ant/ext_hal (inv) and per_brev/long (ev) all ride `ankle_pf`/`ankle_df` | (see actions above, "inv"/"ev" tags) [FETCHED] | John 2012 **excluded subtalar motion**; in/everter contributions unmeasured, likely < ABD/GAS/SOL on level ground; subtalar RoM in level walking <5° [FETCHED p.2442] | **No subtalar group exists in `muscle_map.py`** (`GROUPS` tuple, `muscle_map.py:19-21` — 10 groups, no subtalar). The MJCF model has subtalar joints with ligament-surrogate springs [REPO-MEMORY AGENTS 09-10/11]. A future frontal adductor/abductor pair (§3) is the natural precedent if subtalar pools are ever split out |
| **Trunk** | ercspn (trunk_ext); intobl, extobl (trunk_flex) | not in the fetched leg-muscle PDF; John 2012 groups them as "Back+abs" [FETCHED Table 1] | **Contralateral back+abs** among the next-largest **MEDIAL** contributors at 14–30 % (free/fast), increasing with speed [FETCHED John p.2440, Fig.4B] | groups `trunk_ext`/`trunk_flex`; W_PF_MN E1/E2 trunk_ext 0.20, F2 trunk_flex 0.05; POSTURE trunk_ext 0.30; BAL_TRK_EXT/BAL_TRK_FLX (`build_network.py:892-899`); VEST; trunk Fmax fix (1 N → stock 2500/900) applied 2026-09-13 [REPO-MEMORY AGENTS] |

**Reading of the table for the PF program** [synthesis of FETCHED + REPO]:
the groups that John 2012 identifies as the ML-GRF workhorses split cleanly into a
**medial pair (hip abductors)** and a **lateral set (adductors with vasti/gastroc/soleus)**.
Our PF drive today is strong on the sagittal set (knee_flex F1 1.80, hip_flex F1 0.75,
hip_ext E1/E2 0.45/0.50) but nearly silent on the frontal set: hip_abd gets 0.05 in E1/E2,
hip_add gets **0.0** everywhere. The frontal plane is currently handled only by posture
tone (0.05) + the reactive BAL_LAT feedback — no CPG pattern. That is the architectural
hole Ben's assignment closes, and John 2012 supplies the timing (ABD = all-stance medial,
ADD = late-stance lateral).

---

## 3. The frontal-plane PF pair — Ben's assignment made concrete

> Ben's request: *"assign hip adductors to one PF agonist pair and hip abductors to that
> PF layer's antagonist MNs"* — design the frontal-plane PF pair (adductor group = agonist
> HC, abductor group = antagonist HC, cross-inhibition, Renshaw per pool), consistent with
> existing `build_network.py` pool naming.

### 3.1 What exists to build on [REPO]

- **PF cells** (default phase mode): `PF_{E1,E2,F1,F2}_{side}`, each driven from
  `RG_E_/RG_F_{side}` via `rg_to_pf`, cross-inhibited through the laminated INs
  `PF_IN_E_{side}` / `PF_IN_F_{side}` — **no direct PF↔PF synapses** (house rule, comment
  at `build_network.py:508-511`).
- **Joint-layer mode** (`G["joint_pf"]`, 2026-09-18): `JPF_HCS = (HIP-E, HIP-F, KNEE-E,
  KNEE-F, ANK-E, ANK-F)`; per-muscle weights in `joint_pf_weights.json`; today
  `JPF_GROUP2HC` maps `hip_add → (HIP-E, HIP-F)` and `hip_abd → (HIP-E,)`
  (`build_network.py:163-167`) — i.e., in the joint-layer skeleton the adductors are
  currently split across the sagittal hip HCs and the abductors ride HIP-E only.
- **Antagonist bookkeeping already exists**: `ANTAGONIST["hip_add"] = ("hip_abd",)` and
  vice versa (`build_network.py:190-196`) — this drives Ia reciprocal (IaIN), IIIN, and
  IBIN↔IBIN mutual-inhibition wiring between adductor and abductor pools.
- **Renshaw per pool already exists**: when `G["renshaw"] > 0` every muscle gets `RC_{act}`
  with MN→RC exc (1.0), RC→MN inh (g_r), and **RC↔RC mutual inhibition for every ordered
  pair within a side** (`build_network.py:740-756`) — adductor-pool RCs vs abductor-pool
  RCs are wired by that generic loop with no new code.
- **MN membership per side** (`muscle_map.py:50-61, 25-31`):
  - agonist (hip_add primary, 6): `add_long, add_brev, add_mag1, add_mag2, pect, grac`
  - antagonist (hip_abd, 6): `glut_med1, glut_med2, glut_med3, glut_min1, glut_min2, glut_min3`
  - **excluded**: `add_mag3` (primary `hip_ext` — keep it on the sagittal path; its
    late-stance extensor burst should not acquire a frontal component), `sar`/`tfl`/`grac`
    ride-along secondaries (removed 2026-09-11 for the splay problem).

### 3.2 Design (recommended: default phase-cell mode, conditional topology, default OFF)

**New populations per side (+4; +8 network total, 410 → 418 when enabled):**

| population | class/tau | drive | role |
|---|---|---|---|
| `PF_ADD_{side}` | same PF cell as the others (`self._add(..., TAU["pf"]*tau_m)`), E2-shaped `(tau_m, adapt)` | `RG_E_{side}` via `rg_to_pf` | **agonist HC — adductor window** (late stance; John 2012 Fig. 4C: ADD lateral-GRF peak at 40–56 %) [FETCHED timing] |
| `PF_ABD_{side}` | PF cell, E1-shaped | `RG_E_{side}` via `rg_to_pf` | **antagonist HC — abductor window** (loading; John 2012 Fig. 4A + "at all speeds, in all periods" [FETCHED p.2440]) |
| `PF_IN_ADD_{side}` | IN, `TAU["pf"]` | ← `PF_ADD_{side}` exc (`pf_recip_inh`) | laminated suppressor of the abductor HC |
| `PF_IN_ABD_{side}` | IN, `TAU["pf"]` | ← `PF_ABD_{side}` exc (`pf_recip_inh`) | laminated suppressor of the adductor HC |

**Wiring (mirrors `_build_pf`, `build_network.py:493-553`):**
1. `RG_E_{side} → PF_ADD_{side}` and `RG_E_{side} → PF_ABD_{side}`, gain `G["rg_to_pf"]`.
2. Laminated cross-inhibition only: `PF_ADD → PF_IN_ADD (exc) → PF_ABD (inh)` and
   `PF_ABD → PF_IN_ABD (exc) → PF_ADD (inh)`, gain `G["pf_recip_inh"]`. **No direct
   ADD↔ABD synapse** — same rule as the E/F lamination.
3. **Both HCs hang off RG_E deliberately.** This is a staggered-stance pair, not an
   E/F pair: John 2012 places BOTH groups in stance (abductors all three windows;
   adductors early + late stance) [FETCHED]. Alternation comes from the PF_SHAPE-style
   `(tau, adapt)` stagger plus the cross-inhibition (exactly the E1/E2 pattern, which are
   also both RG_E-driven), and it removes frontal co-contraction — historically our
   splay/splay-back fights were ADD/ABD co-drive (`muscle_map.py` comments; adduction
   ±8° residual at drive 2.5 [REPO-MEMORY AGENTS 09-10/11 night]).
   *If Ben prefers a true bistable half-center (persistent-Na, self-alternating), build the
   pair from `NonSpikingNeuronWithPersistentSodiumChannel` (ships in sns-toolbox 1.5.2,
   Tutorial 8 — verified present) instead; but first re-check its tau_h semantics per the
   quenching finding in DESIGN.md. Default recommendation stays feed-forward PF cells.*
4. **MN routing** in `_wire_muscle` (new branch beside the existing PF block,
   `build_network.py:637-645`): when the pair exists —
   - `mi.groups[0] == "hip_add"` → `PF_ADD_{mi.side} → MN_{act}`, gain `G["pf_to_mn"] * W_FRONT_ADD` (per-muscle weights, e.g. all 1.0 initially; add_mag3 excluded);
   - `mi.groups[0] == "hip_abd"` → `PF_ABD_{mi.side} → MN_{act}`, same form.
5. **De-duplicate the old routes while the pair is on:** zero the `hip_abd` entries in
   W_PF_MN E1/E2 (0.05) and keep F1 `hip_add = 0.0` — single-source the frontal drive into
   the pair (same spirit as the E-sharpening edit of 2026-09-10 in `params.py:322-326`).
   Posture tone (0.05/0.05) stays — it is the standing-solve channel, not the gait pattern.
6. **Renshaw per pool:** nothing new to wire — `RC_{act}` pools and RC↔RC mutual
   inhibition already cover the 12 frontal MN pools when `G["renshaw"] > 0`
   (`build_network.py:743-756`); RC→IaIN recurrent disinhibition (`:778-780`) and the
   ANTAGONIST Ia/II/IB antagonist loops (`:764-810`) apply to adductor↔abductor pairs for
   free. Optionally later: let the IaIN phase gate for these pools use `PF_ABD_{side}`
   instead of `PF_F1_{side}` (`:772-776` defines the gate).
7. **Interaction with BAL_LAT:** BAL_LAT_R/L keep injecting into abductor MNs reactively
   (`build_network.py:889-891`); the pair supplies the feed-forward pattern and BAL_LAT
   the error correction — they sum at the MN, exactly like the trunk pattern/BAL_TRK pair.
   Expect the s3k BAL_LAT optimum to shift once the pair is tuned; retune, don't remove.

**Gating + contracts (house rules):**
- New top-level `G["front_pf"]` (default **0.0**) gates **population creation** — at 0 the
  four populations must not be built (conditional-topology rule: zero-g synapses alone
  change BLAS summation order and shift chaos-sensitive evals; v5 lesson, byte-identity at
  0 is the regression contract).
- **JSON RULE** (AGENTS.md): `front_pf`, any W_FRONT weights, and the PF_SHAPE entries for
  ADD/ABD must be saved into study jsons AND get `if key in best` branches in the
  `--best*` loader — the v7 renshaw omission cost a silent 0.14 kine delta; make the
  loader branch before the first study runs.
- Runner: no new input ports needed for the pair itself (PF cells take no external
  current); the sensory extension in §5 would add ports later, again default-off.

### 3.3 Joint-layer (`joint_pf`) variant

If the joint-layer mode becomes production first, add the pair there instead:
extend `JPF_HCS` with `("ADD-E", "ABD-E")` (both `endswith("E")` → both driven from RG_E
by the existing suffix logic at `build_network.py:562-567`), point
`JPF_GROUP2HC` to `hip_add → ("ADD-E",)`, `hip_abd → ("ABD-E",)` (replacing
`("HIP-E","HIP-F")` / `("HIP-E",)`), add the two laminated INs, and put the per-muscle
weights into `joint_pf_weights.json`. Same MN membership, same Renshaw coverage.
Open item either way: the s3k/s3b winners live on phase cells — port only after the
phase-mode pair is validated, or accept re-tuning on the joint-layer skeleton.

---

## 4. Why the literature backs this pair (one paragraph)

John 2012's induced-accelerations result is, in our vocabulary, a frontal **agonist /
antagonist decomposition**: the **abductors** are the dominant **medial**-GRF producers
throughout stance (all three windows, all speeds), while the **adductors** join
vasti/gastroc/soleus as **lateral**-GRF producers in early and late stance — and the
double-support weight transfer is executed precisely by modulating the leading vs trailing
**abductor** contributions against those lateral forces [FETCHED, §1.3 items 3–5]. A PF
layer that alternates an ADD window (late stance) against an ABD window (loading, all
stance) therefore encodes the documented ML weight-transfer strategy, gives the CPG a
feed-forward handle on the frontal plane that BAL_LAT (reactive-only) lacks, and targets
exactly the pools that are currently at/near zero in W_PF_MN. The 92 %-muscle share
finding also says our contact/balance gap will not be closed by passive dynamics tuning —
it needs muscle-pattern structure, which is this pair.

---

## 5. John 2012 → a 3-DoF foot force sensor

### 5.1 What a foot force sensor must measure, and who uses it when

Per foot (heel and toe regions kept separate — the existing region machinery already
splits `calcn` vs `toes`), three signed components in WORLD coordinates:

| component | gait-phase signature | who consumes it / for what |
|---|---|---|
| **Vertical (Fz)** | double peak at LR and TSt–PSw, mid-stance trough ~0.7–0.8 BW [LIT-STD]; repo anchor: subject01 replay measured peak vGRF **1.09/1.10 BW** and vertical impulse ≈ BW×duration [REPO `DESIGN.md:569-570` (audit_ik_ground)] | already the carrier: `LOAD_c` (Ib stance gating), `contact_onset` load/unload edges, per-side contact-reset phase machine (load > 20 N), duty metric (`runner.py:1377-1407, 1431-1435`) |
| **Mediolateral (Fy, signed)** | [FETCHED John 2012, §1.3 item 1] **early stance 0–6 %: LATERAL** peak on the leading foot; **early single support 14–30 %: MEDIAL** peak on the stance foot; **late stance 40–56 %: MEDIAL** on the trailing foot (→ flips lateral just before toe-off, Fig. 5); magnitude ~single-digit %BW | new: frontal-phase signal for the §3 pair (medial phase ≈ ABD window, late-stance medial→lateral flip ≈ ADD window / toe-off cue); support for BAL_LAT (frontal COM control is exactly the ML balance problem, MacKinnon–Winter per John 2012 p.2441); stance/swing symmetry diagnostics |
| **Fore-aft (Fx, signed)** | braking (negative) from IC through MSt, propulsion (positive) in TSt–PSw, zero crossing ≈ MSt/TSt boundary [LIT-STD — not in the fetched PDFs; John 2012 cites Liu 2008 for support/progression] | new: push-off detection (a cleaner stance-exit cue than pure unload), CoP×Fx lever arm as an ankle push-off estimator, speed-regulation feedback (net impulse ↔ velocity) |

Note the sensor does NOT need force-plate accuracy on ML/Fx: signals are ~10 % of Fz, but
the *sign and timing* carry the information (which window, which transition), and the
normalization constants should reflect the scale difference (below).

### 5.2 Extending the existing heel/toe mechanosensors

Today (`runner.py:1348-1425`): for every contact whose geom body maps to a foot region,
`mj_contactForce` fills a 6-vector and only **`con_force[0]` (normal)** is summed per
region → heel/toe normal loads; normalized by BW fractions (heel 0.35·BW, toe 0.50·BW,
load 0.60·BW, clamp 1.5), fed to network inputs `HEEL_c_{s}` / `TOE_c_{s}` / `LOAD_c_{s}`
(created at `build_network.py:435-447` as `HEEL_/TOE_/LBIN_{side}` IN populations), with
edge transients (`contact_onset`) and the contralateral swing kick (`contra_swing`).

Extension (design):
1. **Keep HEEL_c/TOE_c/LOAD_c exactly as they are** (normal force only) — they are tuned,
   regression-gated, and their npz schema (`contact` 2-channel log, `runner.py:1155-1160`)
   is consumed by figures/objectives. Old runs must stay comparable.
2. **Add signed world-frame channels per side**: `ML_c_r/l` and `FA_c_r/l` (foot total, and
   optionally per heel/toe region), normalized by much smaller constants (first guess:
   ~0.10·BW for ML, ~0.20·BW for Fx, same clamp-and-saturate pattern) — order-of-magnitude
   scale per §5.1 [FETCHED %BW framing for ML; LIT-STD for Fx].
3. **Network entry, conditional topology as always**: new IN populations (e.g.
   `MLIN_{side}`, `FAIN_{side}` analog of `HEEL_{side}`), default gains 0.0 → not built →
   bit-identical. Natural first targets: (a) `MLIN → PF_ABD/PF_ADD` of the §3 pair
   (sensory phase support matching John 2012's group timing); (b) `FAIN` negative-positive
   zero-crossing as a stance-exit cue into the existing contact-reset phase machine
   (`pm_*`, `runner.py:1426-1435`); (c) a ML asymmetry term in the balance objective
   (stage-4 "balance" curriculum already has a `contact_onset`-adjacent symmetry term).
4. **Sign audit before trusting any of it** (house practice: `audit_signs.py` /
   `_muscle_direction_test.py` precedent): verify in-world ML is positive-medial or
   positive-lateral by checking the 14–30 % window of a real walk is medial-directed on the
   stance foot [FETCHED John 2012], and that Fx brakes early / propels late. Our rig: world
   z = up; **mujoco −y = opensim +z (medio-lateral)** per the BAL_LAT comment
   (`build_network.py:875-878`); x = fore/aft.

### 5.3 MuJoCo-side recipe (DESIGN, not code)

Given MuJoCo 2.3.7 (our pinned version):
1. Loop contacts as today: `for ci in range(data.ncon)`, `c = data.contact[ci]`, region via
   `model.geom_bodyid[c.geom1 / c.geom2]` against the calcn/toes body sets
   (`runner.py:1356-1369`). Sum per region AND per side (linear — forces add).
2. `mujoco.mj_contactForce(model, data, ci, buf)` gives the force in the **contact frame**:
   `buf[0]` = normal, `buf[1]`, `buf[2]` = the two tangential (friction) components. Today
   only `buf[0]` is used.
3. The contact frame orientation is `c.frame` (9 floats, row-major; rows are the
   contact-frame axes expressed in world coordinates: `frame[0:3]` = normal,
   `frame[3:6]` = tangent-1, `frame[6:9]` = tangent-2). World force of that contact:
   `F_world = buf[0]*frame[0:3] + buf[1]*frame[3:6] + buf[2]*frame[6:9]`.
4. Accumulate `F_world` per foot; then: vertical = `F_world[2]`, fore/aft = `F_world[0]`,
   mediolateral = `F_world[1]` (world-frame roles per our rig, §5.2 item 4).
5. `buf[0]` is the constraint normal push (non-negative in the separating direction; the
   current `max(buf[0], 0)` stays harmless). Tangential signs are friction signs — do NOT
   clip them; rotate first, then interpret. (A one-off zero-ctrl press test — push the
   stance foot fore/aft against high friction and check the channel — is the cheap sign gate.)
6. **Capture-completeness check:** the XML contact-set surgery (`runner.py:422-437`)
   prunes to foot↔floor pairs, so per-contact summation over foot-geoms is the full GRF
   only if every floor contact involves calcn/toes bodies — assert `ncon` outside the foot
   sets stays zero in ground mode (rig contact lives on pelvis springs; those are wanted
   separately, not in the foot channels).
7. **Cross-checks:** (a) per-side |F_world| total should match the existing 2-channel
   `log_contact` (same numbers, more channels); (b) whole-body option as a one-off
   validator, not a runtime channel: `data.cfrc_int` per foot body (com-based internal
   force) or `data.efc_force` projected on contact rows — use only to reconcile, because
   they mix passive/inertial terms.
8. **npz schema:** extend, don't mutate — add e.g. `contact3` (nsteps × 6: r/l ×
   [Fz, Fx, Fy]); keep `contact` untouched so v6/s3k repro numbers and figures keep loading.
9. Cost: O(ncon) vector ops per step at 2 ms — same order as the existing loop; no solver
   changes; works with the 2.3.7/implicitfast stack and does not interact with
   `--contact-damp` variants beyond the forces they produce.

---

## 6. Sources

### Fetched this session (2026-09-24)

1. **OpenSim Confluence — "Gait 2392 and 2354 Models"**
   `https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53086215/Gait+2392+and+2354+Models`
   — fetched 3 ways (WebFetch summary; raw HTML via curl → `D:\temp\gait2392_page.html`;
   web-reader markdown). Section refs: Overview and Authors; Kinematics → Joint
   geometry / Muscle geometry; Dynamics → Actuators (Peak isometric force).
2. **Attachment "Gait 2392 vs. Gait 2354.pdf"** (att53093029; the page's "what muscles are
   included in each of the model" PDF)
   `https://opensimconfluence.atlassian.net/wiki/rest/api/content/53086215/child/attachment/att53093029/download`
   — 2 pages; the 92-muscle abbreviation list by model. Local: `D:\temp\gait2392_vs_2354.pdf`.
   (Quirks noted: "Gluteus Medius 2, Right glut_med3_r" label typo row; "fixme gem" row.)
3. **Attachment "MuscleIsometricForces.pdf"** (att53092274; the max-isometric-forces +
   scale-factor PDF; columns = joint/type, muscle, Gait2392 [Carhart 2000], Delp 1990,
   scale factors)
   `https://opensimconfluence.atlassian.net/wiki/rest/api/content/53086215/child/attachment/att53092274/download`
   — 1 page. Local: `D:\temp\muscle_isometric_forces.pdf`.
4. **John et al. 2012**, `https://nmbl.stanford.edu/publications/pdf/John2012.pdf`
   — J Biomech 45:2438–2443; 6 pages, full text extracted to `D:\temp\john2012_text.txt`.
   Locator style used above: "p.2439" etc. = journal page numbers (2438 = p.1).
5. **Confluence attachment-listing API** (to resolve the two inline PDFs):
   `https://opensimconfluence.atlassian.net/wiki/rest/api/content/53086215/child/attachment?limit=50`.
   Page inventory (NOT fetched): Anderson2001.pdf, Anderson&Pandy1999.pdf, Yamaguchi1989.pdf,
   Delp1990.pdf, Gait2392ComparisonResultsCMC.pdf, "MuscleIsometricForces 2.pdf"
   (duplicate of #3), 3 figure PNGs (Delp-1990 body frames / knee geometry /
   ankle-subtalar-MTP axes). Model-distribution links (NOT fetched):
   `https://simtk.org/frs/download.php?file_id=3857`, `https://simtk.org/frs/index.php?group_id=91`.

### Repo files read this session

- `Code\MuJoCo_SNS\spinal\muscle_map.py` (whole file) — group definitions, per-muscle
  membership, EXTENSOR_STANCE_GROUPS, splay-removal comments.
- `Code\MuJoCo_SNS\spinal\build_network.py` — PF build `:493-553`; joint-layer build
  `:555-597` + `JPF_*` tables `:156-196`; MN wiring `:618-684`; Renshaw `:740-756`; IaIN
  `:758-780`; heel/toe ports `:420-459`; BAL routing `:871-899`; VEST `:915-926`.
- `Code\MuJoCo_SNS\spinal\runner.py` — PRUNE_MUSCLES `:70-71`; heel/toe sensor block
  `:1348-1425`; phase machine `:1426-1435`; contact log `:1155-1160`; model path `:55`;
  contact-set surgery `:422-437`.
- `Code\MuJoCo_SNS\spinal\params.py` — W_PF_MN `:327-334`; W_POSTURE `:341-343`;
  POSTURE_OVERRIDE `:348-353`; PF_SHAPE `:358`.
- `Code\MuJoCo_SNS\spinal\DESIGN.md` — Architecture table `:1928-1947`; joint-layer /
  IK-audit section `:564-623` (peak vGRF 1.09/1.10 BW at `:569-570`).
- MJCF actuator count check: `Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml`
  → 92 `class="muscle"` general actuators (command in §1.2).

### Repo memory relied on but not re-verified here [REPO-MEMORY]

- AGENTS.md: prune list provenance (quad_fem/gem/peri = Ben's list), trunk Fmax fix
  (ercspn 2500 / intobl 900 / extobl 900 / ext_hal 162), 09-10/11 splay/adduction history,
  subtalar ligament-surrogate springs, JSON RULE, conditional-topology/bit-identity
  contract, world-remap conventions (mujoco −y = opensim +z).
- DESIGN.md 2026-09-16 toolbox correction (persistent-Na class exists; tau_h caveat) —
  cited in §3.2 item 3 as the caution for a NaP-based variant of the pair.

### Standard-literature claims used without a fetched source [LIT-STD]

- Perry-window phase boundaries and swing-phase roles (hip flexors PSw–ISw, hamstrings
  TSw+IC–LR, tib_ant swing+IC–LR, vertical/fore-aft GRF double-peak + braking/propulsion
  shapes). Sources: Perry 1992; Liu 2008 — both cited *within* John 2012 but their
  specific timing tables were not fetched this session. Every such cell in §2/§5 is tagged.

### Not done / open

- The four journal PDFs on the page (Delp 1990, Yamaguchi 1989, Anderson & Pandy 1999,
  Anderson 2001) and Gait2392ComparisonResultsCMC.pdf were inventoried but not fetched;
  per-muscle *moment-arm* timing (vs. GRF-contribution timing used here) would come from
  Delp 1990 / Arnold-Asakawa-style work — fetch on request.
- John 2012's supplementary figures (S1–S4) live behind the DOI, not fetched.
- The §3 pair is a design only — no build_network/params/runner edits were made, no code
  was run against the network, and no study was launched (per the ask: DESIGN, not code).
