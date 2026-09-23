# Shevtsova replication draft — FULL EXTRACTION 2026-09-23 (ZCode)

Source: Shevtsova, Lockhart, Rybak, Magnuson & Danner, "Linking spinal
circuit reorganization to recovery after thoracic spinal cord injury",
eLife 14:RP107480 (DOI 10.7554/eLife.107480; the local cache prints
"eLife 2025;14:RP107480"; the project refers to it as Shevtsova 2026).
Local full text: `spinal\shevtsova_2026_fulltext.txt` (123,725 chars,
from Ben's Zotero storage 6LM33WX6 .zotero-ft-cache). Everything below
with a Table/Figure number was read from that full text THIS session;
the first-pass keyword-count section that used to live here is superseded.

FIRST PASS by ZCode — Ben edits before build. Machine-readable edge list:
`shevtsova_rules.json` (w2l-style nodes/synapses schema).

## 1. What the model IS (verified)

The RAT adaptation of the laminar quadruped locomotor circuit lineage
(Danner 2017 eLife 6:e31050; Zhang 2022), subjected to lateral hemisection
and midline contusion (paper Abstract + "Computational model" sections).
It is a **gait-expression model at the RG/CIN/LPN level ONLY**: persistent-
Na RG half-centers + commissural + long-propriospinal populations, brainstem
drive ramps, no PF layer, no motoneurons, no muscles, no afferents (checked:
the full text's afferent/sensory mentions all refer to experiment, not the
model; "muscle" appears only in in-vivo data discussion). For our bipedal
replication the LUMBAR half is the relevant circuit; cervical copies,
Sh2/dLPNi/dV0D/dV2a/dV0V/aV3/aSh2/dSh2 inter-girdle LPNs are
quadruped-specific (adaptation decision for Ben, flagged in the JSON).

## 2. Neuron populations (Fig 2 schematic + Model description)

Per limb (cervical ×2 + lumbar ×2), one hemisphere each:
RG-F, RG-E (NaP half-centers); IniF, IniE (ipsilateral mutual-inhibition
INs); V2a; V0V; V0D; V3-E; V3-F (lumbar) / V3-E (cervical); Ini (lumbar
and cervical); InE1 (lumbar); Sh2; plus inter-girdle dLPNi, dV0D, dV2a,
dV0V, aV3, aSh2, dSh2. Population sizes: 200 neurons per RG center, 50
per other population in the lineage models (that number is from the 2015
review's Fig 7 caption for Shevtsova 2015 — NOT re-verified for 2026).

## 3. Full connection table (Table 1 "Connection weights in the intact
## model", exact values; i- ipsilateral, c- contralateral, f- fore, h- hind)

**Within cervical and lumbar circuits** (the shared core — our biped
reduction uses the lumbar instance):

| source | target(s) | weight |
|---|---|---|
| RG-F | i-InF | 0.4 |
| RG-F | i-V0D | 0.7 |
| RG-F | i-V2a | 1 |
| RG-E | i-InE | 0.4 |
| RG-E | i-V3-E | 0.35 |
| RG-E | i-Sh2 | 0.5 |
| IniF | i-RG-E | −1 |
| IniE | i-RG-F | −0.1 |
| V2a | i-V0V | 1 |
| V0V | c-Ini | 0.6 |
| V0D | c-RG-F | −0.07 |
| V3-E | c-RG-E | 0.02 |

**Within lumbar circuits** (added on top of the shared core):

| source | target | weight |
|---|---|---|
| RG-F | i-V3-F | 0.4 |
| RG-F | i-aV3 | 0.3 |
| Ini | i-RG-F | −0.075 |
| V3-F | c-RG-F | 0.03 |
| V3-E | c-InE1 | 1 |
| InE1 | c-RG-E | −0.045 |

**Within cervical circuits** (quad-specific): RG-F → i-dLPNi (0.7),
i-dV0D (0.5), i-dV2a (0.5); Ini → i-RG-F (−0.0375); dV2a → i-dV0V (0.9).

**Between cervical and lumbar** (quad-specific): dSh2 → ih-RG-F (0.005);
aSh2 → if-RG-F (0.04); dLPNi → ih-RG-F (−0.01); dV0D → ch-RG-F (−0.075);
dV0V → ch-RG-F (0.02); aV3 → cf-RG-F (0.065).

**Injury reconfiguration** (Table 2 hemisection): crossed LPN weights
×0.4 (e.g. hl-aV3 0.065→0.026, fr-dSh2 0.005→0.002, fl-dV0D −0.075→−0.03);
ipsilesional lumbar RG-F drive 0.1→0.009; both V0V inhibitory drives
halved. (Table 3 contusion analogous, ×0.05 scale-downs + bias currents
on h-V0D/h-V0V.) These are the paper's plasticity levers, not connectivity
rewiring — relevant if Ben ever wants an SCI-robustness experiment.

## 4. Dynamics notes (Methods "Model parameters" section, exact)

- Activity-based populations: C·dV/dt = −INaP − IL − ISynE − ISynI − INoise
  (eq 1, RG centers only have INaP); other populations eq 2 (no INaP).
- Parameters: C = 10 pF; gL = 4.5 nS (RG) / 2.8 nS (others); ḡNaP = 4.5 nS;
  gSynE = gSynI = 10 nS; EL = −62.5 mV (RG) / −60 mV; ENa = 50 mV;
  ESynE = −10 mV; ESynI = −75 mV; Vthr = −50 mV; Vmax = 0 mV;
  V1/2,m = −40 mV, km = −6 mV; V1/2,h = −45 mV, kh = 4 mV;
  τh(V) = τ0 + (τmax−τ0)/cosh((V−V1/2,τ)/kτ) with τmax = 400 ms,
  τ0 = 150 ms, V1/2,τ = −35 mV (**kτ value not present in my extraction —
  not verified; check the PDF**). Noise = OU process (eq 14, σNoise 1.1 pA
  for variability runs); integration Cash-Karp RK, dt = 1 ms.
- Output f(V) piecewise-linear Vthr→Vmax (eq 12) = our clip(V/5 mV,0,1)
  family; weights w>0 exc / w<0 inh rectified by S(x) (eq 11).
- Drive: D = d·α + b per population (eq 13). Intact values: extensor
  centers dE=0, bE=0.1; **flexor centers dE=0.1, bE=0** (the FLEXOR center
  carries the speed drive — frequency rises by SHORTENING the extensor
  phase; contrast with our DRIVE→E bias, DESIGN.md 2026-09-12/13 note);
  homologous V0D dI=0.75, diagonal V0D dI=1.5, V0V dI=0.25 (fore) / 0.15
  (hind); inhibitory drives to V0D/V0V/dV0D are the gait-switch levers.
- **tau_h TRAP (project-known)**: their τh(V) bell keeps 150–400 ms — an
  order of magnitude above the sns_toolbox schedule that collapses to
  ~0.1 ms and quenches half-centers. Our RG runs the toolbox NaP class
  with FIXED τh = 350 ms (`params.py:34`) — inside their range; keep fixed.

## 5. Mapping onto our SNS classes

| Shevtsova 2026 element | our element | status |
|---|---|---|
| RG-F/RG-E NaP centers | RG_E/RG_F toolbox NaP, fixed τh | EXISTS |
| IniF/IniE laminated RG inhibition | RG-E→InE→RG-F + mirror (g 4.0/4.0) | EXISTS (weights differ: theirs 0.4 → −1/−0.1 asymmetric) |
| RG-F→i-V0D (0.7) → c-RG-F (−0.07) | RG-F→CIN_F (4.0) → contra RG-F (4.0) | EXISTS (our crossed gain is ~50× theirs in relative terms; they tune weak crossed inhibition) |
| V2a→V0V→c-Ini→i-RG-F (1 / 0.6 / −0.075) chain | — no V0V/V2a chain | NEW (disynaptic crossed excitation→inhibition onto the ipsi RG-F via contra-driven Ini; note OUR CIN_E chain is the V3-E flavor) |
| V3-E→c-RG-E (+0.02) and V3-E→c-InE1 (1)→c-RG-E (−0.045) | RG-E→CIN_E (1.2)→contra InE (1.2) | EXISTS for the InE1 variant (v3_gain>0, default 0); the DIRECT V3-E→c-RG-E edge we deliberately did not build (Ben's figure reading, DESIGN.md 2026-09-16) |
| V3-F→c-RG-F (+0.03) with RG-F→i-V3-F (0.4) | — | **NEW — the crossed flexor-excitatory edge** (also Rybak 2015 Fig 7A CINe-F) |
| Ini (driven by c-V0V 0.6) inhibiting i-RG-F (−0.075) | — | NEW (no interneuron-targeting crossed input except CIN_E→contra InE) |
| Sh2/LPN machinery | — | N/A for biped (excluded in the JSON with citation) |
| brainstem α ramp | DRIVE schedule | EXISTS (runner input schedule) |

## 6. EXISTS vs NEW (cross-checked against the compiled net THIS session)

Ground truth: `_net_edges.py` default run (82 edges, exit 0) + all-
conditionals probe (150 edges, 68 neurons, exit 0; myo env). See
rybak_draft.md §6 for the commands.

- EXISTS: NaP RG; InE/InF lamination; CIN_F crossed flexor inhibition;
  CIN_E→contra InE (= their V3-E→InE1 chain) conditional on v3_gain;
  CIN_E→IBEXC crossed extensor support (v3_to_ibexc, default 0) — the
  paper has no IBEXC analog (no motor level at all).
- NEW (not in our net at any gain): (1) **V3-F crossed flexor EXCITATION**
  (RG-F→i-V3-F 0.4, V3-F→c-RG-F 0.03); (2) the **V2a/V0V chain** — crossed
  signal that lands on an INI inhibiting the ipsilateral RG-F (their
  alternation-at-speed mechanism); (3) per-population d/b drive-slope
  allocation (our DRIVE is two scalars; theirs assigns slopes to
  flexor centers and inhibitory drives to CIN classes — a tuning-frame
  difference, not an edge).

## 7. Which Shevtsova circuits address our known gaps

- **Frozen left leg / crossed flexor-excitatory edge**: V3-F→c-RG-F (+0.03)
  is the paper's own weak crossed flexor excitation; combined with V0D
  (−0.07) it forms their speed-dependent alternation pair. For us the
  candidate is a CIN_F_E-class excitatory crossed edge (or re-using
  v3_to_ibexc's pattern flexor-side). Sign/weight: they keep excitation
  ~2.3× weaker than the matched inhibition.
- **Sensory phase reset into the RG**: NOT addressed by this paper (no
  afferents in the model) — see rybak_draft.md §7 and shinohara_draft.md.
- **IaIN↔IaIN / Ib-IN mutual inhibition**: NOT in this paper (no motor
  level) — Rybak 2006a Table 2 / Di Russo rule 3 respectively.

## 8. Replication plan (kept from the first pass, now grounded)

1. EXACT SOURCE: eLife publishes code — pull the model from the paper's
   Data availability (link NOT yet extracted from the fulltext — next
   pass; the Danner 2017 / Zhang 2022 model files are the fallback).
2. Biped reduction: lumbar instance only; drop cervical + inter-girdle
   LPNs (decision for Ben, default in shevtsova_rules.json = dropped,
   each dropped edge carries its citation).
3. Build via the connectome editor path (CONNECTOME.md) after Ben's edits.
4. SUCCESS CRITERIA (from the paper's own reported behavior):
   - [ ] intact: speed-dependent stepping with walk→trot transition under
     the α ramp (their Fig 3-4 family; exact figure numbers not yet
     extracted — not verified);
   - [ ] half-center rhythm survives deafferentation trivially (no
     afferents exist — the original README criterion is vacuous for this
     paper; replaced by the α-ramp criterion);
   - [ ] hemisection reconfiguration (Table 2 levers) reproduces their
     asymmetry direction.

## Verified / not verified

- Verified (read from `shevtsova_2026_fulltext.txt` this session):
  Table 1/2/3 weights; Methods neuron equations 1-14; parameter list;
  drive allocation; Fig 2 caption; population inventory; absence of
  PF/MN/afferent machinery.
- Not verified: kτ of the τh bell; the model-code link (Data availability
  section not yet read); which exact result figures show the intact
  walk→trot (not extracted); neuron counts for the 2026 model (the
  200/50 numbers are the 2015 lineage's); anything about "Poirazi/Smith"
  authorship in the old DESIGN.md citation line — the fulltext header
  lists Shevtsova/Lockhart/Rybak/Magnuson/Danner only.
