# BW Tier 2 — deviation table (bilateralrg2 source-sign / conductance deviations)

Ben, 2026-10-24. The prompt's Tier-2 step: "list the template's documented
source-sign and conductance-mapping deviations and close each: either it
carries a Ben ruling (the sign stays) or it is a transcription bug, fixed in
the template now." Sources: `w2l_cpg/make_bilateralrg2.py` docstring (the
7-point GUID-level delta), `w2l_cpg/build_w2l_net.py` W2L_GAINS comments,
DESIGN.md 2026-10-20 open issues 1-2, `reports_20261020/w2l_contact2/README.md`.

| # | Deviation | Status / disposition | Evidence |
|---|---|---|---|
| 1 | "R Hip MN flx RE → R Hip MN ext RE" edge carries the V3-Commissural-EXCITE synapse TYPE (GUI clone-type artifact) — animates as excitation; Ben's 2026-10-03 export rules it **inh / rc_mutual** | **CLOSED as RULED — sign stays INH** (not a transcription bug). Remains flagged for Ben's GUI confirmation (DESIGN open issue 1; the .aproj's own type is his to retype in a GUI pass — out of scope here per the batch text). | make_bilateralrg2.py docstring item 7; DESIGN.md 2026-10-20 open issue 1 |
| 2 | Heel-contact cascade + the crossed heel→contra-RG-flx edge read at the **WEAK class** (dedicated Heel_* types, SynAmp 0.01 / 1e-8 S) rather than the strong generic class | **CLOSED as RULED — type-level effective strength** (the .pl-chain doctrine: the synapse TYPE's conductance IS the effective strength; Ben built dedicated types precisely to set them). Effective-strength reading = type-level, confirmed. | make_bilateralrg2.py items 2-3; DESIGN.md open issue 2; build_w2l_net.py W2L_GAINS `contact_hc`/`contact_contra_flex` comments |
| 3 | R-knee-extensor Renshaw node carried a **stale "L" label** (GUID a3828afc actually serves the R side; Ben's 2026-10-03 export renamed it RC_8971) | **CLOSED as FIXED (label-only)** — the miner relabels it "R Knee MN ext RE"; every edge keeps its GUID-resolved endpoint. | make_bilateralrg2.py item 6 |
| 4 | Toe-contact cascade trimmed to toe → ONLY ipsilateral ankle-MN ext (old strong cascade removed); heel → ipsi Ank MN ext kept at the GENERIC 5e-7 class | **CLOSED as TRANSCRIPTION** — faithful to the corrected aproj (items 1, 4 of the miner delta); classes mapped contact_ank=1.0 (generic) vs contact_hc=0.08 (weak), order preserved. | make_bilateralrg2.py items 1, 4; W2L_GAINS comment block "2026-10-18 corrected aproj" |
| 5 | AnimatLab contact "gain 20" is the ADAPTER current gain (nA/contact into the contact neuron), not a synapse conductance | **CLOSED as DOCUMENTED MAPPING** — adapter value applied runner-side (CONTACT_CURRENT_NA / ctn_amp knob); synapses carry class conductances. | W2L_GAINS `contact` comment; test_w2l_ground contact-encoder doc |
| 6 | Core rhythm/reflex classes (rg_laminate 4.0, rg_direct_exc 0.5, pf_drive 2.4, pf_cross 4.0, pf_to_mn 2.0, Renshaw pair 1.0/0.5, Ia chain 0.6/0.4/0.5, aff_HC_excite 0.08) were HARDCODED (no knobs) | **CLOSED this batch** — exposed as synapse-class SCALE knobs `s_*` (default 1.0 = anchor = bit-identical build, proven by signature comparison 2026-10-24); magnitudes now searchable, signs/topology fixed. | test_w2l_ground.py KV s_* + gains.update block; this folder's campaign |
| 7 | Commissural classes c1 3.0 / V3 0.08 (Shinohara ratio) | **CLOSED as DOCUMENTED** (c1 period-locking measured 2026-09-24; anti-latch V3 weakness) — already knobbed via `comm` (3.0·comm / 0.08·comm), carried into this study. | W2L_GAINS c1_inh/v3_weak comments |

No open transcription bugs found: the miner was GUID-level verified against
the aproj (delta report at generation time), and both sign-level deviations
(items 1-2) carry explicit Ben rulings. Out-of-scope flags (unchanged, for
Ben's eventual GUI pass): the .aproj's V3-clone synapse typing on the item-1
edge, and the weak-class effective-strength confirmation (DESIGN open
issues 1-2).
