# Vision review — spinal circuit figures (2026-09-23)

`/circuit-vision` Mode B run (audit our own figures), triggered by Ben.
Scope: the circuit diagrams in `spinal\figures\`. No figure, tex, or circuit was
edited — findings only (any regeneration needs the M2 supervisor gate + Ben's OK).

## Method + tooling incident

- `mcp__zai-mcp-server__understand_technical_diagram` TIMED OUT (60 s) on every
  attempt, including a 44 KB image — the server is effectively unavailable today.
  Fallback used: 4.5v vision model (`mcp__4_5v__analyze_image`) on downscaled PNGs
  (max side 2200 px, originals untouched; downscaled copies in
  `%TEMP%\cv_figs`). Every vision claim was then adjudicated against ground
  truth: the figure generators' source + `_net_edges.py` + the alt-text files.
- Vision models misread small type at this scale; several dramatic "errors" were
  hallucinations. Only source-verified findings are reported as real below.

Images analyzed:
- `spinal\figures\circuit_rules.png` (drawn 2026-09-22, `draw_circuit_rules.py`)
- `spinal\figures\circuit_dengstyle.png` (2026-09-17, `draw_circuit.py` deng mode)
- `spinal\figures\circuit_literature.png` (2026-09-17, `draw_literature_circuit.py`)
- `spinal\figures\sns_diagram_panels.png` (+ `sns_layer_{rg,pf,motor}.png`,
  toolbox-native renderer, 2026-09-17)
- `curr3i_channels.png` / gait figs are trace/data plots, not circuits — only the
  alt-text consistency of `curr3i_channels` was checked (see A3).

## Ground truth

1. `_net_edges.py` dump (representative 4-muscle build, v5/v6 gains forced):
   52 neurons, 82 connections. Key groups: DRIVE→RG-E/RG-F exc; POSTURE→RG-E exc
   + POSTURE→MN exc; RG-E→InE exc, InE→RG-F inh (laminated, mirrored for F);
   RG-E/RG-F→PF exc; PF→PF mixed x8; PF→MN exc; Ia-aff→MN exc(g=0.6) +
   Ia-aff→MN inh(g=0.4, direct reciprocal at ia_in=0); II-aff→MN exc(g=0.4);
   Ib-aff→MN inh(g=0.35) + Ib-aff→IBEXC_IN exc; RG-E→IBEXC_IN exc(=1, stance
   gate); IBEXC_IN→MN exc; KINH_IN→MN inh(g=0.59) gated by PF→KINH_IN exc;
   RG-F→CIN_F exc + CIN_F→contra RG-F inh.
   **Caveat (F8):** at default gains the dump does NOT contain conditional
   topology — CIN_E/V3 crossed-excitatory relay, IaIN, LBIN, Renshaw, HEEL/TOE
   chains. Production (`connectome_gains.json`, c1_gain/v3_gain etc. > 0 in all
   s3g+ winners) enables them; `draw_circuit.py` deng mode forces the tuned
   gains and asserts every drawn group. Cross-checks for those classes rely on
   the generator contract, not this default dump.
2. `draw_circuit_rules.py` (complete source read): panels A–F, markers 1–3,
   4-line captions, 7.5×10 in, Arial 10–12 pt, no italics, Tol palette —
   complies with the 2026-09-21 project figure standards (it was built after
   them, v4 layout, "three visual-judge passes").
3. Alt texts: `circuit_rules_alt.txt`, `curr3i_channels.alt.txt`.

## Adjudicated vision claims (dismissed)

| # | Fig | Vision claim | Ground truth (source) | Verdict |
|---|-----|--------------|------------------------|---------|
| 1 | rules A | POSTURE→RG-F edge drawn | No POSTURE wire at all — deliberately omitted (source line 169) | vision error (but see F2) |
| 2 | rules A | RG-E→InF / RG-F→InE (crossed) | Own-side: RG-E→InE, RG-F→InF (lines 164–165; = inventory) | misread |
| 3 | rules A | II→MN + Ib→MN drawn excitatory | No per-afferent MN edges in A; one muscle→afferent-row edge (line 183) | hallucinated |
| 4 | rules A | c1→V3 inhibitory | c1→contra-RG inh; V3→contra-RG exc (lines 188–191; = CIN_F inventory) | misread |
| 5 | rules B | IaIN ag→MN ant edge missing | Drawn — long diagonal (line 217) | misread |
| 6 | rules B/C/F | marker "4" referenced, no markers on some panels | Markers: B②(line 213), C①②(246,250), D①(267), F③(324); no marker 4 exists | misread (but see F4) |
| 7 | rules D | Ib→IBEXC sign wrong (should be inh?) | Excitatory is CORRECT (line 269; inventory Ib-aff→IBEXC exc g=0.5) | vision error |
| 8 | rules D | direct Ib→RG-E line | Arrow is the RG-E→IBEXC stance gate, direction read backwards (line 272; inventory g=1) | misread (but see F5) |
| 9 | literature | duplicate formula "a=clip(II/3,0,1)" = copy-paste bug | Actual label "a = clip(V/5 mV, 0, 1)" at 5.1 pt, one per MN column BY DESIGN (lines 453, 457) | misread (but see F1) |
| 10 | deng | nodes "YEEK/TOE", "VRS", "TK cell" | HEEL/TOE chains, V3 (+V0D/V0c strip), tag text — all present in source, edges asserted vs compiled net | misreads |
| 11 | literature S1 | IbIN-E line colored green though inhibitory | Line color = source tint; SIGN is carried by the glyph (filled dot), per Ben's convention | by-design, legibility note |

## Genuine findings (source-verified)

- **F1 — pre-standards type sizes + italics (deng + literature figures).**
  `draw_literature_circuit.py` renders nearly everything at 4.6–7.3 pt and uses
  `style="italic"` (line 445); `draw_circuit.py` deng mode tags at 5.0–6.4 pt
  (HEEL/TOE fs=5.4, KINH fs=5.2, footer 6.4). Project standard (2026-09-21):
  ≥10 pt Arial, no italics. Both figures are 09-17 vintage — they predate the
  rule; `circuit_rules` complies. The 5.1 pt conversion-map label is exactly
  what caused vision misread #9.
- **F2 — rules panel A: POSTURE port drawn with no wire, no in-figure
  explanation.** Source comment says "omitted for clarity (feeds all MNs; see
  caption)" but the caption does NOT mention POSTURE; only the alt text does.
  A cold reader sees a dead-ended input port.
- **F3 — rules panel B: PF F1 gate drawn but unexplained in-figure.** The
  element is labeled, but its role (phase-gating the IaIN pathway) appears only
  in the alt text, not the 4-line caption.
- **F4 — marker key is color-only; marker ③ caption text is about something
  else.** Captions render marker notes in magenta but never repeat the numbered
  circle glyph (linkage by color only). In panel F the drawn mark ③ sits on the
  contra heel-IN→KINH pathway (alt text: "which excites KINH (marker 3)") while
  the caption's ③ line describes the RUNNER-SIDE phase machine (pm_*). Not
  wrong, but the marker points at one thing and explains another.
- **F5 — rules panel D vs alt text: "Ib additionally projects directly to the
  extensor centers" is not drawn in panel D.** The direct Ib→RG-E edge exists
  only as a dashed element in the deng figure (`draw_circuit.py` lines
  1208–1210); panel D routes everything through IBEXC/LBIN.
- **F6 — dengstyle label collisions (needs full-res human eye).** Two
  independent vision passes reported edge-adjacent label overlaps and possible
  clipping near the canvas edge (MLR band, footer strip) on the 5482×5559
  canvas. The structure-assert contract guarantees every edge GROUP is drawn,
  but says nothing about text collisions — worth one visual pass at full res.
- **F7 — sns_diagram_panels has no legend**; polarity is arrowhead-only
  (graphviz render). Fine for internal docs; add a legend if it ever ships.
- **F8 — inventory caveat** (see Ground truth): default-gain `_net_edges.py`
  omits conditional topology; use a tuned-gain build for IaIN/LBIN/RC/V3/c1/
  HEEL/TOE cross-checks.

## Alt-text consistency

- `circuit_rules_alt.txt` matches the drawn panels A–F edge-for-edge (verified
  against source), except the F4/F5 caveats above.
- `curr3i_channels.alt.txt` correctly documents the frozen-left-leg limitation —
  consistent with walker status at s3i (2026-09-21). Note s3j/s3k ran since;
  if the frozen-left diagnosis changed, that figure+alt text are s3i-era.

## Open questions for Ben

1. Regenerate `circuit_dengstyle` / `circuit_literature` to the 10 pt / no-italic
   standard, or are they internal-only? (Regeneration = M2 gate + your OK.)
2. Add numbered key glyphs to the `circuit_rules` caption lines (fix F4)?
3. Should `_net_edges.py` grow a `--tuned` flag that forces the connectome
   gains so the dump always shows production topology (F8)?
