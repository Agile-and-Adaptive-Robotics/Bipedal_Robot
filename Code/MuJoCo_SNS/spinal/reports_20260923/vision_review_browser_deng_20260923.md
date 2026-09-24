# Circuit-vision review — Ben's browser, Deng partial connectome (2026-09-23)

Mode: live screen audit. Screenshot of the virtual screen (4740x1440), browser window
located on the left monitor, cropped + upscaled for label reading. A Zotero popup
occluded the right-lower canvas (the Ib_ext_R region); that zone is unverified visually.
Ground truth for content: `Neuromechanical_Models\Mujoco_SNS_models\Deng_partial.json`
(26 nodes / 46 edges) — Ben exported it from this very canvas, so canvas == JSON by
construction (serialize() reads the live state).

## What the screen shows

connectome_block_editor.html open with the Deng partial: top band RG (HC-RG-E_72 /
HC-RG-F_71 with IN-InE_73 / IN-InF_74 above), afferent cluster on the right (SN-Ia "from
HIp", SN-Ib "from Ankle", two SN-II hip), PF band (HC-PF-F_55 / HC-PF-E_56 + IN-PF_57/58),
MN band (MN_ext_R / MN_flex_R) with IaIN/IbIN rows, muscles flexor_R / extensor_R, Renshaw
cells at the bottom. Right-leg only (the earlier low-res pass that suggested "bilateral
with commissurals" was a misread — corrected by the zoomed pass). Consistent with Ben's
description: RG layer + one PF layer + one MN layer.

## Wiring flags (from the JSON, for Ben's review — he owns the connectome)

1. **Direct HC-RG-F↔HC-RG-E EXCITATORY edges (g=0.5 both ways)** alongside the correct
   laminated path (HC-F→InF→HC-E inh, HC-E→InE→HC-F inh). Deng has NO direct half-center→
   half-center excitation; mutual excitation between antagonist half-centers pushes toward
   synchrony/latching. Suspect leftover or duplicate edges — delete or confirm intent.
2. **IbIN sign asymmetry**: IbIN_ext_R→MN_ext_R is INHIBITORY (autogenic, tag ib_auto)
   but IN-IbIN_85→MN_flex_R is EXCITATORY (Deng's Ib excitatory reversal). One side
   classic autogenic inhibition, the other Deng-style. Decide the doctrine per side.
3. **Crossed Renshaw wiring**: MN_flex→RC_ext_R (exc) + RC_flex_R→MN_ext_R (inh), and the
   mirror pair — antagonist-coupled RC loops beyond the canonical own-pool (rc_own),
   mutual (rc_mutual), and RC→IaIN disinhibition (rc_disinhib) routes that are also
   present. Direct RC→antagonist-MN inhibition is non-canonical; RC→IaIN is the usual route.
4. **Label hygiene**: "SN-Ia_83 from HIp" (capital I typo); two near-duplicate hip II
   afferents ("SN-II_88 from hip" vs "SN-II_94 from Hip") differing only by case; mixed
   suffix styles (_R vs bare node numbers). Renaming now avoids key collisions when the
   left leg / full model is assembled.
5. **No afferent input to RG**: all four named proximal afferents (hip Ia, ankle Ib,
   hip II x2) target PF interneurons only. Known gap in a partial model (Deng routes
   hip afferents into the RG for phase reset) — fine as "partial", flagged so it isn't
   forgotten.

## Layout / spacing (his standing complaint)

Good y-banding (7 distinct rows). The crowded zone is the MN/IaIN/RC/muscle cluster:
14 nodes and ~24 edges in y 407..708, with three long diagonals (RC_ext→IaIN_ext,
RC_flex→IaIN_flex, and the crossed MN→RC wires) cutting across the muscle afferent fan.
Tight pairs: IN-InF_74/IN-InE_73 (91 px apart, same row — labels nearly touch) and
SN-Ia_83/SN-Ib_84 (45 px vertical). Ib_ext_R sits far right (x=1125) under where the
Zotero popup was. Suggest: pull the RC row into clear corners, keep muscles flanking,
and re-capture with the popup moved for a final visual sign-check of the RG band glyphs
(the zoomed pass left the HC↔HC wire glyph ambiguous vs the JSON's 'exc').

## Follow-ups offered

- Cross-check the full edge list against Deng Table A6 (M3 gate before any
  "literature-faithful" claim).
- Add a visible Save button to the editor (requested earlier, still open).
