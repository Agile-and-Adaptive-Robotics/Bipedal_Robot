"""Slides 15-27: Ch. 4 isometric knee torque (identification, validation, redesign ladder).

All numbers verified against ProofFinal/chapters/30-torque.tex and Appendix C
(94-AppendixC.tex); anchors recorded in evidence_b.md.
"""
import os
from pptx.enum.text import PP_ALIGN
from deck_lib import *

def s15_stands(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 isometric knee torque",
        "Two test stands: a pinned knee exposes the losses; the biomimetic four-bar validates them",
        cite="Dissertation Ch. 4 (Isometric Knee Torque), robot leg model and test stands; four-bar knee: Steele et al. 2018",
        notes="Bone geometry from scans of a 6 ft (1.83 m) person. Both joints share the same center of rotation in the femur frame at 0 deg, and the knee axis is co-axial with gravity so measurements stay isometric. Q: why two stands? A: the pinned stand isolates one DoF for identification; the four-bar validates transfer to the moving-ICR joint the robot uses.")
    pic(s, os.path.join(FIG, "testJigs_knee.png"), box=(0.6, 1.9, 4.2, 4.8))
    tx(s, 0.6, 6.72, 4.2, 0.35, "pinned-joint stand (CAD + test photo)", size=11.5, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    pic(s, os.path.join(FIG, "b_kneeicr_trim.png"), box=(4.9, 1.8, 4.95, 5.02), valign="top")
    tx(s, 4.9, 6.85, 4.95, 0.35, "biomimetic knee: migrating ICR", size=11.5, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    card(s, 9.55, 2.1, 3.2, 4.4, fill=CARD)
    rich(s, 9.78, 2.28, 2.75, 4.1, [
        {"t": "FOUR-BAR KNEE", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "Crossed links give a migrating instantaneous center of rotation - posterior and distal with flexion, like the human knee.", "size": 12.5, "color": TEXT, "spacing": 1.1, "space_after": 8},
        {"t": "Torque depends on the ICR: getting it right matters for both prediction and experiment.", "size": 12.5, "color": TEXT, "spacing": 1.1, "space_after": 8},
        {"t": "Both stands tested with \u230010 mm BPAs; the biomimetic knee also with a \u230020 mm flexor, its routes mirrored about Z.", "size": 12.5, "color": MUTED, "spacing": 1.1},
    ])
    pageno(s, n)

def s16_hybrid(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 hybrid method",
        "The hybrid method shows the error is kinematic, not force",
        cite="Ch. 4, Fig. FlxPin_group A-B (pre-optimized pinned flexor); manuscript bolen_isometric_2025 (in prep.)",
        notes="Yellow = hybrid (measured Lm and rk + model force). It tracks measurement, so the validated force model is not the problem; the kinematics were. Q: is the hybrid a fit? A: no; it substitutes measured Lm and rk into the same M = rF physics, isolating which model term was wrong.")
    pic(s, os.path.join(FIG, "b_flxpin_preopt.png"), box=(0.5, 2.05, 7.55, 3.9))
    rich(s, 0.5, 6.1, 7.55, 0.8, [
        {"t": "yellow = hybrid; it sits on the blue measured dots where they agree - that agreement is the result", "size": 10.5, "color": FAINT, "italic": True, "spacing": 1.05, "space_after": 2},
        {"t": "blue = measured \u00b7 black dashed = original model \u00b7 dashed red = improved BPA characterization", "size": 10.5, "color": FAINT, "italic": True, "spacing": 1.05},
    ])
    card(s, 8.2, 2.1, 4.5, 3.95, fill=CARD)
    rich(s, 8.5, 2.32, 3.95, 3.6, [
        {"t": "HYBRID TORQUE", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "M = r \u00d7 F with measured muscle length and moment arm, model force", "size": 14.5, "color": TEXT, "spacing": 1.1, "space_after": 10},
        {"t": "Hybrid \u2248 measured  \u2713", "size": 16, "color": GREEN, "bold": True, "space_after": 2},
        {"t": "the validated force model is not the problem", "size": 12.5, "color": MUTED, "space_after": 10},
        {"t": "Rigid prediction \u2260 measured  \u2717", "size": 16, "color": RED, "bold": True, "space_after": 2},
        {"t": "the kinematics (L\u2098, r\u2096) are what the model gets wrong", "size": 13, "color": MUTED, "spacing": 1.05},
    ])
    pageno(s, n)

def s17_corrections(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 improved model",
        "Three physical effects close the gap: length offset, series compliance, wrap loss",
        cite="Ch. 4, improved torque model (Eqs. lm, K_br, delta-l); bracket frames: Fig. bktFrame; two-bracket projection: Fig. xiProjection",
        notes="chi1/chi2 are EFFECTIVE system stiffness terms (bracket, fixtures, winch, and tendon lumped together), not literal beam stiffness; every stiffness value quoted is an optimizer output. Q: why not model the bracket by FEA? A: printed anisotropic parts, braided sleeve, and construction-dependent tendon stiffness are unmeasured; the lumped terms are identified from the measured torque itself.")
    cw = 3.95; gap = 0.18
    # chi0
    x = MX
    card(s, x, 1.95, cw, 2.75, fill=CARD)
    tx(s, x + 0.25, 2.10, cw - 0.5, 0.45, "\u03c7\u2080 \u00b7 constant length offset", size=15, color=BROWN, bold=True)
    eq(s, "lm", box=(x + 0.25, 2.58, cw - 0.5, 0.52))
    tx(s, x + 0.25, 3.22, cw - 0.5, 1.4, "Measured muscle length is consistently shorter than predicted at every angle - including at \u03b5 = 1, F = 0. That is a constant offset.", size=12.5, color=TEXT, spacing=1.1)
    # chi1 chi2
    x = MX + cw + gap
    card(s, x, 1.95, cw, 2.75, fill=CARD)
    tx(s, x + 0.25, 2.10, cw - 0.5, 0.45, "\u03c7\u2081, \u03c7\u2082 \u00b7 bracket compliance", size=15, color=BROWN, bold=True)
    eq(s, "kbr", box=(x + 0.25, 2.58, cw - 0.5, 0.52))
    tx(s, x + 0.25, 3.22, cw - 0.5, 1.4, "Hooke law in bracket principal directions; compliance projected on the force direction \u00fb; artificial tendon added in series.", size=12.5, color=TEXT, spacing=1.1)
    # chi3
    x = MX + 2 * (cw + gap)
    card(s, x, 1.95, cw, 2.75, fill=CARD)
    tx(s, x + 0.25, 2.10, cw - 0.5, 0.45, "\u03c7\u2083 \u00b7 wrapping loss (\u03b8\u2096 \u2264 \u03b8_wrap)", size=15, color=BROWN, bold=True)
    eq(s, "wrap", box=(x + 0.25, 2.58, cw - 0.5, 0.52))
    tx(s, x + 0.25, 3.22, cw - 0.5, 1.4, "A BPA bent over the joint loses usable arc length: the wall on the inside of the bend carries less tension. Loss grows with the complement of strain, \u03b6 = 1 \u2212 \u03b5*.", size=12.5, color=TEXT, spacing=1.1)
    tx(s, MX, 4.86, 6.2, 0.26, "bracket frames, pinned configurations (A, B)", size=10.5, color=FAINT, italic=True)
    pic(s, os.path.join(FIG, "b_bktframe_pin.png"), box=(MX, 5.14, 6.2, 1.78), valign="top")
    card(s, 6.8, 4.95, 5.95, 1.85, fill=CARD2)
    rich(s, 7.05, 5.14, 5.5, 1.6, [
        {"t": "Solved by nonlinear force balance", "size": 12, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "F\u2098\u2090\u2093 \u00b7 F*(\u03b5*(\u03b4), P*) \u2212 k\u2091\u2094 \u03b4 = 0  \u2192  \u03b4 (bracketed fzero), then length, path, moment arm, and torque recomputed at the deflected geometry.", "size": 12.5, "color": TEXT, "spacing": 1.08},
    ])
    pageno(s, n)

def s18_identification(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 identification protocol",
        "Terms identified by multiobjective optimization, judged on held-out tests",
        cite="Ch. 4, identification (gamultiobj / NSGA-II, Deb 2002); bounds: Appendix C, Table C.3 (app_opt_bounds)",
        notes="gamultiobj: flexor two-bracket fit evolves 150 for up to 600 generations, extensor refit 50 for 150; parameters searched normalized (cm, log stiffness). The encoder-corrected 47.9 cm test is never trained on. Q: why hold one test out? A: the pinned-flexor data admit a family of stiffness pairs of comparable fit, so the adopted set is selected by held-out prediction, not residual minimization.")
    node(s, MX, 2.0, 3.3, 1.0, "5 pinned-flexor tests\n(48.5 / 45.7 / 47.9 / 40.6-tendon / 41.0 cm)", fill=CARD2, size=12.5)
    node(s, MX + 3.9, 2.0, 2.6, 1.0, "train on 4\nhold out 47.9 cm", fill=CARD, size=12.5)
    node(s, MX + 7.1, 2.0, 2.3, 1.0, "Pareto front\n{RMSE, FVU, max res.}", fill=CARD, size=12.5)
    node(s, MX + 9.85, 2.0, 2.35, 1.0, "pick: held-out fit +\ncross-configuration", fill=CARD, size=12.5)
    arrow(s, MX + 3.35, 2.5, MX + 3.85, 2.5); arrow(s, MX + 6.55, 2.5, MX + 7.05, 2.5); arrow(s, MX + 9.45, 2.5, MX + 9.8, 2.5)
    node(s, MX, 3.55, 5.6, 0.95, "flexor fit solves (\u03c7\u2080, \u03c7\u2081, \u03c7\u2082) with TWO brackets\n(origin + insertion), shared stiffness pair", fill=CARD, size=12.5)
    node(s, MX + 6.6, 3.55, 5.6, 0.95, "extensor refit LOCKS (\u03c7\u2081, \u03c7\u2082) to the flexor-identified values,\nsolves only (\u03c7\u2080, \u03c7\u2083); \u03c7\u2080 allowed negative", fill=CARD, size=12.5)
    card(s, MX, 4.9, 12.2, 1.65, fill=CARD2)
    rich(s, MX + 0.3, 5.05, 11.6, 1.4, [
        {"t": "Bounds of record (Appendix C)", "size": 12, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "flexor:  \u03c7\u2080 \u2208 [0, 20] mm \u00b7 \u03c7\u2081 \u2208 [3\u00d710\u2074, 10\u2076] N/m \u00b7 \u03c7\u2082 \u2208 [5\u00d710\u00b3, 2\u00d710\u2074] N/m", "size": 14, "color": TEXT, "space_after": 2},
        {"t": "extensor:  \u03c7\u2080 \u2208 [\u221220, +5] mm \u00b7 \u03c7\u2083 \u2208 [0, 1] \u00b7 (\u03c7\u2081, \u03c7\u2082) = flexor pair, never searched", "size": 14, "color": TEXT},
    ])
    pageno(s, n)

def s19_flexor_opt(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 pinned flexor",
        "With corrections, prediction tracks measurement - including the dip at small angles",
        cite="Ch. 4, Fig. FlxPin_group C-D; Table TorqueFit rows 1-2 (initial single-bracket identification)",
        notes="Panels C (left) and D (right) of Fig. FlxPin_group. Panel C: fit on the 48.5 cm BPA (the panel title reads 48.0 cm; the chapter caption and text use 48.5 cm). Panel D: 45.7 cm, never used in the fit. GoF: 48.5 cm RMSE 1.616, FVU 0.0958, max 3.2527; 45.7 cm validation RMSE 2.1214, FVU 0.0878, max 5.1690. Q: why does the corrected model capture the dip? A: series compliance relieves contraction where force is high and contraction is low, exactly the small-angle regime.")
    pic(s, os.path.join(FIG, "b_flxpin_opt.png"), box=(4.55, 2.15, 7.7, 3.6))
    tx(s, 4.55, 5.7, 7.7, 0.3, "blue = measured \u00b7 dashed orange = BPA characterization alone \u00b7 magenta = optimized prediction", size=10.5, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    rich(s, MX, 2.1, 3.9, 4.4, [
        {"t": "THE MISSING DIP", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "The rigid model completely missed the torque reduction at small flexion angles.", "size": 14.5, "color": TEXT, "spacing": 1.12, "space_after": 12},
        {"t": "The corrected model captures it.", "size": 15.5, "color": GREEN, "bold": True, "space_after": 14},
        {"t": "VALIDATION", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "45.7 cm BPA: different length, different working range, never in the fit.", "size": 14.5, "color": TEXT, "spacing": 1.12},
    ])
    pageno(s, n)

def s20_flexor_numbers(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 pinned flexor",
        "Corrections cut flexor training-mean RMSE from 8.12 to 1.56 N\u00b7m; held-out validates at 2.29",
        cite="Ch. 4, Table FlxPinGof (corrected model, adopted terms); Appendix C, Table app_gof",
        notes="Zero-correction per-test RMSE spans 6.37-9.69 N-m (improved BPA characterization alone); corrected FVU stays at or below 0.17 on all five tests, the 47.9 cm test being pure held-out validation. Q: why is held-out slightly worse than the training mean? A: 2.29 vs 1.56 N-m is expected; it is still 3.3x better than its own zero-correction RMSE of 7.67.")
    bar_chart(s, 0.7, 2.0, 8.1, 4.5,
              ["48.5 cm", "45.7 cm", "47.9 cm*", "40.6 cm-tendon", "41.0 cm"], [2.165, 1.368, 2.287, 1.611, 1.084],
              [GREEN, GREEN, TEAL, GREEN, GREEN],
              title="per-test torque RMSE at adopted \u03c7 values (N\u00b7m)", num_fmt="0.000", max_scale=2.6)
    card(s, 9.2, 2.2, 3.55, 4.2, fill=CARD)
    rich(s, 9.45, 2.42, 3.1, 3.9, [
        {"t": "FVU (all five tests)", "size": 12, "color": TEAL, "bold": True, "space_after": 2},
        {"t": "0.17 / 0.05 / 0.16 / 0.04 / 0.03", "size": 14.5, "color": TEXT, "spacing": 1.1, "space_after": 12},
        {"t": "* 47.9 cm test", "size": 12, "color": TEAL, "bold": True, "space_after": 2},
        {"t": "encoder-corrected, never trained on: pure held-out validation.", "size": 13.5, "color": TEXT, "spacing": 1.12, "space_after": 12},
        {"t": "Every test predicts.", "size": 14, "color": BROWN, "bold": True, "spacing": 1.1},
    ])
    pageno(s, n)

def s21_extensor(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 pinned extensor",
        "The extensor wraps the joint - and \u03c7\u2083 recovers the length it loses there",
        cite="Ch. 4, Table ExtPinGof + follow-up identification (adopted \u03c7 values); wrap loss: Eq. \u0394l = \u03c7\u2083 R |\u0394\u03b8| \u03b6\u00b2",
        notes="Validation test 45.7 cm: RMSE 2.178 to 0.490 N-m, FVU 1.059 to 0.054. The stiffness pair is locked to the flexor-identified values and never searched. Initial single-bracket extensor values were chi0 = -10.9 mm, chi3 = 0.123 (Table OptCoeff). Q: why does chi0 flip sign? A: wrapping consumes usable length, so the extensor behaves as if longer; its zero point shifts opposite to the flexor's.")
    pic(s, os.path.join(FIG, "ExtPin_group.png"), box=(0.6, 1.9, 7.2, 4.85))
    stats = [
        ("\u03c7\u2080 = \u22126.38 mm", "extensor length offset: OPPOSITE sign to the flexor's +3.94 mm", RED),
        ("\u03c7\u2083 = 0.158", "fractional wrap loss acting on arc length \u00d7 \u03b6\u00b2", BROWN),
        ("0.49\u20131.53 N\u00b7m", "corrected RMSE, improved in every pool test (6/6)", GREEN),
        ("2.60 \u2192 1.16 N\u00b7m", "fitting-set mean RMSE; FVU < 1 on all but the 43.2 cm tendon test", TEAL),
    ]
    y = 2.0
    for v, lab, c in stats:
        card(s, 8.1, y, 4.6, 1.08, fill=CARD if c is not RED else CARD2)
        rich(s, 8.35, y + 0.12, 4.1, 0.95, [
            {"t": v, "size": 16.5, "color": c, "bold": True, "space_after": 1},
            {"t": lab, "size": 11.5, "color": MUTED, "spacing": 1.0},
        ])
        y += 1.2
    pageno(s, n)

def s22_xi_record(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 values of record",
        "One stiffness pair of record spans the flexor and extensor configurations",
        cite="Appendix C, Table app_xi_values; flexor front row 77 of 212, extensor flxr77 row 32",
        notes="Mats of record: minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat row 77, and minimizeExt10mmX3_results_20260920_flxr77.mat row 32 (pair locked to row 77). Initial single-bracket values (Table OptCoeff): flexor 11.9 mm / 1.02e5 / 1.06e4 N/m; extensor -10.9 mm / 0.123. Q: how sensitive is the redesign to the pair choice? A: the adoption campaign screened the front for a pair meeting both targets simultaneously: +5.00% flexor and +5.043% extensor at adoption, refined to +5.01% and +5.02% (Appendix C).")
    hdr = ["Term", "Value", "Source", "Row"]
    rows = [
        ("Flexor \u03c7\u2080", "+3.94 mm", "App. C", "77"),
        ("Flexor \u03c7\u2081", "3.998\u00d710\u2074 N/m", "App. C", "77"),
        ("Flexor \u03c7\u2082", "1.473\u00d710\u2074 N/m", "App. C", "77"),
        ("Extensor \u03c7\u2080", "\u22126.38 mm", "App. C", "32"),
        ("Extensor \u03c7\u2081 (locked)", "3.998\u00d710\u2074 N/m", "App. C", "77"),
        ("Extensor \u03c7\u2082 (locked)", "1.473\u00d710\u2074 N/m", "App. C", "77"),
        ("\u03c7\u2083 (extensor)", "0.158", "App. C", "32"),
    ]
    cw = [2.6, 2.4, 1.9, 0.8]
    xx = 0.9
    for i, h_ in enumerate(hdr):
        tx(s, xx, 1.95, cw[i], 0.4, h_, size=13, color=TEAL, bold=True, align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
        xx += cw[i]
    y = 2.4
    for r in rows:
        xx = 0.9
        for i, val in enumerate(r):
            bold = i <= 1
            colr = TEXT
            if "Flexor" in r[0] and i <= 1: colr = GREEN
            if "Extensor" in r[0] and i <= 1: colr = BLUE
            if r[0] == "\u03c7\u2083 (extensor)" and i <= 1: colr = BROWN
            tx(s, xx, y, cw[i], 0.4, val, size=14 if i != 1 else 13.5, color=colr, bold=bold,
               align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
            xx += cw[i]
        y += 0.5
    card(s, 8.85, 1.95, 3.9, 4.35, fill=CARD2)
    rich(s, 9.1, 2.15, 3.45, 4.05, [
        {"t": "READ THIS TABLE LIKE THIS", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "The pinned-flexor data admit a family of stiffness pairs of comparable fit quality.", "size": 13.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "The adopted pair was chosen by a requirement external to the fit: both redesign configurations must clear the required 5% torque-overshoot target.", "size": 13.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "The extensor never searches stiffness; it inherits it.", "size": 13.5, "color": BROWN, "bold": True, "spacing": 1.1},
    ])
    pageno(s, n)

def s23_bioflexor(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 biomimetic knee",
        "The biomimetic knee meets or exceeds human flexor torque over most of the range of motion",
        cite="Ch. 4, Fig. FlxGrp + Table TorqueFit; human curve: Gait2392 maximum-isometric torque (Delp 2007; Seth 2018)",
        notes="Pinned-flexor chi values transferred unchanged. Panel A: \u230010 mm, l0 = 41.5 cm, 21 mm tendon, 604 kPa, RMSE 0.3007, FVU 0.0102, max 0.4647. Panel B: \u230020 mm, l0 = 42.0 cm, 15 mm tendon, 620 kPa, RMSE 1.8024, FVU 0.0499, max 3.5678. Note: the chapter text gives 42.0 cm for the \u230020 test while its summary table row lists 41.5 cm; the slide follows the text. Q: is the 60-100 deg shortfall fixable? A: yes by route redesign: the two-BPA flexor route exceeds the human envelope at every angle (next slides).")
    pic(s, os.path.join(FIG, "b_flxgrp_trim.png"), box=(0.5, 1.95, 10.05, 4.7))
    tx(s, 0.5, 6.7, 10.05, 0.32, "blue = measured \u00b7 magenta = corrected prediction \u00b7 dashed orange = BPA characterization alone \u00b7 black dotted = human maximum-isometric torque", size=10, color=FAINT, italic=True)
    card(s, 10.7, 2.0, 2.35, 4.7, fill=CARD)
    rich(s, 10.86, 2.18, 2.05, 4.4, [
        {"t": "SHORTFALL", "size": 11, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "60\u2013100\u00b0\n(\u00b120%: 78\u2013103\u00b0)", "size": 13, "color": RED, "bold": True, "spacing": 1.05, "space_after": 10},
        {"t": "\u230020 mm, l\u2080 = 42.0 cm\nRMSE 1.80 \u00b7 FVU 0.05\nmax resid 3.57 N\u00b7m", "size": 12, "color": TEXT, "spacing": 1.05, "space_after": 10},
        {"t": "\u230010 mm panel:\nRMSE 0.30 \u00b7 FVU 0.01", "size": 12, "color": TEXT, "spacing": 1.05, "space_after": 10},
        {"t": "At 620 kPa the BPA exceeded the test femur's ultimate strength; remaining points taken at 560/421/325/281 kPa.", "size": 10.5, "color": MUTED, "italic": True, "spacing": 1.05},
    ])
    pageno(s, n)

def s24_bioext(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 biomimetic knee",
        "Biomimetic extensor: corrections cut the error, but prediction still sits above measurement",
        cite="Ch. 4, Fig. Ext10mm_52cm + Table TorqueFit; adopted values: Appendix C, adopted correction-term values",
        notes="Measured torque is significantly smaller than both the original and the hybrid predictions; the corrected prediction more closely matches but remains higher than measurement, so level validation is not claimed. With the pinned-extensor chi values (the figure): RMSE 0.8920, FVU 0.3396, max 1.4188. At the values of record: max residual 2.21 N-m against a rigid baseline of 5.59. Q: why does extensor validation lag the flexor? A: the route wraps the joint through sequential contacts and passes between the femoral condyles near the ICR; the path model itself is the approximation.")
    pic(s, os.path.join(FIG, "b_ext10mm_trim.png"), box=(0.4, 1.95, 8.35, 3.75))
    tx(s, 0.4, 5.74, 8.35, 0.3, "\u230010 mm extensor on the four-bar knee, l\u2080 = 51.8 cm", size=12, color=MUTED, align=PP_ALIGN.CENTER)
    stats = [("0.89", "RMSE (N\u00b7m) with pinned-extensor \u03c7 (figure)", GREEN),
             ("1.17", "RMSE (N\u00b7m) at values of record (rigid: 3.03)", TEAL),
             ("0.67", "FVU at values of record (rigid: 4.49)", TEAL)]
    x = 0.4
    for v, lab, c in stats:
        card(s, x, 6.02, 2.62, 1.0, fill=CARD)
        rich(s, x + 0.1, 6.14, 2.42, 0.85, [
            {"t": v, "size": 16, "color": c, "bold": True, "align": PP_ALIGN.CENTER, "space_after": 1},
            {"t": lab, "size": 9.5, "color": MUTED, "align": PP_ALIGN.CENTER, "spacing": 0.95},
        ])
        x += 2.72
    pic(s, os.path.join(FIG, "b_passthru.png"), box=(8.55, 1.95, 4.25, 2.75))
    tx(s, 8.55, 4.78, 4.25, 0.55, "the route passes between the femoral condyles near the ICR: an approximation with inherent error", size=10.5, color=MUTED, italic=True, spacing=1.05)
    card(s, 8.55, 5.42, 4.25, 1.6, fill=CARD2)
    rich(s, 8.8, 5.56, 3.8, 1.4, [
        {"t": "CLAIM CONTROL", "size": 11, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "Improvement is unambiguous; the corrected prediction still sits above the measured torque, so level validation is not claimed.", "size": 12.5, "color": TEXT, "spacing": 1.1},
    ])
    pageno(s, n)

def s25_ladder(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 route redesign",
        "The \u03c7 values of record drive a redesign ladder: direct routes, a cam, then pulleys",
        cite="Ch. 4, follow-up identification and route redesign + redesigned-route results; solver: surrogateopt \u2192 patternsearch",
        notes="Direct flexor: predicted knee-axis torque -81.23 to -3.86 N-m vs human-scaled requirement -23.16 to -2.46 N-m (peak 23.16), exceeded at every angle; margin fraction 0.0501 vs required 0.05. Direct extensor: 35.45-109.4 N-m vs 20.25-50.32 N-m, zero shortfall points. Solver budgets: surrogateopt 3,000 evals (direct extensor), 1,000 (flexor), 1,000 + 1,500 refinement (cam), then patternsearch 15,000. Q: why is the flexor margin also about 5%? A: the adoption rule for the stiffness pair required both configurations to clear the +5% envelope, and the flexor model still assumes both actuators share one attachment route; making the mirrored pair explicit is the next refinement.")
    nw = 2.83; gap = 0.30
    y0 = 2.15; nh = 2.45
    node(s, MX, y0, nw, nh,
         "1 \u00b7 DIRECT ROUTES\ntwo \u230020 mm BPAs each\n\nflexor  +5.01% @ \u2212107\u00b0\nextensor  +5.02% @ \u221220.2\u00b0\n(later superseded by cam)",
         fill=CARD, size=12)
    node(s, MX + nw + gap, y0, nw, nh,
         "2 \u00b7 CAM EXTENSOR\ndesign of record\n\n+30.83% vs required envelope\n+39.77% vs plain human\n269\u2013507 kPa to meet it",
         fill=CARD2, line=TEAL, size=12)
    node(s, MX + 2 * (nw + gap), y0, nw, nh,
         "3 \u00b7 MONO-PULLEY\nknee flexor, 1:1 tackle\n\npeak 81.23 \u2192 85.74 N\u00b7m\n+4.51 N\u00b7m (+5.6%)\nno mechanical advantage",
         fill=CARD, size=12)
    node(s, MX + 3 * (nw + gap), y0, nw, nh,
         "4 \u00b7 BI-PULLEY SCREEN\nfive anatomical sites\n\nall five rejected\nclearances 4.26\u201319.14 mm\nvs 25 mm required",
         fill=CARD, size=12)
    ym = y0 + nh / 2
    arrow(s, MX + nw + 0.03, ym, MX + nw + gap - 0.03, ym)
    arrow(s, MX + 2 * nw + gap + 0.03, ym, MX + 2 * (nw + gap) - 0.03, ym)
    arrow(s, MX + 3 * nw + 2 * gap + 0.03, ym, MX + 3 * (nw + gap) - 0.03, ym)
    card(s, MX, 4.85, 12.23, 1.05, fill=CARD2)
    rich(s, MX + 0.3, 4.98, 11.6, 0.85, [
        {"t": "shared plan: two \u230020 mm BPAs per route \u00b7 620 kPa ceiling \u00b7 maximum usable contraction 0.255", "size": 12.5, "color": TEXT, "space_after": 3},
        {"t": "targets: vastus medialis (extensor), biceps femoris short head (flexor) \u00b7 required envelope = human curve + 5%", "size": 12.5, "color": TEXT},
    ])
    card(s, MX, 6.05, 12.23, 0.85, fill=CARD)
    tx(s, MX + 0.3, 6.17, 11.6, 0.65,
       "flexor route: \u221281.23 to \u22123.86 N\u00b7m vs human-scaled \u221223.16 to \u22122.46 N\u00b7m \u00b7 direct extensor: 35.45\u2013109.4 vs 20.25\u201350.32 N\u00b7m \u00b7 envelope exceeded at every angle",
       size=12, color=MUTED, spacing=1.05)
    pageno(s, n)

def s26_cam(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 route redesign",
        "The cam extensor is the design of record: +30.83% margin with hardware-safe geometry",
        cite="Ch. 4, redesigned-route results (cam paragraph); figure: Mesh_Optimization\\Results\\Vas_Pam_20mm_cam_winner_4.png (print-QC pass)",
        notes="The final geometry folds the optimizer winner together with a verified adjustment of the distal route point and tendon length; the winner alone carried +31.30% and +42.34% at the same checks. Envelope pressures sit comfortably above the 150 kPa regulated floor of the VTUG valve manifold that supplies the hardware. Hardware: the first and third points sit at essentially the same knee angle and span a twentyfold range of held torque; the full sweep and the flexor measurements are pending. Q: is the cam margin sensitive to the adopted chi pair? A: margins are evaluated under the values of record, whose adoption rule already required both routes to clear +5%.")
    pic(s, os.path.join(FIG, "b_cam_moment_arm.png"), box=(0.55, 1.9, 6.15, 3.9))
    tx(s, 0.55, 5.85, 6.15, 0.35, "extensor cam, design of record: moment arm vs knee angle", size=11.5, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    stat(s, 7.0, 1.95, 2.8, "+30.83%", "margin vs required envelope (human + 5%)", color=GREEN, vsize=26, lsize=11)
    stat(s, 10.0, 1.95, 2.8, "+39.77%", "margin vs plain human target", color=GREEN, vsize=26, lsize=11)
    stat(s, 7.0, 3.15, 2.8, "269\u2013507 kPa", "pressure to meet envelope (VTUG floor 150 kPa)", color=TEAL, vsize=24, lsize=11)
    stat(s, 10.0, 3.15, 2.8, "96.9 mm", "tibial-wrap bend radius", color=BROWN, vsize=26, lsize=11)
    stat(s, 7.0, 4.35, 2.8, "43.1 mm", "condyle clearance at 10\u00b0 extension", color=BROWN, vsize=26, lsize=11)
    stat(s, 10.0, 4.35, 2.8, "72.4 mm", "mean moment arm (span 59.4\u2013102.8 mm)", color=BROWN, vsize=26, lsize=11)
    card(s, MX, 6.3, 12.23, 0.8, fill=CARD2)
    rich(s, MX + 0.3, 6.4, 11.6, 0.65, [
        {"t": "HARDWARE BEGUN \u00b7 rebuilt \u230020 mm extensor route", "size": 11, "color": TEAL, "bold": True, "space_after": 2},
        {"t": "first isometric points: 1.163 / 7.849 / 23.437 N\u00b7m at 10\u00b0 / 11\u00b0 / 10\u00b0, hung loads 3.7 / 25.2 / 76 N at 388.5 / 388 / 592 kPa; full sweep pending", "size": 11.5, "color": TEXT},
    ])
    pageno(s, n)

def s27_pulleys(prs, n):
    s = content_slide(prs, "Ch. 4 \u00b7 route redesign",
        "Pulleys buy reach at the knee, and screen the anatomy: geometry binds, not torque",
        cite="Ch. 4, mono- and bi-pulley tendon routing + Table bipulley_batch",
        notes="Pulley run: collision-feasible, zero bound violations, pattern-search exit flag 1; objective values are not commensurable across families (0.08465 pulley vs 74.47 direct). Insertion-side measure of the pulley peak: 85.56 N-m. Bi-pulley: the medial-gastrocnemius run records no finite clearance (worst constraint 10, exit flag -2); rectus femoris worst constraint 0.2423. These are optimizer-stage results, not yet propagated into the OpenSim or MuJoCo model trees. Q: why 25 mm clearance? A: two 10 mm actuator radii plus a 5 mm bundle margin. Q: does the bi-pulley failure kill the concept? A: no; it redirects design toward clearance-relaxed attachment sites.")
    card(s, MX, 1.95, 5.9, 4.35, fill=CARD)
    rich(s, MX + 0.3, 2.15, 5.3, 4.05, [
        {"t": "MONO-PULLEY \u00b7 KNEE FLEXOR", "size": 12, "color": TEAL, "bold": True, "space_after": 10},
        {"t": "81.23 \u2192 85.74 N\u00b7m", "size": 24, "color": GREEN, "bold": True, "space_after": 1},
        {"t": "peak predicted flexion torque (knee axis), single tackle-line part at 1:1 gain", "size": 11.5, "color": MUTED, "spacing": 1.0, "space_after": 12},
        {"t": "+4.51 N\u00b7m (+5.6%)", "size": 19, "color": BROWN, "bold": True, "space_after": 1},
        {"t": "from the redirect alone, no mechanical advantage", "size": 11.5, "color": MUTED, "spacing": 1.0, "space_after": 12},
        {"t": "\u2264 0.803 N\u00b7m off-axis", "size": 19, "color": BROWN, "bold": True, "space_after": 1},
        {"t": "under 1% of peak \u00b7 human flexor scale 23.16 N\u00b7m exceeded \u00d73.51 (direct) / \u00d73.70 (pulley)", "size": 11.5, "color": MUTED, "spacing": 1.0},
    ])
    tx(s, 6.85, 2.0, 6.0, 0.35, "BI-PULLEY SCREEN \u00b7 FIVE OPENSIM MUSCLE SITES", size=12, color=TEAL, bold=True)
    hdr = ["Muscle site", "Objective", "Outcome"]
    rows = [
        ("soleus", "0.997", "nonpositive margin"),
        ("tibialis anterior", "0.929", "nonpositive margin"),
        ("rectus femoris", "3.42", "worst constraint 0.242"),
        ("medial gastrocnemius", "101*", "sentinel \u00b7 exit flag \u22122"),
        ("biceps femoris sh. head", "101*", "sentinel \u00b7 exit flag \u22122"),
    ]
    cw = [2.55, 1.15, 2.3]
    xx = 6.85
    for i, h_ in enumerate(hdr):
        tx(s, xx, 2.45, cw[i], 0.35, h_, size=11.5, color=TEAL, bold=True, align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
        xx += cw[i]
    y = 2.85
    for r in rows:
        xx = 6.85
        for i, val in enumerate(r):
            tx(s, xx, y, cw[i], 0.35, val, size=12.5, color=TEXT, bold=(i == 1),
               align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
            xx += cw[i]
        y += 0.42
    tx(s, 6.85, y + 0.02, 6.0, 0.3, "* infeasibility sentinel objective (worst constraint 10)", size=10, color=FAINT, italic=True)
    card(s, 6.85, 5.15, 5.95, 1.15, fill=CARD2)
    rich(s, 7.1, 5.3, 5.45, 0.95, [
        {"t": "4.26\u201319.14 mm vs 25 mm required", "size": 15, "color": RED, "bold": True, "space_after": 2},
        {"t": "minimum route-to-route clearances (two \u230010 mm radii + 5 mm bundle margin): geometry binds, not torque", "size": 11.5, "color": MUTED, "spacing": 1.0},
    ])
    tx(s, MX, 6.5, 12.2, 0.4, "optimizer-stage results; not yet propagated into the OpenSim or MuJoCo model trees", size=11, color=FAINT, italic=True)
    pageno(s, n)

def build(prs, start=15):
    s15_stands(prs, start)
    s16_hybrid(prs, start + 1)
    s17_corrections(prs, start + 2)
    s18_identification(prs, start + 3)
    s19_flexor_opt(prs, start + 4)
    s20_flexor_numbers(prs, start + 5)
    s21_extensor(prs, start + 6)
    s22_xi_record(prs, start + 7)
    s23_bioflexor(prs, start + 8)
    s24_bioext(prs, start + 9)
    s25_ladder(prs, start + 10)
    s26_cam(prs, start + 11)
    s27_pulleys(prs, start + 12)
    return start + 12
