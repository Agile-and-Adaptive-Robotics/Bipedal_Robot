"""Slides 15-25: Aim 2 (knee torque, Xi corrections, biomimetic validation)."""
import os
from pptx.enum.text import PP_ALIGN
from deck_lib import *

def s15_stands(prs, n):
    s = content_slide(prs, "Aim 2 \u00b7 joint torque",
        "Two test stands: a pinned knee exposes the losses; the biomimetic four-bar validates them",
        cite="Bolen et al., in preparation 2025 (Ch. 3-4); four-bar knee: Steele 2018",
        notes="Same CoR in femur frame at 0 deg. Co-axial with gravity so measurements stay isometric.")
    pic(s, os.path.join(FIG, "testJigs_knee.png"), box=(0.6, 1.9, 4.5, 4.8))
    tx(s, 0.6, 6.72, 4.5, 0.35, "pinned-joint stand (CAD + test photo)", size=11.5, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    pic(s, os.path.join(FIG, "KneeICR.png"), box=(5.5, 1.9, 3.6, 4.8))
    tx(s, 5.5, 6.72, 3.6, 0.35, "biomimetic knee: migrating ICR", size=11.5, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    card(s, 9.4, 2.1, 3.35, 4.4, fill=CARD)
    rich(s, 9.65, 2.3, 2.9, 4.1, [
        {"t": "FOUR-BAR KNEE", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "Crossed links give a migrating instantaneous center of rotation - posterior and distal with flexion, like the human knee.", "size": 13.5, "color": TEXT, "spacing": 1.1, "space_after": 8},
        {"t": "Torque depends on the ICR: getting it right matters for both prediction and experiment.", "size": 13.5, "color": TEXT, "spacing": 1.1, "space_after": 8},
        {"t": "Each joint: one flexor + one extensor BPA, 10 mm (and 20 mm on the biomimetic flexor).", "size": 13.5, "color": MUTED, "spacing": 1.1},
    ])
    pageno(s, n)

def s16_hybrid(prs, n):
    s = content_slide(prs, "Aim 2 \u00b7 joint torque",
        "The hybrid method shows the error is kinematic, not force",
        cite="Bolen et al., in preparation 2025, Fig. 6 A-B (panels: pre-optimized pinned flexor)",
        notes="Yellow = hybrid (measured Lm, rk + model force). It tracks measurement: so the force model is fine; Lm and rk predictions were wrong.")
    pic(s, os.path.join(FIG, "flxpin_AB.png"), box=(0.6, 1.9, 7.3, 4.85))
    card(s, 8.2, 2.1, 4.5, 3.95, fill=CARD)
    rich(s, 8.5, 2.32, 3.95, 3.6, [
        {"t": "HYBRID TORQUE", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "M = r \u00d7 F with measured muscle length and angle, model force", "size": 14.5, "color": TEXT, "spacing": 1.1, "space_after": 10},
        {"t": "Hybrid \u2248 measured  \u2713", "size": 16, "color": GREEN, "bold": True, "space_after": 2},
        {"t": "the validated force model is not the problem", "size": 12.5, "color": MUTED, "space_after": 10},
        {"t": "Rigid prediction \u2260 measured  \u2717", "size": 16, "color": RED, "bold": True, "space_after": 2},
        {"t": "the kinematics (L\u2098, r\u2096) are what the model gets wrong", "size": 13, "color": MUTED, "spacing": 1.05},
    ])
    tx(s, 8.2, 6.25, 4.5, 0.6, "yellow = hybrid \u00b7 blue = measured \u00b7 black dashed = original model", size=11.5, color=FAINT, italic=True, spacing=1.05)
    pageno(s, n)

def s17_corrections(prs, n):
    s = content_slide(prs, "Aim 2 \u00b7 improved torque model",
        "Three physical effects close the gap: length offset, series compliance, wrap loss",
        cite="Bolen et al., in preparation 2025, Eqs. 10-16; bracket frames: Fig. 3.4",
        notes="chi1/chi2 are EFFECTIVE system stiffness (bracket + fixtures + winch), not literal beam stiffness - be ready for that question.")
    cw = 3.95; gap = 0.18
    # chi0
    x = MX
    card(s, x, 1.95, cw, 3.0, fill=CARD)
    tx(s, x + 0.25, 2.12, cw - 0.5, 0.45, "\u03c7\u2080 \u00b7 constant length offset", size=15, color=BROWN, bold=True)
    eq(s, "lm", box=(x + 0.25, 2.62, cw - 0.5, 0.52))
    tx(s, x + 0.25, 3.3, cw - 0.5, 1.55, "Measured muscle length is consistently shorter than predicted at every angle - including at \u03b5 = 1, F = 0. That is a constant offset.", size=13, color=TEXT, spacing=1.1)
    # chi1 chi2
    x = MX + cw + gap
    card(s, x, 1.95, cw, 3.0, fill=CARD)
    tx(s, x + 0.25, 2.12, cw - 0.5, 0.45, "\u03c7\u2081, \u03c7\u2082 \u00b7 bracket compliance", size=15, color=BROWN, bold=True)
    eq(s, "kbr", box=(x + 0.25, 2.62, cw - 0.5, 0.52))
    tx(s, x + 0.25, 3.3, cw - 0.5, 1.55, "Hooke law in bracket principal directions; compliance projected on the force direction \u00fb; artificial tendon added in series.", size=13, color=TEXT, spacing=1.1)
    # chi3
    x = MX + 2 * (cw + gap)
    card(s, x, 1.95, cw, 3.0, fill=CARD)
    tx(s, x + 0.25, 2.12, cw - 0.5, 0.45, "\u03c7\u2083 \u00b7 wrapping loss (\u03b8\u2096 \u2264 \u03b8_wrap)", size=15, color=BROWN, bold=True)
    eq(s, "wrap", box=(x + 0.25, 2.62, cw - 0.5, 0.52))
    tx(s, x + 0.25, 3.3, cw - 0.5, 1.55, "A BPA bent over the joint loses usable arc length: inner side carries less tension. Loss grows with the complement of strain, \u03b6 = 1 \u2212 \u03b5*.", size=13, color=TEXT, spacing=1.1)
    pic(s, os.path.join(FIG, "bktFrame.png"), box=(MX, 5.2, 5.9, 1.55))
    tx(s, MX, 6.78, 5.9, 0.35, "bracket frames per configuration", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    card(s, 6.8, 5.2, 5.95, 1.55, fill=CARD2)
    rich(s, 7.05, 5.35, 5.5, 1.3, [
        {"t": "Solved by nonlinear force balance", "size": 12, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "F\u2086\u2082\u2080(l\u2080) \u00b7 F*(\u03b5*(\u03b4), P*) \u2212 k\u2091\u2094 \u03b4 = 0  \u2192  \u03b4 (Newton / bracketed solve), then length, path, moment arm, torque recomputed.", "size": 12.5, "color": TEXT, "spacing": 1.08},
    ])
    pageno(s, n)

def s18_identification(prs, n):
    s = content_slide(prs, "Aim 2 \u00b7 identification protocol",
        "Terms identified by multiobjective optimization, judged on held-out tests",
        cite="gamultiobj / NSGA-II (Deb 2002); sources: minimizeFlxPin10mm_2brk.m, minimizeExt10mmX3.m",
        notes="Key protocol choices: encoder-corrected 47 cm is NEVER trained on by the flexor fit; extensor locks stiffness to the flexor front.")
    node(s, MX, 2.0, 3.3, 1.0, "5 pinned-flexor tests\n(48 / 46 / 47 / 40-tendon / 41 cm)", fill=CARD2, size=12.5)
    node(s, MX + 3.9, 2.0, 2.6, 1.0, "train on 4\nhold out 47 cm", fill=CARD, size=12.5)
    node(s, MX + 7.1, 2.0, 2.3, 1.0, "Pareto front\n{RMSE, FVU, max res.}", fill=CARD, size=12.5)
    node(s, MX + 9.85, 2.0, 2.35, 1.0, "pick: held-out fit +\ncross-configuration", fill=CARD, size=12.5)
    arrow(s, MX + 3.35, 2.5, MX + 3.85, 2.5); arrow(s, MX + 6.55, 2.5, MX + 7.05, 2.5); arrow(s, MX + 9.45, 2.5, MX + 9.8, 2.5)
    node(s, MX, 3.55, 5.6, 0.95, "flexor fit solves (\u03c7\u2080, \u03c7\u2081, \u03c7\u2082) with TWO brackets\n(origin + insertion), shared stiffness pair", fill=CARD, size=12.5)
    node(s, MX + 6.6, 3.55, 5.6, 0.95, "extensor refit LOCKS (\u03c7\u2081, \u03c7\u2082) to the flexor front,\nsolves only (\u03c7\u2080, \u03c7\u2083); \u03c7\u2080 allowed negative", fill=CARD, size=12.5)
    card(s, MX, 4.9, 12.2, 1.65, fill=CARD2)
    rich(s, MX + 0.3, 5.05, 11.6, 1.4, [
        {"t": "Bounds of record (Appendix C)", "size": 12, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "flexor:  \u03c7\u2080 \u2208 [0, 20] mm \u00b7 \u03c7\u2081 \u2208 [3\u00d710\u2074, 10\u2076] N/m \u00b7 \u03c7\u2082 \u2208 [5\u00d710\u00b3, 2\u00d710\u2074] N/m", "size": 14, "color": TEXT, "space_after": 2},
        {"t": "extensor:  \u03c7\u2080 \u2208 [\u221220, 0] mm \u00b7 \u03c7\u2083 \u2208 [0, 1] \u00b7 (\u03c7\u2081, \u03c7\u2082) = flexor front, never searched", "size": 14, "color": TEXT},
    ])
    pageno(s, n)

def s19_flexor_opt(prs, n):
    s = content_slide(prs, "Aim 2 \u00b7 pinned flexor",
        "With corrections, prediction tracks measurement - including the dip at small angles",
        cite="Bolen et al., in preparation 2025, Fig. 6 C-D",
        notes="Panel C: fit on 48 cm. Panel D: 45.7 cm never used in the fit - different muscle over a different ROM.")
    pic(s, os.path.join(FIG, "flxpin_CD.png"), box=(4.75, 1.9, 7.4, 4.85))
    rich(s, MX, 2.1, 3.9, 4.4, [
        {"t": "THE MISSING DIP", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "The rigid model completely missed the torque reduction at small flexion angles.", "size": 14.5, "color": TEXT, "spacing": 1.12, "space_after": 12},
        {"t": "The corrected model captures it.", "size": 15.5, "color": GREEN, "bold": True, "space_after": 14},
        {"t": "VALIDATION", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "45.7 cm BPA: different length, different working range, never in the fit.", "size": 14.5, "color": TEXT, "spacing": 1.12},
    ])
    pageno(s, n)

def s20_flexor_numbers(prs, n):
    s = content_slide(prs, "Aim 2 \u00b7 pinned flexor",
        "Held-out flexor errors stay near 1-2 N\u00b7m across all five tests",
        cite="minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat, front row 107 of 212 (adopted); FVU at right",
        notes="47 cm is pure hold-out with the encoder correction (+5.3 deg on experimental angles). FVU < 0.17 everywhere.")
    bar_chart(s, 0.7, 2.0, 8.1, 4.5,
              ["48 cm", "46 cm", "47 cm*", "40 cm-tendon", "41 cm"], [1.94, 1.51, 2.29, 1.63, 1.21],
              [GREEN, GREEN, TEAL, GREEN, GREEN],
              title="per-test torque RMSE at adopted \u03c7 values (N\u00b7m)", num_fmt="0.00", max_scale=2.6)
    card(s, 9.2, 2.2, 3.55, 4.2, fill=CARD)
    rich(s, 9.45, 2.42, 3.1, 3.9, [
        {"t": "FVU", "size": 12, "color": TEAL, "bold": True, "space_after": 2},
        {"t": "0.14 / 0.06 / 0.17 / 0.04 / 0.03", "size": 14.5, "color": TEXT, "spacing": 1.1, "space_after": 12},
        {"t": "* 47 cm test", "size": 12, "color": TEAL, "bold": True, "space_after": 2},
        {"t": "encoder read 5.3\u00b0 low; corrected on experimental angles - and never trained on.", "size": 13.5, "color": TEXT, "spacing": 1.12, "space_after": 12},
        {"t": "No test is garbage; every test predicts.", "size": 14, "color": BROWN, "bold": True, "spacing": 1.1},
    ])
    pageno(s, n)

def s21_extensor(prs, n):
    s = content_slide(prs, "Aim 2 \u00b7 pinned extensor",
        "The extensor wraps the joint - and \u03c7\u2083 recovers the length it loses there",
        cite="minimizeExt10mmX3_results_20260910_noT3.mat, driver pick 1 (six-test validation pool)",
        notes="chi0 flips sign relative to flexor: -10.1 vs +8.9 mm. Held-out RMSE 0.56-1.88 vs baseline 0.95-3.62: improved in every pool test.")
    pic(s, os.path.join(FIG, "ExtPin_group.png"), box=(0.6, 1.9, 7.2, 4.85))
    stats = [
        ("\u03c7\u2080 = \u221210.1 mm", "extensor offset: same magnitude scale, OPPOSITE sign to flexor (+8.9 mm)", RED),
        ("\u03c7\u2083 = 0.621", "fractional wrap loss acting on arc length \u00d7 \u03b6\u00b2", BROWN),
        ("0.56-1.88 N\u00b7m", "held-out RMSE with corrections", GREEN),
        ("vs 0.95-3.62 N\u00b7m", "uncorrected rigid baseline - improved in every pool test (6/6)", TEAL),
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
    s = content_slide(prs, "Aim 2 \u00b7 values of record",
        "One stiffness pair of record spans the flexor and extensor configurations",
        cite="Dissertation App. C, Table C.2; flexor mat: ...20260908_2brkt_2trans_noT3, extensor mat: ...20260910_noT3",
        notes="Advisor requirement satisfied: one (chi1, chi2) across configurations. Note both rows are operating points on ONE flexor Pareto front - the data admit a family; this is the honest statement.")
    hdr = ["Term", "Value", "Source", "Row"]
    rows = [
        ("Flexor \u03c7\u2080", "+8.9 mm", "flexor mat", "107"),
        ("Flexor \u03c7\u2081", "5.62\u00d710\u2074 N/m", "flexor mat", "107"),
        ("Flexor \u03c7\u2082", "1.85\u00d710\u2074 N/m", "flexor mat", "107"),
        ("Extensor \u03c7\u2080", "\u221210.1 mm", "extensor mat", "1"),
        ("Extensor \u03c7\u2081 (locked)", "4.354\u00d710\u2074 N/m", "flexor mat", "1"),
        ("Extensor \u03c7\u2082 (locked)", "1.701\u00d710\u2074 N/m", "flexor mat", "1"),
        ("\u03c7\u2083", "0.621", "extensor mat", "1"),
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
            if r[0] == "\u03c7\u2083" and i <= 1: colr = BROWN
            tx(s, xx, y, cw[i], 0.4, val, size=14 if i != 1 else 13.5, color=colr, bold=bold,
               align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
            xx += cw[i]
        y += 0.5
    card(s, 8.85, 1.95, 3.9, 4.35, fill=CARD2)
    rich(s, 9.1, 2.15, 3.45, 4.05, [
        {"t": "READ THIS TABLE LIKE THIS", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "Rows 107 and 1 are two operating points on one and the same flexor Pareto front.", "size": 13.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "The pinned-flexor data admit a family of stiffness pairs of comparable fit; the adopted pair is the one that also performs across configurations.", "size": 13.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "The extensor never searches stiffness - it inherits it.", "size": 13.5, "color": BROWN, "bold": True, "spacing": 1.1},
    ])
    pageno(s, n)

def s23_bioflexor(prs, n):
    s = content_slide(prs, "Aim 2 \u00b7 biomimetic knee",
        "The biomimetic knee meets or exceeds human flexor torque over most of the range of motion",
        cite="Bolen et al., in preparation 2025, Fig. 9; human curve: Gait2392 (Delp 2007; Seth 2018)",
        notes="Panel B is the headline. Be upfront about the 60-100 deg shortfall (78-103 at +/-20%). The femur note usually gets a laugh: at 620 kPa we exceeded the ultimate strength of the test femur, so higher-pressure points were measured at reduced P.")
    pic(s, os.path.join(FIG, "FlxGrp.png"), box=(0.55, 1.95, 10.4, 4.9))
    card(s, 11.15, 2.0, 1.95, 4.7, fill=CARD)
    rich(s, 11.32, 2.18, 1.65, 4.4, [
        {"t": "SHORTFALL", "size": 11, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "60-100\u00b0\n(\u00b120%: 78-103\u00b0)", "size": 13, "color": RED, "bold": True, "spacing": 1.05, "space_after": 10},
        {"t": "\u03c620 mm, 41.5 cm\nRMSE 1.80\nFVU 0.05", "size": 12.5, "color": TEXT, "spacing": 1.05, "space_after": 10},
        {"t": "At 620 kPa the BPA exceeded the test femur's ultimate strength - remaining points taken at reduced P.", "size": 11.5, "color": MUTED, "italic": True, "spacing": 1.08},
    ])
    pageno(s, n)

def s24_bioext(prs, n):
    s = content_slide(prs, "Aim 2 \u00b7 biomimetic knee",
        "Biomimetic extensor: the corrections halve the error, but variance is still under-predicted",
        cite="l\u2080 = 51.8 cm test; values of record: RMSE 2.20 / FVU 2.35 / max resid 4.73 (baseline 3.03 / 4.49 / 5.59)",
        notes="Honesty slide. Improvement is unambiguous; validation is NOT claimed for the bio-extensor. FVU > 1 means the corrected model under-predicts the measured torque variance.")
    pic(s, os.path.join(FIG, "Ext10mm_52cm.png"), box=(0.6, 1.95, 7.6, 3.4))
    tx(s, 0.6, 5.4, 7.6, 0.35, "\u03c610 mm extensor on the four-bar knee, l\u2080 = 51.8 cm", size=12, color=MUTED, align=PP_ALIGN.CENTER)
    stats = [("2.20", "RMSE (N\u00b7m)\nvs 3.03 baseline", GREEN), ("2.35", "FVU\nvs 4.49 baseline", RED),
             ("4.73", "max resid (N\u00b7m)\nvs 5.59 baseline", GREEN)]
    x = 8.55
    for v, lab, c in stats:
        card(s, x, 2.0, 1.45, 2.3, fill=CARD)
        rich(s, x + 0.08, 2.25, 1.3, 1.9, [
            {"t": v, "size": 20, "color": c, "bold": True, "align": PP_ALIGN.CENTER, "space_after": 2},
            {"t": lab, "size": 10.5, "color": MUTED, "align": PP_ALIGN.CENTER, "spacing": 1.0},
        ])
        x += 1.55
    card(s, 8.55, 4.55, 4.2, 1.6, fill=CARD2)
    rich(s, 8.8, 4.7, 3.75, 1.4, [
        {"t": "CLAIM CONTROL", "size": 11, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "FVU > 1: no validated bio-extensor prediction is claimed. Improvement over baseline is unambiguous; closure is not.", "size": 12.5, "color": TEXT, "spacing": 1.1},
    ])
    pageno(s, n)

def s25_redesign(prs, n):
    s = content_slide(prs, "Aim 2 \u2192 robot redesign",
        "The values of record now drive the full leg redesign optimizer",
        cite="Opt_run.m / Opt_run_Ext.m; latest: Vas_Pam_20mm_Result_20260910_0528.mat (feasible)",
        notes="Bridge slide to current work: same Xi loaded directly from the results mats by the context builders.")
    card(s, MX, 1.95, 6.0, 4.0, fill=CARD)
    rich(s, MX + 0.3, 2.15, 5.4, 3.7, [
        {"t": "COST (route redesign)", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "J = 10\u2075\u00b7J_worst + 10\u00b3\u00b7J_torque + aux + regularizers", "size": 15, "color": TEXT, "spacing": 1.1, "space_after": 8},
        {"t": "Targets", "size": 12, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "human knee torque from Gait2392 - flexor vs biceps femoris short head, extensor vs vastus medialis", "size": 13.5, "color": TEXT, "spacing": 1.08, "space_after": 8},
        {"t": "Solver", "size": 12, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "surrogateopt (\u22487,000 evals) \u2192 patternsearch (\u224815,000), parallel", "size": 13.5, "color": TEXT, "spacing": 1.08},
    ])
    card(s, 6.95, 1.95, 5.8, 4.0, fill=CARD2)
    rich(s, 7.25, 2.15, 5.2, 3.7, [
        {"t": "LATEST RESULT - EXTENSOR REDESIGN", "size": 12, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "FEASIBLE with the \u03c7 values of record", "size": 16, "color": GREEN, "bold": True, "space_after": 10},
        {"t": "thinnest torque margin  +0.063%  at \u221218.9\u00b0", "size": 14, "color": TEXT, "space_after": 5},
        {"t": "contraction  0.69 \u00b7 KMAX  (within limit)", "size": 14, "color": TEXT, "space_after": 5},
        {"t": "max path length  0.832 m", "size": 14, "color": TEXT, "space_after": 5},
        {"t": "binding constraint \u22120.000442 (near-active at optimum)", "size": 12.5, "color": MUTED, "spacing": 1.05},
    ])
    card(s, MX, 6.15, 12.2, 0.75, fill=CARD)
    tx(s, MX + 0.3, 6.28, 11.6, 0.55,
       "flexor target carries a +1% torque margin above the human curve \u00b7 extensor target is the raw curve \u00b7 flexor routes mirrored about the mid-sagittal plane \u00b7 Xi loaded directly from the results mats",
       size=12.5, color=TEXT, spacing=1.05)
    pageno(s, n)

def build(prs, start=15):
    s15_stands(prs, 15)
    s16_hybrid(prs, 16)
    s17_corrections(prs, 17)
    s18_identification(prs, 18)
    s19_flexor_opt(prs, 19)
    s20_flexor_numbers(prs, 20)
    s21_extensor(prs, 21)
    s22_xi_record(prs, 22)
    s23_bioflexor(prs, 23)
    s24_bioext(prs, 24)
    s25_redesign(prs, 25)
    return 25
