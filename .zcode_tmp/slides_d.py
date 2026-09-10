"""Appendix slides A1-A10 (deck slides 38-47)."""
import os
from pptx.enum.text import PP_ALIGN
from deck_lib import *

def a1_force_model(prs, n):
    s = content_slide(prs, "Appendix A \u00b7 model",
        "The complete actuator force model",
        cite="Dissertation Eqs. 3.1-3.5; implementation: Functions/f_festo.m, maxBPAforce.m",
        notes="One-look reference card.")
    rows = [
        ("strain measures", "strain", 0.52),
        ("normalized surface", "fstar", 0.52),
        ("maximum force (general P)", "fmax3d", 0.5),
        ("isometric force", "force", 0.5),
    ]
    y = 2.0
    for lab, name, h in rows:
        card(s, MX, y, 12.2, h + 0.55, fill=CARD if name != "fstar" else CARD2)
        tx(s, MX + 0.25, y + 0.12, 3.0, 0.35, lab.upper(), size=11, color=TEAL, bold=True)
        eq(s, name, MX + 0.4, y + 0.42, h=h)
        y += h + 0.75
    pageno(s, n)

def a2_coefficients(prs, n):
    s = content_slide(prs, "Appendix B \u00b7 coefficients",
        "Coefficient values with 95% confidence intervals",
        cite="Bolen et al. 2026, Actuators, Tables 1-2",
        notes="Tight CIs relative to coefficient magnitudes; fitted with fittype NLLS + LAR robustness, lower bounds [0,0,0].")
    tx(s, MX, 2.0, 6.0, 0.4, "F\u2086\u2082\u2080 arctan fits", size=14, color=TEAL, bold=True)
    hdr = ["", "b\u2081 (N)", "b\u2082 (1/m)", "adj. R\u00b2", "RMSE (N)"]
    rows = [["\u03c610 mm", "303.5 (300, 308)", "19.03 (17.5, 20.6)", "0.9854", "14.7"],
            ["\u03c620 mm", "922.4 (914, 931)", "15.37 (14.8, 16.0)", "0.9945", "23.8"]]
    cw = [1.15, 2.0, 2.2, 1.0, 1.05]
    y = 2.45
    for r_i, r in enumerate([hdr] + rows):
        xx = MX
        for i, val in enumerate(r):
            st = {"size": 12, "color": TEAL if r_i == 0 else TEXT, "bold": r_i == 0 or i == 0,
                  "align": PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER}
            tx(s, xx, y, cw[i] + (0.3 if i == 0 else 0), 0.45, val, **st)
            xx += cw[i] + (0.3 if i == 0 else 0)
        y += 0.52
    tx(s, MX, 3.9, 6.0, 0.4, "F* surface coefficients", size=14, color=TEAL, bold=True)
    hdr2 = ["", "c\u2080", "c\u2081", "c\u2082", "R\u00b2 fit / val."]
    rows2 = [["\u03c610 mm", "0.5682", "4.254", "0.5597", "0.9998 / 0.9994"],
             ["\u03c620 mm", "0.2579", "6.477", "1.321", "0.992 / 0.9943"]]
    y = 4.35
    for r_i, r in enumerate([hdr2] + rows2):
        xx = MX
        for i, val in enumerate(r):
            st = {"size": 12, "color": TEAL if r_i == 0 else TEXT, "bold": r_i == 0 or i == 0,
                  "align": PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER}
            tx(s, xx, y, cw[i] + (0.3 if i == 0 else 0) + (0.85 if i == 4 else 0), 0.45, val, **st)
            xx += cw[i] + (0.3 if i == 0 else 0) + (0.85 if i == 4 else 0)
        y += 0.52
    card(s, 9.15, 2.45, 3.6, 3.6, fill=CARD)
    rich(s, 9.38, 2.65, 3.15, 3.3, [
        {"t": "ALSO ON RECORD", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "F\u2086\u2082\u2080 arctan fits: RMSE vs the FESTO TOOL: 189.9 N (10 mm), 668.4 N (20 mm) - the tool cannot track resting length.", "size": 12.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "A 40 mm fit exists in code (c\u2080 0.1224, c\u2081 10.47, c\u2082 2.023) but is unused in the dissertation.", "size": 12.5, "color": MUTED, "spacing": 1.12},
    ])
    pageno(s, n)

def a3_model_compare(prs, n):
    s = content_slide(prs, "Appendix C \u00b7 comparison",
        "Full goodness-of-fit against published BPA models",
        cite="Bolen et al. 2026, Actuators, Table 3; RMSE / max error in N; FVU = fraction of variance unexplained",
        notes="Sarosi params exist per length; Martens-Boblan is geometry-sensitive - note the negative-force failure at 450 mm.")
    hdr = ["", "Bolen RMSE", "FVU", "max", "S\u00e1rosi RMSE", "FVU", "max", "M&B RMSE", "FVU", "max"]
    rows = [
        ["10 mm, 257 mm", "22.9", "0.04", "46.0", "143.8", "1.39", "235.8", "73.8", "0.37", "122.9"],
        ["10 mm, 233 mm", "16.8", "0.02", "30.3", "132.7", "1.17", "222.6", "62.6", "0.26", "101.2"],
        ["20 mm, 300 mm", "35.6", "0.01", "114.0", "132.1", "0.20", "313.4", "72.8", "0.06", "137.8"],
        ["20 mm, 450 mm", "67.6", "0.03", "166.6", "177.7", "0.18", "294.8", "862.1", "4.20", "1321.1"],
    ]
    cw = [1.9] + [1.15] * 3 + [1.35] * 6
    cw = [1.9, 1.15, 0.85, 1.0, 1.35, 0.85, 1.0, 1.35, 0.85, 1.0]
    xx = MX; y = 2.15
    for i, h_ in enumerate(hdr):
        col = TEAL if i else MUTED
        tx(s, xx, y, cw[i], 0.4, h_, size=12, color=col, bold=True,
           align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
        xx += cw[i]
    y = 2.62
    for r in rows:
        xx = MX
        for i, val in enumerate(r):
            if i == 0:
                tx(s, xx, y + 0.03, cw[i], 0.4, val, size=12.5, color=TEXT)
            elif 1 <= i <= 3:
                tx(s, xx, y, cw[i], 0.4, val, size=12.5, color=GREEN, bold=(i == 1), align=PP_ALIGN.CENTER)
            else:
                tx(s, xx, y, cw[i], 0.4, val, size=12, color=MUTED, align=PP_ALIGN.CENTER)
            xx += cw[i]
        y += 0.52
    card(s, MX, 5.1, 12.2, 1.3, fill=CARD)
    rich(s, MX + 0.3, 5.28, 11.6, 1.0, [
        {"t": "Takeaway: fixed-parameter models fail when moved off their fit length;", "size": 14, "color": TEXT, "spacing": 1.1},
        {"t": "the normalized fit holds across lengths because each muscle carries its own F\u2086\u2082\u2080 and \u03b5\u2086\u2082\u2080.", "size": 14, "color": TEXT, "spacing": 1.1},
    ])
    pageno(s, n)

def a4_compliance(prs, n):
    s = content_slide(prs, "Appendix D \u00b7 compliance chain",
        "From bracket stiffness matrix to solved deflection",
        cite="Dissertation App. A \u00a7A.3; minimizeFlxPin2brk.m / MonoPam_mult.m",
        notes="For the controls professor: this is a one-dimensional nonlinear root-find per configuration, solved by fzero on a bracketed interval.")
    steps = [
        ("1 \u00b7 project compliance on the force direction", "kbr", CARD),
        ("2 \u00b7 add artificial tendon in series", "series", CARD),
    ]
    y = 2.05
    for lab, name, fill in steps:
        card(s, MX, y, 12.2, 1.15, fill=fill)
        tx(s, MX + 0.25, y + 0.12, 4.4, 0.4, lab, size=11.5, color=TEAL, bold=True)
        eq(s, name, MX + 0.5, y + 0.42, h=0.55)
        y += 1.35
    card(s, MX, y, 12.2, 2.0, fill=CARD2)
    rich(s, MX + 0.3, y + 0.15, 11.6, 1.8, [
        {"t": "3 \u00b7 solve the force balance, then propagate", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "\u03b4 from F\u2086\u2082\u2080\u00b7F*(\u03b5*(\u03b4), P*) \u2212 k\u2091\u2094\u00b7\u03b4 = 0 (fzero, bracket [0, l\u2098 \u2212 l\u2086\u2082\u2080]; \u03b4 \u2261 0 if uncorrected \u03b5* > 1)", "size": 13.5, "color": TEXT, "spacing": 1.1, "space_after": 4},
        {"t": "then vector deflection \u2192 attachment moves \u2192 path length, force direction, moment arm, torque recomputed", "size": 13.5, "color": TEXT, "spacing": 1.1, "space_after": 4},
        {"t": "two-BPA flexor grouping: coupled balance shares one common force, solved on [0, F\u2086\u2082\u2080]", "size": 12.5, "color": MUTED, "spacing": 1.1},
    ])
    pageno(s, n)

def a5_wraploss(prs, n):
    s = content_slide(prs, "Appendix E \u00b7 wrap loss",
        "Wrapping-loss implementation, per contact arc",
        cite="minimizeExtX3.m (Contraction, lines 327-370); MonoPamDataExplicit_balanceX3.m",
        notes="The loss enters only the FORCE-PRODUCING strain; the kinematic (prediction) strain excludes it - worth stating if asked about double-counting.")
    card(s, MX, 2.05, 6.0, 2.5, fill=CARD)
    rich(s, MX + 0.25, 2.25, 5.5, 2.2, [
        {"t": "PINNED EXTENSOR", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "two physical contact arcs:", "size": 13.5, "color": TEXT, "space_after": 3},
        {"t": "R\u2081 = 12 mm from 27\u00b0 to \u221223\u00b0", "size": 14, "color": TEXT, "space_after": 3},
        {"t": "R\u2082 = 40 mm on the second bracket", "size": 14, "color": TEXT, "space_after": 5},
        {"t": "arc-length normalizers: 27.75/(ang\u2081\u2212ang\u2082), 92.46/(ang\u2082+120\u00b0)", "size": 12, "color": MUTED, "spacing": 1.1},
    ])
    card(s, 6.75, 2.05, 6.0, 2.5, fill=CARD)
    rich(s, 7.0, 2.25, 5.5, 2.2, [
        {"t": "BIOMIMETIC CONFIGURATIONS", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "\u0394l = \u03c7\u2083 \u00b7 b(\u03b8\u2096) \u00b7 (1 \u2212 \u03b5*)\u00b2", "size": 15, "color": TEXT, "space_after": 5},
        {"t": "b(\u03b8\u2096) = routed bend measure (arc length swept by the wrap)", "size": 13.5, "color": TEXT, "spacing": 1.1, "space_after": 3},
        {"t": "w Rap = 25 mm bend radius in the route builder", "size": 12, "color": MUTED, "spacing": 1.1},
    ])
    card(s, MX, 4.85, 12.2, 1.35, fill=CARD2)
    rich(s, MX + 0.3, 5.02, 11.6, 1.1, [
        {"t": "Applied to the force-producing strain only", "size": 14.5, "color": BROWN, "bold": True, "space_after": 3},
        {"t": "the kinematic strain that predicts muscle length excludes \u0394l - the loss changes force, not recorded geometry.", "size": 13.5, "color": TEXT, "spacing": 1.1},
    ])
    pageno(s, n)

def a6_protocol(prs, n):
    s = content_slide(prs, "Appendix F \u00b7 protocol",
        "Identification campaign details and data decisions",
        cite="minimizeFlxPin10mm_2brk.m; Collect_ExtPinX3_sweep.m; Dig_foldLeverage.m",
        notes="Two data rulings Ben made: 47 cm encoder fix, and the excluded extensor tests 3/4/9. Fold study: all folds land on the same score plateau.")
    items = [
        ("Encoder correction", "the 47 cm flexor test's encoder read \u22485.3\u00b0 low; +5.3\u00b0 applied to experimental angles only (phiD untouched), at build time"),
        ("Training / validation", "flexor: 4 train, 47 cm pure hold-out \u00b7 extensor: six-test pool, tests 3/4/9 excluded a priori"),
        ("Fold-structure check", "every fold lands on the same score plateau; pooled fronts \u224880% cross-fold duplicates - fold choice does not change conclusions"),
        ("Reproducibility", "re-running the extensor CV with a fresh GA seed reproduced the same 4-member front exactly"),
    ]
    y = 2.05
    for t, sub in items:
        card(s, MX, y, 12.2, 1.06, fill=CARD)
        rich(s, MX + 0.25, y + 0.12, 11.7, 0.9, [
            {"t": t, "size": 14, "color": BROWN, "bold": True, "space_after": 1},
            {"t": sub, "size": 12.5, "color": TEXT, "spacing": 1.02},
        ])
        y += 1.18
    pageno(s, n)

def a7_frame_invariance(prs, n):
    s = content_slide(prs, "Appendix G \u00b7 invariance",
        "Frame convention is irrelevant for y-symmetric stiffness arrays",
        cite="proven numerically 2026-09-08: 1-rotation vs 2-rotation frames agree to ~1e-17",
        notes="For the math professor: the compliance scalar is invariant because the second rotation is about the symmetry axis of K.")
    card(s, MX, 2.1, 7.3, 4.3, fill=CARD)
    rich(s, MX + 0.3, 2.35, 6.7, 3.9, [
        {"t": "CLAIM", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "Let K = diag(a, b, a) - equal stiffness about x and z. For any frame whose second rotation is about y, the projected scalar stiffness", "size": 14, "color": TEXT, "spacing": 1.15, "space_after": 8},
        {"t": "k(\u00fb) = (\u00fb\u1d40 K\u207b\u00b9 \u00fb)\u207b\u00b9", "size": 17, "color": BROWN, "bold": True, "space_after": 8},
        {"t": "is unchanged: a rotation about y mixes only the two coordinates with equal compliance 1/a, and \u00fb\u1d40 K\u207b\u00b9 \u00fb is invariant under orthogonal transforms within that eigenspace.", "size": 14, "color": TEXT, "spacing": 1.15, "space_after": 8},
        {"t": "Verified: two-rotation vs one-rotation implementations agree to \u22481\u00d710\u207b\u00b9\u2077.", "size": 13, "color": GREEN, "bold": True, "spacing": 1.1},
    ])
    card(s, 8.15, 2.1, 4.6, 4.3, fill=CARD2)
    rich(s, 8.4, 2.35, 4.1, 3.9, [
        {"t": "SCOPE", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "Holds for y-symmetric arrays only (K\u2093 = K\u1DA3).", "size": 13.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "The extensor's arch bracket uses K = diag(\u03c7\u2082, \u03c7\u2081, \u03c7\u2082)-type orderings that are NOT symmetric; those axis assignments are a modeling choice carried on physical grounds (buckling argument), not derived from this proof.", "size": 13.5, "color": TEXT, "spacing": 1.12},
    ])
    pageno(s, n)

def a8_bracket_point(prs, n):
    s = content_slide(prs, "Appendix H \u00b7 robustness",
        "Extensor identification is robust to the bracket reference point",
        cite="minimizeExtX3.m Pbr variants, 1trans CVs; both fronts captured",
        notes="Rib midpoint (adopted) vs lower bolt hole - essentially identical outcomes, so the choice does not drive the answer.")
    hdr = ["", "rib midpoint (adopted)", "lower bolt hole"]
    rows = [
        ["bracket point", "[-3.84, -46.44, 62.5] mm", "[-6.26, -29.69, 75.06] mm"],
        ["pick \u03c7\u2083", "0.62", "0.60"],
        ["pick \u03c7\u2080", "\u221210.1 mm", "\u221211.6 mm"],
        ["bio-ext RMSE / FVU / max", "2.195 / 2.35 / 4.73", "2.204 / 2.369 / 4.727"],
        ["filter pass", "72 / 90", "83 / 90"],
    ]
    cw = [3.0, 4.4, 4.4]
    y = 2.2
    for r_i, r in enumerate([hdr] + rows):
        xx = MX
        for i, val in enumerate(r):
            if r_i == 0:
                tx(s, xx, y, cw[i], 0.4, val, size=13, color=TEAL, bold=True,
                   align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
            else:
                bold = i == 0
                colr = TEXT if i == 0 else (GREEN if i == 1 else MUTED)
                tx(s, xx, y, cw[i], 0.45, val, size=14 if i else 13, color=colr, bold=bold,
                   align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
            xx += cw[i]
        y += 0.56
    card(s, MX, 5.5, 12.2, 1.0, fill=CARD2)
    tx(s, MX + 0.3, 5.68, 11.6, 0.7, "Conclusion: the extensor \u03c7 identification does not hinge on where the bracket point is placed.",
       size=14.5, color=TEXT, bold=True)
    pageno(s, n)

def a9_moment_arm(prs, n):
    s = content_slide(prs, "Appendix I \u00b7 kinematics",
        "Moment arm and torque conventions",
        cite="Young 2019 (analyzing moment arms); MonoPamDataExplicit.m lines 156-168",
        notes="Scalar arm is in the tibial frame about z; compressive values discarded at eps < -0.02.")
    card(s, MX, 2.15, 12.2, 1.5, fill=CARD)
    eq(s, "torque", 1.6, 2.62, h=0.6)
    card(s, MX, 3.95, 12.2, 2.4, fill=CARD2)
    rich(s, MX + 0.3, 4.15, 11.6, 2.1, [
        {"t": "CONVENTIONS", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "r = perpendicular foot of the joint origin onto the line of action; identical result to the cross-product projection of Eq. 3.6", "size": 14, "color": TEXT, "spacing": 1.12, "space_after": 5},
        {"t": "scalar moment arm r\u2096 = hypot(in-plane components) in the tibial body frame, about the knee z-axis", "size": 14, "color": TEXT, "spacing": 1.12, "space_after": 5},
        {"t": "torque evaluated only in tension (entries discarded for \u03b5 < \u22120.02)", "size": 14, "color": TEXT, "spacing": 1.12, "space_after": 5},
        {"t": "coupled pathpoints stay equality-constrained (not welded) so moment arms keep their full range", "size": 12.5, "color": MUTED, "spacing": 1.1},
    ])
    pageno(s, n)

def a10_sns_impl(prs, n):
    s = content_slide(prs, "Appendix J \u00b7 implementation",
        "Synthetic nervous system: library blocks and reflex topology",
        cite="Simulink SNS library + KneeReflexDemo; AnimatLab neuron params: App. B of the dissertation",
        notes="Two parallel implementations: AnimatLab (bipedal prototyping) and Simulink/MuJoCo (pipeline). The 7-block library runs in base Simulink.")
    pic(s, os.path.join(FIG, "sns_library.png"), box=(0.6, 1.95, 6.0, 3.3))
    tx(s, 0.6, 5.3, 6.0, 0.35, "seven masked blocks; E/I icons auto-draw from synapse sign", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    pic(s, os.path.join(FIG, "sns_reflex.png"), box=(6.9, 1.95, 5.85, 3.3))
    tx(s, 6.9, 5.3, 5.85, 0.35, "reduced-order knee reflex: Ia/Ib, antagonist MNs, Renshaw-style recursion", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    card(s, 0.6, 5.85, 12.15, 0.85, fill=CARD)
    tx(s, 0.9, 6.02, 11.6, 0.6, "CAD route settled: SolidWorks \u2192 URDF \u2192 smimport (verified); tendon parts modeled, insertion pending",
       size=13, color=TEXT)
    pageno(s, n)

def build(prs):
    a1_force_model(prs, 38)
    a2_coefficients(prs, 39)
    a3_model_compare(prs, 40)
    a4_compliance(prs, 41)
    a5_wraploss(prs, 42)
    a6_protocol(prs, 43)
    a7_frame_invariance(prs, 44)
    a8_bracket_point(prs, 45)
    a9_moment_arm(prs, 46)
    a10_sns_impl(prs, 47)
    return 47
