"""Appendix slides A1-A11 (deck batch d: start .. start+10)."""
import os
from pptx.enum.text import PP_ALIGN
from deck_lib import *

def a1_force_model(prs, n):
    s = content_slide(prs, "Appendix A \u00b7 model",
        "The complete actuator force model",
        cite="Dissertation Ch. 3, Eqs. (3.1), (3.4)\u2013(3.6); reference implementation: App. A \u00a7A.3 (Functions/f_festo.m, maxBPAforce.m)",
        notes="Reference card for the actuator force model of Chapter 3 (muscles). Chain: strain (3.1) normalized by the recorded free-load contraction (3.2)-(3.3); the F* surface (3.4) is fitted per diameter; F620 (3.5) carries the resting-length dependence; the product (3.6) is isometric force. Anticipated question: why can F* be clamped to zero? Answer: the zero-clamp is applied by the callers (e.g. festo4.m line 35), not inside f_festo.m itself - documented in App. A sect. A.3.")
    rows = [
        ("strain measures \u00b7 eq. (3.1)", "strain", 0.52),
        ("normalized surface \u00b7 eq. (3.4)", "fstar", 0.52),
        ("maximum force (general P) \u00b7 eq. (3.5)", "fmax3d", 0.5),
        ("isometric force \u00b7 eq. (3.6)", "force", 0.5),
    ]
    y = 2.0
    for lab, name, h in rows:
        card(s, MX, y, 12.2, h + 0.55, fill=CARD if name != "fstar" else CARD2)
        tx(s, MX + 0.25, y + 0.12, 5.5, 0.35, lab.upper(), size=11, color=TEAL, bold=True)
        eq(s, name, MX + 0.4, y + 0.42, h=h)
        y += h + 0.75
    pageno(s, n)

def a2_coefficients(prs, n):
    s = content_slide(prs, "Appendix B \u00b7 coefficients",
        "Coefficient values with 95% confidence intervals",
        cite="Ch. 3 Tables 3.1\u20133.2 (published as Bolen et al. 2026, Actuators, Tables 1\u20132); fitting: App. A \u00a7A.3",
        notes="Coefficients of record. Fitted with fittype nonlinear least squares + least-absolute-residual robustness, lower bounds [0,0,0] (App. A sect. A.3, bpaFits.m). The 40 mm fit exists in code but is unused in the dissertation. Anticipated question: how tight are the CIs? Answer: b1 for the 10 mm actuator is 303.5 N with CI (300, 308) - roughly 1.3% of the value; c0 for the 10 mm surface is 0.5682 with CI (0.5584, 0.578).")
    tx(s, MX, 2.0, 6.0, 0.4, "F\u2086\u2082\u2080 arctan fits", size=14, color=TEAL, bold=True)
    hdr = ["", "b\u2081 (N)", "b\u2082 (1/m)", "adj. R\u00b2", "RMSE (N)"]
    rows = [["\u230010 mm", "303.5 (300, 308)", "19.03 (17.48, 20.57)", "0.9854", "14.72"],
            ["\u230020 mm", "922.4 (914.2, 930.7)", "15.37 (14.75, 15.98)", "0.9945", "23.83"]]
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
    rows2 = [["\u230010 mm", "0.5682", "4.254", "0.5597", "0.9998 / 0.9994"],
             ["\u230020 mm", "0.2579", "6.477", "1.321", "0.992 / 0.9943"]]
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
        {"t": "F\u2086\u2082\u2080 arctan fits, RMSE vs the Festo tool: 189.9 N (10 mm) and 668.4 N (20 mm); the tool cannot track resting length.", "size": 12.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "A 40 mm fit exists in code (c\u2080 0.1224, c\u2081 10.47, c\u2082 2.023) but is unused in the dissertation.", "size": 12.5, "color": MUTED, "spacing": 1.12},
    ])
    pageno(s, n)

def a3_model_compare(prs, n):
    s = content_slide(prs, "Appendix C \u00b7 comparison",
        "Full goodness-of-fit against published BPA models",
        cite="Ch. 3 Table 3.3; RMSE / max error in N; FVU = fraction of variance unexplained",
        notes="Comparison of record against Sarosi and Martens-and-Boblan. Fixed-parameter models degrade when moved off their fit length; Martens-Boblan starts to predict negative force at the 20 mm, 450 mm point. Anticipated question: why does Sarosi fail here? Answer: its parameters are length-specific (200 mm and 400 mm for the 20 mm actuator), and applying them at our 300 mm and 450 mm data over-predicts force at higher pressures and low contraction (Ch. 3, model-comparison discussion).")
    hdr = ["", "Bolen RMSE", "FVU", "max", "S\u00e1rosi RMSE", "FVU", "max", "M&B RMSE", "FVU", "max"]
    rows = [
        ["\u230010 \u00b7 257 mm", "22.9", "0.04", "46.0", "143.8", "1.39", "235.8", "73.8", "0.37", "122.9"],
        ["\u230010 \u00b7 233 mm", "16.8", "0.02", "30.3", "132.7", "1.17", "222.6", "62.6", "0.26", "101.2"],
        ["\u230020 \u00b7 300 mm", "35.6", "0.01", "114.0", "132.1", "0.20", "313.4", "72.8", "0.06", "137.8"],
        ["\u230020 \u00b7 450 mm", "67.6", "0.03", "166.6", "177.7", "0.18", "294.8", "862.1", "4.20", "1321.1"],
    ]
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
        {"t": "Takeaway: fixed-parameter models fail when moved off their fit length; ", "size": 14, "color": TEXT, "spacing": 1.1},
        {"t": "the normalized fit holds across lengths because each muscle carries its own F\u2086\u2082\u2080 and \u03b5\u2086\u2082\u2080.", "size": 14, "color": TEXT, "spacing": 1.1},
    ])
    pageno(s, n)

def a4_compliance(prs, n):
    s = content_slide(prs, "Appendix D \u00b7 compliance chain",
        "From bracket stiffness matrix to solved deflection",
        cite="Ch. 4 Eqs. (4.18), (4.20), (4.22); implementation: App. A \u00a7A.5; frames and arrays: App. C \u00a7\u00a7C.3\u2013C.4",
        notes="One-dimensional nonlinear root-find per configuration: fzero on [0, l_m - l_620], converged to machine precision by bracket collapse; safeguards return delta = 0 at saturated strain or when no sign change appears (App. A sect. A.5). Anticipated question: is the series reduction exact? Answer: only under the stated conditions - both elements carry the same force, displacements add kinematically along the line of action, u-hat is fixed within one iteration then updated, elements are linear and moment-free - each condition is listed in Ch. 4 (Improved Torque Model).")
    steps = [
        ("1 \u00b7 project compliance on the force direction \u00b7 eq. (4.18)", "kbr", CARD),
        ("2 \u00b7 add artificial tendon in series \u00b7 eq. (4.20)", "series", CARD),
    ]
    y = 2.05
    for lab, name, fill in steps:
        card(s, MX, y, 12.2, 1.15, fill=fill)
        tx(s, MX + 0.25, y + 0.12, 6.4, 0.4, lab, size=11.5, color=TEAL, bold=True)
        eq(s, name, MX + 0.5, y + 0.42, h=0.55)
        y += 1.35
    card(s, MX, y, 12.2, 2.0, fill=CARD2)
    rich(s, MX + 0.3, y + 0.15, 11.6, 1.8, [
        {"t": "3 \u00b7 solve the force balance, then propagate \u00b7 eq. (4.22)", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "\u03b4 from F\u2086\u2082\u2080\u00b7F*(\u03b5*(\u03b4), P*) \u2212 k_eq\u00b7\u03b4 = 0 (fzero, bracket [0, l\u2098 \u2212 l\u2086\u2082\u2080]; \u03b4 \u2261 0 if uncorrected \u03b5* > 1)", "size": 13.5, "color": TEXT, "spacing": 1.1, "space_after": 4},
        {"t": "then vector deflection \u2192 attachment moves \u2192 path length, force direction, moment arm, torque recomputed", "size": 13.5, "color": TEXT, "spacing": 1.1, "space_after": 4},
        {"t": "two-BPA flexor grouping: coupled balance shares one common force, solved on [0, F\u2086\u2082\u2080]", "size": 12.5, "color": MUTED, "spacing": 1.1},
    ])
    pageno(s, n)

def a5_wraploss(prs, n):
    s = content_slide(prs, "Appendix E \u00b7 wrap loss",
        "Wrapping-loss implementation, per contact arc",
        cite="Ch. 4 Eqs. (4.24)\u2013(4.26); implementation: App. A \u00a7A.6 (minimizeExtX3.m, Contraction); constants: App. C \u00a7C.4",
        notes="The loss enters only the force-producing strain; the kinematic (prediction) strain excludes it. Anticipated question: is the loss double-counted? Answer: no - the kinematic strain that predicts muscle length excludes delta-l, so the loss changes force, not recorded geometry (App. A sect. A.6).")
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
        {"t": "route builder assumes a 25 mm minimum bend radius", "size": 12, "color": MUTED, "spacing": 1.1},
    ])
    card(s, MX, 4.85, 12.2, 1.35, fill=CARD2)
    rich(s, MX + 0.3, 5.02, 11.6, 1.1, [
        {"t": "Applied to the force-producing strain only", "size": 14.5, "color": BROWN, "bold": True, "space_after": 3},
        {"t": "the kinematic strain that predicts muscle length excludes \u0394l; the loss changes force, not recorded geometry.", "size": 13.5, "color": TEXT, "spacing": 1.1},
    ])
    pageno(s, n)

def a6_protocol(prs, n):
    s = content_slide(prs, "Appendix F \u00b7 protocol",
        "Identification campaign details and data decisions",
        cite="Ch. 4 \u00a74.1 (identification); test constants: App. C \u00a7C.1; bounds and adopted values: App. C \u00a7\u00a7C.4\u2013C.5",
        notes="Protocol of record. Data rulings: the +5.3 degree encoder correction on the 47.9 cm flexor test (applied to experimental angles at load time), and the three extensor tests not used in the evaluation (3, 4, and 9). The stiffness pair is not unique - the adopted pair is the one for which both redesign configurations meet the required 5 percent overshoot simultaneously. Anticipated question: why gamultiobj? Answer: the Pareto front is required for holdout validation; settings of record - flexor: population 150, up to 600 generations, spread tolerance 4e-3; extensor: population 50, up to 150 generations, 1e-3.")
    items = [
        ("Encoder correction",
         "the 47.9 cm flexor test\u2019s measured angles carry a +5.3\u00b0 correction at load time (App. C, Table C-1); it was not used in the identification, and its corrected-model fit is a post hoc evaluation"),
        ("Training / validation",
         "flexor: four tests cross-validated, 2 training / 2 held-out per fold (adopted fold: 48.5 + 41.0 cm train, 45.7 + 40.6-tendon validate) \u00b7 extensor: training {43, 43-tendon, 46 cm}, validation {40, 40-tendon, 47 cm}; 42, 42-tendon, and 48 cm not used"),
        ("Non-uniqueness handled explicitly",
         "a family of stiffness pairs fits the training tests nearly equally well; the adopted pair (pooled-front row 77; extensor locked refit row 32) is the one for which both redesigns meet the required 5% overshoot simultaneously (+5.01% flexor, +5.02% extensor)"),
        ("Solver settings of record",
         "NSGA-II gamultiobj \u00b7 flexor: population 150, \u2264600 generations, spread tol 4\u00d710\u207b\u00b3 \u00b7 extensor: population 50, \u2264150 generations, 10\u207b\u00b3 \u00b7 offsets searched in cm, stiffnesses on log\u2081\u2080"),
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
        cite="Ch. 4 \u00a74.1 (Improved Torque Model); bracket frames and arrays: App. C \u00a7\u00a7C.3\u2013C.4",
        notes="Frame-convention invariance, proven in Ch. 4: a diagonal stiffness array is invariant under a rotation about y if and only if the two stiffnesses in the rotation plane are equal (k_x = k_z); with the y-symmetric array the pitch-only and two-rotation conventions give identical torque predictions to numerical precision. Anticipated question: where does this not apply? Answer: non-symmetric axis assignments - the extensor arch bracket sets chi-2 in X_br and chi-1 in Y_br/Z_br on physical grounds (the arch extends 75 mm in Z), a modeling choice carried, not derived.")
    card(s, MX, 2.1, 7.3, 4.3, fill=CARD)
    rich(s, MX + 0.3, 2.35, 6.7, 3.9, [
        {"t": "CLAIM (CH. 4, IMPROVED TORQUE MODEL)", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "Let K = diag(a, b, a); equal stiffness about x and z. For any frame whose second rotation is about y, the projected scalar stiffness", "size": 14, "color": TEXT, "spacing": 1.15, "space_after": 8},
        {"t": "k(\u00fb) = (\u00fb\u1d40 K\u207b\u00b9 \u00fb)\u207b\u00b9", "size": 17, "color": BROWN, "bold": True, "space_after": 8},
        {"t": "is unchanged: the rotation about y mixes only the two coordinates with equal compliance, and the two frame conventions produce identical torque predictions to numerical precision.", "size": 14, "color": TEXT, "spacing": 1.15, "space_after": 8},
        {"t": "Invariant if and only if the x- and z-axis stiffnesses are equal; the proof is the rotation identity Ry(\u03b8)\u00b7K\u00b7Ry(\u03b8)\u1d40 = K with Ry a rotation about y.", "size": 13, "color": GREEN, "bold": True, "spacing": 1.1},
    ])
    card(s, 8.15, 2.1, 4.6, 4.3, fill=CARD2)
    rich(s, 8.4, 2.35, 4.1, 3.9, [
        {"t": "SCOPE", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "Holds for y-symmetric arrays only.", "size": 13.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "The extensor\u2019s arch bracket is not y-symmetric: the bracket forms an arch extending 75 mm in the Z direction, so \u03c7\u2082 sits along the bracket x axis and \u03c7\u2081 along y and z; an assignment carried on physical grounds, not derived from this proof.", "size": 13.5, "color": TEXT, "spacing": 1.12},
    ])
    pageno(s, n)

def a8_bracket_point(prs, n):
    s = content_slide(prs, "Appendix H \u00b7 bracket data",
        "Bracket reference points and adopted correction terms; recorded together",
        cite="Reference points: App. C \u00a7C.3 (Table C-5); adopted values and provenance: App. C \u00a7C.5 (Table C-8)",
        notes="Bracket data of record. Moving a reference point changes the identified stiffness values, so Appendix C records each point alongside the values it produced. Provenance: flexor values are pooled-front row 77 of the two-bracket CV (212 members); the extensor refit locked the pair and solved only (chi-0, chi-3) - flxr77 row 32. The pre-freeze lab study that compared the rib midpoint against the lower bolt hole is not part of the dissertation and is not shown. Anticipated question: how sensitive is the identification to the point choice? Answer: within the dissertation the points are fixed CAD interpretation choices recorded with the values they produced; the explicit sensitivity comparison predates the freeze.")
    tx(s, MX, 2.0, 7.0, 0.4, "Bracket reference points (mm)", size=13.5, color=TEAL, bold=True)
    hdr = ["point", "x, y, z", "configuration / meaning"]
    rows = [
        ["Pbri", "\u221248.11, \u2212107.81, 13.80", "pinned flexor \u00b7 insertion bracket (ICR frame)"],
        ["Pbr2", "\u221252.61, 0, 75.06", "pinned flexor \u00b7 origin 2nd bracket (hip frame)"],
        ["Pbr", "\u22123.84, \u221246.44, 62.50", "pinned extensor \u00b7 rib midpoint (hip frame)"],
        ["Pbr", "8.38, 20.75, 25.10", "biomimetic extensor 10 mm \u00b7 cantilever"],
        ["Pbr", "9.48, \u221233.38, 33.90", "biomimetic flexor 20 mm \u00b7 bolt-pattern centroid"],
        ["Pbr", "\u221219.00, 22.00, 27.60", "biomimetic flexor 10 mm \u00b7 cantilever centroid"],
    ]
    cw = [0.8, 2.25, 3.6]
    y = 2.42
    for r_i, r in enumerate([hdr] + rows):
        xx = MX
        for i, val in enumerate(r):
            if r_i == 0:
                tx(s, xx, y, cw[i] + 0.12, 0.35, val, size=11.5, color=TEAL, bold=True,
                   align=PP_ALIGN.LEFT)
            else:
                tx(s, xx, y + 0.02, cw[i] + 0.12, 0.4, val, size=11 if i == 2 else 12,
                   color=TEXT if i < 2 else MUTED, bold=(i == 0), align=PP_ALIGN.LEFT)
            xx += cw[i] + 0.12
        y += 0.44
    tx(s, 7.95, 2.0, 4.85, 0.4, "Adopted correction terms", size=13.5, color=TEAL, bold=True)
    hdr2 = ["", "\u03c7\u2080 (mm)", "\u03c7\u2081 (N/m)", "\u03c7\u2082 (N/m)", "\u03c7\u2083"]
    rows2 = [["Flexor", "+3.94", "3.998e4", "1.473e4", "\u2014"],
             ["Extensor", "\u22126.38", "locked", "locked", "0.158"]]
    cw2 = [1.0, 0.9, 1.05, 1.05, 0.55]
    y = 2.42
    for r_i, r in enumerate([hdr2] + rows2):
        xx = 7.95
        for i, val in enumerate(r):
            st = {"size": 12, "color": TEAL if r_i == 0 else TEXT, "bold": r_i == 0 or i == 0,
                  "align": PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER}
            tx(s, xx, y, cw2[i], 0.4, val, **st)
            xx += cw2[i]
        y += 0.44
    rich(s, 7.95, 3.95, 4.85, 1.2, [
        {"t": "extensor refit locks (\u03c7\u2081, \u03c7\u2082) to the flexor pair and solves only (\u03c7\u2080, \u03c7\u2083); one stiffness pair across all configurations", "size": 11.5, "color": MUTED, "spacing": 1.1},
    ])
    card(s, MX, 5.55, 12.2, 1.15, fill=CARD2)
    rich(s, MX + 0.3, 5.72, 11.6, 0.95, [
        {"t": "Moving a reference point changes the identified stiffness; App. C records each point alongside the values it produced.", "size": 14, "color": TEXT, "spacing": 1.05, "space_after": 3},
        {"t": "Adoption provenance: flexor pooled-front row 77 (212 members) \u00b7 extensor locked refit row 32 (only \u03c7\u2080, \u03c7\u2083 searched).", "size": 12.5, "color": MUTED, "spacing": 1.05},
    ])
    pageno(s, n)

def a9_moment_arm(prs, n):
    s = content_slide(prs, "Appendix I \u00b7 kinematics",
        "Moment arm and torque conventions",
        cite="Ch. 4 Eqs. (4.1)\u2013(4.2); implementation: App. A \u00a7A.4 (MonoPamDataExplicit.m lines 156\u2013168, 308); human targets: App. A \u00a7A.1",
        notes="Moment-arm conventions. The robot side uses the perpendicular-foot construction of App. A sect. A.4; the human targets are computed in OpenSim with the derivative definition r(theta_k) = d(l_MT)/d(theta_k) across the range of motion (App. A sect. A.1). Anticipated question: are the two equivalent? Answer: yes for a straight muscle path crossing a planar revolute joint; the derivative form remains correct for wrapping and moving-path-point geometry, which is why the human targets are computed in OpenSim rather than by the robot-side construction.")
    card(s, MX, 1.95, 12.2, 1.55, fill=CARD)
    eq(s, "torque", 1.6, 2.42, h=0.6)
    card(s, MX, 3.72, 12.2, 2.75, fill=CARD2)
    rich(s, MX + 0.3, 3.95, 11.6, 2.4, [
        {"t": "CONVENTIONS", "size": 11.5, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "r = perpendicular foot of the joint origin onto the line of action; identical result to the cross-product projection of Eq. (4.2)", "size": 14, "color": TEXT, "spacing": 1.12, "space_after": 5},
        {"t": "scalar moment arm r\u2096 = hypot(in-plane components) in the tibial body frame, about the knee z-axis", "size": 14, "color": TEXT, "spacing": 1.12, "space_after": 5},
        {"t": "torque evaluated only in tension (entries discarded for \u03b5 < \u22120.02, the excessive-stretch region)", "size": 14, "color": TEXT, "spacing": 1.12, "space_after": 5},
        {"t": "human targets: OpenSim derivative arms r(\u03b8\u2096) = \u2202l_MT/\u2202\u03b8\u2096 across the RoM, multiplied by maximum isometric force", "size": 12.5, "color": MUTED, "spacing": 1.1},
    ])
    pageno(s, n)

def a10_sns_impl(prs, n):
    s = content_slide(prs, "Appendix J \u00b7 implementation",
        "Synthetic nervous system: Simulink library and spiking twin",
        cite="Ch. 6 \u00a76.2 (Figs. 6.5 and 6.11); spiking-twin numbers: Ch. 6 Table 6.4; neuron and synapse equations: App. B \u00a7\u00a7B.1\u2013B.4",
        notes="Two parallel implementations: AnimatLab (prototyping and the archived reference models) and Simulink/MuJoCo (the pipeline). The Simulink library runs in base Simulink without extra toolboxes: seven core blocks, twelve before the spiking additions, fifteen after (Ch. 6 sect. 6.2). The spiking twin of the knee-reflex demonstration replaces only the sensory interneuron layer: transient identical (15.79 deg at 0.05 s in both), settle band +0.9 deg (Table 6.4). CAD coupling verified via SolidWorks-URDF-Simscape on a trial knee assembly; tendon hardware modeled in CAD, not yet inserted. Anticipated question: why maintain Simulink alongside MuJoCo? Answer: a standard controls-engineering environment for inspection and the hardware-transfer path, with unit-verified spiking blocks - single-spike conductance transient to 6.2e-07 uS, 62/62 spikes against the toolbox equations, and every committed demonstration reproducing exactly (two-neuron units test 6.07e-04 mV).")
    pic(s, os.path.join(FIG, "d_simulink_library_grid.png"), box=(0.6, 1.80, 6.0, 3.45))
    tx(s, 0.6, 5.30, 6.0, 0.55,
       "seven core masked blocks (library grows to 12, then 15 with the spiking additions); runs in base Simulink",
       size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER, spacing=1.0)
    pic(s, os.path.join(FIG, "d_simulink_knee_spiking_twin.png"), box=(6.9, 1.80, 5.85, 3.45))
    tx(s, 6.9, 5.30, 5.85, 0.55,
       "knee-reflex spiking twin at the matched operating point: transient identical (15.79\u00b0 at 0.05 s), settle +0.9\u00b0",
       size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER, spacing=1.0)
    card(s, 0.6, 5.90, 6.0, 1.0, fill=CARD)
    rich(s, 0.85, 6.04, 5.5, 0.8, [
        {"t": "CAD COUPLING (CH. 6 \u00a76.2)", "size": 11.5, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "SolidWorks \u2192 URDF \u2192 Simscape Multibody import verified on a trial knee assembly; tendon hardware modeled in CAD, not yet inserted.", "size": 12, "color": TEXT, "spacing": 1.05},
    ])
    card(s, 6.9, 5.90, 5.85, 1.0, fill=CARD)
    rich(s, 7.15, 6.04, 5.35, 0.8, [
        {"t": "SPIKING BLOCK VERIFICATION (CH. 6 \u00a76.4)", "size": 11.5, "color": TEAL, "bold": True, "space_after": 3},
        {"t": "single-spike conductance to 6.2\u00d710\u207b\u2077 \u00b5S \u00b7 62/62 spikes vs the toolbox equations \u00b7 committed regression exact (6.07\u00d710\u207b\u2074 mV).", "size": 12, "color": TEXT, "spacing": 1.05},
    ])
    pageno(s, n)

def a11_walker_numbers(prs, n):
    s = content_slide(prs, "Appendix K \u00b7 walker numbers",
        "Walker configurations of record: scores, library sweep, and pathway deletions",
        cite="Ch. 6 Table 6.1 (gait sweep) and Table 6.2 (pathway removal); drive perturbations: Ch. 6 \u00a76.3",
        notes="Walker numbers of record for committee questions. Score is the cycle-shape kinematic score: 0 is perfect agreement, -25 the fixed no-cycle fallback. Library means sit one to two orders of magnitude beyond the tuning-reference scores - an honest transfer measurement, not an artifact of the 32 running trials (s3k walking-only mean -8816 vs running-only -8761). Anticipated question: which deletion matters most? Answer: removing the Ia pathway freezes the production configuration outright (no cycles detected), while the six-synergy variant keeps a bilateral gait with every afferent class removed (-190.5, ten cycles per leg, knee minimum -71.2 deg, 0.684 s period); the W2L-layout variant depends on Ia (-273.5) and on the weak rhythm-generator mutual excitation (-279.9).")
    tx(s, MX, 2.0, 5.9, 0.4, "Gait-validation sweep (43-reference library)", size=13, color=TEAL, bold=True)
    hdr = ["config", "tuning-ref", "library mean"]
    rows = [
        ["s3k", "\u2212186.4", "\u22128776"],
        ["s3kpruned", "\u2212176.1", "\u22126924"],
        ["syn6", "\u2212200.4", "\u22129961"],
        ["w2lvar", "\u2212144.5", "\u22125808"],
    ]
    cw = [1.6, 1.5, 1.7]
    y = 2.42
    for r_i, r in enumerate([hdr] + rows):
        xx = MX
        for i, val in enumerate(r):
            st = {"size": 12.5, "color": TEAL if r_i == 0 else TEXT, "bold": r_i == 0 or i == 0,
                  "align": PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER}
            tx(s, xx, y, cw[i], 0.4, val, **st)
            xx += cw[i]
        y += 0.46
    rich(s, MX, 4.75, 5.9, 1.0, [
        {"t": "most negative reference: Case_40_motion for s3k (\u22128841), s3kpruned (\u22126986), and syn6 (\u221210040); subject08_Run_20002 for w2lvar (\u22125896)", "size": 11.5, "color": MUTED, "spacing": 1.1},
    ])
    tx(s, 6.85, 1.92, 5.9, 0.4, "Pathway deletions, s3k (tuning-reference score)", size=13, color=TEAL, bold=True)
    hdr2 = ["removed pathway", "score"]
    rows2 = [
        ["none (full model)", "\u2212186.4"],
        ["Ia afferents", "\u221225 (no cycles)"],
        ["all afferent classes", "\u2212267.8"],
        ["Renshaw cells", "\u2212267.3"],
        ["group-II afferents", "\u2212255.2"],
        ["Ib afferents", "\u2212180.3"],
        ["contact channels", "\u2212183.9"],
        ["interleg commissurals", "\u2212174.9"],
    ]
    cw2 = [3.3, 1.7]
    y = 2.35
    for r_i, r in enumerate([hdr2] + rows2):
        xx = 6.85
        for i, val in enumerate(r):
            st = {"size": 12 if r_i == 0 else 11.5, "color": TEAL if r_i == 0 else TEXT,
                  "bold": r_i == 0 or i == 0,
                  "align": PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER}
            tx(s, xx, y, cw2[i] + (0.3 if i == 0 else 0), 0.38, val, **st)
            xx += cw2[i] + (0.3 if i == 0 else 0)
        y += 0.375
    card(s, MX, 5.85, 12.2, 1.05, fill=CARD2)
    rich(s, MX + 0.3, 5.98, 11.6, 0.9, [
        {"t": "Drive stimulus: s3k survives +0.01, \u22120.01, \u22120.02 (\u2212201.9, \u2212177.0, \u2212171.2), freezes at +0.02, and never fell in 15 push evaluations (max sway 0.123 m).", "size": 12.5, "color": TEXT, "spacing": 1.05, "space_after": 3},
        {"t": "syn6 and w2lvar fall under every \u00b10.01 / \u00b10.02 drive step (pelvis tilt 39\u201340\u00b0); s3kpruned survives \u22120.02 (\u2212168.2).", "size": 12, "color": MUTED, "spacing": 1.05},
    ])
    pageno(s, n)

def build(prs, start):
    a1_force_model(prs, start)
    a2_coefficients(prs, start + 1)
    a3_model_compare(prs, start + 2)
    a4_compliance(prs, start + 3)
    a5_wraploss(prs, start + 4)
    a6_protocol(prs, start + 5)
    a7_frame_invariance(prs, start + 6)
    a8_bracket_point(prs, start + 7)
    a9_moment_arm(prs, start + 8)
    a10_sns_impl(prs, start + 9)
    a11_walker_numbers(prs, start + 10)
    return start + 10
