"""Slides 1-14: title, framing (four objectives / four gaps), Ch. 3 muscles (BPA force)."""
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt
from deck_lib import *

LOGO = os.path.join(ASSETS, "aarl_logo_0.png")

def s01_title(prs, n):
    s = prs.slides.add_slide(prs.slide_layouts[6]); bg(s)
    pic(s, LOGO, box=(4.47, 0.75, 4.4, 1.34))
    tx(s, 0.7, 2.55, 11.93, 1.9, "Movement and Control of\nBiomimetic Humanoid Robots",
       size=44, color=BROWN, bold=True, align=PP_ALIGN.CENTER, spacing=1.0)
    tx(s, 0.7, 4.45, 11.93, 0.45, "Dissertation Defense", size=20, color=MUTED, align=PP_ALIGN.CENTER)
    tx(s, 0.7, 5.15, 11.93, 0.5, "Ben P. Bolen", size=24, color=TEXT, bold=True, align=PP_ALIGN.CENTER)
    tx(s, 0.7, 5.72, 11.93, 0.4, "Advisor: M. Anthony Hunt, Ph.D.", size=15, color=MUTED, align=PP_ALIGN.CENTER)
    tx(s, 0.7, 6.30, 11.93, 0.4, "Agile & Adaptive Robotics Laboratory  ·  Portland State University  ·  November 2026",
       size=13, color=MUTED, align=PP_ALIGN.CENTER)
    s.notes_slide.notes_text_frame.text = (
        "Good morning. Committee members, thank you for reading the dissertation. "
        "Today: four objectives - characterize the actuator, validate and correct the joint "
        "torque prediction, consolidate the afferent literature, and build the verified control "
        "toolchain - connecting one actuator to a walking, neurally controlled robot. "
        "Anticipated question: why does the title say 'movement and control' rather than a "
        "single topic? Answer: the dissertation spans actuator to controller, and the title "
        "follows the program it serves (Ch. 1, sec. 1.2).")

def s02_roadmap(prs, n):
    s = content_slide(prs, "Roadmap",
        "Four objectives connect one actuator to a robot that walks under neural control",
        cite="Dissertation Ch. 1 \u00a71.3 (research objectives); outcomes in Ch. 3-7",
        notes="One sentence per objective, mirroring Ch. 1 sec 1.3: objective 1 (Ch. 3) is published as "
              "Bolen et al. 2026, Actuators; objective 2 (Ch. 4) stands on the isometric torque manuscript; "
              "objective 3 (Ch. 5) is the Sensory Afferent Database; objective 4 (Ch. 6-7) is the simulation "
              "toolchain with its cross-engine verification. Anticipated question: which objectives are "
              "experimental? Answer: objectives 1-2 are hardware experiments on the jig and the knee stands; "
              "objective 3 is literature curation; objective 4 is simulation tooling - the same split the "
              "introduction states (first two gaps addressed experimentally, the fourth resolved by the "
              "database, the third closed by the tools).")
    items = [
        ("01", "Characterize\nBPA isometric force",
         "Festo 10 and 20 mm BPAs across resting length,\ncontraction, and pressure; improved maximum-force\nand normalized force models vs. existing models\n(Ch. 3)",
         "537 tests  \u00b7  dimensionless\nforce surface", GREEN),
        ("02", "Validate knee\ntorque prediction",
         "Pinned and biomimetic four-bar knees; correction\nterms for length offset, series compliance, and\nwrapping; envelopes vs. human torque\n(Ch. 4)",
         "redesign meets / exceeds\nhuman targets", GREEN),
        ("03", "Consolidate the\nafferent literature",
         "Audit and curate the Ia, Ib, and group II afferent\nrecord into the Sensory Afferent Database;\nAI-assisted triage for expansion\n(Ch. 5)",
         "943-paper audited corpus", BLUE),
        ("04", "Build tools for\nneural control",
         "OpenSim to MuJoCo conversion; SNS-Toolbox two-level\nspinal CPG with proprioceptive and heel/toe contact\nfeedback; staged tuning and verification\n(Ch. 6-7)",
         "verified cross-engine\ntoolchain", TEAL),
    ]
    x = MX; w = 2.92; gap = 0.183
    for i, (num, q, m, tag, c) in enumerate(items):
        cx = x + i * (w + gap)
        card(s, cx, 2.0, w, 4.55, fill=CARD if i < 3 else CARD2)
        tx(s, cx + 0.25, 2.18, 1.6, 0.7, num, size=30, color=BROWN, bold=True)
        tx(s, cx + 0.25, 2.92, w - 0.5, 1.0, q, size=15, color=TEXT, bold=True, spacing=1.0)
        tx(s, cx + 0.25, 4.08, w - 0.5, 1.6, m, size=11.5, color=MUTED, spacing=1.05)
        tx(s, cx + 0.25, 5.78, w - 0.5, 0.72, tag, size=12, color=c, bold=True, spacing=1.0)
        if i < 3:
            arrow(s, cx + w + 0.015, 4.25, cx + w + gap - 0.015, 4.25, color=FAINT, width=2.2)
    pageno(s, n)

def s03_project(prs, n):
    s = content_slide(prs, "The project",
        "A humanoid robot from the trunk down: human-like body plan, synthetic nervous system",
        cite="Dissertation Ch. 1 \u00a71.2; Bolen & Hunt 2019; Morrow et al. 2020; Steele 2018; Scharzenberger et al. 2019",
        notes="Ground the committee: this dissertation is the actuator-and-joint foundation of the AARL humanoid "
              "program (Ch. 1, sec. 1.2). The converted model shown carries all 92 muscle-tendon actuators "
              "(Ch. 6). Anticipated question: what did the master's work establish? Answer: Gait2392-matched "
              "kinematics and muscle-routing tools - and the finding that transplanting human attachments onto "
              "a BPA robot does not reproduce human torque profiles, which motivates the placement optimization "
              "completed in Ch. 4.")
    pic(s, os.path.join(FIG, "mujoco_robot.png"), box=(0.45, 1.75, 3.4, 5.1))
    tx(s, 0.45, 6.55, 3.5, 0.4, "converted leg model, 92 muscle paths", size=11, color=FAINT, italic=True)
    rows = [
        ("Guiding hypothesis", "human-like muscle numbers, architecture, routing + human-like torque capability, controlled by a biologically grounded neural architecture  \u2192  biomimetic movement", None),
        ("Built on", "my master's work: Gait2392-matched kinematics, muscle routing tools  \u00b7  Morrow: optimized attachments  \u00b7  Steele: biomimetic 4-bar knee  \u00b7  Scharzenberger: BPA quadruped", None),
        ("This dissertation", "the quantitative basis: characterize the actuator, correct the joint model against experiment, compare with human torque, consolidate the afferent data, and define the control path", GREEN),
    ]
    y = 1.95
    for t, b, c in rows:
        rich(s, 4.35, y, 5.45, 1.5, [
            {"t": t.upper(), "size": 12, "color": TEAL if not c else c, "bold": True, "space_after": 3},
            {"t": b, "size": 14.5, "color": TEXT, "spacing": 1.08},
        ])
        y += 1.62
    pic(s, os.path.join(MEDIA, "image44.jpg"), box=(10.0, 5.45, 2.8, 1.5))
    tx(s, 4.35, 6.55, 5.45, 0.4, "AARL quadruped: BPA actuators + SNS, already on a treadmill", size=11, color=FAINT, italic=True)
    pageno(s, n)

def s04_motivation(prs, n):
    s = content_slide(prs, "Why biomimetic robots",
        "Robots run neuromechanical experiments that are impossible in living subjects",
        cite="Dissertation Ch. 1 \u00a71.1; Ijspeert 2020; Shin et al. 2018; photo: Asano et al. 2016, Humanoids",
        notes="The biology professor on the committee will recognize this loop: model, perturb, understand. "
              "The experiments listed are the ones the introduction names (Ch. 1, sec. 1.1). "
              "Anticipated question: why physical robots rather than simulation-only studies? Answer: "
              "musculoskeletal humanoids and limb systems support experiments on internal force transmission "
              "and hardware realities that simulation assumes away; the two are complementary.")
    node(s, MX, 2.15, 3.5, 1.15, "Biology\nanatomy \u00b7 architecture \u00b7 afferents", fill=CARD2, size=14, bold=False)
    node(s, MX, 4.9, 3.5, 1.15, "Robot\nactuators \u00b7 routing \u00b7 controller", fill=CARD, size=14)
    arrow(s, MX + 1.0, 3.4, MX + 1.0, 4.8, color=GREEN, width=2.5)
    tx(s, MX + 1.15, 3.85, 2.6, 0.5, "design principles", size=12, color=GREEN, bold=True)
    arrow(s, MX + 2.6, 4.8, MX + 2.6, 3.4, color=BROWN, width=2.5)
    tx(s, MX + 2.72, 3.85, 2.6, 0.5, "test hypotheses", size=12, color=BROWN, bold=True)
    pic(s, os.path.join(MEDIA, "image11.png"), box=(8.9, 1.9, 3.9, 4.8))
    tx(s, 8.9, 6.72, 3.9, 0.35, "human-mimetic musculoskeletal humanoid", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    exps = [
        "Remove one muscle from the system",
        "Block one sensory feedback pathway",
        "Stimulate a named neural pathway",
        "Identical copies - no individual variance",
        "Faster, cheaper, ethical vs. human studies",
    ]
    y = 2.1
    for e in exps:
        tx(s, 4.55, y, 4.2, 0.5, "\u2013  " + e, size=16, color=TEXT, spacing=1.0)
        y += 0.62
    tx(s, 4.55, 5.55, 4.2, 1.2, "Loss-of-function and stimulation experiments\nare the point of the platform.", size=15, color=BROWN, bold=True, spacing=1.05)
    pageno(s, n)

def s05_bpa(prs, n):
    s = content_slide(prs, "The actuator",
        "Braided pneumatic actuators: muscle-like strengths, and three modeling headaches",
        cite="Dissertation Ch. 1 \u00a71.1 + Ch. 2; Hunt et al. 2017; S\u00e1rosi et al. 2012-2017; Martens & Boblan 2017; Festo DMSP datasheets",
        notes="Sets up gap 1. Force peaks AT resting length - unlike a biological muscle whose optimum is "
              "mid-range (Ch. 3, biomimetic-analogy discussion). Anticipated question: why BPAs rather than "
              "electric actuators? Answer: low mass, high force-to-weight, intrinsic compliance, tension-only "
              "operation, and a qualitatively muscle-like force-length curve (Ch. 1, sec. 1.1) - the costs are "
              "the modeling issues on the right, which Ch. 3 quantifies. Off-slide detail: contraction "
              "is capped near epsilon_620, and the design chain these models serve is actuator force -> "
              "joint torque -> comparison with human capability.")
    pic(s, os.path.join(MEDIA, "image5.png"), box=(0.5, 1.95, 5.15, 4.55))
    tx(s, 0.5, 6.58, 5.15, 0.35, "Festo DMSP-10/20 fluidic muscle", size=12, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    plus = [
        "High force-to-weight, low mass, intrinsic compliance",
        "Tension-only, like biological muscle",
        "Muscle-like force-length curve",
    ]
    minus = [
        "Maximum force occurs AT resting length",
        "Strongly nonlinear in l\u2080, P, \u03b5",
        "Manufacturer tool mis-predicts short BPAs",
    ]
    rich(s, 6.1, 1.95, 6.7, 0.4, [{"t": "WHY WE USE THEM", "size": 14, "color": GREEN, "bold": True}])
    y = 2.42
    for p_ in plus:
        tx(s, 6.1, y, 6.7, 0.5, "+  " + p_, size=17, color=TEXT, spacing=1.0); y += 0.62
    rich(s, 6.1, 4.55, 6.7, 0.4, [{"t": "WHY THEY ARE HARD", "size": 14, "color": RED, "bold": True}])
    y = 5.02
    for m_ in minus:
        tx(s, 6.1, y, 6.7, 0.5, "\u2013  " + m_, size=17, color=TEXT, spacing=1.0); y += 0.62
    pageno(s, n)

def s06_gaps(prs, n):
    s = content_slide(prs, "Problem statement",
        "Four gaps stood between the existing design and a robot that walks under neural control",
        cite="Dissertation Ch. 1 \u00a71.1 (summary paragraph) + \u00a71.3",
        notes="The four gaps are the introduction's own summary (Ch. 1, sec. 1.1, closing paragraph). "
              "Anticipated question: why is gap 4 in the literature rather than in methods? Answer: the "
              "physiological properties of the Ia, Ib, and group II afferents that should inform feedback "
              "design are scattered across hundreds of primary papers and reviews with no structured "
              "resource - the dissertation resolves it by consolidation (the SADb) and closes gap 3 with "
              "the tools, while gaps 1-2 are addressed experimentally.")
    rows = [
        ("GAP 1", "No actuator force model valid at arbitrary resting length",
         "Festo tool over-predicts short BPAs; published models fit a single length each",
         GREEN),
        ("GAP 2", "No quantified correction for torque lost through brackets, tendons, and wrapping",
         "optimized biomimetic routing never validated against human isometric torque",
         GREEN),
        ("GAP 3", "No closed-loop pipeline from the robot model to a biologically realistic controller",
         "two-level CPG ingredients existed, never assembled for a BPA biped with human-like routing",
         BLUE),
        ("GAP 4", "No consolidated resource for the afferent data a controller requires",
         "Ia, Ib, and group II properties scattered across hundreds of primary papers and reviews",
         TEAL),
    ]
    y = 1.82
    for tag, head, sub, c in rows:
        card(s, MX, y, 12.2, 1.06, fill=CARD)
        tx(s, MX + 0.28, y + 0.15, 1.2, 0.5, tag, size=14, color=c, bold=True)
        rich(s, MX + 1.55, y + 0.10, 10.4, 0.9, [
            {"t": head, "size": 15, "color": TEXT, "bold": True, "space_after": 1},
            {"t": sub, "size": 12.5, "color": MUTED, "spacing": 1.0},
        ])
        y += 1.17
    tx(s, MX, 6.52, 12.2, 0.4,
       "Objectives 1-2 address gaps 1-2 experimentally  \u00b7  the SADb resolves gap 4  \u00b7  the tools close gap 3",
       size=14, color=BROWN, bold=True)
    pageno(s, n)

def s07_jig(prs, n):
    s = content_slide(prs, "Ch. 3 \u00b7 muscles",
        "A custom jig measured isometric force across resting length, pressure, and contraction",
        cite="Dissertation Ch. 3 \u00a73.1 (test jig); Bolen, Elzein, Pang, and Hunt 2026, Actuators",
        notes="Hardware slide - one breath. The point: wide l0 range, both diameters, 537 pairs. "
              "Anticipated question: how is contraction set? Answer: the distance between the jig "
              "crossmembers is controlled to obtain each contraction, then the muscle is inflated and "
              "the contractile force recorded at the pressure-contraction pairs (Ch. 3, sec. 3.1).")
    pic(s, os.path.join(FIG, "testJigs_force.png"), box=(0.7, 1.9, 4.4, 4.75))
    specs = [
        ("Diameters", "10 mm and 20 mm Festo DMSP"),
        ("Resting lengths", "112-518 mm (10 mm)  \u00b7  300-509 mm (20 mm)"),
        ("Pressure", "up to P\u2086\u2082\u2080 = 620 kPa supply"),
        ("Data", "537 force-pressure-contraction pairs (321 + 216)"),
        ("Instrumentation", "S-type load cells + HX711 \u00b7 MPX5700 pressure \u00b7 Arduino \u2192 MATLAB"),
    ]
    y = 2.1
    for k, v in specs:
        rich(s, 5.6, y, 7.1, 0.9, [
            {"t": k.upper(), "size": 11.5, "color": TEAL, "bold": True, "space_after": 1},
            {"t": v, "size": 15.5, "color": TEXT, "spacing": 1.0},
        ])
        y += 0.95
    pageno(s, n)

def s08_fmax_discovery(prs, n):
    s = content_slide(prs, "Ch. 3 \u00b7 muscles",
        "Maximum isometric force depends on resting length - a previously unreported characteristic",
        cite="Dissertation Ch. 3 \u00a73.2 (F\u2086\u2082\u2080 vs l\u2080); Bolen et al. 2026, Actuators, Figs. 2-3",
        notes="The discovery slide. Festo's own tool predicts the OPPOSITE trend - force rising as length shrinks. "
              "Anticipated question: is the length effect truly unreported? Answer: the chapter states it is a "
              "previously unreported characteristic of these artificial muscles; the Festo tool does predict a "
              "length dependence, but with the wrong sign (force increasing at shorter lengths).")
    pic(s, os.path.join(FIG, "MaxForce10.png"), box=(0.7, 1.9, 5.9, 3.95))
    pic(s, os.path.join(FIG, "MaxForce20.png"), box=(6.75, 1.9, 5.9, 3.95))
    tx(s, 0.7, 5.92, 5.9, 0.4, "\u03c6 10 mm:  F\u2086\u2082\u2080 vs l\u2080", size=13, color=MUTED, align=PP_ALIGN.CENTER)
    tx(s, 6.75, 5.92, 5.9, 0.4, "\u03c6 20 mm:  F\u2086\u2082\u2080 vs l\u2080", size=13, color=MUTED, align=PP_ALIGN.CENTER)
    card(s, 0.7, 6.38, 12.0, 0.52, fill=CARD2)
    tx(s, 1.0, 6.48, 11.5, 0.4,
       "solid: our arctan fit (data-driven)  \u00b7  dashed: Festo tool - predicts force INCREASING as l\u2080 shrinks",
       size=13.5, color=TEAL, bold=True, align=PP_ALIGN.CENTER)
    pageno(s, n)

def s09_fmax_model(prs, n):
    s = content_slide(prs, "Ch. 3 \u00b7 muscles",
        "Two coefficients per diameter capture the resting-length effect",
        cite="Dissertation Ch. 3 \u00a7\u00a73.1-3.2 (equations + coefficient table); Bolen et al. 2026, Actuators, Eqs. 4-6; fit: nonlinear least squares, least-absolute-residual",
        notes="Math committee: arctan form, bounded asymptote, physical offset from end-cap contact in CAD. "
              "Anticipated question: why an arctan? Answer: the measured F620 vs l0 saturates, and the arctan "
              "is the bounded curve that matches it; the l0 offsets (7.5 mm and 13 mm) are the lengths at "
              "which solid modeling shows the end caps contacting, where no contractile stroke remains.")
    card(s, MX, 1.95, 7.0, 1.5, fill=CARD)
    eq(s, "fmax", box=(1.0, 2.12, 6.3, 1.16))
    card(s, MX, 3.6, 7.0, 1.5, fill=CARD)
    eq(s, "fmax20", box=(1.0, 3.77, 6.3, 1.16))
    card(s, MX, 5.25, 7.0, 1.45, fill=CARD2)
    eq(s, "fmax3d", box=(1.0, 5.45, 6.3, 1.05))
    rich(s, 8.0, 2.0, 4.75, 4.8, [
        {"t": "LIMITS BEHAVE", "size": 12, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "l\u2080 \u2192 \u221e:   470 N (10 mm) vs Festo 490 N", "size": 15.5, "color": TEXT, "space_after": 1},
        {"t": "within Festo's stated 10% manufacturing tolerance (20 mm: 1460 vs 1570 N)", "size": 12.5, "color": MUTED, "space_after": 10},
        {"t": "THE OFFSET IS PHYSICAL", "size": 12, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "7.5 mm / 13 mm = length at which the end caps contact each other in the CAD model - no contractile stroke remains", "size": 14.5, "color": TEXT, "spacing": 1.1, "space_after": 10},
        {"t": "As l\u2080 \u2192 0, measured force \u2192 0.", "size": 15.5, "color": BROWN, "bold": True},
    ])
    pageno(s, n)

def s10_festo_wrong(prs, n):
    s = content_slide(prs, "Ch. 3 \u00b7 muscles",
        "At short resting lengths the manufacturer's tool mis-predicts force by 42%",
        cite="Dissertation Ch. 3 \u00a73.3 (discussion: effects of resting length); example: \u03c610 mm, l\u2080 = 112 mm, 620 kPa",
        notes="The design-relevant failure: anyone specifying BPAs under 300 mm gets burned by the tool. "
              "Anticipated question: is 42% within Festo's tolerance? Answer: no - Festo states roughly 10% "
              "manufacturing variability; 42.1% is roughly four times that, at a length a knee design "
              "actually uses, and the error direction (over-prediction) is the unsafe one for design.")
    bar_chart(s, 0.65, 2.0, 6.3, 4.6,
              ["Measured", "Our model", "Festo tool"], [350.9, 335.3, 498.6],
              [TEAL, GREEN, RED], title="\u03c610 mm at l\u2080 = 112 mm, 620 kPa:  isometric force (N)",
              max_scale=600)
    tx(s, 0.65, 6.65, 6.3, 0.4, "error:  4.5% (ours)  vs  42.1% (Festo)", size=15, color=TEXT,
       bold=True, align=PP_ALIGN.CENTER)
    card(s, 7.35, 2.1, 5.35, 4.3, fill=CARD)
    rich(s, 7.65, 2.35, 4.75, 3.9, [
        {"t": "WHY THE TOOL FAILS SHORT", "size": 12, "color": TEAL, "bold": True, "space_after": 6},
        {"t": "Ideal cylindrical models ignore end effects; the clamped ends force steep stress gradients.", "size": 14.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "The end-effect zone has roughly fixed length - so as l\u2080 shrinks, it dominates.", "size": 14.5, "color": TEXT, "spacing": 1.12, "space_after": 8},
        {"t": "Measured force \u2192 0 as l\u2080 \u2192 0.  The Festo tool instead predicts exponential growth.", "size": 14.5, "color": RED, "spacing": 1.12, "space_after": 8},
        {"t": "\u201cResearchers using BPA resting lengths under 300 mm should take note.\u201d", "size": 14, "color": BROWN, "bold": True, "italic": True, "spacing": 1.1},
    ])
    pageno(s, n)

def s11_surface(prs, n):
    s = content_slide(prs, "Ch. 3 \u00b7 muscles",
        "All 537 measurements collapse onto one dimensionless force surface",
        cite="Dissertation Ch. 3 \u00a73.2 (normalized surface); \u03b5* = \u03b5/\u03b5\u2086\u2082\u2080,  P* = P/620 kPa",
        notes="Normalization by each muscle's own F620 and eps620 is what makes one surface fit all lengths. "
              "Anticipated question: does eps620 need to be measured per actuator? Answer: yes - the chapter "
              "found no conclusive relationship between eps620 and l0 (linear fit adj. R2 = 0.4124), so each "
              "muscle's maximum contraction is recorded for use on the robot.")
    card(s, MX, 1.9, 12.2, 1.05, fill=CARD2)
    eq(s, "fstar", 1.0, 2.12, h=0.62)
    pic(s, os.path.join(FIG, "FStar10.png"), box=(0.75, 3.15, 5.5, 3.55))
    pic(s, os.path.join(FIG, "FStar20.png"), box=(6.6, 3.15, 5.5, 3.55))
    tx(s, 0.75, 6.75, 5.5, 0.35, "F*(\u03b5*, P*) - \u03c610 mm, pressure isoclines", size=12, color=MUTED, align=PP_ALIGN.CENTER)
    tx(s, 6.6, 6.75, 5.5, 0.35, "F*(\u03b5*, P*) - \u03c620 mm", size=12, color=MUTED, align=PP_ALIGN.CENTER)
    pageno(s, n)

def s12_surface_gof(prs, n):
    s = content_slide(prs, "Ch. 3 \u00b7 muscles",
        "Three coefficients per diameter: force predicted out-of-sample to a few percent",
        cite="Dissertation Ch. 3 \u00a73.2 (F* coefficient + goodness-of-fit table); 80/20 train/validation split",
        notes="Stress: validation is 20% held out, and max error under ~11% everywhere - at ANY length. "
              "Anticipated question: what do the max errors mean? Answer: on the held-out set the 10 mm "
              "surface reaches a 10.6% maximum residual (RMSE 0.0245 normalized) and the 20 mm reaches 5.7% "
              "(RMSE 0.0231); adjusted R2 stays at 0.9994 and 0.9943 respectively.")
    cols = [
        ("\u03c6 10 mm", "0.9998", "0.9994", "10.3%", "10.6%"),
        ("\u03c6 20 mm", "0.992", "0.9943", "6.8%", "5.7%"),
    ]
    hdr = ["", "adj. R\u00b2 (fit)", "adj. R\u00b2 (validation)", "max err (fit)", "max err (val.)"]
    x = MX; cw = [1.8, 2.6, 2.9, 2.4, 2.5]
    xx = x
    for i, h_ in enumerate(hdr):
        tx(s, xx + 0.1, 2.05, cw[i] - 0.2, 0.6, h_, size=13.5, color=TEAL, bold=True,
           align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
        xx += cw[i]
    y = 2.75
    for name, f, v, me, mv in cols:
        xx = x
        vals = [name, f, v, me, mv]
        for i, val in enumerate(vals):
            bold = i == 0 or i == 2
            colr = TEXT if i else BROWN
            tx(s, xx + 0.1, y, cw[i] - 0.2, 0.5, val, size=17 if i == 0 else 16, color=colr,
               bold=bold, align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
            xx += cw[i]
        y += 0.62
    tx(s, MX, 4.15, 12.2, 0.4, "RMSE (normalized):  0.005 (10 mm)  \u00b7  0.023 (20 mm)", size=14, color=MUTED)
    card(s, MX, 4.8, 12.2, 1.7, fill=CARD)
    rich(s, MX + 0.3, 5.0, 11.6, 1.4, [
        {"t": "The whole actuator model a designer needs:", "size": 15, "color": MUTED, "space_after": 4},
        {"t": "F*  (3 coefficients, per diameter)   \u00d7   F\u2086\u2082\u2080(l\u2080)  (2 coefficients, per diameter)", "size": 19, "color": BROWN, "bold": True},
    ])
    pageno(s, n)

def s13_vs_models(prs, n):
    s = content_slide(prs, "Ch. 3 \u00b7 muscles",
        "Lowest error of the compared models at every diameter and length tested",
        cite="Dissertation Ch. 3 \u00a73.3 (comparison with published BPA models); Bolen et al. 2026, Actuators, Fig. 5 + Table 3",
        notes="Sarosi params come from specific lengths - extrapolation is where they fail. Martens is geometry-sensitive. "
              "Anticipated question: why does Martens and Boblan blow up at 450 mm? Answer: geometric sensitivity - "
              "the chapter shows its F_PE dD/DL term grows dominant off its fitted geometry, to the point of "
              "predicting negative force at 20 mm, 450 mm.")
    pic(s, os.path.join(FIG, "Comparison.png"), box=(0.55, 1.9, 5.6, 4.95))
    tx(s, 0.55, 6.85, 5.6, 0.35, "measured vs predicted; y = x is perfect fit", size=11.5, color=FAINT,
       italic=True, align=PP_ALIGN.CENTER)
    rows = [
        ("10 mm, 257 mm", "22.9", "143.8", "73.8"),
        ("10 mm, 233 mm", "16.8", "132.7", "62.6"),
        ("20 mm, 300 mm", "35.6", "132.1", "72.8"),
        ("20 mm, 450 mm", "67.6", "177.7", "862.1"),
    ]
    tx(s, 6.6, 1.95, 6.2, 0.4, "RMSE (N) at matched conditions", size=13.5, color=TEAL, bold=True)
    cw = [2.1, 1.35, 1.35, 1.7]
    hdr = ["l\u2080", "Bolen", "S\u00e1rosi", "M&B"]
    xx = 6.6
    for i, h_ in enumerate(hdr):
        tx(s, xx, 2.4, cw[i], 0.4, h_, size=13.5, color=MUTED, bold=True,
           align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER)
        xx += cw[i]
    y = 2.85
    for r in rows:
        xx = 6.6
        for i, val in enumerate(r):
            if i == 1:
                tx(s, xx, y, cw[i], 0.45, val, size=16, color=GREEN, bold=True, align=PP_ALIGN.CENTER)
            elif i == 0:
                tx(s, xx, y + 0.03, cw[i], 0.45, val, size=14, color=TEXT)
            else:
                tx(s, xx, y, cw[i], 0.45, val, size=15, color=MUTED, align=PP_ALIGN.CENTER)
            xx += cw[i]
        y += 0.58
    card(s, 6.6, 5.35, 6.2, 1.35, fill=CARD2)
    rich(s, 6.9, 5.52, 5.6, 1.1, [
        {"t": "Normalization travels: one fit per diameter works across lengths,", "size": 14, "color": TEXT, "spacing": 1.1},
        {"t": "while fixed-parameter models fail when moved off their fit length.", "size": 14, "color": TEXT, "spacing": 1.1},
    ])
    pageno(s, n)

def s14_realtime(prs, n):
    s = content_slide(prs, "Ch. 3 \u2192 Ch. 4",
        "The force model is two cheap maps - built to sit inside a real-time controller",
        cite="Dissertation Ch. 3 \u00a73.3 (implications for real-time control and simulation); dynamics complement: Elzein et al. 2025",
        notes="Static map = feedforward actuator model inside the controller (lookup table). Valve-airflow "
              "dynamics are the Elzein layer. Anticipated question: can a static isometric map really support "
              "control? Answer: the chapter cites Martens and Boblan - BPA dynamics are dominated by valve air "
              "dynamics and mass flow, so the static force characteristic is the crucial piece for torque "
              "control accuracy; the model is computationally inexpensive online and suits lookup tables.")
    card(s, 2.1, 2.1, 4.1, 1.5, fill=CARD2)
    tx(s, 2.3, 2.28, 3.7, 0.4, "MAP 1", size=12, color=TEAL, bold=True)
    tx(s, 2.3, 2.68, 3.7, 0.8, "F\u2086\u2082\u2080(l\u2080)\nresting length \u2192 max force", size=15.5, color=TEXT, spacing=1.05)
    card(s, 7.1, 2.1, 4.1, 1.5, fill=CARD2)
    tx(s, 7.3, 2.28, 3.7, 0.4, "MAP 2", size=12, color=TEAL, bold=True)
    tx(s, 7.3, 2.68, 3.7, 0.8, "F*(\u03b5*, P*)\nstrain, pressure \u2192 fraction", size=15.5, color=TEXT, spacing=1.05)
    arrow(s, 6.28, 2.85, 7.02, 2.85, color=FAINT)
    eq(s, "force", 3.6, 4.05, h=0.62)
    card(s, 2.1, 5.15, 9.1, 1.35, fill=CARD)
    rich(s, 2.4, 5.32, 8.6, 1.05, [
        {"t": "Closed form or lookup table; evaluated per timestep in the robot model.", "size": 15, "color": TEXT, "space_after": 3},
        {"t": "Next: multiply by geometry (moment arm) \u2192 joint torque.  Does it survive contact with the joint?", "size": 15, "color": BROWN, "bold": True},
    ])
    pageno(s, n)

def build(prs, start=1):
    s01_title(prs, start)
    s02_roadmap(prs, start + 1)
    s03_project(prs, start + 2)
    s04_motivation(prs, start + 3)
    s05_bpa(prs, start + 4)
    s06_gaps(prs, start + 5)
    s07_jig(prs, start + 6)
    s08_fmax_discovery(prs, start + 7)
    s09_fmax_model(prs, start + 8)
    s10_festo_wrong(prs, start + 9)
    s11_surface(prs, start + 10)
    s12_surface_gof(prs, start + 11)
    s13_vs_models(prs, start + 12)
    s14_realtime(prs, start + 13)
    return start + 13
