"""Slides 26-37: Aim 3 (path to neural control), synthesis, close."""
import os
from pptx.enum.text import PP_ALIGN
from deck_lib import *

def s26_pipeline(prs, n):
    s = content_slide(prs, "Aim 3 \u00b7 simulation pipeline",
        "The plant moves into simulation: OpenSim \u2192 MuJoCo, verified end to end",
        cite="MyoConverter (Caggiano 2022); MuJoCo (Todorov 2012); model: modified Gait2392 robotbody",
        notes="Conversion is COMPLETE and verified: 92 actuators, attachment geometry audited against OpenSim source; joint sign conventions audited dynamically from simulation torques.")
    pic(s, os.path.join(FIG, "mujoco_robot.png"), box=(0.5, 1.85, 3.3, 5.0))
    tx(s, 0.5, 6.6, 3.3, 0.35, "converted model at reference pose", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    steps = [("Modified OpenSim robot model", "robot routing, 92 muscles, corrected paths"),
             ("MyoConverter", "OpenSim \u2192 MuJoCo, muscle kinematics + kinetics"),
             ("MuJoCo plant", "92 muscle-tendon actuators \u00b7 audited geometry \u00b7 dynamically audited joint signs"),
             ("SNS-Toolbox coupling", "0.1 ms neural step \u00b7 proprioceptive + contact channels in")]
    y = 2.0
    for i, (t, sub) in enumerate(steps):
        node(s, 4.4, y, 5.3, 0.88, t, fill=CARD if i < 3 else CARD2, size=14, bold=True)
        tx(s, 9.95, y + 0.06, 3.1, 0.85, sub, size=12, color=MUTED, spacing=1.05)
        if i < 3:
            arrow(s, 7.0, y + 0.9, 7.0, y + 1.16, color=FAINT)
        y += 1.22
    tx(s, 4.4, 6.72, 8.35, 0.4, "custom BPA force model ported 1:1 to Python - including the \u03c7 force balance", size=13, color=BROWN, bold=True)
    pageno(s, n)

def s27_knee_loop(prs, n):
    s = content_slide(prs, "Aim 3 \u00b7 proof of interface",
        "A synthetic antagonist-BPA knee already closes the loop in MuJoCo",
        cite="SNS-Toolbox (Nourse 2023); minimal SNS with reciprocal inhibition drives flexor + extensor BPAs",
        notes="This is the controller-actuator-joint interface demonstration - NOT yet SNS control of the full converted robot.")
    pic(s, os.path.join(FIG, "mujoco_knee.png"), box=(0.6, 1.9, 5.4, 4.8))
    tx(s, 0.6, 6.75, 5.4, 0.35, "synthetic knee at \u221235\u00b0; cyan/gold = named flexor/extensor routes", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    loop = [("SNS network", "0.1 ms step"), ("MN voltage \u2192 activation", "sigmoid"), ("MuJoCo mj_step", "tendon forces"),
            ("force \u2192 Ib afferent", "feedback"), ("back to SNS", "closed loop")]
    y = 2.1
    for i, (t, sub) in enumerate(loop):
        node(s, 6.6, y, 3.6, 0.72, t, fill=CARD2 if i in (0, 4) else CARD, size=13, bold=(i in (0, 4)))
        tx(s, 10.4, y + 0.12, 2.4, 0.6, sub, size=11.5, color=MUTED)
        if i < 4:
            arrow(s, 8.4, y + 0.74, 8.4, y + 0.94, color=FAINT)
        y += 0.94
    arrow(s, 6.5, 6.15, 6.5, 2.5, color=TEAL, width=2.0)
    tx(s, 4.9, 4.2, 1.5, 0.6, "loop", size=12, color=TEAL, bold=True)
    pageno(s, n)

def s28_rgpf(prs, n):
    s = content_slide(prs, "Aim 3 \u00b7 spinal controller",
        "The spinal controller: two-level rhythm generation + pattern formation, per leg",
        cite="Rybak 2006; McCrea & Rybak 2008; left-right coordination: Molkov 2015, Danner 2016",
        notes="406-neuron network assembled on the converted model. RG verified in ISOLATION: ~0.885 s period, antiphase -0.88. Full-body closed loop not yet running (leg NaNs under investigation).")
    pic(s, os.path.join(FIG, "walker_arch.png"), box=(0.55, 1.95, 7.2, 4.55))
    tx(s, 0.55, 6.55, 7.2, 0.35, "per-side subsystems: RG \u00b7 hip \u00b7 knee-ankle pattern formation", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    card(s, 8.15, 2.0, 4.55, 4.5, fill=CARD)
    rich(s, 8.4, 2.2, 4.05, 4.2, [
        {"t": "BUILT AND VERIFIED", "size": 12, "color": TEAL, "bold": True, "space_after": 5},
        {"t": "406-neuron spinal network on the converted model", "size": 14, "color": TEXT, "spacing": 1.1, "space_after": 6},
        {"t": "per-leg RG half-centers + 4 phase-shifted PF groups", "size": 14, "color": TEXT, "spacing": 1.1, "space_after": 6},
        {"t": "MN pools for all actuators; Ia / Ib / II afferent populations", "size": 14, "color": TEXT, "spacing": 1.1, "space_after": 6},
        {"t": "RG alone oscillates: \u22480.89 s cycle, legs in antiphase (\u22120.88)", "size": 14, "color": GREEN, "bold": True, "spacing": 1.1, "space_after": 6},
        {"t": "Open: full-body closed-loop co-simulation", "size": 12.5, "color": MUTED, "italic": True},
    ])
    pageno(s, n)

def s29_walker(prs, n):
    s = content_slide(prs, "Aim 3 \u00b7 the suspended-walker lesson",
        "Feedforward rhythm produces stepping; only tuned feedback can support the body",
        cite="Li et al. 2017 (biomech. synthesis); Deng et al. 2019 (two-layer CPG); AnimatLab: Cofer 2010",
        notes="This is the honest state of the neural side: rhythm generation is demonstrated; weight-bearing is the critical path, and that is a FEEDBACK tuning problem.")
    pic(s, os.path.join(FIG, "walker_body.png"), box=(0.55, 1.95, 3.4, 4.3))
    tx(s, 0.55, 6.3, 3.4, 0.5, "9-body planar biped,\nantagonist Hill pairs + passive toe", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    pic(s, os.path.join(FIG, "animatlab_phase1_joint_angles.png"), box=(4.15, 1.95, 4.4, 4.3))
    tx(s, 4.15, 6.3, 4.4, 0.5, "phase-1 instrumented run: 10 s,\nhip / knee / ankle, both sides", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    card(s, 8.85, 1.95, 3.85, 4.6, fill=CARD)
    rich(s, 9.1, 2.15, 3.4, 4.3, [
        {"t": "WHAT IT TAUGHT", "size": 12, "color": TEAL, "bold": True, "space_after": 4},
        {"t": "Suspended: alternating stepping in hip, knee, ankle - rhythm does not need sensory feedback (as the run at right shows).", "size": 13.5, "color": GREEN, "bold": True, "spacing": 1.1, "space_after": 8},
        {"t": "On the ground: untuned afferents could not support the pattern.", "size": 13.5, "color": RED, "bold": True, "spacing": 1.1, "space_after": 8},
        {"t": "\u2192 proprioceptive + contact tuning is the critical path to weight-bearing walking.", "size": 13.5, "color": BROWN, "bold": True, "spacing": 1.1, "space_after": 8},
        {"t": "Current: RG latches under tonic drive in the bilateral model - conductance tuning in progress.", "size": 12, "color": MUTED, "italic": True, "spacing": 1.08},
    ])
    pageno(s, n)

def s30_sensory(prs, n):
    s = content_slide(prs, "Aim 3 \u00b7 closing the loop",
        "Four afferent channels close the loop, organized by the Sensory Afferent Database",
        cite="Sensory Afferent Database: Bolen 2023; AI-assisted expansion: Bolen 2024",
        notes="The database (Ben led both papers) is what makes the gain choices literature-grounded rather than arbitrary.")
    pic(s, os.path.join(FIG, "sensorimotor_network.png"), box=(0.55, 1.95, 6.1, 4.6))
    tx(s, 0.55, 6.6, 6.1, 0.35, "sensorimotor organization: MN, IN, Ia / Ib pathways", size=11, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    rows = [
        ("Ia", "muscle length + velocity", "stretch reflex; phase-dependent PF modulation", GREEN),
        ("Ib", "muscle force", "load-dependent regulation of stance extensor activity", GREEN),
        ("II", "muscle length", "position sense; interjoint coordination", GREEN),
        ("heel / toe", "contact", "discrete stance-swing boundaries; phase transitions", TEAL),
    ]
    y = 2.1
    for k, sig, eff, c in rows:
        card(s, 7.0, y, 5.75, 1.02, fill=CARD)
        tx(s, 7.2, y + 0.26, 1.15, 0.5, k, size=15, color=c, bold=True)
        rich(s, 8.45, y + 0.12, 4.15, 0.9, [
            {"t": sig, "size": 13, "color": MUTED, "space_after": 1},
            {"t": eff, "size": 12.5, "color": TEXT, "spacing": 1.0, "bold": True},
        ])
        y += 1.14
    pageno(s, n)

def s31_verification(prs, n):
    s = content_slide(prs, "Aim 3 \u00b7 verification-first",
        "Verify each circuit against tonic-stimulation physiology before any gait tuning",
        cite="air-stepping analogue: Selionov 2009; deletion/stimulation: Rybak 2006; Markin 2010",
        notes="Methodological discipline slide: no walking-speed tuning until sign, gain, and phase structure replicate the literature at each level.")
    levels = [("motoneuron pools", "tonic stim \u2192 force / activation"), ("pattern-formation layer", "tonic stim \u2192 phase-appropriate activation"),
              ("rhythm generator", "tonic stim \u2192 locomotor rhythm")]
    y = 4.6
    for i, (t, sub) in enumerate(levels):
        node(s, MX, y, 6.4, 0.68, t, fill=CARD if i < 2 else CARD2, size=14, bold=True)
        tx(s, 7.15, y + 0.14, 3.2, 0.5, sub, size=11.5, color=MUTED)
        if i < 2:
            arrow(s, MX + 3.2, y - 0.28, MX + 3.2, y - 0.02, color=FAINT)
        y -= 1.28
    tx(s, MX, 6.75, 6.4, 0.4, "stimulate each level; compare with literature", size=12.5, color=TEAL, bold=True)
    card(s, 10.55, 2.35, 2.2, 4.3, fill=CARD)
    rich(s, 10.7, 2.55, 1.95, 4.0, [
        {"t": "VERIFY", "size": 13, "color": BROWN, "bold": True, "align": PP_ALIGN.CENTER, "space_after": 6},
        {"t": "sign", "size": 15, "color": TEXT, "align": PP_ALIGN.CENTER, "space_after": 4},
        {"t": "gain", "size": 15, "color": TEXT, "align": PP_ALIGN.CENTER, "space_after": 4},
        {"t": "phase", "size": 15, "color": TEXT, "align": PP_ALIGN.CENTER, "space_after": 10},
        {"t": "then - and only then - tune walking speed", "size": 12.5, "color": MUTED, "align": PP_ALIGN.CENTER, "spacing": 1.1},
    ])
    pageno(s, n)

def s32_hardware(prs, n):
    s = content_slide(prs, "Aim 3 \u00b7 physical realization",
        "The destination: embedded neural control on a treadmill biped",
        cite="quadruped precedent: Lutz 2025 (MuJoCo-SNS co-sim \u2192 hardware); stack: Teensy + Jetson (SNS-Toolbox backend)",
        notes="The division of labor: Teensy = hard real-time per-joint; Jetson = network evaluation. Static maps of this dissertation = feedforward lookup; Elzein model = pulse-modulated pressure dynamics.")
    pic(s, os.path.join(MEDIA, "image44.jpg"), box=(0.6, 1.95, 5.4, 4.05))
    tx(s, 0.6, 6.05, 5.4, 0.4, "AARL BPA quadruped: SNS control, treadmill hardware", size=11.5, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    rows = [
        ("Teensy per-joint loops", "valve control, pressure + contact sensing (hard real-time)"),
        ("Jetson-class companion", "SNS network evaluation \u2192 activations"),
        ("Dissertation maps on board", "F(l\u2080, P, \u03b5) and \u03c7-corrected torque as feedforward / lookup"),
        ("Contraction lever", "reverse-pulley block-and-tackle: force for travel where routes run long"),
    ]
    y = 1.95
    for t, sub in rows:
        card(s, 6.55, y, 6.2, 0.92, fill=CARD)
        rich(s, 6.8, y + 0.1, 5.7, 0.8, [
            {"t": t, "size": 14, "color": TEXT, "bold": True, "space_after": 1},
            {"t": sub, "size": 12, "color": MUTED, "spacing": 1.0},
        ])
        y += 1.04
    card(s, 6.55, y, 6.2, 0.75, fill=CARD2)
    tx(s, 6.8, y + 0.12, 5.7, 0.55, "aspiration: dynamic whole-body maneuvers - a kickflip exercises every layer",
       size=13, color=BROWN, bold=True)
    pageno(s, n)

def s33_chain(prs, n):
    s = content_slide(prs, "Synthesis",
        "Three contributions, one validated design chain",
        cite="Dissertation Ch. 2-6",
        notes="Each stage feeds the next. The committee read these numbers; this is the one-slide summary.")
    cards_ = [
        ("ACTUATOR", "F(l\u2080, P, \u03b5)", "537 tests \u00b7 dimensionless surface\nadj. R\u00b2 > 0.99 in validation", GREEN),
        ("CORRECTIONS", "\u03c7\u2080 \u03c7\u2081 \u03c7\u2082 \u03c7\u2083", "identified on pinned knee,\none stiffness pair across configs\nbeats baseline in every pool test", GREEN),
        ("JOINT", "meets / exceeds human", "biomimetic knee vs Gait2392 torque\nover most of range of motion", GREEN),
        ("CONTROLLER - NEXT", "RG + PF spinal network", "rhythm verified in isolation;\nfeedback tuning = critical path", TEAL),
    ]
    x = MX; cw = 2.95; gap = 0.15
    for i, (k, big, sub, c) in enumerate(cards_):
        card(s, x, 2.2, cw, 3.6, fill=CARD2 if i == 3 else CARD)
        tx(s, x + 0.2, 2.4, cw - 0.4, 0.4, k, size=11.5, color=c, bold=True)
        tx(s, x + 0.2, 2.9, cw - 0.4, 0.9, big, size=17, color=BROWN, bold=True, spacing=1.0)
        tx(s, x + 0.2, 3.95, cw - 0.4, 1.7, sub, size=12.5, color=MUTED, spacing=1.12)
        if i < 3:
            arrow(s, x + cw + 0.01, 4.0, x + cw + gap - 0.01, 4.0, color=FAINT, width=2.2)
        x += cw + gap
    tx(s, MX, 6.2, 12.2, 0.5, "Each stage feeds the next: force model \u2192 torque model \u2192 plant for controller design \u2192 validated architecture for the robot.",
       size=14.5, color=TEXT, align=PP_ALIGN.CENTER)
    pageno(s, n)

def s34_limits(prs, n):
    s = content_slide(prs, "Honest boundaries",
        "What these results do not claim",
        cite="Dissertation Ch. 5 (Discussion)",
        notes="Pre-empt the committee's sharpest questions by stating them first.")
    items = [
        ("Isometric only", "no dynamics: hysteresis, airflow, damping; dynamic characterization is future work"),
        ("\u03b5\u2086\u2082\u2080 not predictable a priori", "record max contraction per individual muscle; variance across muscles is real"),
        ("Stiffness pair not unique", "the pinned-flexor data admit a family; one consistent pair adopted, ties broken cross-configuration"),
        ("Bio-extensor variance", "FVU 2.35 > 1: improvement over baseline is real, validated prediction is not claimed"),
        ("Static map", "manufacturing tolerance \u00b110%; aging and lifecycle behavior unmeasured"),
    ]
    x = MX; y = 2.1; cw = 6.0
    for i, (t, sub) in enumerate(items):
        cx = x + (i % 2) * (cw + 0.2)
        cy = y + (i // 2) * 1.55
        card(s, cx, cy, cw, 1.38, fill=CARD)
        rich(s, cx + 0.25, cy + 0.14, cw - 0.5, 1.15, [
            {"t": t, "size": 15, "color": BROWN, "bold": True, "space_after": 2},
            {"t": sub, "size": 12.5, "color": TEXT, "spacing": 1.05},
        ])
    pageno(s, n)

def s35_future(prs, n):
    s = content_slide(prs, "Future work",
        "From here: two-muscle tests now, treadmill biped next",
        cite="Dissertation Ch. 6",
        notes="In progress NOW: two-BPA parallel flexor tests + surrogateopt revisit. Then dynamic characterization, closed-loop co-sim, sensorium extensions, hardware.")
    stages = [
        ("NOW", "2-BPA flexor tests;\nsurrogateopt revisit", CARD),
        ("Dynamics", "isokinetic / isobaric /\nisotonic + Elzein model", CARD),
        ("Closed loop", "RG oscillation tuning \u2192\nco-sim \u2192 tonic-stim verification", CARD2),
        ("Sensorium", "vestibular (IMU) \u00b7 ocular drift\ncorrection \u00b7 cerebellar gating", CARD),
        ("Hardware", "treadmill biped:\nTeensy + Jetson + SNS", CARD2),
    ]
    x = MX; cw = 2.32; gap = 0.12
    for i, (t, sub, fill) in enumerate(stages):
        card(s, x, 2.3, cw, 2.5, fill=fill)
        tx(s, x + 0.18, 2.5, cw - 0.36, 0.5, t.upper(), size=13, color=TEAL if i > 0 else GREEN, bold=True)
        tx(s, x + 0.18, 3.05, cw - 0.36, 1.6, sub, size=12.5, color=TEXT, spacing=1.12)
        if i < 4:
            arrow(s, x + cw + 0.005, 3.55, x + cw + gap - 0.005, 3.55, color=FAINT, width=2.0)
        x += cw + gap
    tx(s, MX, 5.3, 12.2, 0.5, "Long-term: dynamic whole-body maneuvers - the kickflip aspiration.",
       size=15, color=BROWN, bold=True, align=PP_ALIGN.CENTER)
    tx(s, MX, 5.95, 12.2, 0.7, "A robot whose nervous system - unlike a living animal's - can be instrumented, stimulated, and lesioned at will.",
       size=14, color=MUTED, italic=True, align=PP_ALIGN.CENTER)
    pageno(s, n)

def s36_ack(prs, n):
    s = content_slide(prs, None,
        "Acknowledgments",
        assert_size=34,
        notes="Keep brief; thank by name.")
    people = [
        ("Advisor", "M. Anthony Hunt"),
        ("Committee", "thank you for reading, commenting, and being here"),
        ("Collaborators", "M. Elzein \u00b7 L. Pang \u00b7 L. Burgess \u00b7 C. Morrow \u00b7 J. Scharzenberger \u00b7 C. Steele \u00b7 J. Chung"),
        ("Lab & support", "Agile & Adaptive Robotics Laboratory \u00b7 Portland State University \u00b7 NeuroNex community"),
    ]
    y = 2.3
    for k, v in people:
        rich(s, 1.6, y, 10.2, 1.0, [
            {"t": k.upper(), "size": 12, "color": TEAL, "bold": True, "space_after": 2},
            {"t": v, "size": 17, "color": TEXT},
        ])
        y += 1.12
    pic(s, os.path.join(ASSETS, "aarl_logo_0.png"), box=(4.87, 6.55, 3.6, 0.8))
    pageno(s, n)

def s37_thanks(prs, n):
    s = prs.slides.add_slide(prs.slide_layouts[6]); bg(s)
    pic(s, os.path.join(ASSETS, "aarl_logo_0.png"), box=(4.47, 1.05, 4.4, 1.34))
    tx(s, 0.7, 3.2, 11.93, 1.1, "Thank you - questions?", size=44, color=BROWN, bold=True, align=PP_ALIGN.CENTER)
    tx(s, 0.7, 4.5, 11.93, 0.5, "Movement and Control of Biomimetic Humanoid Robots", size=16, color=MUTED, align=PP_ALIGN.CENTER)
    tx(s, 0.7, 5.05, 11.93, 0.5, "Ben P. Bolen  \u00b7  bbolen83@gmail.com", size=14, color=MUTED, align=PP_ALIGN.CENTER)
    tx(s, 0.7, 6.6, 11.93, 0.4, "Appendix: model details, protocol, robustness checks", size=12, color=FAINT, italic=True, align=PP_ALIGN.CENTER)
    s.notes_slide.notes_text_frame.text = "Anticipated questions route to the appendix: uniqueness of stiffness pair (G), bracket point (H), frame conventions (G), wrap-loss implementation (E)."

def build(prs, start=26):
    s26_pipeline(prs, 26)
    s27_knee_loop(prs, 27)
    s28_rgpf(prs, 28)
    s29_walker(prs, 29)
    s30_sensory(prs, 30)
    s31_verification(prs, 31)
    s32_hardware(prs, 32)
    s33_chain(prs, 33)
    s34_limits(prs, 34)
    s35_future(prs, 35)
    s36_ack(prs, 36)
    s37_thanks(prs, 37)
    return 37
