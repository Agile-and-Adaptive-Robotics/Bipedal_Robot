"""PowerPoint figure factory for the dissertation schematics.

Builds every flowchart/schematic as NATIVE editable PowerPoint shapes
(Ben's directive: no TikZ; PowerPoint is the source of truth), then a
companion COM step exports each to PDF for the dissertation tree.

Rules baked in: Arial, nothing below 10 pt, explicit box sizes with generous
padding, labels kept short enough to fit with clearance, one arrow style,
light gray fills, dashed outlines for future-work items.
"""
import copy
from pptx import Presentation
from pptx.util import Cm, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn

FIG = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures"
INK = RGBColor(0x1F, 0x1F, 0x1F)
FILL = RGBColor(0xF2, 0xF2, 0xF2)
GRAY = RGBColor(0x8A, 0x8A, 0x8A)

def new_slide(prs, w_cm, h_cm):
    prs.slide_width = Cm(w_cm)
    prs.slide_height = Cm(h_cm)
    return prs.slides.add_slide(prs.slide_layouts[6])

def _ln(shape, dash=None, w_pt=1.0):
    ln = shape.line._get_or_add_ln()
    ln.set('w', str(int(w_pt * 12700)))
    if dash:
        d = ln.find(qn('a:prstDash'))
        if d is None:
            d = ln.makeelement(qn('a:prstDash'), {})
            ln.append(d)
        d.set('val', dash)

def box(slide, x, y, w, h, text, dashed=False, fill=FILL, bold_first=False,
        size=10):
    sp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Cm(x), Cm(y),
                                Cm(w), Cm(h))
    sp.adjustments[0] = 0.10
    sp.fill.solid(); sp.fill.fore_color.rgb = fill
    sp.line.color.rgb = INK; sp.line.width = Pt(1.0)
    _ln(sp, 'dash' if dashed else None)
    tf = sp.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Cm(0.12)
    tf.margin_top = tf.margin_bottom = Cm(0.06)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    lines = text.split('\n')
    for i, line in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = PP_ALIGN.CENTER
        r = p.add_run(); r.text = line
        r.font.size = Pt(size); r.font.name = 'Arial'
        r.font.color.rgb = INK
        if bold_first and i < bold_first:
            r.font.bold = True
    return sp

def label(slide, x, y, w, text, size=10, bold=False, align=PP_ALIGN.CENTER,
          color=INK):
    tb = slide.shapes.add_textbox(Cm(x), Cm(y), Cm(w), Cm(0.75))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Cm(0.03)
    tf.margin_top = tf.margin_bottom = Cm(0.01)
    for i, line in enumerate(text.split('\n')):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = line
        r.font.size = Pt(size); r.font.name = 'Arial'
        r.font.bold = bold; r.font.color.rgb = color
    return tb

def _arrowhead(conn, head=False, tail=True):
    ln = conn.line._get_or_add_ln()
    for tag, on in (('a:headEnd', head), ('a:tailEnd', tail)):
        if on:
            e = ln.find(qn(tag))
            if e is None:
                e = ln.makeelement(qn(tag), {})
                ln.append(e)
            e.set('type', 'triangle'); e.set('w', 'med'); e.set('len', 'med')

def arrow(slide, x1, y1, x2, y2, both=False):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Cm(x1), Cm(y1),
                                   Cm(x2), Cm(y2))
    c.line.color.rgb = INK; c.line.width = Pt(1.0)
    _arrowhead(c, head=both)
    return c

def seg(slide, x1, y1, x2, y2, dash=None):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Cm(x1), Cm(y1),
                                   Cm(x2), Cm(y2))
    c.line.color.rgb = INK; c.line.width = Pt(1.0)
    _ln(c, dash)
    return c

def leader(slide, x1, y1, x2, y2, text, size=10):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Cm(x1), Cm(y1),
                                   Cm(x2), Cm(y2))
    c.line.color.rgb = INK; c.line.width = Pt(0.75)
    d = c.line._get_or_add_ln().find(qn('a:prstDash'))
    label(slide, x2 - (4.0 if x2 > x1 else 0.0), y2 - 0.38, 4.0, text,
          size=size, align=PP_ALIGN.LEFT if x2 > x1 else PP_ALIGN.RIGHT)
    return c

# ------------------------------------------------------------- 1 roadmap
prs = Presentation()
s = new_slide(prs, 16.5, 6.9)
W, H = 4.35, 1.55
row1, row2 = 1.7, 4.9
xs = [0.3, 6.1, 11.9]
top = [("Human biomechanics\nbenchmark (Ch. 2)", False),
       ("BPA force models,\nvalidated (Ch. 3)", False),
       ("Torque correction\nfactors and routes (Ch. 4)", False)]
bot = [("Sensory Afferent\nDatabase (Ch. 5)", False),
       ("Neuromechanical\nsimulation tools (Ch. 6)", False),
       ("Physical BPA robot\nwith SNS control", True)]
for i, (t, d) in enumerate(top):
    box(s, xs[i], row1, W, H, t, dashed=d, bold_first=1)
for i, (t, d) in enumerate(bot):
    box(s, xs[i], row2, W, H, t, dashed=d, bold_first=1)
ymid1 = row1 + H / 2; ymid2 = row2 + H / 2
arrow(s, xs[0] + W, ymid1, xs[1] - 0.05, ymid1)
label(s, xs[0] + W, row1 - 0.65, xs[1] - xs[0] - W, "targets")
arrow(s, xs[1] + W, ymid1, xs[2] - 0.05, ymid1)
label(s, xs[1] + W, row1 - 0.65, xs[2] - xs[1] - W, "actuator model")
arrow(s, xs[0] + W, ymid2, xs[1] - 0.05, ymid2)
label(s, xs[0] + W, row2 - 0.65, xs[1] - xs[0] - W, "sensory channels")
arrow(s, xs[1] + W, ymid2, xs[2] - 0.05, ymid2)
label(s, xs[1] + W, row2 - 0.65, xs[2] - xs[1] - W, "controller and\nbody model")
# torque -> tools elbow
seg(s, xs[2] + W / 2, row1 + H, xs[2] + W / 2, row1 + H + 0.55)
seg(s, xs[2] + W / 2, row1 + H + 0.55, xs[1] + W / 2, row1 + H + 0.55)
arrow(s, xs[1] + W / 2, row1 + H + 0.55, xs[1] + W / 2, row2 - 0.05)
label(s, xs[2] + W / 2 + 0.15, row1 + H + 0.10, 3.4, "routes and\ntorque model",
      align=PP_ALIGN.LEFT)
# feedback return along the bottom into the benchmark's left side
yb = row2 + H + 0.75
seg(s, xs[2] + W / 2, row2 + H, xs[2] + W / 2, yb)
seg(s, xs[2] + W / 2, yb, 0.0 - 0.85 + 0.3 - 0.3, yb)
seg(s, 0.0 - 0.85 + 0.3 - 0.3, yb, -0.55, yb)
seg(s, -0.55, yb, -0.55, ymid1)
arrow(s, -0.55, ymid1, xs[0] - 0.05, ymid1)
label(s, 4.5, yb - 0.62, 6.0, "biological insight")
# key
seg(s, 0.3, 0.85, 0.9, 0.85)
label(s, 1.0, 0.55, 5.4, "validated in this dissertation", align=PP_ALIGN.LEFT)
seg(s, 6.6, 0.85, 7.2, 0.85, dash='dash')
label(s, 7.3, 0.55, 2.6, "future work", align=PP_ALIGN.LEFT)
prs.save(FIG + r"\10-introduction\components\roadmap.pptx")

# ------------------------------------------------------------- 2 pipeline
prs = Presentation()
s = new_slide(prs, 17.2, 7.6)
def band(x, w, title):
    sp = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Cm(x), Cm(1.35),
                            Cm(w), Cm(5.85))
    sp.fill.background()
    sp.line.color.rgb = GRAY; sp.line.width = Pt(0.75)
    _ln(sp, 'dash')
    label(s, x, 0.55, w, title, bold=True)
band(0.15, 4.05, "Model conversion")
band(4.75, 6.05, "Closed-loop simulation")
band(11.35, 5.7, "Hardware deployment")
box(s, 0.45, 1.65, 3.45, 1.15, "OpenSim robot-body model\n(.osim, Gait2392 routing)")
box(s, 0.45, 3.25, 3.45, 0.85, "MyoConverter")
box(s, 0.45, 4.95, 3.45, 1.15, "MuJoCo MJCF model\n(activation-state actuators)")
arrow(s, 2.17, 2.80, 2.17, 3.20)
arrow(s, 2.17, 4.10, 2.17, 4.90)
box(s, 5.05, 1.80, 2.55, 1.35, "MuJoCo physics with\nBPA force model\n(Ch. 3)")
box(s, 5.05, 5.05, 2.55, 1.15, "SNS-Toolbox\nnetwork (Python)")
arrow(s, 5.55, 5.00, 5.55, 3.20)
label(s, 3.55, 3.72, 1.95, "per-muscle\nactivation", align=PP_ALIGN.RIGHT)
arrow(s, 7.30, 3.20, 7.30, 5.00)
label(s, 7.40, 3.72, 3.30, "Ia, Ib, and II afferents;\nheel and toe contact",
      align=PP_ALIGN.LEFT)
box(s, 11.65, 1.65, 4.15, 1.05, "Jetson Orin:\nSNS network evaluation")
box(s, 11.65, 3.55, 4.15, 1.45, "Teensy-class MCUs:\nper-joint valve, pressure,\nand contact loops")
box(s, 11.65, 5.85, 4.15, 1.0, "Treadmill biped")
arrow(s, 13.7, 2.70, 13.7, 3.50)
label(s, 13.85, 2.85, 3.3, "commands", align=PP_ALIGN.LEFT)
arrow(s, 13.7, 5.00, 13.7, 5.80)
label(s, 13.85, 5.15, 3.3, "valve drive", align=PP_ALIGN.LEFT)
arrow(s, 12.15, 5.80, 12.15, 5.02)
label(s, 9.30, 5.15, 2.75, "sensors", align=PP_ALIGN.RIGHT)
# cross-phase connectors
seg(s, 3.90, 5.52, 4.45, 5.52)
seg(s, 4.45, 5.52, 4.45, 2.48)
arrow(s, 4.45, 2.48, 5.00, 2.48)
seg(s, 7.60, 5.62, 10.95, 5.62)
seg(s, 10.95, 5.62, 10.95, 2.18)
arrow(s, 10.95, 2.18, 11.60, 2.18)
label(s, 2.2, 7.30, 12.5, "Computer-in-the-loop: the embedded stack first executes against the simulated body, then against the physical one.")
prs.save(FIG + r"\60-sns\components\pipeline.pptx")

# --------------------------------------------------------- 3 curriculum
prs = Presentation()
s = new_slide(prs, 16.5, 10.4)
LX = {"stage": 1.45, "cap": 5.35, "score": 10.55}
label(s, LX["stage"], 0.30, 3.4, "Stage", bold=True, align=PP_ALIGN.LEFT)
label(s, LX["cap"], 0.30, 4.6, "Capability added", bold=True, align=PP_ALIGN.LEFT)
label(s, LX["score"], 0.30, 5.4, "Stage score", bold=True, align=PP_ALIGN.LEFT)
rows = [
    ("1", "Deafferented\nair stepping", "half-center-free oscillation,\nlegs suspended, feedback off",
     "rhythm and swing depth\n(rhythm-gated)", "D9E8F4"),
    ("2", "Afferented\nair stepping", "heel, toe, and load\nfeedback active",
     "rhythm and swing depth\n(rhythm-gated)", "BFD7EA"),
    ("3", "Standing\nbalance", "posture held under\npush perturbations",
     "posture maintenance\nunder push", "9BC0DD"),
    ("4", "Ground\nwalking", "ground contact\nand forward gait",
     "kinematic score vs human\nreference cycle", "6FA3CE"),
    ("5", "Contact-driven\nground walking", "contact-driven pathway\nvariant of the walker",
     "kinematic score vs human\nreference cycle", "4B8ABE"),
]
y = 1.05
for num, name, cap, score, chip in rows:
    sp = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Cm(0.3), Cm(y),
                            Cm(15.9), Cm(1.35))
    sp.fill.solid(); sp.fill.fore_color.rgb = FILL
    sp.line.color.rgb = INK; sp.line.width = Pt(1.0)
    ch = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Cm(0.45), Cm(y + 0.35),
                            Cm(0.65), Cm(0.65))
    ch.fill.solid()
    ch.fill.fore_color.rgb = RGBColor.from_string(chip)
    ch.line.color.rgb = INK; ch.line.width = Pt(0.75)
    tf = ch.text_frame; tf.word_wrap = False
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = num
    r.font.size = Pt(11); r.font.bold = True; r.font.name = 'Arial'
    r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF) if num in "45" else INK
    label(s, LX["stage"], y + 0.28, 3.6, name, bold=True, align=PP_ALIGN.LEFT)
    label(s, LX["cap"], y + 0.28, 5.0, cap, align=PP_ALIGN.LEFT)
    label(s, LX["score"], y + 0.28, 5.2, score, align=PP_ALIGN.LEFT)
    y += 1.62
for i in range(4):
    arrow(s, 0.775, 1.05 + 1.35 + i * 1.62 + 0.02, 0.775, 1.05 + (i + 1) * 1.62 - 0.02)
label(s, 0.3, 9.35, 15.9, "Each stage's winner seeds the next; only the newly activated pathway gains are searched.",
      align=PP_ALIGN.LEFT)
label(s, 0.3, 9.75, 15.9, "Finalists at every stage must pass the perturbation and basin gates. Darker chips indicate less harness support.",
      align=PP_ALIGN.LEFT)
prs.save(FIG + r"\60-sns\components\sns_curriculum.pptx")

# ------------------------------------------------------ 4 design workflow
prs = Presentation()
s = new_slide(prs, 16.5, 7.3)
box(s, 0.2, 0.6, 3.3, 1.45, "Human torque and\nrange-of-motion targets\n(OpenSim benchmark)")
box(s, 5.0, 0.6, 3.3, 1.45, "Select actuator set,\ndiameter, and\nresting length")
box(s, 0.2, 4.0, 3.3, 1.2, "BPA force model\n(Ch. 3)")
box(s, 5.0, 4.0, 3.3, 1.45, "Design routes under the\ncorrected torque model\n(Ch. 4)")
dia = s.shapes.add_shape(MSO_SHAPE.DIAMOND, Cm(9.4), Cm(3.85), Cm(3.1), Cm(1.75))
dia.fill.solid(); dia.fill.fore_color.rgb = FILL
dia.line.color.rgb = INK; dia.line.width = Pt(1.0)
tf = dia.text_frame; tf.word_wrap = True
p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
r = p.add_run(); r.text = "Meets\ntargets?"
r.font.size = Pt(10); r.font.name = 'Arial'; r.font.color.rgb = INK
box(s, 13.2, 4.0, 3.0, 1.2, "Commit to\nhardware")
box(s, 11.2, 6.0, 4.6, 1.1, "Verify in closed-loop\nsimulation (Ch. 6)", dashed=True)
arrow(s, 3.55, 1.32, 4.95, 1.32)
label(s, 3.55, 0.70, 1.4, "targets")
arrow(s, 3.55, 4.60, 4.95, 4.60)
label(s, 3.55, 4.95, 1.4, "force surface")
arrow(s, 6.65, 2.10, 6.65, 3.95)
label(s, 6.85, 2.70, 2.4, "actuator set", align=PP_ALIGN.LEFT)
arrow(s, 8.35, 4.72, 9.35, 4.72)
label(s, 7.3, 3.30, 4.4, "predicted torque and\nrange of motion")
arrow(s, 12.55, 4.72, 13.15, 4.72)
label(s, 12.35, 4.15, 0.9, "yes")
seg(s, 10.95, 5.62, 10.95, 6.35)
arrow(s, 10.95, 6.35, 11.15, 6.35)
seg(s, 10.95, 5.62, 6.65, 5.62)
seg(s, 6.65, 5.62, 6.65, 5.50)
arrow(s, 6.65, 6.35, 6.65, 5.50)
label(s, 6.9, 5.75, 7.6, "no: redesign routes or reselect actuators", align=PP_ALIGN.LEFT)
seg(s, 14.7, 5.20, 14.7, 5.85)
arrow(s, 14.7, 5.85, 14.15, 5.88)
prs.save(FIG + r"\70-overall\components\design_workflow.pptx")

# -------------------------------------------------- 5 sensorium hierarchy
prs = Presentation()
s = new_slide(prs, 15.6, 8.6)
label(s, 0.2, 0.25, 8.0, "Supraspinal additions (future work)", bold=True, align=PP_ALIGN.LEFT)
label(s, 0.2, 3.15, 8.0, "Spinal network (implemented, Ch. 6)", bold=True, align=PP_ALIGN.LEFT)
label(s, 0.2, 6.05, 8.0, "Body (demonstrated in hardware)", bold=True, align=PP_ALIGN.LEFT)
box(s, 0.2, 0.85, 4.35, 1.55, "Basal ganglia\naction and mode\nselection", dashed=True)
box(s, 5.6, 0.85, 4.35, 1.55, "Cerebellar adaptive\ngateway: gait gating\nand phase tuning", dashed=True)
box(s, 11.0, 0.85, 4.4, 1.55, "Vestibular and ocular\nsensing: IMU orientation,\nvisual drift correction", dashed=True)
box(s, 0.2, 3.75, 4.35, 1.45, "Two-level CPG: rhythm\ngenerator and pattern\nformation")
box(s, 5.6, 3.75, 4.75, 1.45, "Motor-neuron pools;\nIa, Ib, and II afferents;\nheel and toe contact")
box(s, 0.2, 6.65, 13.0, 1.2, "BPA actuators and robot body")
arrow(s, 9.55, 2.45, 9.55, 3.70)
label(s, 9.75, 2.75, 3.4, "descending\ncommands", align=PP_ALIGN.LEFT)
arrow(s, 6.20, 5.25, 6.20, 6.60)
label(s, 3.60, 5.55, 2.45, "activation", align=PP_ALIGN.RIGHT)
arrow(s, 8.55, 6.60, 8.55, 5.28)
label(s, 8.75, 5.55, 3.8, "proprioceptive and\ncontact feedback", align=PP_ALIGN.LEFT)
arrow(s, 13.20, 6.60, 13.20, 2.45)
label(s, 13.40, 4.35, 2.2, "body orientation\nand vision\n(future)", align=PP_ALIGN.LEFT)
prs.save(FIG + r"\70-overall\components\sensorium_hierarchy.pptx")

# ------------------------------------------------- 6 BPA construction FBD
prs = Presentation()
s = new_slide(prs, 16.5, 10.6)
label(s, 0.2, 0.15, 10.0, "(A) Longitudinal section with loads and constraints",
      bold=True, align=PP_ALIGN.LEFT)
# geometry (cm on slide): actuator axis horizontal at y=2.6
x0, x1, ym = 3.2, 12.4, 2.75
ro, ri = 0.95, 0.55
# bladder
bl = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Cm(x0 + 0.75), Cm(ym - ri * 0.62),
                        Cm(x1 - x0 - 1.5), Cm(2 * ri * 0.62))
bl.fill.solid(); bl.fill.fore_color.rgb = RGBColor(0xE3, 0xE3, 0xE3)
bl.line.color.rgb = GRAY; bl.line.width = Pt(0.75)
bl.adjustments[0] = 0.5
# sleeve
sv = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Cm(x0 + 0.45), Cm(ym - ro),
                        Cm(x1 - x0 - 0.9), Cm(2 * ro))
sv.fill.background()
sv.line.color.rgb = INK; sv.line.width = Pt(1.75)
sv.adjustments[0] = 0.5
# braid hatch: crossing thin lines
import math
n = 12
span = (x1 - x0 - 1.3)
for k in range(n):
    xa = x0 + 0.75 + k * (span / n)
    xb = xa + span / n
    seg(s, xa, ym - ro * 0.85, xb, ym + ro * 0.85)
    seg(s, xa, ym + ro * 0.85, xb, ym - ro * 0.85)
# end caps with fixed-support hatching
for xc in (x0, x1):
    cap = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Cm(xc - 0.5), Cm(ym - 1.15),
                             Cm(1.0), Cm(2.3))
    cap.fill.solid(); cap.fill.fore_color.rgb = RGBColor(0xC9, 0xC9, 0xC9)
    cap.line.color.rgb = INK; cap.line.width = Pt(1.0)
    for k in range(-3, 4):  # hatch marks beyond the fixed wall
        seg(s, xc - 0.5 - 0.35, ym + k * 0.4 - 0.15, xc - 0.5, ym + k * 0.4 + 0.15)
# internal pressure arrows (inside -> out)
for k in range(7):
    xa = x0 + 1.4 + k * ((x1 - x0 - 2.8) / 6)
    arrow(s, xa, ym - 0.08, xa, ym - 0.62)
    arrow(s, xa, ym + 0.08, xa, ym + 0.62)
# axial reaction and applied load
arrow(s, x1 + 0.35, ym, x1 + 1.85, ym, both=False)
label(s, x1 + 0.55, ym - 0.72, 2.6, "axial force F\n(ends fixed)")
arrow(s, x0 - 1.85, ym, x0 - 0.35, ym, both=False)
label(s, x0 - 2.6, ym - 0.72, 2.1, "reaction R")
# contact pair symbols between layers
arrow(s, 7.2, ym - ri * 0.62 - 0.35, 7.2, ym - ri * 0.62 + 0.05, both=True)
arrow(s, 8.9, ym + ri * 0.62 + 0.35, 8.9, ym + ri * 0.62 - 0.05, both=True)
# symmetry centerline, dash-dot
cl = seg(s, 1.3, ym, 14.9, ym, dash='sysDot')
label(s, 12.6, ym + 0.12, 2.4, "symmetry axis")
# callouts
leader(s, 5.4, ym - ro, 3.6, 0.75, "braided aramid sleeve\n(two handed weaves)")
leader(s, 10.6, ym + ri * 0.5, 11.4, 0.75, "elastomer bladder")
leader(s, x0 - 0.2, ym + 1.15, 1.2, 4.05, "end cap, clamped\n(Onyx or aluminum)")
leader(s, 7.2, ym - ri * 0.62 - 0.4, 4.9, 4.05, "layer contact pairs\n(braid--bladder, braid--cap)")
label(s, 0.2, 5.0, 15.0, "Pressure acts inside-out on the bladder; the braided sleeve carries the tension and converts radial expansion to axial contraction.", size=10)

# panel B: exploded view along the axis
label(s, 0.2, 5.75, 10.0, "(B) Exploded view", bold=True, align=PP_ALIGN.LEFT)
yy = 6.35
parts = [
    ("end cap", 2.2, 1.6, RGBColor(0xC9, 0xC9, 0xC9)),
    ("braided sleeve", 3.9, 4.6, None),
    ("elastomer bladder", 8.3, 3.4, RGBColor(0xE3, 0xE3, 0xE3)),
    ("end cap", 12.6, 1.6, RGBColor(0xC9, 0xC9, 0xC9)),
]
# alignment axis
seg(s, 1.2, yy + 0.5, 15.3, yy + 0.5, dash='sysDot')
ycap1 = yy - 0.55
cap1 = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Cm(2.2), Cm(ycap1), Cm(1.0), Cm(2.1))
cap1.fill.solid(); cap1.fill.fore_color.rgb = RGBColor(0xC9, 0xC9, 0xC9)
cap1.line.color.rgb = INK; cap1.line.width = Pt(1.0)
slv = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Cm(3.9), Cm(ycap1 - 0.15),
                         Cm(4.6), Cm(2.4))
slv.fill.background(); slv.line.color.rgb = INK; slv.line.width = Pt(1.75)
slv.adjustments[0] = 0.5
for k in range(8):
    seg(s, 4.15 + k * 0.52, ycap1 - 0.05, 4.15 + (k + 0.5) * 0.52, ycap1 + 2.15)
    seg(s, 4.15 + k * 0.52, ycap1 + 2.15, 4.15 + (k + 0.5) * 0.52, ycap1 - 0.05)
bld = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Cm(8.3), Cm(ycap1 + 0.15),
                         Cm(3.4), Cm(1.7))
bld.fill.solid(); bld.fill.fore_color.rgb = RGBColor(0xE3, 0xE3, 0xE3)
bld.line.color.rgb = GRAY; bld.line.width = Pt(0.75)
bld.adjustments[0] = 0.5
cap2 = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Cm(12.6), Cm(ycap1), Cm(1.0), Cm(2.1))
cap2.fill.solid(); cap2.fill.fore_color.rgb = RGBColor(0xC9, 0xC9, 0xC9)
cap2.line.color.rgb = INK; cap2.line.width = Pt(1.0)
label(s, 1.7, ycap1 + 2.35, 2.2, "end cap")
label(s, 4.6, ycap1 - 0.85, 3.2, "braided sleeve")
label(s, 8.6, ycap1 - 0.85, 3.0, "elastomer bladder")
label(s, 12.4, ycap1 + 2.35, 2.0, "end cap")
label(s, 0.2, 10.05, 15.5, "Assembly axis shown as the dash-dot centerline; components separate along it.", size=10)
prs.save(FIG + r"\70-overall\components\bpa_construction.pptx")

# ------------------------------------------- 7 SNS library block overview
prs = Presentation()
s = new_slide(prs, 16.2, 6.6)
blocks = [
    ("Graded neuron", "membrane voltage is\nthe output", 0.3, 1.3),
    ("Graded synapse", "sign from the\nconnection", 5.85, 1.3),
    ("Spiking neuron", "leaky integrate-\nand-fire", 11.4, 1.3),
    ("Ia encoder", "length and\nvelocity", 0.3, 3.9),
    ("Ib encoder", "muscle force", 5.85, 3.9),
    ("Activation\ndynamics", "muscle activation\nfilter", 11.4, 3.9),
]
for name, sub, x, y in blocks:
    box(s, x, y, 4.5, 1.75, name + "\n" + sub, bold_first=1)
box(s, 0.3, 5.55, 15.6, 0.95,
    "plus the synapse-sum junction, three diameter-specific BPA force blocks, and the biological muscle block")
prs.save(FIG + r"\60-sns\components\sns_library_overview.pptx")

print("all pptx written")
