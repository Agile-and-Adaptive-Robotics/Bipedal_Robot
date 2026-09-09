"""Assemble a figure review with native application exports.

This is a review packet, not a compilation of the dissertation's LaTeX.
Run using the bundled Python with reportlab and pypdf installed.
"""
from io import BytesIO
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph

DISS = Path(__file__).resolve().parents[1]
FIG = DISS / "ProofFinal/figs/Preliminary"
OUT = DISS / "Notes/neuromechanical_figure_review.pdf"
writer = PdfWriter()
style = ParagraphStyle("caption", fontName="Helvetica", fontSize=10, leading=14,
                       textColor=colors.HexColor("#263747"))


def page(title, caption, assets, landscape=False):
    """assets = (path relative to figure directory, x, y, width, height)."""
    width, height = (792, 612) if landscape else (612, 792)
    stream = BytesIO()
    cv = canvas.Canvas(stream, pagesize=(width, height))
    cv.setTitle("Neuromechanical figures - review")
    cv.setFillColor(colors.HexColor("#203e55"))
    cv.setFont("Helvetica-Bold", 17)
    cv.drawString(45, height - 47, title)
    cv.setFillColor(colors.HexColor("#667887"))
    cv.setFont("Helvetica", 9)
    cv.drawString(45, height - 65, "Figure review | dissertation placement: Methods unless identified as Results")
    para = Paragraph(caption, style)
    _, ph = para.wrap(width - 90, 110)
    assert ph <= 96, (title, ph)
    para.drawOn(cv, 45, 55)
    cv.setFont("Helvetica", 8)
    cv.setFillColor(colors.HexColor("#667887"))
    cv.drawString(45, 26, "Layout preview; final LaTeX pagination and compilation remain to be checked in Overleaf.")
    cv.drawRightString(width - 45, 26, str(len(writer.pages) + 1))
    vector_assets = []
    for name, x, y, bw, bh in assets:
        path = FIG / name
        assert path.is_file(), path
        assert y > 55 + ph + 8, (title, name, y, ph)
        if path.suffix.lower() == ".pdf":
            vector_assets.append((path, x, y, bw, bh))
        else:
            image = ImageReader(str(path))
            iw, ih = image.getSize()
            scale = min(bw / iw, bh / ih)
            dw, dh = iw * scale, ih * scale
            cv.drawImage(image, x + (bw-dw)/2, y + (bh-dh)/2,
                         width=dw, height=dh, mask="auto")
    cv.save()
    result = PdfReader(stream).pages[0]
    for path, x, y, bw, bh in vector_assets:
        source = PdfReader(path).pages[0]
        sx, sy = float(source.mediabox.left), float(source.mediabox.bottom)
        iw, ih = float(source.mediabox.width), float(source.mediabox.height)
        scale = min(bw / iw, bh / ih)
        tx, ty = x + (bw-iw*scale)/2, y + (bh-ih*scale)/2
        transform = Transformation().translate(-sx, -sy).scale(scale).translate(tx, ty)
        result.merge_transformed_page(source, transform)
    writer.add_page(result)


page("AnimatLab: body and controller hierarchy",
     "The supplied body and subsystem illustrations show the segmented walker, its joints, "
     "and left/right neural hierarchy. These are architecture illustrations, not images from the new phase-1 trajectory.",
     [("AnimatLab_reference/OurWalker.png", 80, 420, 452, 260),
      ("AnimatLab_reference/2layerImprove.png", 45, 165, 522, 230)])

page("AnimatLab: neural circuit organization",
     "Reference schematics: rhythm generation (upper left), pattern formation (upper right), "
     "and sensorimotor pathways (below). The source labels and synaptic symbols are preserved. "
     "These show controller organization, not newly validated simulation results.",
     [("AnimatLab_reference/RGnetwork.PNG", 60, 500, 230, 186),
      ("AnimatLab_reference/PFnetwork.PNG", 320, 500, 230, 186),
      ("AnimatLab_reference/sensoryMotor_Network.PNG", 130, 155, 352, 328)])

page("AnimatLab: working-project interface",
     "Actual application captures of Biped_2xCPG_wSubs: body at initialization (top), "
     "and left-side controller subsystem (bottom). Opening and viewing the project did not run a new simulation.",
     [("animatlab_body_gui.png", 45, 430, 522, 270),
      ("animatlab_network_gui.png", 45, 145, 522, 270)])

page("MuJoCo: converted robot and BPA knee",
     "Native MuJoCo renders. Left: converted Gait2392 robot at keyframe 0 with muscle paths. "
     "Right: the independent synthetic knee at an illustrative static -35-degree pose, showing cyan flexor "
     "and gold extensor routes. Different scales. No claim of SNS control of the converted robot or validated torque.",
     [("mujoco_converted_robot.png", 45, 170, 253, 520),
      ("mujoco_synthetic_knee.png", 314, 170, 253, 520)])

if (FIG / "mujoco_robot_gui.png").is_file():
    page("MuJoCo: native viewer",
         "The converted robot displayed in MuJoCo's native viewer at its saved reference pose. "
         "This is an interface capture with no simulation stepping or SNS control of this robot.",
         [("mujoco_robot_gui.png", 45, 175, 522, 510)])

page("Simulink: SNS library and neuron",
     "Native high-resolution exports of all seven saved masked blocks and the non-spiking neuron implementation. "
     "The library is arranged for display. The BPA force block and plant parameters remain provisional; "
     "these diagrams do not establish quantitative regulation or a CAD-coupled Simscape Multibody plant.",
     [("simulink_sns_library_grid.png", 45, 350, 522, 335),
      ("simulink_neuron_detail.png", 45, 155, 522, 180)])

page("Simulink: complete knee-reflex topology",
     "Native high-resolution export of the complete reduced-order reflex model. The display was automatically arranged "
     "and every input port's source block and port were checked against the original. Source models were not saved. "
     "Mechanical and BPA force parameters remain provisional.",
     [("simulink_knee_reflex_arranged.png", 45, 140, 702, 390)], landscape=True)

page("Preliminary Results: AnimatLab joint motion",
     "Bilateral hip, knee, and ankle angles from the new phase-1 run. All 50,001 samples from 0 through 10 s "
     "are retained without filtering; nine trailing zero-padding rows were excluded. The model applied an "
     "external horizontal force during the first second. Stable walking and effective ground-contact gating are not established.",
     [("animatlab_phase1_joint_angles.pdf", 45, 165, 522, 520)])

writer.add_metadata({"/Title": "Neuromechanical figure review", "/Author": "Ben Bolen",
                     "/Subject": "Preliminary Methods and Results figures; not a compiled dissertation"})
with OUT.open("wb") as stream:
    writer.write(stream)
print(f"Wrote {len(writer.pages)} pages: {OUT}")
