"""Prepare Chapter 6 workflow and copy verified portable vector evidence.

Run from any directory with Python (matplotlib, pypdf, pypdfium2).
No experimental recordings or numerical results are changed.
"""
from pathlib import Path
import hashlib
import json
import shutil
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject
import pypdfium2 as pdfium

ROOT = Path(__file__).resolve().parents[6]
DISS = ROOT / 'Documentation/Reports and Papers/Dissertation'
HANDOFF = DISS / 'Defense/Session1_Handoff_20261104'
OUT = DISS / 'Figures/60-sns'
REPORT = DISS / 'Notes/ch6_revision_20261116'
PREVIEW = ROOT / 'tmp/ch6_20261116/previews'
OUT.mkdir(parents=True, exist_ok=True)
PREVIEW.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({'font.family': 'Arial', 'font.size': 10.4,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(6.0, 6.9))
fig.subplots_adjust(0.02, 0.02, 0.98, 0.98)
ax.set(xlim=(0, 1), ylim=(0, 1))
ax.axis('off')
ink, blue, green, grey = '#202A31', '#4477AA', '#228833', '#667078'

def box(x, y, w, h, text, fill='#EAF1F7', edge=blue, size=10.4):
    patch = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.008,rounding_size=0.012',
                          facecolor=fill, edgecolor=edge, linewidth=1.2)
    ax.add_patch(patch)
    ax.text(x+w/2, y+h/2, text, ha='center', va='center', color=ink,
            fontsize=size, linespacing=1.35)

def arrow(points, color=ink, style='-'):
    for a, b in zip(points[:-2], points[1:-1]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=1.25, ls=style)
    ax.add_patch(FancyArrowPatch(points[-2], points[-1], arrowstyle='-|>',
                                mutation_scale=10, color=color, lw=1.25,
                                linestyle=style))

ax.text(0.02, .975, 'PREPARE AND CHECK THE MECHANICS', fontsize=10.8,
        fontweight='bold', va='top', color=ink)
box(.025, .795, .445, .115, 'OpenSim\nConversion and repair\nMuJoCo model')
box(.53, .795, .445, .115, 'SOLIDWORKS CAD\nImport or separate assembly\nSimscape / MuJoCo model')

ax.text(.02, .689, 'DEFINE AND EXECUTE THE COUPLED LOOP', fontsize=10.8,
        fontweight='bold', va='top', color=ink)
box(.025, .485, .265, .13, 'SNS controller\nGraphical blocks\nor scripted network')
box(.365, .485, .265, .13, 'Actuator interface\nActivation / pressure\nForce calculation')
box(.705, .485, .265, .13, 'Mechanical model\nSimscape or\nMuJoCo')
arrow([(.29, .55), (.365, .55)])
arrow([(.63, .55), (.705, .55)])
ax.text(.325, .635, 'command', fontsize=9.7, ha='center', color=grey)
ax.text(.675, .635, 'tension', fontsize=9.7, ha='center', color=grey)
arrow([(.247, .785), (.247, .72), (.84, .72), (.84, .625)], color=blue)
arrow([(.75, .785), (.75, .72)], color=blue)
box(.32, .285, .37, .105, 'Sensory encoding\nLength, velocity, force, contact',
    fill='#EAF5EE', edge=green, size=10.2)
arrow([(.84, .48), (.84, .337), (.70, .337)], color=green)
arrow([(.31, .337), (.16, .337), (.16, .475)], color=green)

ax.text(.02, .225, 'CALIBRATE, TEST, AND INSPECT', fontsize=10.8,
        fontweight='bold', va='top', color=ink)
box(.025, .035, .945, .13,
    'Functional targets and parameter calibration\nPhysical curriculum, interventions, and fresh replays\nAligned neural and mechanical recordings',
    fill='#F4F4F0', edge=grey)
fig.savefig(OUT / 'simulation_workflow_20261116.pdf')
fig.savefig(OUT / 'simulation_workflow_20261116.png', dpi=220)
plt.close(fig)
(OUT / 'simulation_workflow_20261116.alt.txt').write_text(
    'Two checked mechanical preparation routes feed a coupled SNS controller, '
    'actuator-force interface, and mechanical model. Mechanical measurements return '
    'through sensory encoders. Functional calibration, physical-stage tests, '
    'interventions, and replay provide the experimental protocol. The diagram '
    'does not assert automatic Simscape-to-MuJoCo translation.', encoding='utf-8')

items = [
    (HANDOFF / 'experiments_20261104/Gait2392_BiW2LvarX_pose_stage2_verification.pdf',
     OUT / 'biw2lvar_air_stage2_20261116.pdf',
     'Recorded afferented air-stage Gait2392 experiment with BiW2Lvar+; aligned '
     'neural and mechanical panels, suspension active, ground disabled and foot '
     'forces zero. Compatibility selector hunt does not attribute the architecture.'),
    (HANDOFF / 'media/CAD_reference_tension_77p5pctBW.png',
     OUT / 'cad_supported_20261116.png',
     'Sagittal and oblique recorded views of the guided CAD reference-tension '
     'model at 22.9 seconds in the selected 28-second run. Final 16–28 seconds use '
     '77.5 percent body-weight vertical assistance, with lateral and angular '
     'guidance retained. SNS timing plus engineering force allocation; not '
     'dynamically calibrated BPA actuation.'),
]
manifest = []
for src, dst, alt in items:
    shutil.copy2(src, dst)
    dst.with_suffix('.alt.txt').write_text(alt, encoding='utf-8')
    manifest.append({'source': str(src.relative_to(ROOT)),
                     'destination': str(dst.relative_to(ROOT)),
                     'source_sha256': hashlib.sha256(src.read_bytes()).hexdigest(),
                     'destination_sha256': hashlib.sha256(dst.read_bytes()).hexdigest(),
                     'pages': len(PdfReader(dst).pages) if dst.suffix=='.pdf' else None})

# Keep the complete two-page source unchanged in the handoff. Rearrange exact
# vector strips for the main chapter without digitizing/recalculating traces.
src = items[0][0]
reader = PdfReader(src)
strips = [(0, (0, 544, 470, 662)), (0, (0, 166, 470, 282)),
          (1, (0, 411, 470, 536)), (1, (0, 0, 470, 157))]
height = sum(b[3]-b[1] for _, b in strips) + 10*(len(strips)-1)
writer = PdfWriter()
page = writer.add_blank_page(width=470, height=height)
cursor = height
for index, rect in strips:
    clipped = copy.deepcopy(reader.pages[index])
    clipped.cropbox = RectangleObject(rect)
    clipped.mediabox = RectangleObject(rect)
    cursor -= rect[3]-rect[1]
    page.merge_transformed_page(clipped,
        Transformation().translate(tx=-rect[0], ty=cursor-rect[1]))
    cursor -= 10
writer.write(items[0][1])
manifest[0].update({'destination_sha256': hashlib.sha256(items[0][1].read_bytes()).hexdigest(),
                    'pages': 1, 'operation': 'Exact vector strip rearrangement; no trace recalculation',
                    'source_page_and_cropbox_points': strips})
for path in [OUT / 'biw2lvar_air_stage2_20261116.pdf']:
    doc = pdfium.PdfDocument(str(path))
    for n in range(len(doc)):
        doc[n].render(scale=1.8).to_pil().save(PREVIEW / f'{path.stem}_{n+1}.png')
    doc.close()
(REPORT / 'CHAPTER_ASSET_MANIFEST.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print(json.dumps(manifest, indent=2))
