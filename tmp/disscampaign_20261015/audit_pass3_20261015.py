# -*- coding: utf-8 -*-
"""Pass 3: full Output-written line, warning census, resolved-figure dump,
pre-restructure (Overleaf pull) reference check for now-orphaned moved files."""
import os, re
from pathlib import Path

BASE = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation")
PF = BASE / "ProofFinal"
OL = Path(r"D:\Github\Bipedal_Robot\tmp\disscampaign_20261015\overleaf_pull")

log = (PF / "build_wf" / "main.log").read_text(encoding="utf-8", errors="replace")
for l in log.split("\n"):
    if "Output written" in l:
        print("OUTPUT LINE:", l.strip())

warns = [l.strip() for l in log.split("\n") if "LaTeX Warning" in l or "Package .* Warning" in l]
print("\n'LaTeX Warning' lines: %d" % len(warns))
for w in warns[:20]:
    print("  ", w[:150])
cit = [l.strip() for l in log.split("\n") if "Citation" in l and "undefined" in l.lower()]
print("undefined-citation lines: %d" % len(cit))

# resolved figure names per chapter (reuse pass-2 logic quickly)
def strip_comments(text):
    out = []
    for line in text.split("\n"):
        l = line.replace("\\%", "\x00")
        i = l.find("%")
        if i >= 0:
            l = l[:i]
        out.append(l.replace("\x00", "\\%"))
    return "\n".join(out)

inc_re = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]*)\}")
main_code = strip_comments((PF / "main.tex").read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n"))
includes = re.findall(r"\\include\{([^}]*)\}", main_code)
print("\nRESOLVED \\includegraphics names per chapter:")
for inc in includes:
    fp = PF / (inc + ".tex")
    tc = strip_comments(fp.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n"))
    names = [m.group(1) for m in inc_re.finditer(tc)]
    if names:
        print("  %s:" % inc)
        for n in names:
            print("     {%s}" % n)

# pre-restate: which now-orphaned files were referenced in the Overleaf pull (pre-restructure)?
ol_tex = sorted(OL.rglob("*.tex"))
print("\nOverleaf-pull .tex files: %d" % len(ol_tex))
probes = ["CPG_Output_MN", "footMechano", "fullNetwork", "jointAngles",
          "00_GraphicalAbstract", "muscle_force_compare.csv",
          "animatlab_spiking_conversion", "circuit_core", "circuit_full",
          "route_geometry_ext20mm_20260925", "xiFrameTransform.png", "xiProjection.png",
          "xiProjection.svg", "xiProjection.pptx", "mujoco_render_provenance",
          "spiking_mirror_rhythm.png", "spiking_mirror_ground_gait.png",
          "simulink_knee_spiking_twin.png", "simulink_knee_reflex_arranged.pdf",
          "simulink_neuron_detail.pdf", "simulink_sns_library_grid.pdf",
          "bone_mesh_static_extensor20mm", "muscle_force_compare.png"]
for probe in probes:
    hits = []
    for t in ol_tex:
        txt = t.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
        if probe in txt:
            # active (non-commented)?
            active = any(not ln.lstrip().startswith("%") and probe in ln for ln in txt.split("\n"))
            hits.append((t.name, "ACTIVE" if active else "commented"))
    print("  %-38s -> %s" % (probe, hits if hits else "NOT REFERENCED pre-restructure"))
