# -*- coding: utf-8 -*-
"""Phase-A dissertation restructure, 2026-10-15.
Splits ProofFinal monograph chapters into advisor-approved topic chapters.
Sections move VERBATIM (exact line slices). Zero prose rewritten except the
deduplicated publication notes (authorized by the ask) and skeleton scaffolding.

Log: prints a full provenance + conservation report (captured to restructure_log.txt).
"""
import os, re, shutil, sys
from collections import Counter

ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
PF   = os.path.join(ROOT, "ProofFinal")
CH   = os.path.join(PF, "chapters")
FIG  = os.path.join(ROOT, "Figures")

def rd(name):
    with open(os.path.join(CH, name), "r", encoding="utf-8", newline="") as f:
        return f.read().splitlines(keepends=True)

def wr(path, lines):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("".join(lines))

report = []
def say(s=""):
    report.append(s)
    print(s)

say("HOSTNAME: %s" % os.environ.get("COMPUTERNAME", "?"))
say("Python: %s" % sys.version.split()[0])
say("")

# ---------------------------------------------------------------- donors
M  = rd("20-methods.tex")    # 20-methods
R  = rd("30-results.tex")    # 30-results
D  = rd("40-discussion.tex") # 40-discussion
F  = rd("50-futurework.tex") # 50-futurework
C  = rd("60-conclusion.tex") # 60-conclusion
B  = rd("15-background.tex") # background (gains prior_walker)

def sl(lines, a, b):
    """1-based inclusive slice."""
    return lines[a-1:b]

say("donor lengths (lines): methods=%d results=%d discussion=%d futurework=%d conclusion=%d background=%d"
    % (len(M), len(R), len(D), len(F), len(C), len(B)))

# ---------------------------------------------------------------- publication notes (deduplicated)
NOTE_MUSCLES = (
    "{\\doublespacing\\normalsize\\upshape\n"
    "\\noindent Publication note: Sections~\\ref{sec:force_model}--\\ref{sec:force_experiment}, "
    "Section~\\ref{sec:results_force}, and Sections~\\ref{sec:disc_effects}--\\ref{sec:disc_bio}, "
    "with their associated figures and tables, appeared in \\citet{bolen_2026} (Bolen, Elzein, Pang, "
    "and Hunt, \\emph{Actuators}, 2026). The knee-torque sections of Chapter~\\ref{ch:torque} are also "
    "the subject of an unpublished manuscript in preparation by Bolen, Burgess, Morrow, Elzein, Pang, "
    "Scharzenberger, and Hunt \\citep{bolen_isometric_2025}; that manuscript has not been submitted "
    "for publication. For both works, B.B.\\ performed the conceptualization, methodology, software, "
    "validation, formal analysis, investigation, data curation, visualization, and original draft "
    "preparation; M.E., L.P., L.B., and C.M.\\ contributed methodology, software, investigation, and "
    "data curation; and A.H.\\ supervised the work, administered the projects, and acquired funding. "
    "All co-authors contributed to writing, review, and editing.\\par}\n"
)
NOTE_TORQUE = (
    "{\\doublespacing\\normalsize\\upshape\n"
    "\\noindent Publication note: Sections~\\ref{sec:robot_model}--\\ref{sec:improved_model}, "
    "Section~\\ref{sec:results_torque}, and Section~\\ref{sec:disc_error} are also the subject of an "
    "unpublished manuscript in preparation by Bolen, Burgess, Morrow, Elzein, Pang, Scharzenberger, "
    "and Hunt \\citep{bolen_isometric_2025}. This manuscript has not been submitted for publication. "
    "Author contributions for both works are stated in Chapter~\\ref{ch:muscles}.\\par}\n"
)

PROV = "% [ZCODE 2026-10-15 Phase-A restructure: sections moved verbatim from the former monograph chapters; provenance in Notes/revision_campaign_20261015/migration_map.md]\n"

def block(*chunks):
    out = []
    for c in chunks:
        out.extend(c)
        out.append("\n")
    return out

# ---------------------------------------------------------------- 20-muscles.tex
muscles = (
    ["\\chapter{Artificial Muscles}\\label{ch:muscles}\\label{ch:methods}\\label{ch:discussion}\n",
     "\\graphicspath{{../Figures/20-muscles/}}\n", "\n", PROV, "\n",
     NOTE_MUSCLES, "\n"]
    + block(sl(M, 13, 68),  sl(M, 70, 87),  sl(R, 11, 230),
            sl(D, 8, 19),   sl(D, 22, 30),  sl(D, 32, 93),
            sl(D, 96, 101), sl(D, 104, 112))
)
wr(os.path.join(CH, "20-muscles.tex"), muscles)

# ---------------------------------------------------------------- 30-torque.tex
torque = (
    ["\\chapter{Isometric Knee Torque}\\label{ch:torque}\\label{ch:results}\n",
     "\\graphicspath{{../Figures/30-torque/}}\n", "\n", PROV, "\n",
     NOTE_TORQUE, "\n"]
    + block(sl(M, 91, 125),  sl(M, 128, 137), sl(M, 140, 155),
            sl(M, 159, 367), sl(M, 369, 385),  sl(R, 233, 330),
            sl(R, 332, 376), sl(D, 117, 135),  sl(D, 162, 169))
)
wr(os.path.join(CH, "30-torque.tex"), torque)

# ---------------------------------------------------------------- 40-ankle.tex
ankle = (
    ["\\chapter{Inverted-Pendulum Balance Platform}\\label{ch:ankle}\n",
     "\\graphicspath{{../Figures/40-ankle/}}\n", "\n", PROV, "\n"]
    + block(sl(M, 387, 404), sl(R, 378, 384))
)
wr(os.path.join(CH, "40-ankle.tex"), ankle)

# ---------------------------------------------------------------- 50-database.tex (skeleton)
database = [
    "\\chapter{The Sensory Afferent Database}\\label{ch:database}\n",
    "\\graphicspath{{../Figures/50-database/}}\n", "\n",
    "% [ZCODE 2026-10-15 Phase-A restructure: SKELETON ONLY -- \\chapter + \\graphicspath + section\n",
    "% headings; Phase B writers fill the content. Existing SADb prose lives in\n",
    "% Section~\\ref{sec:bg_sns} (background, unchanged this phase) and in the foundation summary\n",
    "% of Chapter~\\ref{ch:overall}; Phase B expands that material here.]\n", "\n",
    "\\section{A Structured, Model-Oriented Afferent Resource}\n", "\n",
    "% Phase B: content.\n", "\n",
    "\\section{Audited Corpus and Curation Coverage}\n", "\n",
    "% Phase B: content (audited corpus, curation notes, afferent classifications).\n", "\n",
    "\\section{AI-Assisted Literature Triage}\n", "\n",
    "% Phase B: content.\n", "\n",
    "\\section{Toward Robot and Simulation Use}\n", "\n",
    "% Phase B: content (pathway-translation notes, citation graph, explorer).\n", "\n",
]
wr(os.path.join(CH, "50-database.tex"), database)

# ---------------------------------------------------------------- 60-sns.tex
sns = (
    ["\\chapter{Synthetic Nervous System Simulation}\\label{ch:sns}\n",
     "\\graphicspath{{../Figures/60-sns/}{../Figures/60-sns/AnimatLab_reference/}}\n", "\n", PROV, "\n"]
    + block(sl(M, 414, len(M)), sl(R, 386, 473), sl(R, 475, len(R)), sl(D, 146, 160))
)
wr(os.path.join(CH, "60-sns.tex"), sns)

# ---------------------------------------------------------------- 70-overall.tex (skeleton + moved)
overall = (
    ["\\chapter{Overall System Integration and Future Work}\\label{ch:overall}\\label{ch:futurework}\n",
     "\\graphicspath{{../Figures/70-overall/}}\n", "\n",
     "% [ZCODE 2026-10-15 Phase-A restructure: skeleton + verbatim-moved sections (methods Overview,\n",
     "% foundation studies, synthesis, and the former Future Work chapter). Phase B writers\n",
     "% reorganize and fill; nothing was rewritten this phase.]\n", "\n"]
    + block(sl(M, 7, 11), sl(R, 7, 9), sl(D, 138, 144), sl(F, 4, len(F)))
)
wr(os.path.join(CH, "70-overall.tex"), overall)

# ---------------------------------------------------------------- 80-conclusion.tex
concl = (
    ["\\chapter{Conclusion}\\label{ch:conclusion}\n",
     "\\graphicspath{{../Figures/80-conclusion/}}\n", "\n", PROV, "\n"]
    + block(sl(C, 4, len(C)))
)
wr(os.path.join(CH, "80-conclusion.tex"), concl)

# ---------------------------------------------------------------- 15-background gains prior_walker
prior_walker = sl(M, 406, 412)
bg_new = B[:150] + ["\n"] + prior_walker + ["\n"] + B[150:]
wr(os.path.join(CH, "15-background.tex"), bg_new)

# ---------------------------------------------------------------- conservation check
# Donor lines dropped on purpose: chapter headers, graphicspath, pub notes, the two
# "%%%%%%%%%%" divider comment lines, and 40-discussion blank separators.
dropped = set()
dropped |= {( "20-methods", i) for i in [1, 2, 4, 5, 89]}
dropped |= {( "30-results", i) for i in [1, 2, 4, 5, 231]}
dropped |= {("40-discussion", i) for i in [1, 2, 4, 5, 6, 7]}
dropped |= {("50-futurework", i) for i in [1, 2, 3]}
dropped |= {("60-conclusion", i) for i in [1, 2, 3]}

donor_pool = Counter()
for name, lines in [("20-methods", M), ("30-results", R), ("40-discussion", D),
                    ("50-futurework", F), ("60-conclusion", C)]:
    for i, ln in enumerate(lines, 1):
        if ("dummy", i) is None:
            pass
        if (name, i) in dropped:
            continue
        s = ln.strip()
        if s:
            donor_pool[s] += 1

target_pool = Counter()
new_names = ["20-muscles.tex", "30-torque.tex", "40-ankle.tex", "50-database.tex",
             "60-sns.tex", "70-overall.tex", "80-conclusion.tex"]
for name in new_names:
    for ln in rd(name):
        s = ln.strip()
        if s:
            target_pool[s] += 1
# 15-background: new file = old content + prior_walker lines (already counted in donor pool once;
# the copy in bg_new is the same text). Count bg_new only for the ADDED copy? The prior_walker text
# exists exactly once in donors and once in targets (inside 15-background) -- no double counting
# needed because we did NOT put prior_walker in any other new file.

missing = donor_pool - target_pool   # donor content that failed to land anywhere
added   = target_pool - donor_pool   # lines that are new (skeleton, notes, headers)

say("")
say("=== CONSERVATION CHECK (non-blank line multiset) ===")
say("donor content lines: %d ; target lines: %d" % (sum(donor_pool.values()), sum(target_pool.values())))
say("MISSING (donor lines that landed nowhere): %d" % sum(missing.values()))
for ln in sorted(missing.elements()):
    say("  MISSING> " + ln[:160])
say("ADDED (new skeleton/header/note lines): %d" % sum(added.values()))

# ---------------------------------------------------------------- label conservation
lab = re.compile(r"\\label\{([^}]+)\}")
def labels_of(lines):
    return Counter(lab.search(ln).group(1) for ln in lines if lab.search(ln))
donor_labels = labels_of(M) + labels_of(R) + labels_of(D) + labels_of(F) + labels_of(C) + labels_of(B)
target_labels = Counter()
for name in new_names:
    target_labels += labels_of(rd(name))
target_labels += labels_of(bg_new)

say("")
say("=== LABEL CHECK ===")
lab_missing = donor_labels - target_labels
lab_added   = target_labels - donor_labels
say("labels lost: %d %s" % (sum(lab_missing.values()), sorted(lab_missing.elements())))
say("labels added: %d %s" % (sum(lab_added.values()), sorted(lab_added.elements())))

# ---------------------------------------------------------------- delete emptied donors
for name in ["20-methods.tex", "30-results.tex", "40-discussion.tex",
             "50-futurework.tex", "60-conclusion.tex"]:
    os.remove(os.path.join(CH, name))
    say("deleted chapters/%s" % name)

# ---------------------------------------------------------------- main.tex include order
mp = os.path.join(PF, "main.tex")
main = open(mp, "r", encoding="utf-8", newline="").read()
old_inc = ("\\include{chapters/20-methods}\n\\include{chapters/30-results}\n"
           "\\include{chapters/40-discussion}\n\\include{chapters/50-futurework}\n"
           "\\include{chapters/60-conclusion}\n")
new_inc = ("\\include{chapters/20-muscles}\n\\include{chapters/30-torque}\n"
           "\\include{chapters/40-ankle}\n\\include{chapters/50-database}\n"
           "\\include{chapters/60-sns}\n\\include{chapters/70-overall}\n"
           "\\include{chapters/80-conclusion}\n")
assert old_inc in main, "main.tex include block not found"
main = main.replace(old_inc, new_inc)
old_cmt = ("% Monograph structure: Introduction, Background, Materials and Methods,\n"
           "% Results, Discussion, Future Work, Conclusion, Bibliography")
new_cmt = ("% Topic-chapter structure (advisor-approved, 2026-10-15): Introduction, Background,\n"
           "% Artificial Muscles, Isometric Knee Torque, Inverted-Pendulum Balance Platform,\n"
           "% Sensory Afferent Database, Synthetic Nervous System Simulation, Overall System\n"
           "% Integration and Future Work, Conclusion, Bibliography")
assert old_cmt in main, "main.tex structure comment not found"
main = main.replace(old_cmt, new_cmt)
wr(mp, main.splitlines(keepends=True))
say("main.tex include order updated")

# ---------------------------------------------------------------- figure moves
def mv(src_rel, dst_dir, name=None):
    src = os.path.join(FIG, src_rel)
    dst = os.path.join(FIG, dst_dir, name or os.path.basename(src))
    assert os.path.isfile(src), "missing figure: " + src
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.move(src, dst)
    say("fig-move: Figures/%s -> Figures/%s" % (src_rel.replace("\\", "/"), os.path.relpath(dst, FIG).replace("\\", "/")))

say("")
say("=== FIGURE MOVES ===")
# muscles figures
for f in ["02a_MaxForce10.eps", "02a_MaxForce10-eps-converted-to.pdf",
          "02b_MaxForce20.eps", "02b_MaxForce20-eps-converted-to.pdf",
          "03a_MaxForce10_all.eps", "03a_MaxForce10_all-eps-converted-to.pdf",
          "03b_kmax10.eps", "03b_kmax10-eps-converted-to.pdf",
          "04a_FStar10.eps", "04a_FStar10-eps-converted-to.pdf",
          "04b_FStar20.eps", "04b_FStar20-eps-converted-to.pdf"]:
    mv("30-results/" + f, "20-muscles")
for f in ["05_Comparison.eps", "05_Comparison-eps-converted-to.pdf"]:
    mv("40-discussion/" + f, "20-muscles")
mv("20-methods/01_testJigs1.pdf", "20-muscles")
for f in ["02_Fmax_fit.eps", "03_F-P_contraction.eps", "04_GeneralizedForceFit.eps"]:
    mv("30-results/Deprecated/" + f, "20-muscles/Deprecated")

# torque figures
for f in ["testJigs_reviewed.pdf", "KneeICR.eps", "KneeICR-eps-converted-to.pdf",
          "bktFrame.pdf", "xiFrameTransform.pdf", "xiFrameTransform.png",
          "xiProjection.pdf", "xiProjection.png", "xiProjection.svg", "xiProjection.pptx"]:
    mv("20-methods/" + f, "30-torque")
for f in ["FlxPin_group.eps", "FlxPin_group-eps-converted-to.pdf",
          "ExtPin_group.eps", "ExtPin_group-eps-converted-to.pdf",
          "Ext10mm_52cm.eps", "Ext10mm_52cm-eps-converted-to.pdf",
          "FlxGrp.eps", "FlxGrp-eps-converted-to.pdf", "KneeExt_10mm_passthru.png"]:
    mv("30-results/" + f, "30-torque")
for f in ["testJigs.pdf", "testJigs_alternate1.pdf"]:
    mv("20-methods/Deprecated/" + f, "30-torque/Deprecated")
mv("30-results/Deprecated/paired_coefficient_tables_2026-09-29.jpg", "30-torque/Deprecated")

# ankle figures
mv("20-methods/inverted_pendulum_photo.jpg", "40-ankle")

# sns figures
for f in ["animatlab_body_gui.png", "animatlab_network_gui.png",
          "mujoco_converted_robot.png", "mujoco_synthetic_knee.png",
          "simulink_knee_reflex_arranged.png", "simulink_knee_reflex_arranged.pdf",
          "simulink_neuron_detail.png", "simulink_neuron_detail.pdf",
          "simulink_sns_library_grid.png", "simulink_sns_library_grid.pdf",
          "mujoco_render_provenance.json"]:
    mv("20-methods/" + f, "60-sns")
shutil.move(os.path.join(FIG, "20-methods", "AnimatLab_reference"),
            os.path.join(FIG, "60-sns", "AnimatLab_reference"))
say("fig-move: Figures/20-methods/AnimatLab_reference -> Figures/60-sns/AnimatLab_reference (whole subfolder)")
for f in ["simulink_knee_reflex.pdf", "simulink_knee_reflex.png",
          "simulink_sns_library.pdf", "simulink_sns_library.png"]:
    mv("20-methods/Deprecated/" + f, "60-sns/Deprecated")
mv("30-results/Deprecated/animatlab_phase1_preliminary.pdf", "60-sns/Deprecated")
for f in ["animatlab_phase1_joint_angles.pdf", "sns_w2l_traces.png",
          "spiking_mirror_rhythm.pdf", "spiking_mirror_rhythm.png", "spiking_mirror_rhythm_alt.txt",
          "spiking_mirror_ground_gait.pdf", "spiking_mirror_ground_gait.png",
          "spiking_mirror_ground_gait_alt.txt",
          "simulink_knee_spiking_twin.pdf", "simulink_knee_spiking_twin.png",
          "simulink_knee_spiking_twin_alt.txt"]:
    mv("30-results/" + f, "60-sns")
shutil.move(os.path.join(FIG, "30-results", "CPG_airstepping_figs"),
            os.path.join(FIG, "60-sns", "CPG_airstepping_figs"))
say("fig-move: Figures/30-results/CPG_airstepping_figs -> Figures/60-sns/CPG_airstepping_figs (whole subfolder, unreferenced by ProofFinal)")

# empty target folders for skeleton chapters
for d in ["50-database", "70-overall", "80-conclusion"]:
    os.makedirs(os.path.join(FIG, d), exist_ok=True)
    say("fig-dir created: Figures/%s (empty; skeleton chapter)" % d)

# remove emptied old figure folders
for d in ["20-methods", "30-results", "40-discussion", "50-futurework", "60-conclusion"]:
    p = os.path.join(FIG, d)
    leftover = [os.path.join(r, f) for r, _, fs in os.walk(p) for f in fs] if os.path.isdir(p) else []
    if leftover:
        say("NOT removed (still has files): Figures/%s -> %s" % (d, leftover))
    else:
        shutil.rmtree(p)
        say("removed emptied folder: Figures/%s" % d)

say("")
say("=== DONE ===")
with open(os.path.join(ROOT, "Notes", "revision_campaign_20261015", "restructure_log.txt"), "w",
          encoding="utf-8", newline="\n") as f:
    f.write("\n".join(report) + "\n")
print("log written to Notes/revision_campaign_20261015/restructure_log.txt")
