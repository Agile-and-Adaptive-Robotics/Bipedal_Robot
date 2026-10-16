# -*- coding: utf-8 -*-
"""Round-2 audit verification: reproduce the anchor-occurrence counts first-hand.
Read-only. Counts substring occurrences (matching the audit's grep semantics) across
all .tex under ProofFinal (chapters/ + main.tex)."""
import os, glob

CH = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal"
files = sorted(glob.glob(os.path.join(CH, "chapters", "*.tex"))) + [os.path.join(CH, "main.tex")]
texts = {}
for p in files:
    with open(p, "r", encoding="utf-8", errors="replace") as f:
        texts[os.path.relpath(p, CH)] = f.read().replace("\r\n", "\n")

def count(anchor, live_only=False):
    total, locs = 0, []
    for name, txt in texts.items():
        if live_only:
            # strip full-line comments for the 'live text' variant
            txt = "\n".join(ln for ln in txt.split("\n") if not ln.lstrip().startswith("%"))
        c = txt.count(anchor)
        if c:
            total += c
            first = txt.find(anchor)
            line = txt[:first].count("\n") + 1
            locs.append("%s:%d (x%d)" % (name, line, c))
    return total, locs

CHECKS = [
    # (label, anchor, expected-this-round)
    ("C01 'Artificial'", "Artificial", "multi (degenerate)"),
    ("C02 'Second'", "Second", "multi (degenerate)"),
    ("C04 'begin'", "begin", "multi (degenerate substring)"),
    ("C05 'others'", "others", "2 (degenerate)"),
    ("C08 'Neuromechanical'", "Neuromechanical", "multi (degenerate)"),
    ("C09 'possible'", "possible", "multi (degenerate)"),
    ("C18 'ported'", "ported", "multi (degenerate substring)"),
    ("T02 'requirement'", "requirement", "multi (degenerate)"),
    ("T04 'these tools '", "these tools ", "2"),
    ("T05 'also '", "also ", "multi (degenerate)"),
    ("T06 ' walking'", " walking", "multi (degenerate)"),
    ("T01 surviving: sec Motivation", "\\section{Motivation}", "1"),
    ("T01 surviving: Motivation para", "Biomimetic robots provide a platform for investigating how mechanical structure and control interact", "1"),
    ("T03 'these tools were '", "these tools were ", "0"),
    ("T07 'increasing foot compliance...'", "increasing foot compliance and conserving energy", "0"),
    ("T08 concat fragment", "01_testJigs1.pdf\\includegraphics", "0"),
    ("T08 kept: bktFrame.pdf", "{bktFrame.pdf}", "1 (30-torque)"),
    ("T09 old header", "\\chapter{Results}", "0"),
    ("T09 Overleaf path", "figs/Preliminary/", "0"),
    ("T10 frag1", "Where our measured data differs significantly from the Festo model", "0"),
    ("T10 frag2", "end effects of the pressure vessel", "0"),
    ("T10 frag3", "starts to cantilever", "0"),
    ("T10 frag4", "slip off of the anterior bolt head", "0"),
    ("T11 deleted", "CPG acting as a feedforward controller", "0"),
    ("T11 kept1", "regularize the behavior of different locomotor CPG architectures", "1 (70-overall)"),
    ("T11 kept2", "Jetson-class companion computer", "1 (70-overall)"),
    ("T11 kept3", "stimulate a named neural pathway", "1 (70-overall)"),
    ("T12 deleted (live)", "predicts with the wrong trend", "0 live / 1 commented"),
    ("T13 'Fcompressive'", "Fcompressive", "0"),
    ("T13 kept", "as they exist in the project's MATLAB code", "1 (94-AppendixC:4)"),
    ("T15 'Creatingwith'", "Creatingwith", "0"),
    ("T14 hold", "Regenerating the results files updates this table.", "exactly 1"),
    ("T02 post-TC", "first improvement of the torque tools", "1 (15-background:38)"),
    ("T03/T04 post-TC", "I used these tools to convert", "1 (15-background:167)"),
    ("T05 post-TC", "animals also carry real moral cost", "1 (15-background)"),
    ("T06 post-TC", "encountered by a humanoid", "1 (40-ankle:8)"),
    ("T07 post-TC", "simplifying the system", "1 (15-background)"),
    ("C01 anchor para", "Artificial muscles are a key enabling technology", "1 (10-introduction:8)"),
    ("C02 anchor para", "Second, actuator force alone does not determine", "1 (10-introduction:22)"),
]

for label, anchor, expect in CHECKS:
    total, locs = count(anchor)
    live_total, _ = (count(anchor, live_only=True) if label.startswith("T12") else (None, None))
    extra = (" | live-only: %d" % live_total) if live_total is not None else ""
    print("%-42s -> %3d%s  [%s]  %s" % (label, total, extra, expect, "; ".join(locs) if locs else ""))
