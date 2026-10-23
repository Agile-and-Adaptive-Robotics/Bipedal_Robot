# Split full-chapter patches into discrete, independently applicable hunk patches.
import os
import re
import subprocess

GD = r"C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.7\resources\app\git\cmd\git.exe"
ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022"
PATCH = os.path.join(ROOT, "patches")

PLANS = {
    "10-introduction": [
        ("01-gaps-and-objectives", ["One lies in the literature", "Four gaps thus", "audits at 943", "Research Objectives", "resolves the afferent-data gap"]),
        ("02-project-connection-and-person", ["Steele", "Muscle Mutt", "master's project", "builds on its hardware"]),
        ("03-organization-section-fixes", ["spiking-neuron mirrors", "first-pass robot design", "robot-body model, and compares", "my publication", "It also scopes"]),
    ],
    "15-background": [
        ("01-gait2392-naming-kinematics", ["23 degrees of freedom", "planar kinematics", "musculotendon actuators", "Gait2392"]),
        ("02-masters-paragraph", ["In my master's work", "painstaking", "placement optimization algorithm", "From my master's study"]),
        ("03-compneuro-cpg-move", ["persistent sodium", "half-center", "Computational Neuro", "compneuro", "two-level"]),
        ("04-remove-graphical-abstract", ["GraphicalAbstract", "bpa_graphical"]),
        ("05-first-person-and-opener", ["Work from my own program", "my own included", "my laboratory", "my program first applied"]),
    ],
    "30-torque": [
        ("01-physical-validation-sweep", ["Physical testing of the redesigned", "adjoint-reduced", "cam route"]),
        ("02-publication-note", ["follow the manuscript", "postdate the manuscript"]),
    ],
}


def split_hunks(diff_text):
    header = diff_text[: diff_text.index("@@ -")]
    body = diff_text[diff_text.index("@@ -"):]
    hunks = re.split(r"(?=^@@ -)", body, flags=re.M)
    return header, [h for h in hunks if h.strip()]


for chapter, topics in PLANS.items():
    src = open(os.path.join(PATCH, chapter + ".patch"), encoding="utf-8").read()
    header, hunks = split_hunks(src)
    outdir = os.path.join(PATCH, chapter)
    os.makedirs(outdir, exist_ok=True)
    assigned = {name: [] for name, _ in topics}
    unassigned = []
    for h in hunks:
        changed = "".join(l for l in h.splitlines() if l[:1] in "+-")
        hits = [name for name, kws in topics if any(k in changed for k in kws)]
        if len(hits) == 1:
            assigned[hits[0]].append(h)
        elif len(hits) == 0:
            unassigned.append(h)
        else:
            # multi-topic hunk: give it to the first and warn
            assigned[hits[0]].append(h)
            print(f"WARNING {chapter}: hunk matched {hits}; -> {hits[0]}")
    if unassigned:
        print(f"WARNING {chapter}: {len(unassigned)} unassigned hunks -> 99-other")
        assigned["99-other"] = unassigned
    for name, hs in assigned.items():
        if not hs:
            print(f"NOTE {chapter}/{name}: no hunks")
            continue
        out = header + "".join(hs)
        path = os.path.join(outdir, name + ".patch")
        open(path, "w", encoding="utf-8", newline="").write(out)
        # verify it applies cleanly to a copy of the baseline
        basefile = os.path.join(ROOT, "baseline", "chapters", chapter + ".tex")
        tmpf = os.path.join(os.environ.get("TEMP", "/tmp"), "chk_" + chapter + ".tex")
        with open(basefile, "rb") as f:
            base_bytes = f.read()
        with open(tmpf, "wb") as f:
            f.write(base_bytes)
        # normalize patch paths to the temp file
        p = open(path, encoding="utf-8").read()
        p = p.replace(f'"a/baseline\\chapters\\{chapter}.tex"', f'"a/{tmpf}"').replace(
            f'"b/proposals\\{chapter}.tex"', f'"b/{tmpf}"'
        )
        tmppatch = tmpf + ".patch"
        open(tmppatch, "w", encoding="utf-8", newline="").write(p)
        r = subprocess.run([GD, "apply", "--check", tmppatch], capture_output=True, text=True)
        status = "OK" if r.returncode == 0 else "FAIL " + r.stderr.strip()[:120]
        print(f"{chapter}/{name}.patch: {len(hs)} hunks, apply-check {status}")
