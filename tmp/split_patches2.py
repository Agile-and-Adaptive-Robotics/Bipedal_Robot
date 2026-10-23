# Discrete sub-patches: generate -U0 (atomic) diffs, classify hunks by topic,
# verify each applies cleanly to the unmodified baseline.
import os
import re
import shutil
import subprocess

GD = r"C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.7\resources\app\git\cmd\git.exe"
ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022"
TMP = os.path.join(os.environ["TEMP"], "patchchk")
PATCH = os.path.join(ROOT, "patches")

PLANS = {
    "10-introduction": [
        ("01-gaps-and-objectives", ["One lies in the literature", "Four gaps thus", "audits at 943", "resolves the afferent-data gap", "Research Objectives"]),
        ("02-project-connection-and-person", ["Steele", "Muscle Mutt", "master's project", "hardware and demonstrations"]),
        ("03-organization-section-fixes", ["spiking-neuron mirrors", "first-pass robot design", "robot-body model, and compares", "my publication", "five-stage training curriculum"]),
    ],
    "15-background": [
        ("01-gait2392-naming-kinematics", ["23 degrees of freedom", "planar kinematics", "musculotendon actuators"]),
        ("02-masters-paragraph", ["In my master's work", "From my master's study", "placement optimization algorithm", "moment arm is the quantity", "Steele's biomimetic knee", "overview of the mechanics"]),
        ("03-compneuro-cpg-move", ["persistent sodium", "half-centers", "compneuro", "two-level", "Computational Neuro"]),
        ("04-remove-graphical-abstract", ["GraphicalAbstract", "bpa_graphical"]),
        ("05-first-person-and-opener", ["Work from my own program", "my own included", "my laboratory", "my program first applied"]),
    ],
    "30-torque": [
        ("01-physical-validation-sweep", ["Physical testing of the redesigned", "adjoint-reduced", "hung loads", "cam route"]),
        ("02-publication-note", ["follow the manuscript", "postdate"]),
    ],
}

KEYWORDS_ANY = {  # secondary signal for hunks that only delete baseline text
    "10-introduction": ["There are additional gaps", "hundreds if not thousands", "studies by the author"],
    "15-background": ["Gait2392 model", "the author's master's", "23 corresponds", "first-class", "synaptic saturation"],
    "30-torque": ["under way with", "direct route", "swept"],
}


def classify(changed, topics, chapter):
    hits = [name for name, kws in topics if any(k in changed for k in kws)]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        return "MIXED:" + ",".join(hits)
    return None


for chapter, topics in PLANS.items():
    # fresh temp check tree
    if os.path.isdir(TMP):
        shutil.rmtree(TMP)
    os.makedirs(os.path.join(TMP, "baseline", "chapters"))
    basefile = os.path.join(ROOT, "baseline", "chapters", chapter + ".tex")
    shutil.copy(basefile, os.path.join(TMP, "baseline", "chapters", chapter + ".tex"))

    r = subprocess.run(
        [GD, "diff", "--no-index", "-U0", f"baseline/chapters/{chapter}.tex", f"proposals/{chapter}.tex"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    diff = r.stdout
    if "@@" not in diff:
        print(chapter, "NO HUNKS?!")
        continue
    header = diff[: diff.index("@@ -")]
    body = diff[diff.index("@@ -"):]
    hunks = [h for h in re.split(r"(?=^@@ -)", body, flags=re.M) if h.strip()]

    outdir = os.path.join(PATCH, chapter)
    os.makedirs(outdir, exist_ok=True)
    buckets = {name: [] for name, _ in topics}
    problems = []
    for h in hunks:
        changed = "".join(l for l in h.splitlines() if l[:1] in "+-")
        t = classify(changed, topics, chapter)
        if t and t.startswith("MIXED:"):
            problems.append((t, h[:80]))
            buckets[t.split(":")[1].split(",")[0]].append(h)
        elif t:
            buckets[t].append(h)
        else:
            # deletion-only hunk: try matching against baseline keywords
            alt = None
            for name, kws in topics:
                if any(k in changed for k in KEYWORDS_ANY.get(chapter, [])) and any(
                    k in changed for k in kws
                ):
                    alt = name
                    break
            if alt:
                buckets[alt].append(h)
            else:
                problems.append(("UNASSIGNED", changed[:100]))

    for name, hs in buckets.items():
        if not hs:
            print(f"NOTE {chapter}/{name}: no hunks")
            continue
        # same-file header so `git apply` patches the chapter in place
        same_header = (
            f"diff --git a/chapters/{chapter}.tex b/chapters/{chapter}.tex\n"
            f"--- a/chapters/{chapter}.tex\n"
            f"+++ b/chapters/{chapter}.tex\n"
        )
        out = same_header + "".join(hs)
        path = os.path.join(outdir, name + ".patch")
        open(path, "w", encoding="utf-8", newline="").write(out)
        os.makedirs(os.path.join(TMP, "chapters"), exist_ok=True)
        shutil.copy(basefile, os.path.join(TMP, "chapters", chapter + ".tex"))
        chk = subprocess.run(
            [GD, "apply", "--check", "--unidiff-zero", "--ignore-whitespace", path], cwd=TMP, capture_output=True, text=True
        )
        status = "OK" if chk.returncode == 0 else "FAIL " + chk.stderr.strip()[:150]
        print(f"{chapter}/{name}.patch: {len(hs)} hunks, apply-check {status}")
    for t, info in problems:
        print(f"  {chapter} ISSUE {t}: {info}")
