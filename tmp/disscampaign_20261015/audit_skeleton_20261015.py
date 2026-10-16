# -*- coding: utf-8 -*-
"""Skeleton audit of the 2026-10-15 Phase-A dissertation restructure.
Checks: (1) re-anchor ledger quotes, (2) refs/labels + includegraphics resolution,
(3) orphaned figure files. Read-only: prints findings, writes nothing.
"""
import json, os, re, sys
from pathlib import Path

BASE = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation")
PF = BASE / "ProofFinal"
FIGS = BASE / "Figures"
LEDGER = BASE / "Notes" / "revision_campaign_20261015" / "comment_ledger.json"

def norm(s):
    return s.replace("\r\n", "\n").replace("\r", "\n")

def strip_comments(text):
    """Remove LaTeX comments (unescaped % to EOL), line by line."""
    out = []
    for line in text.split("\n"):
        l = line.replace("\\%", "\x00")
        i = l.find("%")
        if i >= 0:
            l = l[:i]
        out.append(l.replace("\x00", "\\%"))
    return "\n".join(out)

# ---------- gather the new tree ----------
tex_files = sorted(PF.rglob("*.tex"))
raw = {p: norm(p.read_text(encoding="utf-8", errors="replace")) for p in tex_files}
code = {p: strip_comments(t) for p, t in raw.items()}   # comment-stripped
names = {p: p.relative_to(PF).as_posix() for p in tex_files}

print("=== TREE: %d .tex files under ProofFinal ===" % len(tex_files))
for p in tex_files:
    print("  ", names[p])

# ---------- CHECK 1: re-anchor ledger quotes ----------
print("\n" + "=" * 70)
print("CHECK 1: re-anchor — every non-empty quotedHighlight vs the new tree")
print("=" * 70)
entries = json.loads(LEDGER.read_text(encoding="utf-8"))

HOME_MAP = {
    "intro": "10-introduction.tex", "background": "15-background.tex",
    "sns": "60-sns.tex", "ankle": "40-ankle.tex", "torque": "30-torque.tex",
    "muscles": "20-muscles.tex", "overall": "70-overall.tex",
    "backmatter": "80-conclusion/92-93-94 appendices",
}

results1 = []
for e in entries:
    q = e.get("quotedHighlight", "")
    if not q or not q.strip():
        continue
    qn = norm(q)
    hits = {}
    for p in tex_files:
        c = raw[p].count(qn)
        if c:
            hits[names[p]] = c
    total = sum(hits.values())
    if total == 1:
        status = "OK(1)"
    elif total == 0:
        status = "MISSING"
    else:
        status = "MULTI(%d)" % total
    results1.append((e["id"], status, hits, e.get("newHome", "")))
    line = "%-5s %-12s newHome=%-10s %s" % (e["id"], status, e.get("newHome", ""), dict(hits))
    print(line)

print("\n-- CHECK 1 summary: %d quoted entries | OK=%d MULTI=%d MISSING=%d" % (
    len(results1),
    sum(1 for r in results1 if r[1] == "OK(1)"),
    sum(1 for r in results1 if r[1].startswith("MULTI")),
    sum(1 for r in results1 if r[1] == "MISSING")))

# For MISSING multi-line quotes: per-line breakdown (prefix-trim for truncated tails)
print("\n-- MISSING quotes: per-line fragment analysis --")
for e in entries:
    q = e.get("quotedHighlight", "")
    if not q or not q.strip():
        continue
    qn = norm(q)
    if any(raw[p].count(qn) for p in tex_files):
        continue
    lines = [l.strip() for l in qn.split("\n") if l.strip()]
    print("\n[%s] %d non-empty lines; per-line hit counts:" % (e["id"], len(lines)))
    for l in lines:
        c = sum(raw[p].count(l) for p in tex_files)
        if c == 0:
            # try trimming trailing words (export truncation)
            words = l.split(" ")
            trimmed = None
            for k in range(len(words) - 1, 3, -1):
                pref = " ".join(words[:k])
                cc = sum(raw[p].count(pref) for p in tex_files)
                if cc:
                    trimmed = (k, cc)
                    break
            if trimmed:
                print("   0 hits (full line) BUT prefix of %d/%d words hits %dx: %r" %
                      (trimmed[0], len(words), trimmed[1], l[:90]))
            else:
                print("   0 hits any prefix: %r" % l[:90])
        else:
            print("   %d hit(s): %r" % (c, l[:90]))

# ---------- CHECK 2: labels / refs / includegraphics ----------
print("\n" + "=" * 70)
print("CHECK 2: \\ref/\\eqref resolution + \\includegraphics under \\graphicspath")
print("=" * 70)

labels_raw = {}
refs_raw = {}
for p in tex_files:
    for m in re.finditer(r"\\label\{([^}]*)\}", raw[p]):
        for lab in m.group(1).split(","):
            labels_raw.setdefault(lab.strip(), []).append(names[p])
    for m in re.finditer(r"\\(ref|eqref|pageref|autoref|cref|Cref)\{([^}]*)\}", raw[p]):
        for t in m.group(2).split(","):
            refs_raw.setdefault((m.group(1), t.strip()), []).append(names[p])

# comment-stripped versions (the ones that matter for the build)
labels = {}
refs = {}
for p in tex_files:
    for m in re.finditer(r"\\label\{([^}]*)\}", code[p]):
        for lab in m.group(1).split(","):
            labels.setdefault(lab.strip(), []).append(names[p] + ":code")
    for m in re.finditer(r"\\(ref|eqref|pageref|autoref|cref|Cref)\{([^}]*)\}", code[p]):
        for t in m.group(2).split(","):
            refs.setdefault((m.group(1), t.strip()), []).append(names[p])

dup_labels = {k: v for k, v in labels.items() if len(v) > 1}
print("labels (comment-stripped): %d total; duplicated: %d" % (len(labels), len(dup_labels)))
for k, v in dup_labels.items():
    print("   DUP label %r in %s" % (k, v))

unresolved = []
for (kind, t), srcs in sorted(refs.items()):
    if kind in ("ref", "eqref") and t not in labels:
        unresolved.append((kind, t, srcs))
    elif kind in ("pageref", "autoref", "cref", "Cref") and t not in labels:
        unresolved.append((kind + "(extra)", t, srcs))
print("unresolved \\ref/\\eqref targets: %d" %
      len([u for u in unresolved if not u[0].endswith("(extra)")]))
for kind, t, srcs in unresolved:
    print("   UNRESOLVED %s{%s} used in %s" % (kind, t, srcs))
# refs that exist only thanks to commented-out lines
in_comment_only = []
for (kind, t), srcs in sorted(refs_raw.items()):
    if kind in ("ref", "eqref") and t in labels_raw and t not in labels:
        in_comment_only.append((kind, t))
print("refs resolving only via commented-out \\label lines: %d %s" % (len(in_comment_only), in_comment_only))

# --- include order + graphicspath simulation ---
main_tex = raw[PF / "main.tex"]
includes = re.findall(r"\\include\{([^}]*)\}", main_tex)
print("\nmain.tex \\include order: %s" % includes)

def gp_dirs_of(text):
    """Return list of graphicspath dirs declared in a tex text (in order)."""
    dirs = []
    for m in re.finditer(r"\\graphicspath\{((?:\s*\{[^{}]*\}\s*)+)\}", text):
        for d in re.findall(r"\{([^{}]*)\}", m.group(1)):
            dirs.append(d)
    return dirs

EXTS = [".pdf", ".png", ".jpg", ".jpeg", ".eps", ".svg", ".PPTX".lower()]
gp_global = gp_dirs_of(main_tex)
print("main.tex \\graphicspath: %s" % gp_global)

inc_re = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]*)\}")
current = list(gp_global)
resolved, unresolved_g = [], []
per_file = []
for inc in includes:
    fp = PF / (inc + ".tex")
    if not fp.exists():
        print("  !! include target missing on disk: %s" % inc)
        continue
    t = raw[fp]
    own = gp_dirs_of(t)
    if own:
        current = own   # \graphicspath replaces the list
    hits = []
    for m in inc_re.finditer(t):
        name = m.group(1).strip()
        # candidates: as-given relative to compile CWD (ProofFinal), then each gp dir
        cands = []
        bases = [name] if name else []
        root, ext = os.path.splitext(name)
        if not ext:
            bases = [name + e for e in EXTS]
        elif ext.lower() == ".eps":
            bases = [name, root + "-eps-converted-to.pdf"]
        for b in bases:
            cands.append(Path("AS_GIVEN") / b)
            for d in current:
                cands.append(Path(d) / b)
        existing = []
        for c in cands:
            if c.parts[0] == "AS_GIVEN":
                test = PF / Path(*c.parts[1:])
                via = "cwd"
            else:
                test = (PF / c)
                via = "gp"
            if test.exists():
                existing.append((str(c), via, str(test)))
        hits.append((name, existing))
        if existing:
            resolved.append((inc, name, existing))
        else:
            unresolved_g.append((inc, name, [str(c) for c in cands]))
    per_file.append((inc, len(hits)))

print("\n\\includegraphics occurrences per file: %s" % per_file)
print("resolved: %d | UNRESOLVED: %d" % (len(resolved), len(unresolved_g)))
for inc, name, cands in unresolved_g:
    print("   UNRESOLVED in %s: {%s}  tried: %s" % (inc, name, cands[:8]))
# flag includes that resolve ONLY via cwd (not the chapter's declared graphicspath)
cwd_only = [(inc, n, e) for inc, n, e in resolved if not any(v == "gp" for _, v, _ in e)]
print("resolve ONLY via compile-CWD as-given path (not chapter \\graphicspath): %d" % len(cwd_only))
for inc, n, e in cwd_only:
    print("   cwd-only in %s: {%s} -> %s" % (inc, n, e[:2]))

# ---------- CHECK 3: orphans ----------
print("\n" + "=" * 70)
print("CHECK 3: orphaned files under Figures/")
print("=" * 70)
referenced = set()
for inc, name, existing in resolved:
    for _, via, abspath in existing:
        referenced.add(Path(abspath).resolve())
# also: for .eps targets, the converted artifact is used by the build chain
extra_refs = set()
for r in referenced:
    if r.suffix.lower() == ".eps":
        c = r.with_name(r.stem + "-eps-converted-to.pdf")
        if c.exists():
            extra_refs.add(c.resolve())
referenced |= extra_refs

allfigs = sorted(p for p in FIGS.rglob("*") if p.is_file())
orphans = [p for p in allfigs if p.resolve() not in referenced]
print("figure files under Figures/: %d total | referenced by build: %d | ORPHANS: %d"
      % (len(allfigs), len(referenced), len(orphans)))
dep, act = [], []
for p in orphans:
    (dep if "Deprecated" in p.parts else act).append(p)
print("-- orphans OUTSIDE Deprecated/ (active tree): %d" % len(act))
for p in act:
    print("   %s" % p.relative_to(FIGS).as_posix())
print("-- orphans INSIDE Deprecated/: %d" % len(dep))
for p in dep:
    print("   %s" % p.relative_to(FIGS).as_posix())

# cross-check: are any orphans referenced by .tex files elsewhere in Dissertation?
other_tex = [p for p in BASE.rglob("*.tex") if PF not in p.parents]
print("\ncross-check: %d .tex files outside ProofFinal under Dissertation/" % len(other_tex))
for p in act:
    rel = p.relative_to(FIGS).as_posix()
    base = p.name
    stem = p.stem
    hits = []
    for t in other_tex:
        txt = norm(t.read_text(encoding="utf-8", errors="replace"))
        if base in txt or ("Figures/" + rel) in txt or stem in txt:
            hits.append(t.relative_to(BASE).as_posix())
    if hits:
        print("   active orphan %s MENTIONED in: %s" % (rel, hits))
print("\nDONE")
