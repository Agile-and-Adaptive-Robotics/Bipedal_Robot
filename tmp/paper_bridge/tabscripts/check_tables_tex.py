# -*- coding: utf-8 -*-
"""Static structural check of tmp/paper_bridge/sections/tables_main.tex
(the compile itself is pending: the local MiKTeX install is unfinished).
Checks: environment balance, per-row column counts vs the tabular preamble,
caption-above-table ordering, label uniqueness, \\ref closure, and that no
campaign/file artifact strings leaked in."""
import re
import sys

PATH = r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\sections\tables_main.tex"
src = open(PATH, encoding="utf-8").read()
lines = src.splitlines()
fails = []

# 1. Environment balance
for env in ("table", "tabular"):
    n_begins = len(re.findall(r"\\begin\{%s\}" % env, src))
    n_ends = len(re.findall(r"\\end\{%s\}" % env, src))
    print("env %-8s: %d begin / %d end -> %s"
          % (env, n_begins, n_ends, "OK" if n_begins == n_ends else "FAIL"))
    if n_begins != n_ends:
        fails.append(("env-imbalance", env))

# 2. Per-table structure: caption before tabular, label present, one float each
floats = re.findall(
    r"\\begin\{table\}.*?\\end\{table\}", src, flags=re.S)
print("floats found: %d" % len(floats))
if len(floats) != 5:
    fails.append(("float-count", len(floats)))

PRE_COLS = re.compile(r"\{(@\{\})?(.*?)(@\{\})?\}$")
NCOLS_EXPECTED = {  # by label
    "tab:FlxPinGof": 7, "tab:ExtPinGof": 7, "tab:AdoptedXi": 4,
    "tab:Redesign": 6, "tab:BiPulley": 4,
}

labels = re.findall(r"\\label\{(.*?)\}", src)
print("labels: " + repr(labels))
if len(labels) != len(set(labels)):
    fails.append(("duplicate-labels", labels))

# 3. Column counts: parse each tabular preamble, count &'s per body row
for f in floats:
    lab = re.search(r"\\label\{(.*?)\}", f).group(1)
    cap_i, tab_i, end_i = f.find("\\caption"), f.find("\\begin{tabular}"), f.find("\\end{tabular}")
    if not (cap_i < tab_i):
        fails.append(("caption-not-above", lab))
    body = f[f.find("\\begin{tabular}"):]
    body = body[:body.find("\\end{tabular}")]
    # preamble = first {...} after \begin{tabular}
    m = re.search(r"\\begin\{tabular\}\{(.*)\}", body)
    preamble = m.group(1)
    # count effective columns: split on | not inside {}, count l/r/c/p{..}
    cells = re.findall(r"[lrc]|p\{[^}]*\}", preamble)
    ncols = len(cells)
    if ncols != NCOLS_EXPECTED[lab]:
        fails.append(("preamble-cols", lab, ncols, NCOLS_EXPECTED[lab]))
    # rows: split body on \\\\ up to \bottomrule
    body_rows = body[m.end():]
    for row in re.split(r"\\\\", body_rows):
        row = row.strip()
        if not row or row.startswith("\\cmidrule") or row.startswith("\\toprule") \
           or row.startswith("\\midrule") or row.startswith("\\bottomrule"):
            continue
        n_amp = row.count("&")
        # multicolumn{n}{..} spans n cells but contributes 1 '&'-separated cell
        multicols = [int(x) for x in re.findall(r"\\multicolumn\{(\d+)\}", row)]
        # each &-separated field = 1 column unless multicolumn (span-1 extra each)
        expected_amp = ncols - 1 + sum(mc - 1 for mc in multicols)
        if n_amp != expected_amp:
            fails.append(("row-cols", lab, row[:60], n_amp, expected_amp))
    print("  %-14s preamble cols=%d, caption-above=%s"
          % (lab, ncols, "OK" if cap_i < tab_i else "FAIL"))

# 4. \ref closure within the fragment (assembler dedupes later)
refs = re.findall(r"\\ref\{(.*?)\}", src)
for r in set(refs):
    if r not in labels:
        fails.append(("dangling-ref", r))
print("refs: " + repr(sorted(set(refs))))

# 5. banned leakage (campaign dates, filenames, pick/trial artifacts, .mat/.csv)
banned = re.findall(
    r"202609\d\d|20260920|20260925|20260926|\.mat\b|\.csv\b|row 77|pick=|trial \d+|row77|flxr77",
    src)
print("banned-token hits: " + repr(banned))
if banned:
    fails.append(("banned-tokens", banned))

# 6. banned AI-slop / banned words per the editorial contract
slop = re.findall(r"delve|crucial role|it is important to note|moreover|furthermore|plant\b", src, re.I)
print("slop/banned-word hits: " + repr(slop))
if slop:
    fails.append(("slop", slop))

print("\n" + ("FAILURES: " + repr(fails) if fails else "ALL STATIC CHECKS PASSED"))
sys.exit(1 if fails else 0)
