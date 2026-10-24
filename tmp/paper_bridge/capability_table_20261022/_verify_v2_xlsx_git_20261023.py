# _verify_v2_xlsx_git_20261023.py -- xlsx cell-by-cell vs md/tex, and git-state check.
import openpyxl, re, os, subprocess

BASE = os.path.dirname(os.path.abspath(__file__))
problems = []
def bad(m): problems.append(m)

# ---- xlsx ----
wb = openpyxl.load_workbook(os.path.join(BASE, "capability_table_20261023.xlsx"), data_only=True)
print("sheets:", wb.sheetnames)
for ws in wb.worksheets:
    print(f"--- sheet {ws.title!r} ({ws.max_row}x{ws.max_column}) ---")
    for row in ws.iter_rows():
        for c in row:
            if c.value is not None:
                print(f"  {c.coordinate}: {c.value!r}")

# expected headline numbers from md/tex
EXP = {"72.49", "96.64", "65.19", "74.27"}
found = set()
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value, (int, float)):
                s = f"{c.value:.2f}"
                if s in EXP:
                    found.add(s)
print("\nheadline values found in xlsx:", sorted(found))
if found != EXP:
    bad(f"xlsx headline values {sorted(found)} != expected {sorted(EXP)}")

# any percent-looking cells not in the expected set (excluding maybe per-muscle sheets)
for ws in wb.worksheets:
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value, (int, float)) and f"{c.value:.2f}" not in EXP:
                bad(f"unexpected numeric cell {ws.title}!{c.coordinate} = {c.value} (not one of {sorted(EXP)})")

# ---- md / tex cell extraction ----
md = open(os.path.join(BASE, "capability_table_20261023.md"), encoding="utf-8").read()
tex = open(os.path.join(BASE, "capability_table_20261023.tex"), encoding="utf-8").read()
md_rows = re.findall(r"^\| (Hand-tuned|Optimized)[^|]*\| ([\d.]+)% \| ([\d.]+)% \|$", md, re.M)
tex_rows = re.findall(r"^(Hand-tuned|Optimized)[^&]*& ([\d.]+)\\% & ([\d.]+)\\% \\\\$", tex, re.M)
print("md rows:", md_rows)
print("tex rows:", tex_rows)
if len(md_rows) != 2: bad(f"md table rows parsed = {len(md_rows)}, expected 2")
if len(tex_rows) != 2: bad(f"tex table rows parsed = {len(tex_rows)}, expected 2")
if md_rows != tex_rows: bad(f"md rows {md_rows} != tex rows {tex_rows}")
expected_rows = [("Hand-tuned", "72.49", "96.64"), ("Optimized", "65.19", "74.27")]
if md_rows != expected_rows: bad(f"md rows {md_rows} != expected {expected_rows}")

# labels match counts: md claims 27 of 27 graded for both rows, 27 built for optimized
for label, graded, built in (("Hand-tuned", "27 of 27 graded", None),
                             ("Optimized", "27 of 27 graded", "27 built")):
    pat = label + " (" + graded + ("; " + built if built else "") + ")"
    in_md = pat in md
    in_tex = re.search(re.escape(label) + r"[^&\\]*" + re.escape(graded) + (r"[^&\\]*" + re.escape(built) if built else "") + r"[^&\\]*&", tex) is not None
    print(f"label {pat!r}: in md = {in_md}, in tex = {in_tex}")
    if not in_md: bad(f"md missing label {pat!r}")
    if not in_tex: bad(f"tex missing label {pat!r}")

print("\n===== XLSX/MD/TEX PROBLEMS =====" if problems else "\n===== XLSX/MD/TEX: NO PROBLEMS =====")
for m in problems: print("[DISCREPANCY]", m)

# ---- git ----
print("\n===== GIT STATUS (porcelain, whole repo) =====")
r = subprocess.run(["git", "status", "--porcelain"], cwd=r"D:\GitHub\Bipedal_Robot",
                   capture_output=True, text=True, shell=False)
print(r.stdout)
print("stderr:", r.stderr.strip() or "(none)")
r2 = subprocess.run(["git", "log", "--oneline", "-3"], cwd=r"D:\GitHub\Bipedal_Robot",
                    capture_output=True, text=True, shell=False)
print("last commits:", r2.stdout.strip())
