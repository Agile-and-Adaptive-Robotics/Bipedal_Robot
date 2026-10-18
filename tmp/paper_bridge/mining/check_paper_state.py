# Cross-checks for paper_state.md (run once, read-only on the manuscript tree).
# Outputs a report to stdout; run with easteregg2 base Anaconda python.
import re, os, sys

PROOF = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal"
MAIN = os.path.join(PROOF, "Main.tex")
SUPP = os.path.join(PROOF, "frontiers_SupplementaryMaterial.tex")
BIB  = os.path.join(PROOF, "BolenFrontiers22.bib")
FIGS = os.path.join(PROOF, "Figures")

def read(p):
    with open(p, "r", encoding="utf-8", errors="replace") as f:
        return f.read()

def strip_comments_lines(text):
    """Return (active_text, n_commented_lines) — drop full-line comments only."""
    active, commented = [], 0
    for line in text.splitlines():
        ls = line.lstrip()
        if ls.startswith("%"):
            commented += 1
            continue
        active.append(line)
    return "\n".join(active), commented

cite_re = re.compile(r"\\(?:cite|citep|citet|citealt|citealp|citeauthor|citeyear)\s*(?:\[[^\]]*\]\s*){0,2}\{([^}]*)\}")
inc_re  = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\]\s*)?\{([^}]*)\}")
label_re = re.compile(r"\\label\{([^}]*)\}")
ref_re  = re.compile(r"\\(?:ref|eqref|pageref|autoref)\{([^}]*)\}")
todo_re = re.compile(r"(TODO|FIXME|XXX|\btodo\b|\\todo|fillme|TBD|PLACEHOLDER|\?\?\?)", re.IGNORECASE)

def keys_from(text):
    keys = []
    for m in cite_re.finditer(text):
        for k in m.group(1).split(","):
            k = k.strip()
            if k:
                keys.append(k)
    return keys

report = []
report.append("=== CITE KEY CHECK ===")
main_raw = read(MAIN); supp_raw = read(SUPP)
main_active, main_cmt = strip_comments_lines(main_raw)
supp_active, supp_cmt = strip_comments_lines(supp_raw)

main_keys = keys_from(main_active)
main_keys_c = keys_from(main_raw)  # includes commented
supp_keys = keys_from(supp_active)
supp_keys_c = keys_from(supp_raw)

bib_text = read(BIB)
bib_keys = re.findall(r"^@\w+\{([^,\s]+),", bib_text, re.M)
bib_set = set(bib_keys)

report.append(f"Main.tex active cite instances: {len(main_keys)}; unique: {len(set(main_keys))}")
report.append(f"Main.tex commented lines skipped: {main_cmt}")
report.append(f"Supp active cite instances: {len(supp_keys)}; unique: {len(set(supp_keys))}")
report.append(f"Bib entries: {len(bib_keys)} (unique {len(bib_set)})")

dup_bib = [k for k in set(bib_keys) if bib_keys.count(k) > 1]
report.append(f"Bib duplicate keys: {dup_bib if dup_bib else 'none'}")

missing_main = sorted(set(main_keys) - bib_set)
report.append(f"MISSING from bib (Main.tex ACTIVE cites): {missing_main if missing_main else 'NONE'}")
for k in sorted(set(main_keys)):
    report.append(f"  main-active: {k}  ({'OK' if k in bib_set else 'MISSING'}, {main_keys.count(k)}x)")

# keys only in commented-out text of Main (would matter if uncommented)
cmt_only = set(main_keys_c) - set(main_keys)
report.append(f"Main.tex cites inside comments only (currently inert): {sorted(cmt_only)}")

missing_supp = sorted(set(supp_keys) - bib_set)
report.append(f"MISSING from bib (Supp ACTIVE): {missing_supp if missing_supp else 'NONE'}")

unused_bib = sorted(bib_set - set(main_keys_c) - set(supp_keys_c))
report.append(f"Bib entries cited NOWHERE (main+supp incl. comments): {len(unused_bib)}")
for k in unused_bib:
    report.append(f"  unused-bib: {k}")

report.append("")
report.append("=== INCLUDEGRAPHICS CHECK ===")
main_figs = inc_re.findall(main_active)
supp_figs_raw = inc_re.findall(supp_raw)  # include commented for supp
supp_figs = inc_re.findall(supp_active)
all_on_disk = set(os.listdir(FIGS))
report.append(f"Figures/ on disk ({len(all_on_disk)}): {sorted(all_on_disk)}")
used = set()
for tgt in main_figs:
    used.add(tgt)
    exists = os.path.exists(os.path.join(PROOF, tgt))
    report.append(f"  main: {tgt} -> {'EXISTS' if exists else 'MISSING'}")
for tgt in supp_figs:
    used.add(tgt)
    exists = os.path.exists(os.path.join(PROOF, tgt))
    report.append(f"  supp-active: {tgt} -> {'EXISTS' if exists else 'MISSING'}")
cmt_figs = set(supp_figs_raw) - set(supp_figs)
report.append(f"supp figures referenced only in comments: {sorted(cmt_figs)}")
# unused files in Figures/
unused_files = []
for f in sorted(all_on_disk):
    base = "Figures/" + f
    hit = any(t.replace("\\", "/").endswith(f) for t in used)
    if not hit:
        unused_files.append(f)
report.append(f"Files in Figures/ NOT referenced by any ACTIVE includegraphics: {unused_files}")

report.append("")
report.append("=== REF/LABEL CHECK ===")
labels = set(label_re.findall(main_raw)) | set(label_re.findall(supp_raw))
refs_main = ref_re.findall(main_active)
refs_supp = ref_re.findall(supp_active)
report.append(f"labels (main+supp, incl comments): {sorted(labels)}")
undef_main = sorted(set(refs_main) - labels)
undef_supp = sorted(set(refs_supp) - labels)
report.append(f"Main.tex ACTIVE refs with NO label anywhere: {undef_main if undef_main else 'NONE'}")
report.append(f"Supp ACTIVE refs with NO label anywhere: {undef_supp if undef_supp else 'NONE'}")

report.append("")
report.append("=== TODO MARKERS ===")
for name, text in (("Main.tex", main_raw), ("Supp", supp_raw)):
    hits = []
    for i, line in enumerate(text.splitlines(), 1):
        if todo_re.search(line):
            hits.append((i, line.strip()[:120]))
    report.append(f"{name}: {len(hits)} hits")
    for i, l in hits:
        report.append(f"  {name}:{i}: {l}")

report.append("")
report.append("=== SECTION OUTLINE (Main.tex active) ===")
sec_re = re.compile(r"\\(section|subsection|subsubsection)\*?\{")
for i, line in enumerate(main_active.splitlines(), 1):
    m = re.match(r"\s*\\(section|subsection|subsubsection)(\*)?\{(.*)\}", line)
    if m:
        lvl = {"section": "", "subsection": "  ", "subsubsection": "    "}[m.group(1)]
        star = "*" if m.group(2) else ""
        report.append(f"  L{i}: {lvl}{m.group(1)}{star}: {m.group(3)[:80]}")

print("\n".join(report))
