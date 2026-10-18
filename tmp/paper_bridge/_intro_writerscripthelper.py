# Scratch helper for the introduction writer (tmp/paper_bridge is the sanctioned scratch area).
# Dumps key + title (+ year) for candidate citation keys from ProofFinal/BolenFrontiers22.bib
# so one-clause characterizations in the intro are grounded in the actual bibliography.
import re, os, sys

BIB = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\BolenFrontiers22.bib"

keys_wanted = [
    "chou_measurement_1996", "tondu_static_2022", "sarosi_new_2012",
    "martens_modeling_2017", "hunt_modeling_2017", "hunt_development_2017",
    "morrow_optimization_2020", "kurumaya_musculoskeletal_2016",
    "hitzmann_anthropomorphic_2018", "asano_musculoskeletal_2019",
    "asano_design_2017", "liu_robotic_2018", "ijspeert_amphibious_2020",
    "shin_constructive_2018", "schilling_adaptive_2022",
    "goldsmith_investigating_2021", "aschenbeck_design_2006",
    "veale_towards_2016", "madden_artificial_2004", "suzumori_trends_2018",
    "bolen_2026", "bolen_determination_2019", "brand_model_1982",
    "menegaldo_moment_2004", "suderman_moving_2012", "young_analyzing_2019",
    "steele_biomimetic_2018", "steele_development_2017",
    "steele_biomimetic_2018-1", "seth_opensim_2018", "delp_interactive_1990",
]

with open(BIB, encoding="utf-8", errors="replace") as f:
    src = f.read()

entries = {}
# entries look like @article{key, ... } — split on lines starting with @
for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", src):
    etype, key = m.group(1), m.group(2)
    start = m.end()
    # find the matching closing brace at depth 0
    depth = 1
    i = start
    while i < len(src) and depth > 0:
        if src[i] == "{":
            depth += 1
        elif src[i] == "}":
            depth -= 1
        i += 1
    body = src[start:i]
    tm = re.search(r"title\s*=\s*[{\"](.+?)[}\"]\s*,?\s*\n", body, re.S | re.I)
    ym = re.search(r"year\s*=\s*[{\"]?(\d{4})", body, re.I)
    title = " ".join(tm.group(1).split()) if tm else "<no title>"
    year = ym.group(1) if ym else "?"
    entries[key] = (etype, year, title)

print(f"total bib entries parsed: {len(entries)}")
for k in keys_wanted:
    if k in entries:
        et, yr, ti = entries[k]
        print(f"{k:38s} @{et:14s} {yr:5s} {ti[:95]}")
    else:
        print(f"{k:38s} *** MISSING FROM BIB ***")
