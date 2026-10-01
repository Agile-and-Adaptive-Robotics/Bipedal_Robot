"""Probe: per-tag census of the aff template vs the hardcoded expectations
in W2LAffNet._assert_complete (build_w2l_aff_net.py). Run with the myo env."""
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = Path(r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco")
sys.path.insert(0, str(HERE))
import build_w2l_split_net as SPLIT  # noqa: E402
import build_w2l_aff_net as AFF      # noqa: E402

TPL = Path(r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\tmp\li_fetch"
           r"\connectome_templates_ee2.json")

tpl = json.loads(TPL.read_text(encoding="utf-8"))
split = SPLIT.make_split_template(tpl["w2laproj"], comm=1.0)
aff = AFF.make_aff_template(split, heel=True, toe=True, ib=True)
n_types = {n["label"]: n["type"] for n in aff["nodes"]}
expected = {}
for e in aff["edges"]:
    if e["tag"] == "mn_to_muscle":
        continue
    if n_types[e["from"]] == "PORT-load":
        continue
    expected[e["tag"]] = expected.get(e["tag"], 0) + 1
print("aff census total:", sum(expected.values()))
for k in sorted(expected):
    print(f"  {k}: {expected[k]}")
want = {"heel_rge_in": 4, "heel_pf": 4, "toe_in": 2, "toe_df": 2,
        "toe_df_mn": 2, "ib_load": 8, "pf_drive": 8, "rg_laminate": 12,
        "comm_c1": 2, "c1_SynAmp2.749": 2}
print("mismatches vs _assert_complete expectations:")
for k, v in want.items():
    got = expected.get(k, 0)
    if got != v:
        print(f"  {k}: got {got}, want {v}")
# also census the SPLIT template alone (its own gate says 87/166)
se = {}
stypes = {n["label"]: n["type"] for n in split["nodes"]}
for e in split["edges"]:
    if e["tag"] == "mn_to_muscle":
        continue
    if stypes[e["from"]] == "PORT-load":
        continue
    se[e["tag"]] = se.get(e["tag"], 0) + 1
print("split census total:", sum(se.values()))
for k in sorted(se):
    print(f"  {k}: {se[k]}")
