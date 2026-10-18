r"""W2L ROBOT-SIDE RULES AUDIT (overnight campaign 2026-10-16, egg).

Audits spinal\w2l_cpg (bilateralrg contact-driven CPG) + spinal\w2l_mujoco
(split-RG MJCF robot, M5 per-rules afferents) against Ben's rules JSON
(same canonical pair the sibling audits use):
    spinal/ben_rules_20260924.json
    Neuromechanical_Models/Mujoco_SNS_models/
        Circuit_rules_CONNECTOME_md__connectome.json

Checks (all against the COMPILED builders' own template dicts — the
builders assert 1:1 transcription, so the template edges ARE the built
edges; gains come from each net's merged G table):
  1. rules-file identity (both copies byte-identical payload)
  2. heel family: rules 5 edges vs implemented (ipsi InE, contra InF,
     hip/knee PF-E INs; ankle PF-E IN N/A - no ankle PF layer)
  3. toe family: chain gains/signs vs rules (1 / 5 / 2.749) + terminal
     reduction note + the rules' second TOEDF-IN input (ankle PF-E HC)
  4. Ib load family: rules 3 out-edges vs implemented 4 per side
  5. afferent -> RG/PF HALF-CENTER enumeration in every net
     (intended-per-rules vs other)
  6. autogenic structure: per-pool ib_auto / ia_recip chains + the
     active Ib-group signal source
  7. muscle-map TODO scan: who imports w2l_cpg / muscle_map
  8. census + per-side MN/muscle label dump (drive-plan inputs)

Output: stdout (captured to rules_audit.log by the driver).
"""
from __future__ import annotations

import hashlib
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = Path(__file__).resolve().parent          # .../reports_20261016/w2l_audit
SPINAL = HERE.parent.parent                     # .../spinal
REPO = SPINAL.parent.parent.parent              # repo root
sys.path.insert(0, str(SPINAL))
sys.path.insert(0, str(SPINAL / "w2l_cpg"))
sys.path.insert(0, str(SPINAL / "w2l_mujoco"))

RULES_A = SPINAL / "ben_rules_20260924.json"
RULES_B = (REPO / "Neuromechanical_Models" / "Mujoco_SNS_models" /
           "Circuit_rules_CONNECTOME_md__connectome.json")

print("=" * 78)
print("1. RULES FILE IDENTITY")
a = json.loads(RULES_A.read_text(encoding="utf-8"))
b = json.loads(RULES_B.read_text(encoding="utf-8"))
sa = json.dumps(a, sort_keys=True)
sb = json.dumps(b, sort_keys=True)
print(f"   {RULES_A}")
print(f"   {RULES_B}")
print(f"   payload identical: {sa == sb}  "
      f"(sha {hashlib.sha256(sa.encode()).hexdigest()[:16]})")
nodes = {n['label']: n for n in a['nodes']}
redges = a['edges']

print()
print("=" * 78)
print("2-4. RULES FAMILIES (verbatim from the JSON)")
fam = {}
for e in redges:
    if e['from'] == 'heel: heel':
        fam.setdefault('heel', []).append(e)
    if e['from'] == 'toe: toe':
        fam.setdefault('toe', []).append(e)
    if e['from'] == 'toe: toe IN':
        fam.setdefault('toe_in', []).append(e)
    if e['from'] == 'ib_load: Ib grp':
        fam.setdefault('ib_load', []).append(e)
    if e['from'] == 'IN-PF_dorsiflexion_inhibit':
        fam.setdefault('toedf_out', []).append(e)
    if nodes.get(e['from'], {}).get('type', '') in ('SN-Ia', 'SN-II') \
            and nodes.get(e['to'], {}).get('type', '').startswith('HC-'):
        fam.setdefault('sn_ia_ii_to_HC', []).append(e)
for k in ('heel', 'toe', 'toe_in', 'toedf_out', 'ib_load', 'sn_ia_ii_to_HC'):
    print(f"   [{k}] {len(fam.get(k, []))} edges")
    for e in fam.get(k, []):
        dt = nodes.get(e['to'], {}).get('type', '?')
        print(f"      {e['from']!r} -> {e['to']!r} ({dt}) "
              f"sign={e['sign']} gain={e.get('gain')}")
# extra: edges INTO the TOEDF IN (second input family)
toedf_in = [e for e in redges if e['to'] == 'IN-PF_dorsiflexion_inhibit']
print(f"   [edges INTO IN-PF_dorsiflexion_inhibit] {len(toedf_in)}")
for e in toedf_in:
    st = nodes.get(e['from'], {}).get('type', '?')
    print(f"      {e['from']!r} ({st}) -> sign={e['sign']} gain={e.get('gain')}")

print()
print("=" * 78)
print("5-8. BUILT NETS")
import build_w2l_net as CPG          # noqa: E402  (bilateralrg, w2l_cpg)
import build_w2l_split_net as SPLIT  # noqa: E402
import build_w2l_aff_net as AFF      # noqa: E402

AFF_TPL_RAW = json.loads((SPINAL / "connectome_templates.json").
                         read_text(encoding="utf-8"))


def census(name, net):
    print(f"   [{name}] neurons={len(net.idx)} synapses={net.n_synapses} "
          f"inputs={len(net.inputs)} outputs={len(net.muscle_outputs)}")
    print(f"      ports: {net.inputs}")


def edges_of(net):
    return net.template["edges"]


def types_of(net):
    return {n["label"]: n["type"] for n in net.template["nodes"]}


def dump_family(net, tags, title):
    T = types_of(net)
    print(f"   -- {title}")
    for e in edges_of(net):
        if e.get("tag") in tags:
            print(f"      {e['from']!r} -> {e['to']!r}  "
                  f"({T.get(e['to'],'?')}) sign={e['sign']} "
                  f"g={net.G.get('aff_' + e['tag'], net.G.get(e['tag'],'?'))} "
                  f"tag={e['tag']}")


print()
print("--- w2l_cpg / bilateralrg (contact-driven CPG) ---")
cpg = CPG.build()
census("bilateralrg", cpg)
dump_family(cpg, {"aff_HC_excite", "contact_C20"},
            "afferent->HC (aff_HC_excite) + contact (contact_C20)")

print()
print("--- w2l_mujoco split-RG (M4) ---")
split = SPLIT.build()
census("split coupled", split)

print()
print("--- w2l_mujoco aff net (M5, Ben's-rules afferents ON) ---")
aff = AFF.build()
census("aff", aff)
dump_family(aff, {"heel_rge_in", "heel_pf", "toe_in", "toe_df", "toe_df_mn",
                  "ib_load"},
            "M5 afferent families (built-in gains)")

print()
print("=" * 78)
print("6. AFFERENT -> RG/PF HALF-CENTER EDGES (every net, enumerated)")
AFF_TYPES = ("SN-Ia", "SN-Ib", "SN-II", "SN-heel", "SN-toe", "IN-C")
for name, net, extra_aff in (("bilateralrg", cpg, set()),
                             ("split", split, set()),
                             ("aff", aff,
                              {"heel L", "heel R", "toe L", "toe R",
                               "toe IN L", "toe IN R", "TOEDF IN L",
                               "TOEDF IN R", "Ib grp L", "Ib grp R"})):
    T = types_of(net)
    hits = []
    for e in edges_of(net):
        src_t = T.get(e["from"])
        is_aff = (src_t in AFF_TYPES) or (e["from"] in extra_aff)
        if not is_aff:
            continue
        dst_t = T.get(e["to"], "")
        if dst_t.startswith("HC-") or dst_t in ("IN-InE", "IN-InF"):
            hits.append((e["from"], e["to"], dst_t, e["sign"], e["tag"]))
    print(f"   [{name}] {len(hits)} afferent-source edges onto HC/InE/InF:")
    for f, t, dt, s, tag in hits:
        print(f"      {f!r} -> {t!r} ({dt}) sign={s} tag={tag}")

print()
print("=" * 78)
print("7. AUTOGENIC STRUCTURE (template per-pool chains + active Ib source)")
for name, net in (("bilateralrg", cpg), ("aff", aff)):
    T = types_of(net)
    ib_auto = [(e["from"], e["to"]) for e in edges_of(net)
               if e["tag"] == "ib_auto"]
    ia = [(e["from"], e["to"], e["sign"]) for e in edges_of(net)
          if e["tag"] == "ia_recip"]
    print(f"   [{name}] ib_auto (SN-Ib -> homonymous MN) {len(ib_auto)} edges:"
          f" sources={sorted({f for f, _ in ib_auto})}")
    print(f"      ia_recip {len(ia)} edges; SN-Ia sources="
          f"{sorted({f for f, _, _ in ia if T.get(f)=='SN-Ia'})}")

print()
print("=" * 78)
print("8. MUSCLE-MAP TODO SCAN (spinal root, one level of .py)")
hits = []
for p in sorted(SPINAL.glob("*.py")):
    txt = p.read_text(encoding="utf-8", errors="replace")
    for needle in ("w2l_cpg", "muscle_map"):
        if needle in txt:
            hits.append((p.name, needle))
w2lm = sorted((SPINAL / "w2l_mujoco").glob("*.py"))
print(f"   w2l_mujoco .py files: {len(w2lm)}")
for p in sorted(SPINAL.glob("*.py")):
    txt = p.read_text(encoding="utf-8", errors="replace")
    if "w2l_cpg" in txt:
        print(f"   spinal/{p.name}: imports/mentions w2l_cpg")
print("   (runner.py wiring w2l_cpg into the 92-muscle model would appear "
      "above; README says NOT wired)")

print()
print("=" * 78)
print("9. PER-SIDE MN + MUSCLE LABELS (drive-plan entry points)")
for name, net in (("bilateralrg", cpg), ("aff", aff)):
    print(f"   [{name}] mn_to_muscle map (MN -> MUSCLE node label):")
    for e in edges_of(net):
        if e["tag"] == "mn_to_muscle":
            print(f"      {e['from']!r} -> {e['to']!r}")
    print(f"      actuator map (MUSCLE_MAP): "
          f"{getattr(net, 'muscle_outputs', {})}")

print()
print("AUDIT SCRIPT DONE")
