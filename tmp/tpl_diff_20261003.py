"""Diff guard (the _tpl_diff_groups_20260924 pattern): every pre-existing
template entry must be byte-identical after adding spiking_mirror; report
the new entry's note + counts."""
import io, json, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

OLD = r"D:\Github\Bipedal_Robot\tmp\connectome_templates_before_20261003.json"
NEW = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json"

old = json.load(open(OLD, encoding="utf-8"))
new = json.load(open(NEW, encoding="utf-8"))

added = [k for k in new if k not in old]
removed = [k for k in old if k not in new]
changed = []
for k in old:
    if k in new and json.dumps(old[k], sort_keys=True) != \
            json.dumps(new[k], sort_keys=True):
        changed.append(k)

print("added entries:  ", added)
print("removed entries:", removed)
print("CHANGED pre-existing entries:", changed if changed else "NONE")

sm = new.get("spiking_mirror", {})
print("\nspiking_mirror: %d nodes / %d edges" %
      (len(sm.get("nodes", [])), len(sm.get("edges", []))))
print("groups:", sorted(sm.get("groups", {}).keys()))
print("\nnote:\n" + sm.get("_note", ""))

# quick integrity: types all in catalog, no dangling edges, no dup labels
TYPES = set("""SN-Ia SN-II SN-Ib SN-heel SN-toe PORT-load IN-V0D IN-V0V IN-V1
IN-V2a IN-V2b IN-V3 IN-C HC-RG-E HC-RG-F IN-InE IN-InF HC-PF-E HC-PF-F IN-PF
IN-IaIN IN-IbIN IN-IIe IN-IIi IN-Ib+ IN-LBIN IN-KINH RC MN MUSCLE SUB""".split())
labels = [n["label"] for n in sm["nodes"]]
lab = set(labels)
bad_t = [n["type"] for n in sm["nodes"] if n["type"] not in TYPES]
dang = [e for e in sm["edges"] if e["from"] not in lab or e["to"] not in lab]
dups = len(labels) - len(lab)
signs = {e["sign"] for e in sm["edges"]}
spk_tags = sum(1 for e in sm["edges"] if e["tag"].endswith(" spk"))
grd_tags = sum(1 for e in sm["edges"] if e["tag"].endswith(" grd"))
print("\nintegrity: bad types=%s | dangling=%d | dup labels=%d | signs=%s"
      % (sorted(set(bad_t)) if bad_t else "NONE", len(dang), dups, signs))
print("edge coding: %d spike-synapse edges, %d graded edges"
      % (spk_tags, grd_tags))
sys.exit(1 if (changed or removed or bad_t or dang or dups) else 0)
