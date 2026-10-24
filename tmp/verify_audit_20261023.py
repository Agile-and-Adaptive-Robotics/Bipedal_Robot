# One-shot spot-check of AUDIT.md's count-based claims (read-only).
import json, os, re

REPO = r"D:\Github\Bipedal_Robot"
SP = os.path.join(REPO, "Code", "MuJoCo_SNS", "spinal")
AUD = os.path.join(SP, "reports_20261023", "connectome_audit")

# 1. HTML embedded library vs sidecar
html = open(os.path.join(SP, "connectome_block_editor.html"), encoding="utf-8", errors="replace").read()
side = json.load(open(os.path.join(SP, "connectome_templates.json"), encoding="utf-8"))
print("HTML 'bilateralrg2' occurrences:", html.count("bilateralrg2"))
print("HTML '_meta' occurrences (proxy for embedded entries):", html.count('"_meta"'))
print("sidecar top-level keys (%d):" % len(side), sorted(side.keys()))

# 2. inventory.json
inv = json.load(open(os.path.join(AUD, "inventory.json"), encoding="utf-8"))
emb = inv["embedded_templates"]
emb_keys = [e["key"] for e in emb]
print("inventory embedded entries:", len(emb))
print("all matches_sidecar true:", all(e.get("matches_sidecar") is True for e in emb))
print("bilateralrg2 in embedded inventory:", "bilateralrg2" in emb_keys)
print("duplicate embedded keys:", len(emb_keys) - len(set(emb_keys)))
print("inventory other sections:", [k for k in inv.keys() if k != "embedded_templates"])

# 3. W2L .aproj page counts + XML well-formedness
w2l_dir = os.path.join(REPO, "Neuromechanical_Models", "Walker_2_Layer_CPG_BilateralRG")
import xml.dom.minidom as md
for f in ["Walker_2_Layer_CPG_BilateralRG.aproj",
          "Walker_2_Layer_CPG_BilateralRG_control20261020.aproj",
          "Walker_2_Layer_CPG_BilateralRG_tuned20261020.aproj"]:
    p = os.path.join(w2l_dir, f)
    s = open(p, encoding="utf-8", errors="replace").read()
    nsub = s.count("<Subsystem>")
    ndiag = s.count("<DiagramXml")
    ok = "?"
    try:
        md.parse(p)
        ok = "WELL-FORMED"
    except Exception as e:
        ok = "PARSE FAIL: %s" % e
    print("%s | <Subsystem>: %d | <DiagramXml: %d | %s" % (f, nsub, ndiag, ok))

# 4. dirusso draft
p = os.path.join(SP, "replication", "dirusso_connectome_draft.json")
d = json.load(open(p, encoding="utf-8"))
print("dirusso top-level keys:", sorted(d.keys()) if isinstance(d, dict) else type(d))
blob = json.dumps(d)
print("dirusso 'phase' mentions:", len(re.findall(r"phase", blob, re.I)))
print("dirusso 'synergy|primitive|nmf' mentions:", len(re.findall(r"synergy|primitive|nmf", blob, re.I)))

# 5. ben_shinohara_20260924.json
cands = []
for root in [os.path.join(REPO, "Neuromechanical_Models", "Mujoco_SNS_models"),
             SP, os.path.join(SP, "w2l_cpg"), os.path.join(SP, "replication")]:
    if os.path.isdir(root):
        for dirpath, _, files in os.walk(root):
            for fn in files:
                if fn == "ben_shinohara_20260924.json":
                    cands.append(os.path.join(dirpath, fn))
print("shinohara json found:", cands)
if cands:
    sh = json.load(open(cands[0], encoding="utf-8"))
    nodes = sh.get("nodes", [])
    edges = sh.get("edges", [])
    mns = [n for n in nodes if "MN" in str(n.get("type", "")) or "MN" in str(n.get("label", ""))]
    subs = [n for n in nodes if str(n.get("type", "")).lower().startswith("sub") or n.get("subsystem")]
    print("shinohara nodes:", len(nodes), "edges:", len(edges), "MN-ish nodes:", len(mns), "subsystem-ish:", len(subs))
    print("shinohara node types:", sorted({str(n.get('type')) for n in nodes}))
