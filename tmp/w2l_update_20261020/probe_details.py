"""GUID-level probe: syn_types gain table, c1/V3/RG edges, knee-RE cluster,
stimulus targets, muscle map, contact adapters."""
import json

M = json.load(open(r"D:\GitHub\Bipedal_Robot\tmp\w2l_update_20261020\w2l_aproj_mine_20261018.json"))
nodes = {n["guid"]: n for n in M["nodes"]}

print("== synapse TYPES with gains (NonSpikingChemical only) ==")
for name, t in sorted(M["syn_types"].items()):
    if t["class"] != "NonSpikingChemical":
        continue
    eref = t.get("eref")
    sign = "EXC" if (eref is not None and eref > -0.055) else "inh"
    print(f"   {name:<34} eref={eref}  gmax={t.get('gmax')}  syn_g={t.get('syn_g')}  -> {sign}")

print("\n== commissural + RG-type edges (GUID level) ==")
for e in M["edges"]:
    if e["type"] in ("RG Excite", "RG Inhibit", "RG to RG Excite", "RG to PF Excite",
                     "c1 Commissural Inhibit", "V3 Commissural Excite") and (
            "c1" in str(e["src"]) or "c1" in str(e["dst"]) or "V3" in str(e["src"])
            or "V3" in str(e["dst"]) or "RG" in str(e["src"])):
        if "RG" in str(e["src"]) or "c1" in str(e["src"]) or "V3" in str(e["src"]):
            print(f"   {e['src']:<16} -> {e['dst']:<16} type={e['type']:<24} g={e['g']} eref={e['eref']}")

print("\n== knee RE cluster (GUID level, both duplicates) ==")
for e in M["edges"]:
    s, d = str(e["src"]), str(e["dst"])
    if "Knee MN" in s + d and "RE" in s + d:
        print(f"   {e['src_guid'][:8]}:{s:<20} -> {e['dst_guid'][:8]}:{d:<20} type={e['type']} g={e['g']} eref={e['eref']}")

print("\n== the two L Knee MN ext RE node GUIDs ==")
for n in M["nodes"]:
    if (n.get("name") or "") == "L Knee MN ext RE":
        print("   ", n["guid"], n.get("TonicStimulus"), n.get("TimeConstant"))

print("\n== V3-type edge R Hip (the flagged odd edge) ==")
for e in M["edges"]:
    if e["type"] == "V3 Commissural Excite":
        print(f"   {e['src_guid'][:8]}:{e['src']:<20} -> {e['dst_guid'][:8]}:{e['dst']:<20} g={e['g']} eref={e['eref']}")

print("\n== contact neuron nodes + their adapters ==")
for n in M["nodes"]:
    nm = n.get("name") or ""
    if "contact" in nm.lower() and n["class"] == "NonSpiking":
        print("   node:", n["guid"][:8], nm, "tau=", n.get("TimeConstant"), "tonic=", n.get("TonicStimulus"))
for n in M["nodes"]:
    nm = n.get("name") or ""
    if n["class"] == "PhysicalToNodeAdapter" and ("contact" in nm.lower() or nm == ""):
        pass

print("\n== all edges FROM contact neurons (incl gains) ==")
for e in M["edges"]:
    if "contact" in str(e["src"]):
        print(f"   {e['src']:<16} -> {e['dst']:<16} type={e['type']:<30} g={e['g']}")

print("\n== afferent-node tonic/timeconst sample (II relays) ==")
for n in M["nodes"]:
    nm = n.get("name") or ""
    if nm.endswith("flx II") or nm.endswith("ext Ib"):
        print(f"   {nm:<18} tau={n.get('TimeConstant')} tonic={n.get('TonicStimulus')}")

print("\n== Stimulus elements ==")
for n in M["nodes"]:
    if n["class"] == "Stimulus":
        print("   ", n)

print("\n== Muscle nodes (names via Text) ==")
mus = [n for n in M["nodes"] if n["class"] == "Muscle"]
print("   count:", len(mus))
