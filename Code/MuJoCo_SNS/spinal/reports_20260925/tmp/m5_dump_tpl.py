"""M5 pre-work: dump the w2laproj template node labels + edge families so the
afferent extension can reference exact labels."""
import io, json, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

tpl = json.load(open(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json", encoding="utf-8"))
t = tpl["w2laproj"]
print("== nodes ==")
for n in t["nodes"]:
    print(f"  {n['type']:<12} | {n['label']}")
print("\n== edges (from -> to | sign | gain | tag) ==")
for e in t["edges"]:
    print(f"  {e['from']:<28} -> {e['to']:<28} {e.get('sign','?'):<4} "
          f"g={e.get('gain','')}  tag={e.get('tag','')}")
