"""Schema discovery + miner for the 2026-10-18 W2L BilateralRG aproj.

Phase 1 (--schema): census element types, dump one sample of each neural-ish element.
Phase 2 (--mine): full inventory -> w2l_aproj_mine_20261018.json
"""
import sys
import json
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

APROJ = r"D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\Walker_2_Layer_CPG_BilateralRG.aproj"
OUT = r"D:\GitHub\Bipedal_Robot\tmp\w2l_update_20261020\w2l_aproj_mine_20261018.json"

INTEREST = ("Neuron", "Synapse", "Adapter", "Contact", "Stimulus", "Motor", "Muscle",
            "Spiking", "NonSpiking", "Gain", "Receiver", "Organism", "Node")


def sample_el(el, depth=0):
    """Compact dump of an element: tag, key attrs, child tags."""
    attrs = {k: v for k, v in el.attrib.items()}
    kids = defaultdict(int)
    for ch in el:
        kids[ch.tag] += 1
    return {"tag": el.tag, "attrs": attrs, "child_tags": dict(kids)}


def main(phase):
    tree = ET.parse(APROJ)
    root = tree.getroot()

    if phase == "--schema":
        c = Counter(el.tag for el in root.iter())
        print("== element census (neural-ish) ==")
        shown = 0
        for tag, n in c.most_common():
            if any(k.lower() in tag.lower() for k in INTEREST):
                print(f"{n:6d}  {tag}")
        print("\n== one sample each of the main neural elements ==")
        seen = set()
        for el in root.iter():
            t = el.tag
            if t in seen:
                continue
            if any(k.lower() in t.lower() for k in ("Neuron", "Synapse")):
                seen.add(t)
                s = sample_el(el)
                print(json.dumps(s, indent=1)[:1200])
                print("-" * 60)
        # find where synapses live (parent chains)
        print("\n== parent chain of first Synapse-ish element ==")
        parent_map = {c2: p for p in root.iter() for c2 in p}
        for el in root.iter():
            if "Synapse" in el.tag:
                chain = []
                cur = el
                while cur in parent_map and len(chain) < 12:
                    cur = parent_map[cur]
                    chain.append(f"{cur.tag}(id={cur.attrib.get('ID','?')})")
                print(" <- ".join(chain))
                break
        return

    if phase == "--mine":
        items = defaultdict(list)
        for el in root.iter():
            t = el.tag
            if any(k.lower() in t.lower() for k in ("Neuron", "Synapse")):
                d = {"tag": t}
                d.update(el.attrib)
                # numeric children (AnimatLab stores params as child elements)
                pl = {}
                for ch in el:
                    if len(ch) == 0 and (ch.text or "").strip():
                        pl[ch.tag] = ch.text.strip()
                if pl:
                    d["_params"] = pl
                items[t].append(d)
        out = {"source": APROJ, "counts": {k: len(v) for k, v in items.items()}}
        for k, v in items.items():
            out[k] = v
        with open(OUT, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1)
        print("wrote", OUT)
        print(json.dumps(out["counts"], indent=1))
        return


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "--schema")
