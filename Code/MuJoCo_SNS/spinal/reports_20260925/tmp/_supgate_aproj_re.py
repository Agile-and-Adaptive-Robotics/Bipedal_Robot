"""SUPERVISOR GATE: targeted probe of the claimed clone-type bug in Ben's
Walker_2_Layer_CPG_BilateralRG.aproj: link 'R Hip MN flx RE' -> 'R Hip MN
ext RE' claimed typed 'V3 Commissural Excite' (excitatory) where the L
mirror and the 2023 original are 'RE to RE Inhibit'. Read-only.
"""
import io
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

APROJ = (r"D:\GitHub\Bipedal_Robot\Neuromechanical_Models"
         r"\Walker_2_Layer_CPG_BilateralRG"
         r"\Walker_2_Layer_CPG_BilateralRG.aproj")
x = open(APROJ, encoding="utf-8", errors="replace").read()
print(f"aproj bytes: {len(x)}")

# strip CDATA (drawings) so <Node ...> in drawings can't confuse us
x_nc = re.sub(r"<!\[CDATA\[.*?\]\]>", "<![CDATA[]]>", x, flags=re.S)

# synapse types: blocks that carry a Type name + an ID
types = {}
for m in re.finditer(r"<SynapseType>(.*?)</SynapseType>", x_nc, re.S):
    b = m.group(1)
    nm = re.search(r"<Name>([^<]+)</Name>", b)
    i = re.search(r"<ID>([^<]+)</ID>", b)
    if nm and i:
        types[i.group(1)] = nm.group(1)
# fallback: any Type-name-bearing blocks
if len(types) < 3:
    for m in re.finditer(r"<(\w+[^>]*Type)>(.*?)</\1>", x_nc, re.S):
        pass
print(f"synapse types found: {len(types)}")

# neurons: <Node> blocks in the object tree with <Text> label + <ID>
nodes = {}
for m in re.finditer(r"<Node>(.*?)</Node>", x_nc, re.S):
    b = m.group(1)
    t = re.search(r"<Text>([^<]+)</Text>", b)
    i = re.search(r"<ID>([^<]+)</ID>", b)
    if t and i:
        nodes[i.group(1)] = t.group(1)
print(f"nodes found: {len(nodes)}")
by_name = {}
for i, n in nodes.items():
    by_name.setdefault(n, i)

pairs = [("R Hip MN flx RE", "R Hip MN ext RE"),
         ("L Hip MN flx RE", "L Hip MN ext RE")]
links = re.findall(r"<Link>(.*?)</Link>", x_nc, re.S)
print(f"link blocks: {len(links)}")
for src, dst in pairs:
    sid, did = by_name.get(src), by_name.get(dst)
    if not sid or not did:
        print(f"{src} -> {dst}: NODE NOT FOUND ({sid},{did})")
        continue
    hits = []
    for b in links:
        o = re.search(r"<OriginID>([^<]+)</OriginID>", b)
        d = re.search(r"<DestinationID>([^<]+)</DestinationID>", b)
        if o and d and o.group(1) == sid and d.group(1) == did:
            hits.append(b)
    for b in hits:
        tid = re.search(r"<SynapseTypeID>([^<]+)</SynapseTypeID>", b)
        if not tid:
            tid = re.search(r"<TypeID>([^<]+)</TypeID>", b)
        tn = types.get(tid.group(1)) if tid else None
        eq = re.search(r'EquilibriumPotential[^/]*Actual="([^"]+)"', b)
        mc = re.search(r'MaxConductance[^/]*Actual="([^"]+)"', b)
        g = re.search(r"<G>([^<]+)</G>", b)
        print(f"{src} -> {dst}:")
        print(f"    type = {tn!r} (tid={tid.group(1) if tid else None})")
        print(f"    equil={eq.group(1) if eq else '?'} "
              f"maxcond={mc.group(1) if mc else '?'} G={g.group(1) if g else '?'}")
    if not hits:
        print(f"{src} -> {dst}: LINK NOT FOUND")
