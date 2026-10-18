anchors = [
    "The moment arm is the quantity",
    "The human knee is the joint where",
    "The mechanics of that master",
    "fig:masters_routing_and_arms",
    "fig:masters_torque_comparison",
    "Three published circuits",
    "fig:deng_connectome",
    "fig:afferent_connectomes",
    "fig:shinohara_connectome",
    "approximations with known limits",
    "accurate for a purpose",
    "includegraphics",
]
path = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\15-background.tex"
lines = open(path, encoding="utf-8").read().splitlines()
for a in anchors:
    found = False
    for i, l in enumerate(lines, 1):
        if a in l:
            print(f"{i:4d} | {a[:46]}")
            found = True
            break
    if not found:
        print(f"  ?? | {a[:46]} NOT FOUND")
print("--- all includegraphics lines:")
for i, l in enumerate(lines, 1):
    if "\\includegraphics" in l:
        print(f"{i:4d} | {l.strip()[:100]}")
