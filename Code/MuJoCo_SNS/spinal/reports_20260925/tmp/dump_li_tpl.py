import json, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
d = json.load(open(r'D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json'))
t = d['li']
print("=== NODES (%d) ===" % len(t['nodes']))
for n in t['nodes']:
    print(json.dumps(n, ensure_ascii=False))
print("=== EDGES (%d) ===" % len(t['edges']))
for e in t['edges']:
    print(json.dumps(e, ensure_ascii=False))
