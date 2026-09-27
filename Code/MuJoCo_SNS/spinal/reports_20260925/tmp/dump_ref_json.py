import json, sys
sys.path.insert(0, r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
import kine_ref

def conv(x):
    return x.tolist() if hasattr(x, "tolist") else x

ref = kine_ref.load_reference()
out = {}
for k, v in ref.items():
    out[k] = {kk: conv(vv) for kk, vv in v.items()} if isinstance(v, dict) else conv(v)
json.dump(out, open(sys.argv[1], "w"), indent=1)
print("dumped", sys.argv[1])
