"""Which joints drive the equality-coupled pathpoint followers?"""
import collections

import mujoco

import runner

m = mujoco.MjModel.from_xml_path(str(runner.MODEL))
jn = lambda j: mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_JOINT, j)
pairs = []
for i in range(m.neq):
    if m.eq_type[i] != mujoco.mjtEq.mjEQ_JOINT:
        continue
    follower = jn(int(m.eq_obj1id[i])) or "?"
    driver = jn(int(m.eq_obj2id[i])) or "?"
    pairs.append((follower, driver))

print("vastii pathpoint followers:")
for f, d in pairs:
    if "vas" in f:
        print(f"  {f} -> {d}")

print("\nall coupler drivers (counts):")
for d, n in collections.Counter(d for _, d in pairs).most_common():
    print(f"  {d}: {n}")
