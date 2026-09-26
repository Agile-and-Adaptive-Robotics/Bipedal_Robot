# -*- coding: utf-8 -*-
"""Print corrected joint axes + anchors from the dump (for the report table)."""
import json
D = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_source_dump.json"
dump = json.load(open(D))
for b in dump["bodies"]:
    if "joint" in b:
        j = b["joint"]
        ax = j["axis_world_al"]
        print("%-8s body=%-8s axis_world_al=(% .4f, % .4f, % .4f) limits=[%s, %s] deg anchor=(%.4f, %.4f, %.4f)"
              % (j["name"], b["name"], ax[0], ax[1], ax[2], j["lower_deg"], j["upper_deg"],
                 *j["anchor_world_al"]))
