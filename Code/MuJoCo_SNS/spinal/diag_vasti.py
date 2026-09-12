"""Inspect the vastii tendon routes + which body the final sites sit on."""
import re

t = open(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
         r"\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml",
         encoding="utf-8").read()

# site -> owning body: walk every <body ...> open tag, remember the last one
# seen before each <site ... name="..."/> definition
site_body = {}
cur_body, cur_depth, site_depths = None, 0, {}
tokens = re.finditer(r'<(/?)body[^>]*?name="([^"]+)"[^>]*?(/?)>'
                     r'|<(/?)site\s+([^>]*?)/?>', t)
stack = []
for m in tokens:
    if m.group(2) is not None:            # body open/close
        closing, bname, selfclose = m.group(1), m.group(2), m.group(3)
        if closing:
            if stack:
                stack.pop()
        else:
            stack.append(bname)
    else:                                  # site
        sname = re.search(r'name="([^"]+)"', m.group(5))
        if sname:
            site_body[sname.group(1)] = stack[-1] if stack else "WORLD?"

for name in ("vas_med_r", "vas_lat_r", "vas_int_r", "rect_fem_r"):
    m = re.search(rf'<spatial name="{name}_tendon">(.*?)</spatial>', t, re.S)
    if not m:
        print(f"{name}: tendon block not found")
        continue
    print(f"{name}_tendon:")
    for s in re.findall(r'site="([^"]+)"', m.group(1)):
        print(f"   {s:44s} body={site_body.get(s, '???')}")
