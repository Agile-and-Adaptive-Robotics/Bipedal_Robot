import re
p = r"D:/GitHub/Bipedal_Robot/Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim"
s = open(p, encoding='utf-8', errors='replace').read()
for m in re.finditer(r'<Thelen2003Muscle name="([^"]+)".*?</Thelen2003Muscle>', s, re.S):
    if m.group(1) in ('bifemsh_r', 'med_gas_r', 'soleus_r', 'tib_ant_r', 'rect_fem_r'):
        blk = m.group(0)
        fmax = re.search(r'<max_isometric_force>([^<]+)</max_isometric_force>', blk)
        pts = re.findall(
            r'<PathPoint name="([^"]+)">\s*<location>([^<]*)</location>\s*'
            r'<socket_parent_frame>([^<]*)</socket_parent_frame>', blk)
        print(f"== {m.group(1)}  fmax={fmax.group(1) if fmax else '?'}")
        for name, loc, fr in pts:
            print(f"   {name}: frame={fr.split('/')[-1]} loc={loc.strip()}")
