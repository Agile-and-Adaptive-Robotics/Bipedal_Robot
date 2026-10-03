"""List top-level block names in SNS_Library: committed HEAD vs working copy."""
import re
import zipfile


def top_blocks(path):
    with zipfile.ZipFile(path) as z:
        root = z.read("simulink/systems/system_root.xml").decode("utf-8", errors="replace")
        # also try blockdiagram
        bd = z.read("simulink/blockdiagram.xml").decode("utf-8", errors="replace")
    names = re.findall(r'<Block BlockType="(\w+)" Name="([^"]+)"', root)
    return names, bd[:0]


for path, label in [
    (r"D:\GitHub\Bipedal_Robot\tmp\audit_gate2\SNS_Library_HEAD.slx", "HEAD (committed R2025b)"),
    (r"D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\SNS_Library.slx", "working (R2025a + spiking)"),
]:
    names, _ = top_blocks(path)
    print(f"== {label}: {len(names)} top-level blocks")
    for bt, nm in sorted(names, key=lambda x: x[1]):
        print(f"   {nm} ({bt})")
