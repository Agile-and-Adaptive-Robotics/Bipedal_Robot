"""Audit: dump version-bearing entries from .slx files."""
import zipfile

FILES = [
    (r"D:\GitHub\Bipedal_Robot\tmp\audit_gate2\SNS_Library_HEAD.slx", "committed HEAD SNS_Library"),
    (r"D:\GitHub\Bipedal_Robot\tmp\audit_gate2\KneeReflexDemo_HEAD.slx", "committed HEAD KneeReflexDemo"),
    (r"D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\SNS_Library.slx", "working SNS_Library"),
    (r"D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\demos\KneeReflexDemo.slx", "working KneeReflexDemo"),
]
ENTRIES = ["metadata/mwcorePropertiesReleaseInfo.xml", "simulink/blockdiagram.xml"]

for path, label in FILES:
    print("=" * 20, label)
    with zipfile.ZipFile(path) as z:
        for e in ENTRIES:
            try:
                data = z.read(e)
                head = data[:600]
                print(f"--- {e} ({len(data)} bytes) first 600:")
                print(head.decode("utf-8", errors="replace"))
            except KeyError:
                print(f"--- {e}: ABSENT")
