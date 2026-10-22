import os

ROOT = r"D:\Github\Ben Bolen".replace("Ben Bolen", "")  # noop guard
ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022"
raw = open(os.path.join(ROOT, r"baseline\chapters\30-torque.tex"), "rb").read()
line7 = raw.split(b"\n")[6]

patch = open(os.path.join(ROOT, r"patches\30-torque\02-publication-note.patch"), "rb").read()
minus = [l for l in patch.split(b"\n") if l.startswith(b"-This chapter connects")][0]

a = line7.rstrip(b"\r")
b = minus.rstrip(b"\r")
print("len base:", len(a), "len patch:", len(b))
for i in range(min(len(a), len(b))):
    if a[i] != b[i]:
        print("first diff at", i)
        print("base :", a[i - 30 : i + 30])
        print("patch:", b[i - 30 : i + 30])
        break
else:
    print("common prefix identical; tail base:", a[len(b) :][:60], "tail patch:", b[len(a) :][:60])
