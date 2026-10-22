import os

ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022"
base = open(os.path.join(ROOT, r"baseline\chapters\30-torque.tex"), "rb").read().split(b"\r\n") if b"\r\n" in open(os.path.join(ROOT, r"baseline\chapters\30-torque.tex"), "rb").read() else open(os.path.join(ROOT, r"baseline\chapters\30-torque.tex"), "rb").read().split(b"\n")
raw = open(os.path.join(ROOT, r"baseline\chapters\30-torque.tex"), "rb").read()
print("baseline CRLF count:", raw.count(b"\r\n"), "LF count:", raw.count(b"\n"))
line7 = raw.split(b"\n")[6]
print("line7 first 80:", line7[:80])
print("line7 last 20:", line7[-20:])
print("line7 ends with CR:", line7.endswith(b"\r"))

patch = open(os.path.join(ROOT, r"patches\30-torque\02-publication-note.patch"), "rb").read()
print("patch CR count:", patch.count(b"\r"), "LF count:", patch.count(b"\n"))
minus = [l for l in patch.split(b"\n") if l.startswith(b"-This chapter connects")][0]
print("patch minus-line ends with CR:", minus.endswith(b"\r"))
print("minus == line7 minus CR:", minus.rstrip(b"\r") == line7.rstrip(b"\r"))
