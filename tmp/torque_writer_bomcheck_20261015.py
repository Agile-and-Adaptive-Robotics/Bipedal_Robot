import io
p = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex"
b = open(p, "rb").read()
print("first bytes:", b[:3])
c = b.decode("utf-8-sig")
print("non-ascii chars:", sum(1 for ch in c if ord(ch) > 127))
lines = c.splitlines()
print("total lines:", len(lines))
print("=== tail ===")
for ln in lines[-8:]:
    print(repr(ln[:110]))
