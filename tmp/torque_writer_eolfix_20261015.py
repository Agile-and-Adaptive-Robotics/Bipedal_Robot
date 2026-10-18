import glob

for p in [r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\20-muscles.tex",
          r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\40-ankle.tex",
          r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\main.tex"]:
    b = open(p, "rb").read()
    print(p.split("\\")[-1], "| BOM:", b[:3] == b"\xef\xbb\xbf",
          "| CRLF:", b.count(b"\r\n"), "| bare LF:", b.count(b"\n") - b.count(b"\r\n"))

p = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex"
b = open(p, "rb").read()
print("30-torque.tex before fix | BOM:", b[:3] == b"\xef\xbb\xbf",
      "| CRLF:", b.count(b"\r\n"), "| bare LF:", b.count(b"\n") - b.count(b"\r\n"))
if b[:3] == b"\xef\xbb\xbf":
    b = b[3:]
    open(p, "wb").write(b)
    print("BOM stripped.")
b = open(p, "rb").read()
print("30-torque.tex after fix  | BOM:", b[:3] == b"\xef\xbb\xbf",
      "| CRLF:", b.count(b"\r\n"), "| bare LF:", b.count(b"\n") - b.count(b"\r\n"))
