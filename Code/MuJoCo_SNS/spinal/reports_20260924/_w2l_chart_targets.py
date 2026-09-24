"""Dump the DataColumn targets of the L Hip PF chart in the modern W2L asim."""
import re

P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG_Standalone_modern.asim"
x = open(P, encoding="utf-8", errors="replace").read()
i = x.find("<OutputFilename>L Hip PF.txt</OutputFilename>")
assert i > 0
blk = x[i:i + 20000]
end = blk.find("</DataChart>")
blk = blk[:end]
for m in re.findall(r"<ColumnName>(.*?)</ColumnName>|<TargetID>(.*?)</TargetID>|"
                    r"<DataType>(.*?)</DataType>", blk):
    a, b, c = m
    if a:
        print("COLUMN:", a)
    elif b:
        print("  target:", b)
    elif c:
        print("  type:", c)
