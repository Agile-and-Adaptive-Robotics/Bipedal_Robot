import zipfile, re, sys
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

def read_xlsx(path, max_rows=14):
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))
    # first worksheet
    sheet_name = [n for n in z.namelist() if re.match(r"xl/worksheets/sheet1\.xml", n)][0]
    root = ET.fromstring(z.read(sheet_name))
    rows = []
    for row in root.iter("{%s}row" % NS["m"]):
        cells = {}
        for c in row.findall("m:c", NS):
            ref = c.get("r"); col = re.match(r"[A-Z]+", ref).group(0)
            t = c.get("t"); v = c.find("m:v", NS)
            if t == "s" and v is not None:
                val = shared[int(v.text)]
            elif t == "inlineStr":
                val = "".join(x.text or "" for x in c.iter("{%s}t" % NS["m"]))
            else:
                val = v.text if v is not None else ""
            cells[col] = val
        if cells:
            rows.append(cells)
    return rows

base = "D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/Tables/"
for name in ["PinnedFlexor_Preopt_GoF.xlsx", "PinnedFlexor_GoF.xlsx", "PinnedExtensor_GoF.xlsx",
             "AdoptedCorrectionTerms.xlsx", "Route_PinnedBrackets.xlsx", "Route_Biomimetic10mm.xlsx",
             "Route_Redesign20mm.xlsx"]:
    print("=" * 20, name)
    try:
        rows = read_xlsx(base + name)
        for cells in rows[:14]:
            cols = sorted(cells, key=lambda c: (len(c), c))
            line = " | ".join(str(cells[c]) for c in cols)
            if line.strip():
                print(line[:210])
    except Exception as e:
        print("ERROR:", repr(e))
