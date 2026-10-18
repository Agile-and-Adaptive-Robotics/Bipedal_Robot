# Probe sheet names of the archived results workbooks (read-only).
import sys
try:
    import openpyxl
    HAVE_OX = True
except Exception:
    HAVE_OX = False
import os
BASE = r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo"
files = [
    "Results_table_10mm_pinned_LoadCell.xlsx",
    "Results_table_10mm_pinned_FishScale.xlsx",
    "Results_table_10mm.xlsx",
    r"D:\Github\Bipedal_Robot\Testing_Data\2026_06_Festo\Results_table_10mm_pinned.xlsx",
]
print("openpyxl:", HAVE_OX)
if HAVE_OX:
    for f in files:
        p = f if os.path.isabs(f) else os.path.join(BASE, f)
        try:
            wb = openpyxl.load_workbook(p, data_only=True)
            print(p, "->", wb.sheetnames)
        except Exception as e:
            print(p, "ERR", e)
