"""Dump chart definitions from a sheet."""
from openpyxl import load_workbook

path = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Testing_Data\2026_06_Festo\Results_table_20mm.xlsx"
wb = load_workbook(path)
for name in ["FlxTest20mm_42cm (3)", "ExtTest40mm_1"]:
    ws = wb[name]
    print(f"=== {name}: {len(ws._charts)} charts")
    for ch in ws._charts:
        print(" type:", type(ch).__name__)
        print(" title:", ch.title)
        print(" anchor:", ch.anchor if isinstance(ch.anchor, str) else type(ch.anchor).__name__)
        try:
            a = ch.anchor
            if hasattr(a, "_from"):
                print("  from col,row:", a._from.col, a._from.row)
            if hasattr(a, "to") and a.to is not None:
                print("  to col,row:", a.to.col, a.to.row)
        except Exception as e:
            print("  anchor detail err", e)
        for s in ch.series:
            xr = s.xVal.numRef.f if s.xVal and s.xVal.numRef else (s.xVal.strRef.f if s.xVal else None)
            yr = s.yVal.numRef.f if s.yVal and s.yVal.numRef else None
            vr = s.val.numRef.f if s.val and s.val.numRef else None
            print("  series: tx=", s.tx.strRef.f if s.tx and s.tx.strRef else None,
                  "| xVal=", xr, "| yVal=", yr, "| val=", vr)
    print(" sheetview: showGridLines =", ws.sheet_view.showGridLines)
