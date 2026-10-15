"""Dump chart axes/marker details compactly."""
from openpyxl import load_workbook

path = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Testing_Data\2026_06_Festo\Results_table_20mm.xlsx"
wb = load_workbook(path)
ws = wb["FlxTest20mm_42cm (3)"]
ch = ws._charts[0]
print("scatterStyle:", ch.scatterStyle)
print("x_axis: delete=", ch.x_axis.delete, " title=", ch.x_axis.title is not None)
if ch.x_axis.title is not None:
    try:
        print("  x title text:", ch.x_axis.title.tx.rich.p[0].r[0].t)
    except Exception:
        print("  x title: <complex>")
print("y_axis: delete=", ch.y_axis.delete, " title=", ch.y_axis.title is not None)
if ch.y_axis.title is not None:
    try:
        print("  y title text:", ch.y_axis.title.tx.rich.p[0].r[0].t)
    except Exception:
        print("  y title: <complex>")
for i, s in enumerate(ch.series):
    print(f"series {i}: marker=", s.marker.symbol if s.marker else None,
          " line-none=", (s.graphicalProperties.line.noFill if s.graphicalProperties and s.graphicalProperties.line else None),
          " smooth=", s.smooth)
print("legend pos:", ch.legend.position if ch.legend else None)
print("style:", ch.style)
