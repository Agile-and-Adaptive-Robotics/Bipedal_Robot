"""Add ExtTest20mm_1 and FlxTest20mm_51cm tabs to Results_table_20mm.xlsx.

Both tabs are stylized after the existing 'FlxTest20mm_42cm (3)' tab:
same label skeleton, same purple data-column styling, same row-16 torque
formula pattern -- except the extensor tab drops the leading negative
(extension torque is positive). Data cells are left blank: these are
data-entry templates for the upcoming 20 mm tests.

Idempotent: refuses to run if either sheet already exists.
"""
from copy import copy
from openpyxl import load_workbook
from openpyxl.chart import ScatterChart, Reference, Series
from openpyxl.comments import Comment
from openpyxl.chart.marker import Marker
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.drawing.line import LineProperties
from openpyxl.drawing.spreadsheet_drawing import TwoCellAnchor, AnchorMarker
from openpyxl.utils import get_column_letter

PATH = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Testing_Data\2026_06_Festo\Results_table_20mm.xlsx"
TPL = "FlxTest20mm_42cm (3)"

# data columns C..X == test slots 1..22 (same span as the template tab)
FIRST_COL, LAST_COL = 3, 24

NEW_TABS = [
    # (sheet name, A1 title, chart title, negative sign on row 16, insert index)
    ("ExtTest20mm_1", "Extensor Test 20mm", "20 mm Extensor Torque", False, 2),
    ("FlxTest20mm_51cm", "Flexor Test 20mm", "20 mm Flexor Torque", True, None),
]

# Label cells copied verbatim from the template (same coordinates).
LABEL_CELLS = [
    "A2", "B2", "B3", "D3", "B4", "D4", "B5", "B6", "B7", "B8",
    "A9", "B9", "B10", "B11", "A12", "B12", "B13", "B14", "B15",
    "B16", "B17", "B21", "C21",
]

wb = load_workbook(PATH)
tpl = wb[TPL]

for name in [t[0] for t in NEW_TABS]:
    if name in wb.sheetnames:
        raise SystemExit(f"Sheet '{name}' already exists -- aborting without changes.")

for sheet_name, title, chart_title, negative, insert_at in NEW_TABS:
    ws = wb.create_sheet(sheet_name, insert_at)

    # Column widths from the template
    for col_letter, dim in tpl.column_dimensions.items():
        if dim.width is not None:
            ws.column_dimensions[col_letter].width = dim.width

    # Title + labels (clone template style + text)
    for coord in ["A1"] + LABEL_CELLS:
        src = tpl[coord]
        dst = ws[coord]
        dst._style = copy(src._style)
        dst.value = title if coord == "A1" else copy(src.value)

    # Test-number header + blank data rows, uniform purple (template col C style)
    data_rows = [5, 6, 7, 8, 10, 11, 13, 14, 16, 17]
    for r in data_rows:
        style_src = tpl.cell(row=r, column=3)  # column C of this row
        for cidx in range(FIRST_COL, LAST_COL + 1):
            dst = ws.cell(row=r, column=cidx)
            dst._style = copy(style_src._style)
            if r == 5:
                dst.value = cidx - FIRST_COL + 1  # Test # 1..22
            elif r == 16:
                col = get_column_letter(cidx)
                sign = "-" if negative else ""
                dst.value = (
                    f"={sign}{col}6*COS(RADIANS({col}10-2.83))*{col}13/1000"
                )

    # Workflow note on the Torque label (row 16/17 semantics)
    sign_note = "flexion (negative)" if negative else "extension (positive)"
    ws["B16"].comment = Comment(
        "Row 16 'Torque' = planar approximation "
        "-Load*cos(LC angle - 2.83 deg)*Tibia origin/1000 "
        f"({sign_note} sign convention; update the 2.83 deg offset per rig "
        "measurement).\n"
        "Row 17 'Torque actual' = measured torque via the Adjoint transform, "
        "printed by " + ("Knee_Flexor_Data_20mm.m" if negative else "Knee_Extensor_20mm.m") +
        " (FlxTest20mm_51cm / ExtTest20mm_1 section) for paste here.",
        "ZCode")

    # Chart: markers-only scatter, knee angle (row 7) vs Torque (row 16)
    ch = ScatterChart()
    ch.scatterStyle = "lineMarker"
    ch.title = chart_title
    xref = Reference(ws, min_col=FIRST_COL, max_col=LAST_COL, min_row=7, max_row=7)
    yref = Reference(ws, min_col=FIRST_COL, max_col=LAST_COL, min_row=16, max_row=16)
    s = Series(yref, xref, title=None)
    s.marker = Marker(symbol="circle", size=7)
    s.graphicalProperties = GraphicalProperties(ln=LineProperties(noFill=True))
    s.smooth = False
    ch.series.append(s)
    ch.x_axis.delete = False
    ch.y_axis.delete = False
    _from = AnchorMarker(col=9, colOff=0, row=22, rowOff=0)
    to = AnchorMarker(col=16, colOff=0, row=37, rowOff=0)
    ch.anchor = TwoCellAnchor(editAs="twoCell", _from=_from, to=to)
    ws.add_chart(ch)

    print(f"added sheet '{sheet_name}'")

wb.save(PATH)
print("saved", PATH)
print("sheet order:", wb.sheetnames)
