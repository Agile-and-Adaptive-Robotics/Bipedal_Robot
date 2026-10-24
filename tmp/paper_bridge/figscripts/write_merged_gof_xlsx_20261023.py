# write_merged_gof_xlsx_20261023.py -- companion workbook for the merged
# pinned-flexor GoF table (Tables/PinnedFlexor_GoF_merged_20261023.xlsx).
# Mirrors the Main.tex merged sidewaystable: four flexor tests with
# Training/Validation groups and in-group means, four prediction families
# (Hunt equation, improved characterization, back-calculated hybrid,
# corrected model), GoF = RMSE / FVU / Max. Hunt and hybrid group means are
# arithmetic means of the two in-group rows, rounded half-up to 4 decimals.
from decimal import Decimal, ROUND_HALF_UP
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment

def q4(x):
    return float(Decimal(str(x)).quantize(Decimal('0.0001'), rounding=ROUND_HALF_UP))

def q3(x):
    return float(Decimal(str(x)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP))

# rows: (l0, lten, hunt(rmse,fvu,max), char(...), hybrid(...), corr(...))
rows = [
    ("Training", 48.5, None,
     (5.779, 1.1761, 9.798), (7.617, 2.1278, 25.531), (1.635, 0.0991, 3.842), (2.165, 0.1718, 5.467)),
    ("Training", 41.0, None,
     (6.633, 3.2313, 11.409), (9.694, 2.1149, 21.076), (2.187, 0.1184, 4.243), (1.084, 0.0277, 1.766)),
    ("Training", None, None, "Mean", None, None, None),
    ("Validation", 45.7, None,
     (5.315, 0.9408, 9.991), (8.803, 1.9380, 26.264), (1.246, 0.0409, 2.275), (1.368, 0.0468, 3.614)),
    ("Validation", 40.6, 35,
     (3.945, 0.6750, 7.527), (6.373, 0.5935, 16.859), (2.085, 0.0755, 3.073), (1.611, 0.0427, 5.891)),
    ("Validation", None, None, "Mean", None, None, None),
]
# corrected/char means of record (from the former Table 2, full-precision run)
mean_char = {"Training": (8.656, 2.1213, 23.303), "Validation": (7.588, 1.2657, 21.562)}
mean_corr = {"Training": (1.624, 0.0998, 3.617), "Validation": (1.490, 0.0447, 4.752)}

headers = ["Group", "l0_cm", "lten_mm", "RowLabel",
           "Hunt_RMSE", "Hunt_FVU", "Hunt_Max",
           "CharAlone_RMSE", "CharAlone_FVU", "CharAlone_Max",
           "Hybrid_RMSE", "Hybrid_FVU", "Hybrid_Max",
           "Corrected_RMSE", "Corrected_FVU", "Corrected_Max"]

wb = Workbook()
ws = wb.active
ws.title = "FlxPin_GoF_merged"
ws.append(headers)
for c in range(1, len(headers) + 1):
    ws.cell(row=1, column=c).font = Font(bold=True)
    ws.cell(row=1, column=c).alignment = Alignment(horizontal="center", wrap_text=True)

def mean3(a, b):
    # RMSE and Max means to 3 dp, FVU means to 4 dp, half-up, matching the
    # LaTeX table cells (the of-record char/corrected means carry the same
    # precision: RMSE/Max 3 dp, FVU 4 dp)
    return (q3((a[0] + b[0]) / 2), q4((a[1] + b[1]) / 2), q3((a[2] + b[2]) / 2))

for r in rows:
    grp, l0, lt, hunt, char, hyb, corr = r
    if hunt == "Mean":
        tr = [x for x in rows if x[0] == grp and isinstance(x[3], tuple)]
        a, b = tr[0], tr[1]
        hm = mean3(a[3], b[3])
        qm = mean3(a[5], b[5])
        ws.append([grp, None, None, "Mean",
                   hm[0], hm[1], hm[2],
                   mean_char[grp][0], mean_char[grp][1], mean_char[grp][2],
                   qm[0], qm[1], qm[2],
                   mean_corr[grp][0], mean_corr[grp][1], mean_corr[grp][2]])
    else:
        ws.append([grp, l0, lt, "",
                   *hunt, *char, *hyb, *corr])

out = r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/Tables/PinnedFlexor_GoF_merged_20261023.xlsx"
wb.save(out)
print("WROTE", out)
for row in ws.iter_rows(values_only=True):
    print(row)
