# -*- coding: utf-8 -*-
"""Rewrite Tables/PinnedFlexor_Preopt_GoF.xlsx with the corrected three-group
common-support numbers matching tab:FlxPinPreGof in Main.tex."""
from openpyxl import Workbook
from openpyxl.styles import Font

path = r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/Tables/PinnedFlexor_Preopt_GoF.xlsx"
wb = Workbook()
ws = wb.active
ws.title = "Preopt GoF"

hdr = ["Test", "Hunt RMSE", "Hunt FVU", "Hunt Max",
       "Improved char alone RMSE", "Improved char alone FVU", "Improved char alone Max",
       "Back-calculated (hybrid) RMSE", "Back-calculated (hybrid) FVU", "Back-calculated (hybrid) Max"]
rows = [
    ["48.5 cm (training)", 5.779, 1.1761, 9.798, 4.871, 0.8354, 9.865, 1.635, 0.0991, 3.842],
    ["41.0 cm (training)", 6.633, 3.2313, 11.409, 6.074, 2.7090, 9.262, 2.187, 0.1184, 4.243],
    ["45.7 cm (validation)", 5.315, 0.9408, 9.991, 4.465, 0.6640, 9.536, 1.246, 0.0409, 2.275],
    ["40.6 cm tendon (validation)", 3.945, 0.6750, 7.527, 3.226, 0.4512, 7.376, 2.085, 0.0755, 3.073],
    ["Mean (training)", 6.206, 2.2037, 10.604, 5.472, 1.7722, 9.563, 1.911, 0.1088, 4.043],
    ["Mean (validation)", 4.630, 0.8079, 8.759, 3.845, 0.5576, 8.456, 1.666, 0.0582, 2.674],
]
ws.append(hdr)
for r in rows:
    ws.append(r)
ws.append([])
ws.append(["Basis: Hunt-equation and improved-characterization columns are scored head-to-head "
           "on the measured points where the Hunt prediction is defined (75/90, 17/20, 14/19, "
           "15/22 points in panel order); the back-calculated (hybrid) column is scored on all "
           "measured points. Torque residuals in N*m. Matches tab:FlxPinPreGof, Main.tex 2026-10-19."])
for c in ws[1]:
    c.font = Font(bold=True)
wb.save(path)
print("WROTE", path)
