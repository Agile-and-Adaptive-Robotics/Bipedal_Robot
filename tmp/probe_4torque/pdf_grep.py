# Extract text of the old knee-torque draft PDF and print lines around
# RMSE/FVU/table mentions (read-only).
import fitz
P = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\Torque_Values_About_Robotic_Knee_Using_Braided_Pneumatic_Actuators.pdf"
doc = fitz.open(P)
print("pages:", doc.page_count)
import re
for i, page in enumerate(doc):
    t = page.get_text()
    if re.search(r"RMSE|FVU|Table", t):
        print(f"----- page {i+1} -----")
        print(t)
