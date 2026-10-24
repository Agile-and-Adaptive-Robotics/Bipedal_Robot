import fitz

doc = fitz.open("Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.pdf")
for pno in range(29, 34):
    text = doc[pno].get_text()
    if "gives the aggregate" in text or "Table ??" in text:
        i = max(text.find("gives the aggregate"), text.find("Table ??"))
        print(f"page {pno+1}: ...{text[max(0, i-60):i+40]!r}...")
        break
else:
    print("reference sentence not found on pages 30-33")
print("total pages:", doc.page_count)
