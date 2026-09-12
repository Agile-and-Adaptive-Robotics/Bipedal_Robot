import os, sys
sys.path.insert(0, os.path.join(os.environ['TEMP'], 'pdfx'))
import pypdf

src = r"C:\Users\Ben Bolen\Zotero\storage\47QYMID7\Thomas Brown et al. - 1911 - The intrinsic factors in the act of progression in.pdf"
out = r"D:\Github\Bipedal_Robot\SAD_audit\pilot_ground_brown1911_pdf.txt"
r = pypdf.PdfReader(src)
chunks = []
for i, p in enumerate(r.pages):
    t = p.extract_text() or ''
    chunks.append(f"--- page {i+1} ---\n{t}")
text = "\n".join(chunks)
with open(out, 'w', encoding='utf-8') as f:
    f.write(text)
print("pages:", len(r.pages), "chars:", len(text))
