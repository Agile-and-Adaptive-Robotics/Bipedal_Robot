from pathlib import Path

from pypdf import PdfReader


pdf_path = Path(
    r"D:\Github\Bipedal_Robot\output\pdf\Bolen_Dissertation_proposed_edits_yellow.pdf"
)
reader = PdfReader(pdf_path)
print(f"pages\t{len(reader.pages)}")
for page_number, page in enumerate(reader.pages, 1):
    for annotation_ref in page.get("/Annots") or []:
        annotation = annotation_ref.get_object()
        contents = annotation.get("/Contents")
        if contents:
            normalized = str(contents).replace("\n", " | ")
            print(f"{page_number}\t{annotation.get('/Subtype')}\t{normalized}")
