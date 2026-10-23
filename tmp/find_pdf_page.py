import re
import sys

import fitz  # PyMuPDF, if available

pdf = fitz.open(sys.argv[1])
needle = sys.argv[2]
for i, page in enumerate(pdf):
    t = page.get_text()
    if needle in t:
        print("PDF page", i + 1, "::", needle[:60])
