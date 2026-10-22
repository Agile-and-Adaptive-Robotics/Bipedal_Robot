import re

import fitz

doc = fitz.open(r"D:\temp\diss_build_live_20261022\main.pdf")
t = " ".join(doc[i].get_text() for i in range(130, 155))
t = re.sub(r"\s+", " ", t)
for phrase in (
    "unique identifier wherever a record carries one",
    "primary form is the packaged explorer",
    "the papers the prune archived from the working base",
    "943 papers that the explorer carries in full",
):
    print(phrase in t, "::", phrase)
