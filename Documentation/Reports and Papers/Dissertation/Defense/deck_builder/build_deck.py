"""Build the Bolen dissertation defense deck."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image
from deck_lib import new_deck, FIG
import slides_a, slides_b, slides_c, slides_d

# crop FlxPin_group 4-panel figure into A/B (left) and C/D (right) halves
src = os.path.join(FIG, "FlxPin_group.png")
im = Image.open(src)
w, h = im.size
half = w // 2
im.crop((0, 0, half - 60, h)).save(os.path.join(FIG, "flxpin_AB.png"))
im.crop((half + 60, 0, w, h)).save(os.path.join(FIG, "flxpin_CD.png"))

prs = new_deck()
n1 = slides_a.build(prs)
n2 = slides_b.build(prs, n1 + 1)
n3 = slides_c.build(prs, n2 + 1)
n4 = slides_d.build(prs, n3 + 1)

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Bolen_Dissertation_Defense.pptx")
prs.core_properties.title = "Movement and Control of Biomimetic Humanoid Robots"
prs.core_properties.author = "Ben P. Bolen"
prs.save(OUT)

import json
ranges = {"a": [1, n1], "b": [n1 + 1, n2], "c": [n2 + 1, n3], "d": [n3 + 1, n4]}
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "deck_ranges.json"), "w") as f:
    json.dump(ranges, f)
print("slides:", n4, "->", OUT)
