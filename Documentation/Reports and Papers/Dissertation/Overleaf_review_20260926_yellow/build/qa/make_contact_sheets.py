from pathlib import Path
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parent
GROUPS = {
    "qa_frontmatter_org_steele.png": [
        "pages_3_4-003.png", "pages_3_4-004.png",
        "pages_22_24-022.png", "pages_22_24-023.png", "pages_22_24-024.png",
        "pages_31_32-031.png", "pages_31_32-032.png",
    ],
    "qa_methods.png": [
        *[f"pages_50_58-{page:03d}.png" for page in range(50, 59)],
        "pages_73_73-073.png", "pages_96_96-096.png",
    ],
    "qa_discussion_future_conclusion.png": [
        *[f"pages_112_122-{page:03d}.png" for page in range(112, 123)],
    ],
}

THUMB = (408, 528)
LABEL = 30
COLS = 3

for output_name, names in GROUPS.items():
    rows = (len(names) + COLS - 1) // COLS
    sheet = Image.new("RGB", (COLS * THUMB[0], rows * (THUMB[1] + LABEL)), "white")
    draw = ImageDraw.Draw(sheet)
    for index, name in enumerate(names):
        image = Image.open(ROOT / name).convert("RGB")
        image.thumbnail(THUMB)
        x = (index % COLS) * THUMB[0] + (THUMB[0] - image.width) // 2
        y0 = (index // COLS) * (THUMB[1] + LABEL)
        y = y0 + LABEL + (THUMB[1] - image.height) // 2
        sheet.paste(image, (x, y))
        draw.text((x + 6, y0 + 7), name.rsplit("-", 1)[-1].removesuffix(".png"), fill="black")
    sheet.save(ROOT / output_name)
