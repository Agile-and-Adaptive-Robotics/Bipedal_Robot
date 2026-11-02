from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

W = Path(__file__).resolve().parent
source = W / 'build/render'
out = W / 'build/contact'
out.mkdir(parents=True, exist_ok=True)
slides = sorted(source.glob('slide_*.png'))
cols, rows = 4, 3
thumb_w, thumb_h = 400, 225
label_h = 28
font = ImageFont.truetype('C:/Windows/Fonts/arial.ttf', 18)
for page, start in enumerate(range(0, len(slides), cols * rows), 1):
    batch = slides[start:start + cols * rows]
    sheet = Image.new('RGB', (cols * thumb_w, rows * (thumb_h + label_h)), 'white')
    draw = ImageDraw.Draw(sheet)
    for i, path in enumerate(batch):
        image = Image.open(path).convert('RGB').resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        x = (i % cols) * thumb_w
        y = (i // cols) * (thumb_h + label_h)
        sheet.paste(image, (x, y))
        draw.rectangle((x, y + thumb_h, x + thumb_w, y + thumb_h + label_h), fill='#173D24')
        draw.text((x + 8, y + thumb_h + 4), path.stem.replace('slide_', 'Slide '), font=font, fill='white')
    sheet.save(out / f'contact_{page:02d}.jpg', quality=92)
print(f'{len(slides)} slides -> {page} contact sheets')
