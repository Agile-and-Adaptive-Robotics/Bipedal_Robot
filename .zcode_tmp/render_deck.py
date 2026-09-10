import os, sys
import win32com.client

PPTX = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Bolen_Dissertation_Defense_20260910.pptx"
OUT = r"D:\GitHub\Bipedal_Robot\.zcode_tmp\deck_render"
os.makedirs(OUT, exist_ok=True)
for f in os.listdir(OUT):
    os.remove(os.path.join(OUT, f))

app = win32com.client.Dispatch("PowerPoint.Application")
pres = app.Presentations.Open(PPTX, ReadOnly=True, Untitled=False, WithWindow=False)
n = pres.Slides.Count
for i in range(1, n + 1):
    pres.Slides(i).Export(os.path.join(OUT, f"slide_{i:02d}.png"), "PNG", 1600, 900)
pres.Close()
app.Quit()
sizes = sorted(os.listdir(OUT))
print("exported", len(sizes), "slides")
small = [f for f in sizes if os.path.getsize(os.path.join(OUT, f)) < 15000]
print("suspiciously small:", small if small else "none")
