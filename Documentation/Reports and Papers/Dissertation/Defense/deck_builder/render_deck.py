import os, sys
import win32com.client

PPTX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Bolen_Dissertation_Defense.pptx")
if len(sys.argv) > 1:
    PPTX = os.path.abspath(sys.argv[1])  # optional override, e.g. while PowerPoint holds the default open
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "deck_render")
os.makedirs(OUT, exist_ok=True)
for f in os.listdir(OUT):
    os.remove(os.path.join(OUT, f))

# attach to a running PowerPoint if there is one (never Quit() an app the user has
# open, e.g. Ben studying the official A-E template); only spawn-and-quit when none
try:
    app = win32com.client.GetActiveObject("PowerPoint.Application")
    spawned = False
except Exception:
    app = win32com.client.Dispatch("PowerPoint.Application")
    spawned = True
pres = app.Presentations.Open(PPTX, ReadOnly=True, Untitled=False, WithWindow=False)
n = pres.Slides.Count
for i in range(1, n + 1):
    pres.Slides(i).Export(os.path.join(OUT, f"slide_{i:02d}.png"), "PNG", 1600, 900)
pres.Close()
if spawned:
    app.Quit()
sizes = sorted(os.listdir(OUT))
print("exported", len(sizes), "slides")
small = [f for f in sizes if os.path.getsize(os.path.join(OUT, f)) < 15000]
print("suspiciously small:", small if small else "none")
