#!/usr/bin/env python3
"""Drive SOLIDWORKS menu: Simscape Multibody Link > Export as > Simscape Multibody,
then fill the export dialog (filename + STEP) and OK it. PID-scoped pywinauto."""
import os
import time

from pywinauto import Application, keyboard

PID = 26280
EXPORT_DIR = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair\cad\export"
XML = os.path.join(EXPORT_DIR, "ah001_full.xml")
os.makedirs(EXPORT_DIR, exist_ok=True)

app = Application(backend="uia").connect(process=PID)
main = app.top_window()
main.set_focus()
time.sleep(0.5)
print("focused:", main.element_info.name)


def find_menu_item(win, texts, timeout=10):
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            it = win.descendants(title_text=texts, control_type="MenuItem")
            if it:
                return it[0]
        except Exception:
            pass
        time.sleep(0.4)
    return None


# open the add-in menu
m = find_menu_item(main, ["Simscape Multibody Link"])
if m is None:
    # try via menu bar traversal
    mbs = main.descendants(control_type="MenuBar")
    print("menu bars:", [mb.element_info.name for mb in mbs])
    raise SystemExit("menu 'Simscape Multibody Link' not found")
print("found menu:", m.element_info.name)
m.select()          # opens the dropdown
time.sleep(0.8)

sub = find_menu_item(main, ["Export as"])
print("Export as:", sub)
if sub is not None:
    sub.select()
    time.sleep(0.8)
    leaf = find_menu_item(main, ["Simscape Multibody"])
    print("leaf:", leaf)
    if leaf is not None:
        leaf.select()

time.sleep(1.5)
# enumerate fresh top-level windows: the export dialog should be #32770 or similar
app2 = Application(backend="uia").connect(process=PID)
wins = app2.windows()
print("windows now:")
for w in wins:
    print("  ", repr(w.element_info.name), "|", w.element_info.class_name)
