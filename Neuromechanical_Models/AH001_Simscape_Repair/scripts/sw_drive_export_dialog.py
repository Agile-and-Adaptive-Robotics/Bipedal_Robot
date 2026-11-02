#!/usr/bin/env python3
"""Drive the Simscape Multibody Link export dialog inside SLDWORKS (PID-scoped)."""
import time

import pywinauto
from pywinauto import Application, Desktop

PID = 26280

app = Application(backend="uia").connect(process=PID)
wins = [w for w in app.windows()]
for w in wins:
    print("WIN:", w.element_info.name, "|", w.element_info.className)

# find the export dialog
target = None
for w in wins:
    nm = (w.element_info.name or "").lower()
    cls = (w.element_info.className or "").lower()
    if "simscape" in nm or "multibody" in nm or "export" in nm or "3d" in cls:
        target = w
        break

if target is None:
    # fallback: top-level dialog that is not the main frame
    for w in wins:
        if w.element_info.className not in ("Qt5152QWindowIcon", "Afx:00400000:0", "wndclass_desked_gsk"):
            print("fallback candidate:", w.element_info.name, w.element_info.className)
            target = w
            break

print("target dialog:", target.element_info.name if target else None)
if target is None:
    raise SystemExit("no dialog found")

target.set_focus()
time.sleep(0.5)


def dump(el, depth=0, maxdepth=4):
    try:
        ei = el.element_info
        print("  " * depth + f"[{ei.control_type}] {ei.name!r} auto={ei.automation_id}")
    except Exception:
        return
    if depth >= maxdepth:
        return
    try:
        for c in el.children():
            dump(c, depth + 1, maxdepth)
    except Exception:
        pass


dump(target.wrapper_object(), 0, 3)
