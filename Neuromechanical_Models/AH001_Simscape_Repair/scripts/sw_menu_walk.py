#!/usr/bin/env python3
"""Alt-walk the SOLIDWORKS menu bar, dumping each popup menu's items."""
import time

from pywinauto import Application
from pywinauto.keyboard import send_keys

PID = 26280
app = Application(backend="uia").connect(process=PID)
main = app.top_window()
main.set_focus()
time.sleep(0.4)


def menu_snapshot():
    """Return list of visible MenuItem texts across the process windows."""
    items = []
    try:
        app2 = Application(backend="uia").connect(process=PID)
        for w in app2.windows():
            try:
                for it in w.descendants(control_type="MenuItem"):
                    t = it.element_info.name
                    if t:
                        items.append(t)
            except Exception:
                pass
    except Exception:
        pass
    return items


send_keys("%")           # Alt: activates the menu bar
time.sleep(0.5)

for step in range(24):
    items = menu_snapshot()
    label = " | ".join(items[:14])
    print(f"step {step:02d}: {label[:200]}")
    if any("simscape" in t.lower() or "multibody" in t.lower() for t in items):
        print(">>> Multibody Link menu is open")
        for t in items:
            print("   item:", t)
        break
    send_keys("{RIGHT}")
    time.sleep(0.35)

send_keys("{ESC}", pause=0.1)
send_keys("{ESC}", pause=0.1)
