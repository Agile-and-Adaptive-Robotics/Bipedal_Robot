#!/usr/bin/env python3
"""Find the Simscape Multibody Link SW add-in GUID + DLL path."""
import subprocess
import winreg

root = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"SOFTWARE\SolidWorks\AddIns")
i = 0
link_guid = None
while True:
    try:
        guid = winreg.EnumKey(root, i)
    except OSError:
        break
    i += 1
    try:
        with winreg.OpenKey(root, guid) as k:
            title, _ = winreg.QueryValueEx(k, "Title")
            desc = ""
            try:
                desc, _ = winreg.QueryValueEx(k, "Description")
            except OSError:
                pass
            print(guid, "|", title, "|", desc)
            if "multibody link" in str(title).lower() or "multibody link" in str(desc).lower():
                link_guid = guid
    except OSError as e:
        print(guid, "ERR", e)

print("\nlink GUID:", link_guid)
if link_guid:
    for hive, hw in ((winreg.HKEY_CLASSES_ROOT, None),):
        try:
            with winreg.OpenKey(hive, rf"CLSID\{link_guid}\InprocServer32") as k:
                dll, _ = winreg.QueryValueEx(k, None)
                print("DLL:", dll)
        except OSError as e:
            print("InprocServer32 ERR", e)
