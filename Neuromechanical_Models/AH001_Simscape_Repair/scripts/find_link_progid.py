#!/usr/bin/env python3
import winreg

G = "{2666BDBF-5207-4731-9976-13172BEB124F}"
for sub in (None, "ProgID", "VersionIndependentProgID", "InprocServer32"):
    try:
        with winreg.OpenKey(winreg.HKEY_CLASSES_ROOT, rf"CLSID\{G}" + (f"\\{sub}" if sub else "")) as k:
            print(sub or "(Default)", "=", winreg.QueryValueEx(k, None)[0] if not sub else winreg.QueryValueEx(k, "")[0])
    except OSError as e:
        print(sub, "ERR", e)

# enumerate any subkeys
with winreg.OpenKey(winreg.HKEY_CLASSES_ROOT, rf"CLSID\{G}") as k:
    i = 0
    while True:
        try:
            print("subkey:", winreg.EnumKey(k, i))
        except OSError:
            break
        i += 1
