"""Download/locate a Zotero attachment file path + copy the ft-cache text.

Usage: _zot_file.py <attachment_key> <dest_path>
Uses the Zotero storage directory directly (PDF not synced -> .zotero-ft-cache).
"""
import os
import shutil
import sys

STOR = r"C:\Users\Ben Bolen\Zotero\storage"


def main():
    key, dest = sys.argv[1], sys.argv[2]
    folder = os.path.join(STOR, key)
    if not os.path.isdir(folder):
        print(f"no storage folder for {key}")
        return 1
    files = os.listdir(folder)
    print("files:", files)
    pdf = [f for f in files if f.lower().endswith(".pdf")]
    if pdf:
        shutil.copy(os.path.join(folder, pdf[0]), dest)
        print("copied PDF ->", dest)
        return 0
    cache = os.path.join(folder, ".zotero-ft-cache")
    if os.path.isfile(cache):
        shutil.copy(cache, dest)
        print("copied ft-cache ->", dest)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
