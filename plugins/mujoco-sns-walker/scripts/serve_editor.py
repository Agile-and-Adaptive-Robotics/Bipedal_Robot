"""Serve the spinal folder over http://localhost so both connectome editors
open without file:// quirks (the block editor's w2l template does a sync XHR
to a relative path - reliable as a same-origin GET, browser-dependent on
file://). localStorage works the same or better on localhost.

Usage:
  serve_editor.py                 # port 8765, opens both editors, serves forever
  serve_editor.py --port 8971
  serve_editor.py --no-open       # just print URLs (agent runs me in background)
Stop: Ctrl+C (or kill the background task / close the terminal).
"""
from __future__ import annotations

import argparse
import functools
import http.server
import socket
import sys
import threading
import webbrowser
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import wcommon  # noqa: E402


def main() -> int:
    wcommon.utf8_stdio()
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--no-open", action="store_true")
    args = ap.parse_args()

    spinal = wcommon.resolve_spinal()
    handler = functools.partial(
        http.server.SimpleHTTPRequestHandler, directory=str(spinal))
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", args.port), handler)
    url_rules = f"http://127.0.0.1:{args.port}/connectome_editor.html"
    url_blocks = f"http://127.0.0.1:{args.port}/connectome_block_editor.html"
    wcommon.banner(f"serving {spinal}")
    print(f"rules editor  : {url_rules}   (exports connectome_gains.json - what runner reads)")
    print(f"block editor  : {url_blocks}   (design tool; export is NOT consumed by runner)")
    print("stop: Ctrl+C / kill this process")

    if not args.no_open:
        threading.Timer(0.6, lambda: webbrowser.open(url_rules)).start()
        threading.Timer(1.1, lambda: webbrowser.open(url_blocks)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
