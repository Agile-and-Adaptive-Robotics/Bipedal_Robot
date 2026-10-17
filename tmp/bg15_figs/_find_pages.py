"""Map figure captions to PDF page numbers by per-page txtwrite extraction."""
import re
import subprocess
import sys
import tempfile

MGS = r"D:\MiKTeX\miktex\bin\x64\mgs.exe"


def page_text(pdf, page):
    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
        out = f.name
    subprocess.run([MGS, "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=txtwrite",
                    f"-dFirstPage={page}", f"-dLastPage={page}",
                    f"-sOutputFile={out}", pdf],
                   check=True, capture_output=True)
    t = open(out, encoding="utf-8", errors="replace").read()
    return t


def npages(pdf):
    r = subprocess.run([MGS, "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=txtwrite",
                        "-dFirstPage=99999", "-dLastPage=99999",
                        "-sOutputFile=NUL", pdf],
                       capture_output=True, text=True)
    m = re.search(r"greater than the number of pages in the file: (\d+)",
                  r.stderr + r.stdout)
    return int(m.group(1)) if m else 0


def run(pdf, patterns):
    n = npages(pdf)
    print(f"=== {pdf} ({n} pages) ===")
    for p in range(1, n + 1):
        t = page_text(pdf, p)
        flat = re.sub(r"\s+", " ", t)
        for pat in patterns:
            for m in re.finditer(pat, flat):
                print(f"  p{p}: {m.group(0)[:130]}")


if __name__ == "__main__":
    which = sys.argv[1]
    if which == "deng":
        run(r"D:\Github\Bipedal_Robot\tmp\bg15_figs\deng_2019.pdf",
            [r"Figure \d+\. ?\(?[A-Z]?\)?[^F]{0,110}"])
    elif which == "shinohara":
        run(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
            r"\shinohara_2025_biophiv687930.pdf",
            [r"Figure \d+: ?[^F]{0,110}"])
    elif which == "dirusso":
        run(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\_dirusso_2023.pdf",
            [r"[Ff]igure\s*\d+\s*[.:][^F]{0,110}"])
    elif which == "masters":
        run(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers"
            r"\Bolen_MSME\Bolen Final report.pdf",
            [r"Figure \d+[:.] ?[^F]{0,110}"])
