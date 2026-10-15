"""Extract code cells from SNS-Toolbox tutorial notebooks and write them
as runnable scripts (usage: _extract_tutorials.py [numbers...]).
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).parent
TUT = HERE / "sns_tutorials"

HEADER = (
    "# Auto-extracted from SNS-Toolbox tutorial notebook (cells in order).\n"
    "# matplotlib .show() calls removed; render calls save to file instead.\n"
    "import matplotlib\nmatplotlib.use('Agg')\n"
    "import matplotlib.pyplot as plt\n"
    "def _save_current_fig(name):\n"
    "    plt.gcf().savefig(r'%s' / name, dpi=130,\n"
    "                      bbox_inches='tight')\n" % str(HERE / "sns_tutorials")
)


def extract(num: int):
    src = None
    for p in TUT.glob(f"Tutorial_{num}__*.ipynb"):
        src = p
    if src is None:
        print(f"tutorial {num}: NOT FOUND")
        return
    nb = json.loads(src.read_text(encoding="utf-8"))
    cells = []
    for c in nb["cells"]:
        if c["cell_type"] != "code":
            continue
        lines = []
        for ln in c["source"]:
            if "%matplotlib" in ln or "plt.show()" in ln:
                continue
            if ".show()" in ln and "plt" in ln:
                continue
            if "renderer.render(" in ln and "save" not in ln and "view" in ln:
                continue
            lines.append(ln)
        cells.append("".join(lines).rstrip())
    body = "\n\n# --- cell ---\n".join(cells)
    out = TUT / f"tutorial{num}_run.py"
    out.write_text(HEADER + "\n\n# --- cell ---\n" + body + "\n",
                   encoding="utf-8")
    print(f"tutorial {num}: {len(cells)} code cells -> {out.name}")


if __name__ == "__main__":
    nums = [int(a) for a in sys.argv[1:]] or [1, 2, 4, 8]
    for n in nums:
        extract(n)
