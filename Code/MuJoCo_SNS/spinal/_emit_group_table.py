"""Emit the muscle -> functional-group table (LaTeX + CSV) LIVE from
muscle_map.py, and append it to the dissertation draft section.
"""
import csv
import io
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

from muscle_map import GROUPS, _GROUPS_BY_NAME

HERE = Path(__file__).parent
TEX = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers"
           r"\Dissertation\CPG_spinal_section_draft.tex")
GROUP_LABEL = {"hip_ext": "HIP-E", "hip_flex": "HIP-F", "hip_abd": "HIP-AB",
               "hip_add": "HIP-AD", "knee_ext": "KNEE-E",
               "knee_flex": "KNEE-F", "ankle_pf": "ANK-PF",
               "ankle_df": "ANK-DF", "trunk_ext": "TRK-E",
               "trunk_flex": "TRK-F"}
ORDER = ("hip_ext", "hip_flex", "hip_abd", "hip_add", "knee_ext",
         "knee_flex", "ankle_pf", "ankle_df", "trunk_ext", "trunk_flex")


def main():
    prim = {g: sorted(b for b, gs in _GROUPS_BY_NAME.items() if gs[0] == g)
            for g in GROUPS}
    sec = {g: sorted(b for b, gs in _GROUPS_BY_NAME.items() if g in gs[1:])
           for g in GROUPS}
    # CSV (machine-readable companion)
    with open(HERE / "muscle_groups.csv", "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["group", "role", "muscles"])
        for g in ORDER:
            w.writerow([GROUP_LABEL[g], "primary", ", ".join(prim[g])])
            if sec[g]:
                w.writerow([GROUP_LABEL[g], "secondary", ", ".join(sec[g])])
    # LaTeX
    rows = []
    for g in ORDER:
        p = ", ".join(f"\\textit{{{b}}}" for b in prim[g])
        s = (", ".join(sec[g])) if sec[g] else "---"
        rows.append(f"    {GROUP_LABEL[g]} & {len(prim[g])} & {p} & {s} \\\\")
    tex_block = (
        "\n\\begin{table}[t]\n  \\centering\n"
        "  \\caption{Functional muscle groups used by the spinal circuit "
        "(live from \\texttt{muscle\\_map.py}; $46$ muscles per side, $92$ "
        "total). Primary membership receives the full PF/posture weight; "
        "secondary (biarticular) membership receives half.}\n"
        "  \\label{tab:muscle-groups}\n  \\small\n"
        "  \\begin{tabular}{lc p{4.6cm} p{3.4cm}}\n    \\toprule\n"
        "    Group & $n$ & Primary muscles & Secondary muscles \\\\\n"
        "    \\midrule\n" + "\n".join(rows) + "\n    \\bottomrule\n"
        "  \\end{tabular}\n\\end{table}\n")
    txt = TEX.read_text(encoding="utf-8")
    anchor = "\\end{table}\n\n\\subsection{Sensorimotor Refinement"
    if anchor in txt:
        txt = txt.replace(anchor, "\\end{table}\n" + tex_block +
                          "\n\\subsection{Sensorimotor Refinement", 1)
        note = ("Table~\\ref{tab:muscle-groups} lists the functional "
                "muscle groups and their member muscles.\n\n\\subsection"
                "{Sensorimotor Refinement")
        txt = txt.replace("\n\\subsection{Sensorimotor Refinement",
                          "\n" + note, 1)
    else:
        txt += tex_block
    TEX.write_text(txt, encoding="utf-8")
    n_muscles = len(_GROUPS_BY_NAME)
    print(f"table written into {TEX.name}; {n_muscles} muscles, "
          f"muscle_groups.csv written")


if __name__ == "__main__":
    main()
