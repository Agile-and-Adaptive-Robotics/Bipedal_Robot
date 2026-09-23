---
description: Vision analysis of connectomes - reference paper figures, our draw_circuit outputs, and Ben's correction markups - using the zai-mcp-server vision tools
---

# /circuit-vision

Use the MCP vision tools (`mcp__zai-mcp-server__understand_technical_diagram`,
`analyze_image`, `ui_diff_check`) to actually SEE circuits, in three modes. This is the
standard way to compare a published connectome against ours, audit our own figures for
missing elements, and read Ben's marked-up corrections.

## Mode A - reference paper figure

1. If the figure is in a PDF, render the page to PNG first (vision tools take images):
   ```bat
   D:\MiKTeX\miktex\bin\x64\mgs.exe -dNOPAUSE -dBATCH -dSAFER -sDEVICE=png16m -r200 ^
     -dFirstPage=N -dLastPage=N -sOutputFile=<out>.png <paper>.pdf
   ```
   (MiKTeX path on EB475WS4; easteregg2 uses `D:\Anaconda\envs\gs\Library\bin\gswin64c.exe`.
   Full texts also live in `spinal\` as `*_fulltext.txt` / `lit_*.txt` for the text side.)
2. `understand_technical_diagram` on the PNG with a prompt that asks for an explicit
   EDGE LIST: every connection as `source -> target, sign (exc/inh), where drawn`,
   plus neuron classes/populations and any labels.
3. Cross-check against our compiled edge inventory:
   ```bat
   "%CLAUDE_PLUGIN_ROOT%\scripts\walker.cmd" status.py   (any invocation loads spinal)
   ... and run  <spinal>\_net_edges.py  under the same python for the compiled edge dump
   ```
4. Diff the two lists: pathways present in the paper but absent in ours (or drawn but
   not implemented) are the deliverable - cite figure panel + paper table for each.

Reference targets in scope (see `spinal\replication\` and DESIGN.md): Di Russo 2023
rules 1-5, Deng/Nourse Table A6, Shevtsova 2026 laminar, Shinohara 2025, W2L walker
(`w2l_equivalent_draft.json` for the crossed contact->contra-flexor pattern that the
frozen-left-leg diagnosis points at).

## Mode B - audit our own figure

`understand_technical_diagram` on a `draw_circuit.py` output (circuit_full_vclasses,
circuit_dengstyle, in `Dissertation\CPG_airstepping_figs\` and `spinal\figures\`).
Ask: list every drawn element + every drawn edge, then compare with `_net_edges.py`
output and the figure's own alt-text file. Report: missing elements, unlabeled
elements, sign errors (open triangle vs filled circle), legend vs drawing mismatches.
This is a second pair of eyes BEFORE Ben's visual pass, not a replacement for it.

## Mode C - Ben's markup

`analyze_image` on his screenshot/photo of a corrected figure. Prompt for: each
annotation transcribed verbatim, which element it points at, and the implied edit.
Then produce an actionable edit list; CROSS-CHECK each item against the compiled
network before proposing changes. Wiring/connectome changes are Ben's call - present
the list, do not implement unilaterally.

## Output

Write findings to `<spinal>\vision_review_<topic>_<YYYYMMDD>.md` with: the image
path(s) used, extracted edge/element lists, the diff table, and open questions.
Cite panel/table for every claimed paper pathway. Follow the project figure standards
when any regenerated figure is produced.
