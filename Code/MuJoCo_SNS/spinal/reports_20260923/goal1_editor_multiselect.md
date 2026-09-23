# Goal 1 — Connectome editor: rubber-band multi-select + group move

**For Ben.** What changed, what you must decide, and what was verified.

## The one thing to read first: WHICH file was changed

The task named `connectome_editor.html`, but that file is the **rules/gains form**
(R1a…R7 checkboxes, hops, gains; exports `connectome_gains.json` — the only spec
`runner.py` consumes). It has **no node canvas, no palette, no bLoad/bClear buttons,
no node positions** — none of the things the task describes exist in it.

Every concrete anchor in the task (draggable nodes, palette, bLoad, bClear, saved
node positions) exists only in **`connectome_block_editor.html`** (the block/canvas
editor served by `/connectome-gui`). The "recently fixed palette" is even commented
in that file ("palette (was dropped in the v2.1 rewrite — restored)").

**Decision made:** the feature was implemented in
`Code/MuJoCo_SNS/spinal/connectome_block_editor.html` (v2.1 → **v2.2**).
`connectome_editor.html` was **not touched** (mtime/size unchanged), so the
runner-consumed export path has zero exposure. If you actually wanted the rules
editor changed instead, say so — but that tool has no nodes to select.

No circuit content changed: `CONNECTOME.md` untouched, `connectome_gains.json`
untouched (it doesn't even exist in `spinal/` right now — only exported when you
click Export in the rules editor), no template, type, gain or sign semantics changed.

## New interactions (the spec)

| Gesture | v2.2 behavior |
|---|---|
| **Drag on empty canvas** | Rubber-band rectangle (magenta dashed). On release, every node whose bounding box intersects the rectangle is selected (replaces the previous selection). One hit = normal single selection; ≥2 hits = group. |
| **Drag a node that is in the selection (group ≥2)** | Moves ALL selected nodes together; relative offsets preserved exactly. The group stays selected on release. |
| **Drag an unselected node** | Exactly today's behavior: only that node moves, and on release it becomes the single selection. |
| **Click empty canvas** | Clears the selection (today's behavior, kept). |
| **SHIFT+click a node** | Toggles that node in/out of the selection. |
| **Esc** | Clears selection (also cancels a rubber band in progress and an armed two-click link). |
| **Highlight** | Every selected node keeps the existing `.sel` thick-stroke highlight (the same style single selection already used). |

How positions persist (unchanged code paths): a group move mutates each node's
`n.x`/`n.y` exactly like the single-node drag always did, so **Export JSON
(`serialize()`), undo snapshots, and the localStorage tab data all see group moves**
— verified in the test suite below.

## Interaction changes you should know about (deliberate)

1. **SHIFT+click now toggles selection.** In v2.1 SHIFT+click on a node *armed the
   two-click synapse flow*. The task explicitly requires SHIFT+click = toggle, so
   the two-click link is now **ALT+click A, then click B**. SHIFT/ALT-**drag**
   node-to-node still draws a synapse directly.
2. **The synapse gestures were genuinely broken in the on-disk v2.1 and are now
   fixed** (measured, see verification): a one-gesture SHIFT-drag created **no edge
   at all** (the link branch never attached its mouse listeners; completion relied
   on a stale `linkDrag` being consumed by your *next* gesture), and the documented
   "SHIFT-click A then click B" flow created a **duplicate edge (2× A→B)**. v2.2
   makes the link gesture self-contained: one drag = one edge; two-click flow =
   exactly one edge. This also removes the accidental-edge hazard when a neighbor
   node sits within the 24 px release radius on a tight layout — which is exactly
   your use case.
3. Small extras forced by the above: a no-move SHIFT/ALT press is treated as a
   click (no junk undo entries); after a group move the follow-up click event is
   suppressed so it cannot collapse the selection.

## Known non-goals (not implemented, on purpose)

- No group **delete** — Del still deletes the primary selection only (today's
  behavior; a band-selected group has no primary, so Del is inert until you click
  a node).
- Rubber band **replaces** the selection; there is no additive shift+band.
- Group selection is nodes-only; clicking a synapse still selects just that edge.

## Inherited quirk worth knowing (pre-existing, NOT changed)

`persist()` writes the `models` array, but the active tab's `data` is only
refreshed from the live canvas when you **switch tabs** (`switchTab`). So
localStorage always lags current-tab edits until a tab switch — for single drags
too, in v2.1 and v2.2 alike. Export JSON and undo are always current. If you want
autosave tightened, that's a separate 2-line fix (`models[activeTab].data =
serialize()` inside `snapshot()`); I did not change save semantics without your OK.

## What was verified (all run in this session on EB475WS4)

- `node -v` → `v22.23.2` (node available).
- **Syntax check (the task's named check):** extracted the single inline
  `<script>` to `%TEMP%\cbe_script.js` and ran `node --check` → **exit 0, OK**.
- **Behavioral suite (beyond the required check):** ran the extracted script under
  node against a minimal DOM stub (`%TEMP%\cbe_combined.js`) — **37/37 PASS,
  exit 0**, covering: band start/draw/hit-test (circle/rect/ellipse bounds),
  follow-up-click guard, group drag (+30/+40 move, offsets preserved,
  non-selected untouched, snapshot records moved positions, selection kept on
  release), single-node drag unchanged (incl. select-on-release), SHIFT+click
  add/remove/promote-primary, Esc, ALT+click arming + completion (exactly one
  edge), one-gesture SHIFT-drag synapse, no accidental edge near a 20 px neighbor,
  undo walking back through group and single drags, **palette populated (30 types)
  + placement**, **bLoad opens file picker**, **bClear empties tab**, w2l template
  (18 nodes / 16 synapses), and tab-switch `persist()` storing moved x/y.
- **v2.1 baseline probe** (same harness, pristine pre-edit script):
  one-gesture shift+drag → `edges=0, linkDrag_still_set=true`;
  SHIFT-click A then click B → `edges=2 list=A->B,A->B` — the evidence for
  interaction change #2 above.
- **Diff re-read vs the pre-edit original** (full 255-line diff reviewed): zero
  hunks touch the palette builder, `bSave`/`bLoad`/`fileIn`, `bClear`, templates,
  `serialize`/`persist`, or tab logic.
- **Do-not-touch files:** `connectome_editor.html` (10,438 B, 09-22 12:38 AM) and
  `CONNECTOME.md` (09-22 09:18 AM) mtimes unchanged; `connectome_gains.json` does
  not exist in `spinal/` — nothing written.

## Not verified

- Not opened in a real browser (no GUI session here). The suite exercises the
  exact shipped script through a DOM stub; real-browser rendering/CSS was not
  observed. The `/connectome-gui` command (`serve_editor.py`) serves it unchanged.
- A pre-existing quirk was noticed, not fixed: after a successful synapse drag the
  follow-up canvas click clears the new edge's selection highlight (v2.1 behavior,
  left as-is; the edge itself is created fine).

*Supervisor note: this is tooling only — no literature-derived wiring claims are
made in this change.*
