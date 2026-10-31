# Superseded builders

The editable W2L model is now `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj`. The canonical reference-style/parameter record is in that same folder as `CANONICAL_REFERENCE_STYLE_RECORD.md`.

The following helpers describe earlier generated diagnostics or the rejected first layout. They are **disabled by default**, before file creation, through `historical_build_guard.py`:

- `build_guided_native.py`
- `build_project_authoritative_si.py`
- `build_newer_v3_diagnostic.py`
- `rebuild_subnetwork_style.py`

They can recreate retired model folders or overwrite later canonical work. An explicitly requested historical reconstruction must supply `--historical-regeneration`. That flag identifies a deliberate reconstruction; it is not a routine way to run or update the current model. The positive native and MuJoCo results from the historical variants remain recorded as historical results and are not reassigned to a changed canonical controller.

Read-only inspection, rendering, and verification helpers remain available. New canonical edits must start from the current canonical project and retain its later bilateral RG, contacts, and independent homonymous type-II additions. Do not restore the archived 96-neuron checkpoint over later work.
