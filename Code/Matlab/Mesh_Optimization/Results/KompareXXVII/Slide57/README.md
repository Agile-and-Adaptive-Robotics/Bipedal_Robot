# Defense slide 57 figures

- `robot_routes_2019.png` is the existing right-side robot route view reused from slide 14 of v33a.
- `route_schematic.png` is the right panel embedded on slide 57 of v34. Colors distinguish one direct BPA, parallel BPAs, and parallel BPAs with tackle.
- `route_projection.mat` contains the recorded routes and composed transforms used for that panel. The composition follows the parent folder's `make_table412_figures_20261031.m`. The psoas is drawn separately in its local frame.

From this folder, run `python route_visual.py` with NumPy, SciPy, Matplotlib, and Arial installed. The script reads the parent folder's `route_configurations.csv` and writes `route_schematic.png`. It uses the saved selected configurations and does not run an optimization. The selected 27 sites use 69 BPAs. Tackle line parts describe the modeled transmission ratio. See the parent `REPORT.md` for geometry and packaging limitations.
