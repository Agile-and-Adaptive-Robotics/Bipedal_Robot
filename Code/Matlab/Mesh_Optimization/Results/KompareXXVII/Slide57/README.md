# Defense slide 57 figures

- `robot_routes_2019.png` is the existing right-side robot route view reused from slide 14 of v33a.
- `route_schematic.png` and `route_schematic.pdf` are the replacement right-hand panel over the stock OpenSim Gait2392 skeleton. The erector spinae uses the compliant rerun and is drawn as a single blue route. Psoas and iliacus use the same blue.
- `actuator_site_table.png` and `actuator_site_table.pdf` list the 27 sites. Erector spinae has four BPAs and one line part per cable. The other counts come from the saved survey.
- `route_schematic_original.png` preserves the image embedded on slide 57 of v34.
- `route_schematic_sticks.png` and `route_schematic_sticks.pdf` preserve the previous panel with pivot lines.
- `route_projection.mat` contains the saved survey routes. The psoas is registered using the stock OpenSim neutral torso-to-pelvis transform. Its registered guide coincides with the iliacus guide, and the drawing shows their common distal connection at the saved iliacus attachment.

The erector-spinae rerun starts from the stock OpenSim attachments and compares 1–4 nominal 20 mm BPAs at direct, 2:1, and 3:1 travel ratios. The compliant result uses four direct BPAs. The rigid result uses two BPAs with one moving sheave per cable. Both selected cases maintain tension over the 31-point lumbar-extension check and fall short of the full signed torque target. Drawn sheave positions are schematic.

The replacement panel uses 72 BPAs across 27 sites. Two sites use one, six use two, 18 use three, and one uses four. The saved survey rows use finite bracket stiffness.

The femur and tibia use the rectus femoris survey's saved display transforms. The talus and calcaneus use the soleus survey transforms, and the toes and torso use their stock neutral relative transforms. The routes retain each site's saved display pose, so the overlay serves as a route schematic.

Run `python route_visual.py --schematic-only` with OpenSim, NumPy, SciPy, Matplotlib, and Arial to regenerate the panel. Pass `--geometry` if the OpenSim Geometry folder is outside the default installation. MATLAB checkpoints and the selected records are in `erector_rerun`. See the parent `REPORT.md` for the survey geometry and packaging limitations.
