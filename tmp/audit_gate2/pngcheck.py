"""PNG dimensions + DPI for the goal-2 report figures."""
import struct

for name in ["goal2_knee_spiking.png", "goal2_beer_spiking.png",
             "goal2_rgmn_spiking.png"]:
    p = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_spiking_20261002" + "\\" + name
    with open(p, "rb") as f:
        data = f.read()
    w, h = struct.unpack(">II", data[16:24])
    dpi = None
    # pHYs chunk
    idx = data.find(b"pHYs")
    if idx > 0:
        x, y, unit = struct.unpack(">IIB", data[idx + 4:idx + 13])
        if unit == 1:
            dpi = round(x * 0.0254)
    size_in = (round(w / dpi, 2), round(h / dpi, 2)) if dpi else "?"
    print(f"{name}: {w}x{h} px, dpi={dpi}, physical={size_in} in")
