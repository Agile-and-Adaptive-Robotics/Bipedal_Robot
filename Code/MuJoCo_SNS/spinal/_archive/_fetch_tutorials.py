"""Fetch the SNS-Toolbox tutorial notebooks from GitHub master."""
import urllib.request
from pathlib import Path
from urllib.parse import quote

OUT = Path(__file__).parent / "sns_tutorials"
OUT.mkdir(exist_ok=True)
BASE = "https://raw.githubusercontent.com/williamnourse/SNS-Toolbox/master/tutorials/"
NAMES = [
    "Tutorial 1 - Network Design.ipynb",
    "Tutorial 2 - Simulation.ipynb",
    "Tutorial 3 - Spiking Networks.ipynb",
    "Tutorial 4 - Subnetworks.ipynb",
    "Tutorial 5 - Spiking Transmission Delay.ipynb",
    "Tutorial 6 - Connectivity Patterns.ipynb",
    "Tutorial 7 - Electrical Synapses.ipynb",
    "Tutorial 8 - Voltage-gated Ion Channels.ipynb",
    "Tutorial 9 - Advanced Spiking.ipynb",
]
for nm in NAMES:
    url = BASE + quote(nm)
    dest = OUT / nm.replace(" ", "_").replace("-", "")
    with urllib.request.urlopen(url, timeout=40) as r:
        data = r.read()
    dest.write_bytes(data)
    print(f"{nm}: {len(data)} bytes -> {dest.name}")
