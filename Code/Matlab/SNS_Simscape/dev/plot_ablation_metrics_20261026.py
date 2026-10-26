# Ablation-pair metrics figure from the two saved JSONs (traces were not
# persisted; this plots the saved metrics honestly as metrics).
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results"
a = json.load(open(R + r"\ablation_intact_20261026.json"))
b = json.load(open(R + r"\ablation_ablated_20261026.json"))
fig, ax = plt.subplots(1, 3, figsize=(9, 3.2))
for k, (key, ttl) in enumerate({
        "cycL": "cycles, left leg", "cycR": "cycles, right leg",
        "fwd": "forward progress (m)"}.items()):
    ax[k].bar(["intact", "ablated"], [a[key], b[key]],
              color=["#0072B2", "#D55E00"])
    ax[k].set_title(ttl)
    ax[k].grid(alpha=0.3, axis="y")
    for i, v in enumerate([a[key], b[key]]):
        ax[k].text(i, v, f"{v:.2f}", ha="center", va="bottom")
fig.suptitle("Feedback ablation, frozen tuning (same drives, same "
             "conductances): intact vs contact+Ib ablated, 30 s")
fig.tight_layout()
fig.savefig(R + r"\pictures\v4\ablation_pair_metrics_20261026.png", dpi=220)
fig.savefig(R + r"\pictures\v4\ablation_pair_metrics_20261026.pdf")
print("wrote ablation_pair_metrics_20261026.png/.pdf")
