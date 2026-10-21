import pathlib

src = pathlib.Path(r"D:/GitHub/Bipedal_Robot/tmp/paper_bridge/Opt_run_BiPulley_gain2_20261020.m")
root = src.parent

configs = [
    ("2BPA_G2", 2, 2),
    ("2BPA_G4", 2, 4),
    ("4BPA_G4", 4, 4),
]

text = src.read_text(encoding="utf-8")
assert "tackleLineParts', [2 1]" in text

for tag, nbpa, g in configs:
    t = text
    t = t.replace(
        "struct('tackleLineParts', [2 1], 'routingMode', 'moving_via'))",
        f"struct('nPulleyBPA', {nbpa}, 'tackleLineParts', [{g} 1], 'routingMode', 'moving_via'))",
    )
    t = t.replace(
        "BiPulley_checkpoint_gain2_20261020.mat",
        f"BiPulley_checkpoint_doctrine_{tag}_20261021.mat",
    )
    t = t.replace(
        "BiPulley_batch_summary_gain2_20261020.csv",
        f"BiPulley_batch_summary_doctrine_{tag}_20261021.csv",
    )
    t = t.replace(
        "sprintf('BiPulley_%s_%s_gain2_%s.mat'",
        f"sprintf('BiPulley_%s_%s_{tag}_%s.mat'",
    )
    out = root / f"Opt_run_BiPulley_doctrine_{tag}_20261021.m"
    out.write_text(t, encoding="utf-8", newline="\n")
    print(f"wrote {out.name}: nPulleyBPA={nbpa} tackleLineParts=[{g} 1]")
print("DONE")
