import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
p = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_block_editor.html"
lines = open(p, encoding="utf-8").read().splitlines()
for i, l in enumerate(lines):
    if ("TPL_JSON" in l or "const TYPES" in l or "stampWalkerGrps" in l
            or "TPL-EMBED" in l or 'getElementById("tplSel")' in l
            or "function loadSidecar" in l or "TPL_DATA_EMBED" in l):
        print(i + 1, "|", l.strip()[:260])
