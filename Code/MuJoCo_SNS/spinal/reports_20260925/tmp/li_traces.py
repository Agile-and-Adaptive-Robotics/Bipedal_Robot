import io, sys, glob, os
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

files = sorted(glob.glob(r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\Trace_2020*.txt'))
for f in files:
    print("=" * 100)
    print(os.path.basename(f))
    with open(f, 'r', errors='replace') as fh:
        head = [next(fh) for _ in range(6)]
    for h in head:
        print("   HDR:", h.rstrip()[:200])
    # try to load numeric body
    try:
        d = np.genfromtxt(f, skip_header=5, delimiter=',')
        if d.ndim == 1:
            d = d[:, None]
        print("   shape:", d.shape, " t range:", d[0, 0], "->", d[-1, 0])
        for j in range(1, d.shape[1]):
            col = d[:, j]
            print(f"   col{j}: min={np.nanmin(col):.4g} max={np.nanmax(col):.4g} mean={np.nanmean(col):.4g}")
    except Exception as e:
        print("   parse failed:", e)
