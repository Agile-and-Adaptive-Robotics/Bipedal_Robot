# Read-only dump of the flxr77 extensor mat split fields (train/hold tests)
# and the launcher log header that defines the pool/fold of record.
import scipy.io as sio
import numpy as np
import os
BASE = r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo"
E = sio.loadmat(os.path.join(BASE, "minimizeExt10mmX3_results_20260920_flxr77.mat"),
                squeeze_me=True, struct_as_record=False)
for k in ("allBPA", "numHold", "fold", "holdoutIdx", "trainIdx", "holdTests", "trainTests",
          "pick", "exthold", "numBPA", "tag", "RESULTFILE", "KMODE", "Xi0", "Xi1", "Xi2", "Xi3"):
    if k in E:
        v = E[k]
        print(k, "=", np.asarray(v).flatten() if not np.isscalar(v) else v)
print()
# how many rows in filtered_results belong to each fold value (col 1)?
fr = np.asarray(E["filtered_results"])
print("filtered_results shape:", fr.shape)
u, c = np.unique(fr[:, 0], return_counts=True)
print("col0 value counts:", dict(zip(u.tolist(), c.tolist())))
print("row32 col0..4:", fr[31, :5])
