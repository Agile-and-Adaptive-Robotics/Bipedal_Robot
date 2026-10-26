# Desktop (EB475WS4) replication of the OPTIMIZATION_SPEED_MOD:
# 1) numba smoke (cache=True from a file context), 2) sparse-backend
# validation (must be bit-identical), 3) FastTrial availability check.
import sys

import numpy as np

sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg")

from numba import njit


@njit(cache=True)
def _smoke(x):
    s = 0.0
    for i in range(len(x)):
        s += x[i] * 2
    return s


print("numba smoke:", _smoke(np.array([1.0, 2.0])), "(expect 6.0)")

import fast_eval  # noqa: E402

r = fast_eval.validate()
print("fast_eval.validate ->", r)
