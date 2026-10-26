from numba import njit
import numpy as np

@njit(cache=True)
def test_fn(x):
    s = 0.0
    for i in range(len(x)):
        s += x[i] * 2
    return s

print('numba works:', test_fn(np.array([1.0, 2.0, 3.0])))
print('second call (cached):', test_fn(np.array([4.0, 5.0])))
