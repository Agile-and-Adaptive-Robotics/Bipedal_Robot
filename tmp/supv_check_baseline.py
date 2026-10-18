"""Supervisor spot-check part 3: extensor baselineScores vs Table C9 rigid RMSE column."""
import scipy.io as sio
import numpy as np

ext = sio.loadmat(r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo\minimizeExt10mmX3_results_20260920_flxr77.mat", squeeze_me=True, struct_as_record=False)
bs = ext["baselineScores"]
print("baselineScores (9x3):")
print(np.round(bs, 4))
# Table C9 rigid RMSE (pool tests 1,2,5,6,7,8 = rows 0,1,4,5,6,7 of the 9): 3.09 3.62 2.44 2.88 2.18 0.95
pool = [0, 1, 4, 5, 6, 7]
print("pool rows col0:", np.round(bs[pool, 0], 4))
