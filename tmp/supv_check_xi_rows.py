"""Supervisor spot-check part 2: flexor front row 77 and extensor flxr77 row 32."""
import scipy.io as sio
import numpy as np

flx = sio.loadmat(r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo\minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat", squeeze_me=True, struct_as_record=False)
rs = flx["results_sort"]
print("=== FLEXOR results_sort: header guess — 3 param cols (Xi0 cm? log10 Xi1? log10 Xi2?) + scores ===")
print("first 3 rows (cols 0-7):")
print(np.round(rs[:3, :8], 4))
row77 = rs[76]  # 1-based row 77
print("ROW 77 (0-idx 76), all 13 cols:")
print(np.round(row77, 6))
# candidate interpretations
x0_cm = row77[0]
print("interpretation A: Xi0[cm]=%.6f -> mm %.4f ; 10^c1=%.1f ; 10^c2=%.1f" % (x0_cm, x0_cm*10, 10**row77[1], 10**row77[2]))
print("interpretation B: raw Xi0[m]=%.6g ; Xi1=%.6g ; Xi2=%.6g" % (row77[0], row77[1], row77[2]))
print("xCols:", flx["xCols"], " W:", flx["W"])

ext = sio.loadmat(r"D:\Github\Bipedal_Robot\Testing_Data\2022_02_Festo\minimizeExt10mmX3_results_20260920_flxr77.mat", squeeze_me=True, struct_as_record=False)
xa = ext["results_sort"]
print("\n=== EXTENSOR results_sort (180x15) row 32 (0-idx 31), all 15 cols ===")
row32 = xa[31]
print(np.round(row32, 6))
print("xCols:", ext["xCols"])
print("lb:", np.round(ext["lb"], 4), " ub:", np.round(ext["ub"], 4))
print("interpretation: Xi0[cm]=%.6f -> mm %.4f ; Xi3=%.6f" % (row32[0], row32[0]*10, row32[1] if ext["xCols"].size >= 2 else float("nan")))
print("cols 0-5 of first 3 rows:")
print(np.round(xa[:3, :6], 6))
