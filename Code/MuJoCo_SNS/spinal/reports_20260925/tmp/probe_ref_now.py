import sys, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SPINAL = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, SPINAL)
os.chdir(SPINAL)
import kine_ref as KR
ref = KR.load_reference()
print(f"T_l = {ref['T_l']!r}  T_r = {ref['T_r']!r}")
print(f"duty_l = {ref['duty_l']!r}  knee_min_l = {ref['knee_min_l']!r}")
