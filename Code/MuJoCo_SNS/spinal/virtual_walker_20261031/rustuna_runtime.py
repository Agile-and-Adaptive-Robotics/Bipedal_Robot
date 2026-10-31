"""Load the isolated wheel without changing the scientific environment."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
def load_rustuna():
 if sys.version_info<(3,11):raise RuntimeError('Run Rustuna with D:\\Anaconda\\python.exe (3.11); leave myo3.10 physics in its worker.')
 root=HERE.parents[3]/'tmp/rustuna_runtime_py311'
 if not root.exists():raise RuntimeError('Isolated Rustuna runtime missing: '+str(root))
 sys.path.insert(0,str(root))
 import rustuna
 return rustuna
