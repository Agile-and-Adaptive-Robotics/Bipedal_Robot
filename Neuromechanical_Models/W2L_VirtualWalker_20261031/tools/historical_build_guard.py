"""Prevent accidental restoration of superseded generated models or layouts."""
import sys

def require_historical_regeneration():
    flag = '--historical-regeneration'
    if flag not in sys.argv:
        raise SystemExit(
            'Historical builder disabled by default. The canonical model is '
            r'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj. '
            'This helper may recreate retired copies or overwrite later canonical work. '
            'Use --historical-regeneration only for an explicitly requested historical reconstruction.'
        )
    sys.argv.remove(flag)
