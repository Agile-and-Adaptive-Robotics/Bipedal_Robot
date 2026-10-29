"""Read OpenSim motion tables without importing an audit or changing state."""
import numpy as np
def read_motion(path):
    lines=path.read_text(encoding='utf-8',errors='replace').splitlines()
    end=next(i for i,l in enumerate(lines) if l.strip().lower()=='endheader')
    names=lines[end+1].split()
    values=np.array([[float(x) for x in l.split()[:len(names)]] for l in lines[end+2:] if l.strip()])
    unit=next((l.split('=',1)[1].strip().lower() for l in lines[:end] if l.lower().startswith('indegrees=')),None)
    return values[:,0],names[1:],values[:,1:],unit
