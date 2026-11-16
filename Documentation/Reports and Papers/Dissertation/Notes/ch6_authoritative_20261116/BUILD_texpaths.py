"""Prioritize only installed package directories already read by a build.

The normal installed-tree recursive search remains the final fallback. Package
contents and selected package versions are not changed or reconstructed.
"""
from pathlib import Path
import re

def texinputs(final, reference_log, recursive_fallback=True):
    text=Path(reference_log).read_text(encoding='utf-8',errors='replace')
    directories=set()
    for value in re.findall(r'\((D:[\\/]MiKTeX[\\/]tex[\\/][^()\r\n]*)',text):
        parent=Path(value.strip()).parent
        if parent.is_dir(): directories.add(parent)
    if not directories: raise ValueError('No installed package directories found in reference log')
    base=Path(r'D:\MiKTeX\tex\latex\base')
    priority=[str(final),'.',str(base)]
    priority.extend(str(p) for p in sorted(directories) if p!=base)
    if recursive_fallback: priority.extend([r'D:\MiKTeX\tex//',''])
    result=';'.join(priority)
    if len(result)>30000: raise RuntimeError('Prioritized TEXINPUTS exceeds safe environment length')
    return result, {'reference_log':str(reference_log),'directory_count':len(directories),
                    'environment_length':len(result),'fallback':r'D:\MiKTeX\tex//' if recursive_fallback else None}
