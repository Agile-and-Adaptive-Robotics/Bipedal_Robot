"""Regenerate a small reviewable results table from separate raw campaigns."""
import json
from pathlib import Path
H=Path(__file__).resolve().parent
lines=['# Corrected physical curriculum: measured outcomes','',
'Machine: easteregg2. Native neural outputs, MuJoCo 2.3.7, implicitfast, 2 ms, sparse preallocated Numba SNS backend.',
'','These are stage tests and diagnostic continuations. Later tests do not certify that earlier failed stages graduated.',
'','| Campaign / case | Stage | Pass | Knee cycles R/L | Foot load R/L (BW) | Spring support (BW) | Lateral restraint RMS (BW) |',
'|---|---:|---|---|---|---:|---:|']
for folder in sorted((H/'runs').glob('corrected_curriculum*')):
    for p in sorted(folder.glob('*.result.json')):
        r=json.loads(p.read_text());a=r['assessment'];m=r['measured']
        foot='/'.join(f'{v:.3f}' for v in m['mean_foot_normal_bw'])
        lines.append(f"| {folder.name} / {p.name[:-12]} | {r['stage']} | {a['accepted']} | {a.get('knee_cycles','standing')} | {foot} | {m['mean_spring_support_bw']:.3f} | {m.get('lateral_harness_rms_bw',float('nan')):.3f} |")
lines+=['','Vertical loads in current v2 tests count each foot-ground contact once. The legacy runner assigns each contact to both heel and toe, doubling its force sum. Its original contact metrics remain in raw files; current ground-load verification uses independent world-vertical forces. v1 s1/s2 were air tests; subsequent cases used the corrected force measurement.',
'','The harness wrench is computed from the actual spring and damper forces on pelvis coordinates and projected into native world axes. It is a virtual instrumentation calculation, not a physical single-point six-axis load-cell attachment. Distributed hip/lumbar restrictions remain additional assistance and are explicitly configured.',
'','Native world x is anterior, y is left, and z is superior; anatomical right force is negative native y. Record force and torque with contact timing; a large harness lateral reaction alone cannot identify which foot generated it.',
'','Effort: mean activation squared and excitation squared are recorded. These are dimensionless regularizers, not calibrated metabolic energy, and never override movement, contact, or crossing gates.',
'','Speed-dependent feedback evidence: Hurteau et al. 2017 reports nonlinear cutaneous-reflex modulation in spinal cats, strongest at intermediate speeds, not a universal monotonic decline (https://doi.org/10.1523/JNEUROSCI.3042-16.2017). Mayer et al. 2018 supports a role for spindle feedback in speed-dependent extensor recruitment in mice (https://doi.org/10.1152/jn.00250.2018). These do not establish that Ib becomes unnecessary or that Ia and II should both increase identically.',
'','Separate pathway ablations across achieved speed, assistance, and phase are required before attributing causal importance. A change in signal amplitude or fitted conductance alone is insufficient.']
(H/'PHYSICAL_RESULTS.md').write_text('\n'.join(lines)+'\n')
