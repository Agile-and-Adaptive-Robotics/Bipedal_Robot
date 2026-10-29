# Functional PF five-stage curriculum — Session 1

**Read [the final audit](FINAL_REPORT.md): the original screen omits leg crossing. No valid 3D gait is accepted.**

Machine: easteregg2. Model: articulated human Gait2392. Native SNS outputs; preallocated Numba/NumPy backend.

Each stage seeds the next. Attempted stages are distinguished from accepted stages; failed gates are retained.

| Stage | Completed trials | Best objective | Original screen passed |
|---|---:|---:|---|
| 1 | 100 | 90.128 | False |
| 2 | 100 | 90.827 | False |
| 3 | 100 | 10.643 | False |
| 4 | 100 | 109.905 | True |
| 5 | 100 | 97.510 | True |

Acceptance uses actual bilateral knee/hip excursions and cycle timing, not RG burst counts. Standing has separate tilt/sway/load gates. Ground adds contact-cycle and duty gates. Forward progression and reduced support are not assumed.

The stage-one trials remove muscle afferents, added joint feedback, and interleg wiring. Later stages enable them explicitly. The finite left RG initialization pulse ends at 2.5 s.

Full parameters, metrics, raw trajectories, and logs are stored beside this report. No production winner or database is replaced.
