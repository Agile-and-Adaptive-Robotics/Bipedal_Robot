# Anatomical reference and hinge limits

This replaces the first geometry-only prototype. AnimatLab gives the assembled pose joint coordinate zero. Each native limit is therefore shifted by old anatomical zero minus new anatomical zero. Native world-marker measurements verified a positive relationship between engine rotation and the anatomical hip, knee, and ankle angles. This preserves the original physical ranges without expanding them.

| Joint | Old zero (deg) | New zero (deg) | Relative shift (deg) | New native limits (deg) | Preserved anatomical limits (deg) |
|---|---:|---:|---:|---|---|
| hip_L | -3.0000 | -10.0000 | 7.0000 | -8.0000 to 30.0000 | -18.0000 to 20.0000 |
| knee_L | 3.0000 | 10.0000 | -7.0000 | -7.0000 to 53.0001 | 3.0000 to 63.0001 |
| ankle_L | 8.0002 | 0.0000 | 8.0002 | -11.9998 to 3.0002 | -11.9998 to 3.0002 |
| hip_R | -3.0000 | 20.0000 | -23.0000 | -38.0000 to 0.0000 | -18.0000 to 20.0000 |
| knee_R | 3.0000 | 55.0000 | -52.0000 | -52.0000 to 8.0001 | 3.0000 to 63.0001 |
| ankle_R | 8.0002 | 0.0000 | 8.0002 | -11.9998 to 3.0002 | -11.9998 to 3.0002 |

Both original ankle zeros are approximately +8 degrees; the original physical range is approximately -12 to +3 degrees on each side. The nominal source hip range is -18 to +20 degrees, and knee range +3 to +63 degrees. The source GUI uses degrees and source ASIM uses rounded radians, so their values differ by less than 0.0003 degrees. That existing precision difference is preserved.

Toe relative poses are unchanged, so their relative limits are unchanged. Joint motors, neural parameters, and muscle properties are unchanged. The fixed-root simulation remains a supported motion experiment; success does not establish unsupported walking.

Clearance is calculated from rotated corners of the four solid foot/toe boxes and the actual moving-belt top. Contact sensor boxes are excluded. The final 32-corner native probe test gives minimum sole Y 0.015 m, belt top Y 0.005 m, and 0.010 m initial gap. The earlier 14.291 mm number described a solid toe gap in the first prototype, while its 10 mm number included a lower contact-sensor box. Those first-prototype clearances are superseded.
