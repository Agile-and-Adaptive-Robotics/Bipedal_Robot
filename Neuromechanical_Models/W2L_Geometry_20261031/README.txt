W2L corrected geometry and asymmetric initial pose

Original shared-RG controller

Open W2L_Geometry.aproj. This is a new copy; source and recording folders were not modified. Root remains frozen, with the original moving treadmill. GUI playback is MatchPhysicsStep (normal speed); the ASIM uses PlaybackControlMode=1.

Initial pose, relative to the upright pelvis:
- Left hip: 10 degrees extension; knee: 10 degrees flexion, nearly straight.
- Right hip: 20 degrees flexion; knee: 55 degrees flexion.
- Ankles: 0 degrees anatomically. Feet point down with the shanks (left foot pitch -20 degrees, right -35 degrees). Minimum SOLID foot/toe-box clearance above the actual belt top: 10 mm. Contact sensor boxes are excluded from this clearance definition.

The original tibia center was placed at (0,-0.43,0) in its femur frame without rotating the knee-to-center offset. Its knee anchor therefore landed 11.514 mm anterior to the distal femur center. This was not a nonzero COM field: both limb boxes have COM=(0,0,0). The relevant quantity is each rotated body's center-to-joint vector.

The revised center is calculated as p_child = p_parent_anchor - R_child * q_child_anchor. The tibia's knee anchor is at the proximal box-face center, local (0,0,-0.21) m; the distal femur point is (0,-0.21,0) m. This also removes the old 10 mm hinge-to-box-face discrepancy. Hip, tibia, and foot centers and both GUI LocalMatrix/degree and ASIM radian transforms were updated consistently. Muscle/receptor elements, attachment local coordinates, body IDs, neural wiring, weights, and stimuli are retained byte for byte at element level. The new pose and joint geometry change muscle lengths and mechanics, so these are new experiments, not replays of dissertation results.

Native AnimatSimulator verification: exit 0, ten seconds, 1,001 samples. Both knees start with coincident femur/tibia markers in native world coordinates. Maximum dynamic separation between those markers: 0.549 mm, consistent with finite constraint compliance. Root Y stays 0.982885 m. Native probes at all 32 solid-foot/toe box corners confirm a lowest sole Y of 0.015 m against belt top Y 0.005 m. See GEOMETRY_REPORT.json for full geometry and global hinge coordinates; NATIVE_TEST.json for measurements; PRESERVATION.json for source hashes and neural/muscle preservation checks. All DiagramXml pages parse. Interactive GUI opening was not attempted, to avoid disturbing the current AnimatLab session.

Joint-angle traces are engine coordinates relative to the newly assembled pose; the anatomical initial angles above are geometry angles. Relative hip, knee, and ankle limits were re-referenced to the new assembled zero so the original anatomical ranges are preserved. This corrects the earlier prototype that accidentally trapped the right knee at 55 degrees or more. The current right knee extends from 55 to about 2.2 degrees in the native fixed-root tests (small soft-limit overshoot below the nominal 3-degree limit). Muscle properties and controller gains were not retuned. See ANATOMICAL_LIMITS.md for the full reference-zero and limit table. No claim of stable unsupported walking is made. The separate W2L_VirtualWalker_20261031 copy adds the requested hip-wing/rail support.

Reproduction: original copy uses tools/build_geometry.py. Newer copy can be rebuilt with that script using source Neuromechanical_Models/W2L_Comparison_20261030/recording_projects/bilateral_treadmill and output Neuromechanical_Models/W2L_Geometry_Newer_20261031. Native verification uses the original copy's tools/test_geometry.py.
