# CAD pressure-interface checkpoint

Machine: **easteregg2**. This is opt-in plumbing and a bounded sleeve-placement study, not a mounted pressure-actuated standing result.

The current CAD model is `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\ah001_neural.xml`. It uses x anterior, y superior, z right, and gravity `[0,-9.81,0]`. Do not substitute a legacy z-up frame without an explicit adapter.

## Force and sensory interface

`cad_moving_force.py` provides `CADMovingForce`. The caller must explicitly disable **all** native reference-force gains and biases. The bridge then uses the complete `DynamicRoutes.length_jacobian` and adds `-J_length.T @ positive_tensions` to the caller's force buffer exactly once per physics step. It preserves existing support forces. Do not additionally apply a second moving-site correction to these same external tensions.

`advance_pressure()` returns the actual modeled pressure, positive tension, complete route length, complete route velocity, and generalized force. Neural **Ib force feedback must use these bank tensions**, because native `data.actuator_force` is zero in this all-external-force mode. Ia/II feedback must use the returned lengths/velocities rather than silently reverting to fixed-site derivatives. The caller must record the neural sensor conversion and reference scaling. Human MIF is a reference scale, not a calibrated BPA force capacity.

An alternative implementation could preserve native negative tendon forces and add only the moving-site contribution. That is a separate implementation with its own force/sensor parity gate; `NativeBPAForce` by itself does not implement the CAD moving-site correction. The current all-external bridge avoids that double-counting ambiguity.

`check_cad_moving_force.py` verifies the actual 54-route CAD geometry at three saved poses. Full derivatives match finite differences within `5.78e-10` m per generalized coordinate, and power within `3.4e-16` W. Native fixed-site derivatives alone miss up to `0.04951` m/rad. These prescribed 10–100 N tests establish force-direction/virtual-work plumbing, not pressure dynamics or standing.

## Placement subset and outstanding parameters

The simultaneous candidate subset is **vas_med_r and vas_med_l**, with paired identical sleeve dimensions derived from the larger bilateral route excursion. Sleeve envelopes clear CAD obstacles and each other at 17 saved poses; their minimum mutual surface gap is 63.43 mm. Brackets, fittings, cable routing, Bowden bend radius/friction, attachment loads, added masses, and continuous swept motion remain unverified.

The separately considered bifemsh sleeve pair clears CAD obstacles individually, but its conservative envelopes overlap by 18.07 mm during the sampled motion. Do not count the four body-clear sleeves as a simultaneous four-actuator installation.

The existing route ledger assigns active `vas_med_r/l` to master ID 11 labeled **Vastus Intermedius**. That naming discrepancy is flagged, not silently corrected. The proposal follows the actual active vas_med paths.

No complete hardware assignment is inferred from these placements. Effective series stiffness, exact nonmuscle/Bowden length, pressure response, delay, friction/hysteresis, and component masses remain unassigned. The all-null 54-actuator template deliberately remains non-runnable. The virtual Gait2392 pressure design is a separate assumed design and must not be relabeled as mounted CAD identification.

## Evidence and repeat commands

- Record: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\EXPERIMENT_RECORD.md`
- Full force gate: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\CAD_MOVING_FORCE_PARITY.json`
- Candidate coordinates: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\mount_subset\PAIRED_REMOTE_SLEEVE_MOUNT_PROPOSAL.json`
- Adopted simultaneous subset: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\mount_subset\SLEEVE_PAIR_CLEARANCE.json`
- Inspected geometry figure: `D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\mount_subset\CAD_SLEEVE_MOUNT_PROPOSAL.png` (also PDF/SVG).

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\check_cad_moving_force.py'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\cad_remote_mount_proposal.py'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\draw_cad_mount_proposal.py'
```

These commands regenerate geometry/parity reports and figures. They do not start a whole-CAD neural simulation or assign missing hardware properties.
