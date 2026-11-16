# Independent historical mechanical-evidence review

Bounded follow-up review on easteregg2, 16 November 2026. This review covers the unchanged twelve-muscle force comparison and mechanical-repair claims in Chapter 6 and Appendix D. No canonical manuscript, scientific source, model, or result file was edited. No simulation, forward dynamics, static force replay, parameter fitting, or optimizer ran.

`HISTORICAL_MECHANICAL_audit.py` reads the retained data, parses the primary sources, compiles the current unpatched model, and evaluates the isolated pure XML transformer without importing or running the experiment driver. It does not call `mj_forward` or `mj_step`, and it does not write a patched model into the scientific source tree. Exact input hashes and detailed checks are retained in `HISTORICAL_MECHANICAL_RESULTS.json`.

## Force-comparison verdict

The existing twelve rows match the retained CSV at the reported precision. Exactly eight satisfy both stated conditions: `0.87 <= R2 <= 0.98` and mean MuJoCo/OpenSim force ratio within `[0.9,1.1]`. These are medial gastrocnemius, lateral gastrocnemius, tibialis anterior, vastus lateralis, psoas, gluteus maximus 2, tensor fasciae latae, and peroneus brevis. Soleus and semimembranosus have lower R2. Rectus femoris and biceps femoris short head have negative R2, as reported.

All twelve entries are **right-side muscles**. Adding that scope to the table caption or preceding paragraph would remove a minor ambiguity without changing any result.

| Muscle, right side | Saved R2 | Saved squared Pearson correlation | Saved lag, frames | Saved mean MuJoCo/OpenSim force ratio |
|---|---:|---:|---:|---:|
| Soleus | 0.7423303481 | 0.7997166479 | 1 | 0.9156346470 |
| Medial gastrocnemius | 0.8881319833 | 0.8940968082 | 2 | 0.9095977480 |
| Lateral gastrocnemius | 0.8836081344 | 0.8871205578 | 2 | 0.9272328155 |
| Tibialis anterior | 0.9449255929 | 0.9451627710 | 0 | 1.0063626008 |
| Vastus lateralis | 0.9406639034 | 0.9454483495 | 1 | 0.9265095111 |
| Rectus femoris | −1.1502919379 | 0.3835226163 | 0 | 1.0878738873 |
| Semimembranosus | 0.5650806119 | 0.6206771699 | 0 | 0.7756142078 |
| Biceps femoris short head | −0.2606914436 | 0.8675516761 | 0 | 1.8462983710 |
| Psoas | 0.9268307215 | 0.9339191261 | 1 | 0.9461942468 |
| Gluteus maximus 2 | 0.9259562830 | 0.9430871023 | 1 | 0.9212762572 |
| Tensor fasciae latae | 0.9330872419 | 0.9348311058 | 0 | 1.0180463575 |
| Peroneus brevis | 0.9769681745 | 0.9894936903 | 0 | 0.9169773836 |

The raw local reference inputs support the retained source binding: `bsolve_out.npz` contains 121 poses from 0.5 through 2.5 s, with a median frame interval of 0.017 s. Linear interpolation of the current `zz_bsolve_StaticOptimization_force.sto` onto those saved pose times reproduces the CSV's OpenSim peak force for every one of the twelve muscles with zero difference. The approximately 17 ms frame-lag description is therefore appropriate. The lag is an integer analysis-frame diagnostic, not a physical tendon delay.

The actual generator body, rather than its older opening description, establishes the matched-activation comparison. It reads the OpenSim static-optimization force and activation files. For each saved pose, it projects coupled follower coordinates, calls the MuJoCo forward calculation, assigns interpolated OpenSim activations clipped to `[0,1]`, calls forward again, and takes absolute actuator force. A separate earlier pass evaluates unit activation; those unit-activation results are not the table's R2 values. R2 and squared correlation are calculated on the unshifted matched-activation sequences. The best cross-correlation lag is reported separately and does not realign the sequences before their R2 calculation.

The comparison uses fresh `MjData` without assigning the trajectory velocities and does not integrate activation or fiber dynamics. The manuscript's static-pose/activation scope, warning against dynamic-equivalence claims, and rejection of a tendon-delay interpretation are correct. The comparison script's old docstring refers to MuscleAnalysis tendon-force and unit-activation calculations; the executed body instead deliberately selects the static-optimization force output after documenting zero-control problems in the former tendon-force file. The manuscript correctly follows the executed body.

### CSV schema limitation

The CSV header contains six fields, while every row contains eight values. The primary writer makes their order unambiguous:

1. muscle name;
2. matched-activation R2;
3. matched-activation squared Pearson correlation;
4. lag in frames;
5. mean MuJoCo/OpenSim force ratio;
6. OpenSim peak force;
7. matched-activation MuJoCo peak force;
8. unit-activation R2.

The two peak-force column labels were omitted from the header. This does not change the first five statistics used in the dissertation, but a naive named-column reader would misidentify the later fields. This audit maps them against the primary writer and preserves the original file unchanged.

## Mechanical-repair verification

The current raw converted model compiles with 92 actuator entries. The current `runner.patch_xml` function was isolated with the Python AST and applied only to an in-memory copy. Its assertions and resulting XML establish the following implementation facts; no actuator-dynamics test was rerun.

| Claim | OLD in raw model | NEW in in-memory repair | Finding |
|---|---|---|---|
| Right and left hip-flexion axes | `[0,0,1]` | `[0,0,-1]` | Both negate as stated |
| Right hip-adduction axis | `[1,0,0]` | `[-1,0,0]` | Negated |
| Left hip-adduction axis | `[-1,0,0]` | `[1,0,0]` | Negated |
| Six retained small-muscle actuators | Nonzero source scales | Active and passive scales both zero | `quad_fem`, `gem`, and `peri`, right and left, remain present |
| Bilateral erector-spinae active scale | 1 N | 2500 N | Both agree with source OpenSim maximum isometric force |
| Bilateral internal-oblique active scale | 1 N | 900 N | Both agree with source |
| Bilateral external-oblique active scale | 1 N | 900 N | Both agree with source |
| Bilateral extensor-hallucis-longus active scale | 1 N | 162 N | Both agree with source |
| Passive scales for those eight repaired entries | Source values | Unchanged | Confirmed separately from active correction |
| Anatomical driver-hinge limits | `limited=false` | `limited=true` for 14 hinges | Ranges unchanged |
| Coupled path-point followers | Sampled range attributes | `limited=false`, with retained coupling | Verified separately from anatomical driver limits |
| Explicit ground–body pairs | 19 | Six foot–ground pairs in ground mode, zero in air mode | Confirms contact policy |
| Virtual path-marker collision flags | 32 explicit `2/2` pairs of flags | None retain `2/2`; changed to `0/0` | Confirms virtual geometry is muted |
| Right/left rectus-femoris intermediate tendon site | Own P2 | Corresponding vastus-medialis P4 | Origin P1 and distal P3 retained |

The knee driver axis and range are preserved by this repair; the function explicitly retains the native flexion-negative coordinate. The follower coordinate update in `bsolve_ik.apply_eq_followers` evaluates the joint-coupling polynomial from its driver, supporting the manuscript's warning against posing the knee while ignoring those translations.

The source model currently contains 82 joint-equality follower entries; comments referring to an older count of 36 are not used as quantitative evidence. The dissertation does not assert that older count. The fourteen anatomical driver hinges remain a separate verified count.

The pure transformer also adds the documented harness and passive-restraint settings. This review does not interpret those as intrinsic physiological properties or remove their mechanical influence. The manuscript correctly treats the repaired reference, the BPA-routed counterpart, and the supported controller experiments as distinct configurations.

## Primary files and identity limits

All paths below are relative to `D:\GitHub\Bipedal_Robot`; complete paths and SHA256 values are retained in the JSON audit.

- `Code\MuJoCo_SNS\spinal\_muscle_force_compare.py`
- `Code\MuJoCo_SNS\spinal\muscle_force_compare.csv`
- `Code\MuJoCo_SNS\spinal\runner.py`
- `Code\MuJoCo_SNS\spinal\bsolve_ik.py`
- `Code\MuJoCo_SNS\spinal\bsolve_out.npz`
- `Solid_Models\OpenSim\Gait2392_Robotbody\ResultsBSolve\zz_bsolve_StaticOptimization_force.sto`
- `Solid_Models\OpenSim\Gait2392_Robotbody\ResultsBSolve\zz_bsolve_StaticOptimization_activation.sto`
- `Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_simbody.osim`
- `Solid_Models\OpenSim\Gait2392_Robotbody\subject01_simbody.osim`
- `Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml`

The generator retains the summary CSV and plot but does not save its `F_mj_so` or `F_os` arrays. The historical matched MuJoCo force samples and a compiled-model hash bound to this old comparison were not identified. The pose NPZ likewise does not save that comparison's process environment or chosen `AARL_MODEL` override. The current default points to `gait2392_simbody_cvt3.xml`; the OpenSim static-optimization setup uses `subject01_simbody.osim`. Current source paths, raw reference-force peak agreement, and the retained summary support the reported values, but they do not establish a complete historical compiled-model identity or permit independent R2 recomputation from absent MuJoCo force samples.

No new static force replay was run to fill that gap. The manuscript already calls these retained saved statistics and limits their interpretation to pose-wise comparison. The source implementation of the historical sign repairs is verified; their original dynamic test outcomes were not newly rerun.

## Disposition

No numerical correction is required to the twelve-muscle table or eight-of-twelve statement. The professional interpretation of static agreement, force-scale differences, negative R2, correlation, and lag is appropriate. Add the right-side scope if convenient. Preserve the malformed CSV and old source descriptions as historical evidence rather than silently rewriting them in this dissertation-only pass. The remaining raw-force/model-identity limitations belong in the revision review record.
