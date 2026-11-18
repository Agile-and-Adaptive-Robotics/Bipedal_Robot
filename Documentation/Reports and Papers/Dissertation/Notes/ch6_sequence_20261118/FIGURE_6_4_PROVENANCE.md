# Figure 6.4 — transient knee feedback

The source renderer plotted the injected current from −0.15 to +0.25 s and the knee and motor-neuron responses from −1 to +10 s. The pulse onset consequently occupied different horizontal positions. The revised figure uses one horizontal transform for panels B–D in each column. The left column covers −1 to +10 s and the shaded right column enlarges −20 to +200 ms. Both place zero at stimulus onset. The full-resolution vector shows the original six connections in panel A.

## Inputs and reproduced response

`Code/Matlab/SNS_Simscape/results/knee_transient_reflex_20261028/run_transient.m` and its two `replay_summary` files provide the experiment. A separate MATLAB R2025a run of Appendix D's commands generated `tmp/ch6_reader_20261117/knee/results/probe/reader_summary.mat`. I compared the time vector, all three 50,001 × 9 matrices, and the channel names against the archived arrays. Every stored value matched. The original figure renderer is `Notes/ch6_authoritative_20261116/CAD_render_sns_compact.py`, lines 148–169.

| Source | SHA-256 |
| --- | --- |
| `run_transient.m` | `8ad5b8164153236c2938e723c88cc09a0255a72da3ea56928fd5107d1693fbb8` |
| `replay_summary.mat` | `01a67c749c0a317986cbcc7d55192f0795770c722b15bd4ee52f89cecfd9309c` |
| `replay_summary.json` | `d0ccf9d0f3603949777442a63012453e75ab48d9732d98914463f316c0401e32` |
| Original six-connection schematic vector stream | `5b6a6ef16d83dcf0cb207e256553a9add76da214828e53c0e24be1a71731101c` |

The 2 nA extensor Ia stimulus occupies 30 consecutive 1 ms samples beginning at simulation time 30 s. The next sample is zero. Panel C uses the dissertation's negative-flexion convention, plotting `−rad2deg(th)` in both trials. Relative to the concurrent no-stimulus control, peak extension was 1.0907928305° at 111 ms. From 6.355 s after pulse onset to the end of the trial, the stimulated knee stayed within 0.01° of its baseline and its angular speed stayed at or below 0.0001 rad/s. Panel D shows the stimulated trial's extensor and flexor motor-neuron membrane voltages in mV.

## Vector and print checks

The five left-column paths retain every original sample over the full 50 s trial. An axis clip displays 11,001 samples per signal. Each right-column detail contains all 221 samples in its interval. The pulse has two additional vertices for its vertical rising and falling edges. The angular and voltage traces join adjacent samples. The generator compares the retained panel A stream with the original drawing byte for byte.

`verify_knee_time_pdf.py` reads the clipping rectangles and every plotted curve back from the delivered PDF. The recovered paths place zero at 67.5454545 PDF points in the left column and 274 points in the right, with respective scales of 16.5454545 and 550 points/s in all three rows. The largest discrepancy from any independently calculated coordinate was 0.000498 point. At the chapter's 14 cm width, the figure is 16.65 cm tall. The TeX proof preserves a 10.0006 pt minimum for the circuit and 10.4006 pt for the new text.

`FIGURE_6_4_NUMERICAL_VERIFICATION.json` records the source-array hashes, axis transforms, and sample counts. `FIGURE_6_4_PRINT_VERIFICATION.json` records the standalone TeX placement and recovered PDF coordinates. `FIGURE_6_4_REVIEW.png` is the 175 dpi review image. At printed size, the onset lines align, the pulse is visibly rectangular in the enlarged view, and the peak and recovery annotations identify their associated points.

In the 315-page delivery proof, Figure 6.4 and its complete caption appear on printed page 143 (PDF page 164). I inspected the 144 dpi page in `FIGURE_6_4_DELIVERY_REVIEW.png`. The four-line caption describes the common onset, the seconds and milliseconds windows, and the knee sign. The revised page fits at 14 cm without a float-size warning. `FIGURE_6_4_COMBINED_VERIFICATION.json` records an exact match between the delivered figure and the embedded 5,088,241-byte page-content stream and circuit stream. The delivery PDF SHA-256 is `1fc7380fe1fe61bf0d9be4bc7edf12131912387cdf3ba7bd98e185a909cf8a73`.

## Files and reproduction

The two active files are directly in `Figures/60-sns/`. In the Overleaf project the destination is `figs/60-sns/` with the same names. The chapter includes the PDF at 14 cm. The PDF and its alt text regenerate byte for byte with the commands below. `documentPython` denotes the bundled document runtime with ReportLab and pypdf. The generator uses `D:/Anaconda/envs/myo/python.exe` to load the MATLAB data on easteregg2. On another machine, pass its project environment with `--scientific-python`.

```powershell
$notes = 'Documentation/Reports and Papers/Dissertation/Notes/ch6_sequence_20261118'
$out = 'tmp/ch6_sequence_20261118/fig64_generator_repeat'
$proof = 'Documentation/Reports and Papers/Dissertation/Notes/ch6_sequence_20261118/Bolen_Dissertation_revised.pdf'
& $documentPython "$notes/make_knee_time_figure.py" $out
& $documentPython "$notes/verify_knee_time_pdf.py" "$out/sns_knee_reflex_authoritative.pdf" $proof "$out/print-check.json" --proof-page 164
```

To compare the independent replay, add `--replay 'tmp/ch6_reader_20261117/knee/results/probe/reader_summary.mat'` to the generator call. The verifier accepts any copy of the same delivery PDF as `$proof`. Relative paths above start at the repository root.

| Delivered file | Bytes | SHA-256 |
| --- | ---: | --- |
| `sns_knee_reflex_authoritative.pdf` | 1,506,355 | `80e845185f25adfbf3a69b39f3939982bff086e3fe3b4fb1d53d165fef6a256b` |
| `sns_knee_reflex_authoritative.alt.txt` | 539 | `863579994ed75161c85e29420d870fae7f9491bdaf91d2f2c79f13ca2a832212` |
| `make_knee_time_figure.py` | 15,101 | `a070eccb48e728855e60940d3c73d1c19b01973c0f4f9ecf50dbe6dacc8a8f92` |
| `verify_knee_time_pdf.py` | 6,074 | `145b5c097dd10a74088b80d0863f4264bd4b156a167c441e9d6b245df8a5398f` |

The exact superseded PDF and alt text are in `Figures/60-sns/Deprecated/20261118_knee_time_alignment/` with `manifest.json`. Their SHA-256 values are `771c02ec1777d181ee19d36d69b4dac635ffc2891e2b55be64adf836480a3747` and `25ade8f298cbd15b806cb6ebfa2603ab19ef019e530edecd37b2e92fb80764bf`, respectively. Git ignores the chapter's `Deprecated` directory. For transfer, copy its two files and manifest explicitly, or take the identical two files from the tracked `Notes/ch6_authoritative_20261116/CAD_sns_compact/` archive and restore the manifest. The generator accepts either archived source.
