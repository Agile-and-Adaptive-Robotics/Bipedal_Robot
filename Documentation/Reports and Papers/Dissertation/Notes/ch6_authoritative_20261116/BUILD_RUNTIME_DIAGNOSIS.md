# Canonical build runtime investigation

Host: easteregg2. Repository writes succeeded before this review. The compiler issue below is separate from repository write access.

## Installed runtime

- `D:\MiKTeX\miktex\bin\x64\pdflatex.exe`: MiKTeX-PDFTEX 4.27.0, pdfTeX 1.40.29, MiKTeX 26.5.
- `D:\MiKTeX\miktex\bin\x64\bibtex.exe` is available.
- `D:\Anaconda\python.exe` provides PyMuPDF for compiled-page inspection if compilation succeeds.

The normal minimal invocation was:

```powershell
& 'D:\MiKTeX\miktex\bin\x64\pdflatex.exe' --miktex-disable-maintenance --miktex-disable-installer -interaction=nonstopmode -halt-on-error -file-line-error '-output-directory=D:\GitHub\Bipedal_Robot\tmp\ch6_authoritative_build_probe_20261116' 'D:\GitHub\Bipedal_Robot\tmp\ch6_authoritative_build_probe_20261116\probe.tex'
```

It exited unsuccessfully. The exact runtime errors included:

```text
Access is denied.
Function: RegOpenKeyExW
Result: 5
Data: Software\MiKTeX.org\MiKTeX\2.9\Core
The memory dump file could not be found.
fileName="pdflatex.fmt"
Permission denied: path="C:\Users\Ben Bolen\AppData\Roaming\MiKTeX\2.9\miktex/config\issues.json"
```

## Bounded scratch configuration checks

The investigation created a portable startup/configuration and copied readable existing filename databases and the existing `pdflatex.fmt` into `tmp/ch6_authoritative_build_probe_20261116`. The installed tree, user configuration, registry, and original format were not changed. `BUILD_probe.py` records the bounded final probe and `BUILD_PROBE_RESULT.json` its result.

Explicit `TEXINPUTS` and `TFMFONTS` let one minimal portable trial read the document and reach page output. Missing `pdftex.map` lookup then invoked `miktex-makepk`, whose normal runtime again encountered registry error 5. Existing readable font-map, Type 1, encoding, and generated TFM paths were identified for an additional probe. Automatic format lookup subsequently attempted to start `miktex ... formats build pdflatex --engine pdftex`, again encountering registry error 5. Bounded probe processes were stopped using their own command sessions or exact owned subprocess handles. No existing MATLAB or TeX processes were targeted.

## Successful minimal runtime recovery

The final minimal probe succeeded with exit code 0 in 2.359 s and produced a valid one-page, 14,703-byte PDF with the installed `cmr10.pfb` embedded. `BUILD_PROBE_RESULT.json` records its command and transcript. The successful recipe uses explicit `TEXINPUTS`, `TFMFONTS`, `TEXFONTMAPS`, `T1FONTS`, and `ENCFONTS`, and the existing format copied into the scratch registered format directory. Updating only the scratch format's timestamp prevents MiKTeX from launching automatic regeneration; the original and scratch format bytes remain identical. The original installation, user directories, and registry remain unchanged.

`build_authoritative.py --portable-probe` recreates this runtime inside each fresh build directory, with absolute source/font paths and source/scratch format hashes in the build result. The source-format LaTeX and L3 versions are those of the installed format, rather than a newly generated or downloaded toolchain.

The successful minimal runtime probe is not a successful dissertation build, and no minimal or incomplete PDF is a deliverable. Static source checks are recorded separately in `BUILD_STATIC_CHECK.json`. The canonical build driver records its exact attempted commands, source hashes, elapsed times, and final result in `BUILD_RESULT_<revision>.json`. Bibliography processing and compiled-page inspection are contingent on a successful first manuscript pass.
