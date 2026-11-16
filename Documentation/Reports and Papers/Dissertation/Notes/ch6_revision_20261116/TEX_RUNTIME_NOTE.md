# Scratch-local MiKTeX verification runtime

The normal MiKTeX invocation could not read its registry settings or create its default user configuration directory in this restricted execution environment. A scratch-local configuration was used instead; no installed package, original format, user configuration, registry value, or system setting was changed.

## Verified result

`tmp/ch6_20261116/tex_probe/work/minimal.tex` compiled successfully with exit status 0 to a one-page, 24,354-byte PDF. The transcript records embedded `cmmi10.pfb` and `cmr10.pfb` fonts. This verifies the engine, installed-file lookup, generated format, scratch font map, and PDF output. It does not itself verify the complete dissertation.

## Working environment

Paths below are relative to the repository unless absolute:

- `MIKTEX_USERCONFIG`: `tmp/ch6_20261116/tex_probe/config`
- `MIKTEX_USERDATA`: `tmp/ch6_20261116/tex_probe/data`
- `MIKTEX_USERINSTALL`: `C:/Users/Ben/AppData/Local/Programs/MiKTeX`
- `MIKTEX_USERSTARTUPFILE`: `tmp/ch6_20261116/tex_probe/startup.ini`
- `MIKTEX_COMMONSTARTUPFILE`: the same scratch startup file

All environment values must be expanded to absolute paths before launching an executable. Both `pdflatex` and `bibtex` should receive `--miktex-disable-maintenance --miktex-disable-installer`. Maintenance and automatic package installation are not required for an already installed tree and otherwise trigger operations outside the permitted scratch area.

The scratch startup file selects `Portable` configuration and points its user trees at the above locations. The scratch `miktex.ini` sets `RenewFormatsOnUpdate=false` under `[TeXandFriends]`. Existing filename databases and format assets were copied from the user's readable MiKTeX data directory into scratch. The previously saved format still triggered automatic regeneration, so it was replaced only in scratch by direct initialization:

```text
pdftex.exe --miktex-disable-maintenance --miktex-disable-installer -ini -etex -jobname=pdflatex -progname=pdflatex -interaction=nonstopmode -halt-on-error pdflatex.ini
```

This completed with exit status 0 in `tmp/ch6_20261116/tex_probe/work/`, producing `pdflatex.fmt` for the current engine. Copy that new format into the compilation working directory as `pdflatex.fmt` to use the tested current-directory search behavior. The `--no-mktex=fmt` option was not usable in this MiKTeX build and should not be added.

The new format uses LaTeX2e `2026-06-01` and L3 `2026-08-10`. Because the scratch tree initially lacked `language.dat`, the initializer used the standard English `hyphen.tex`; this is sufficient for the English dissertation but should remain declared when comparing page counts with another installation. Original multilingual configuration files remain unchanged in the user's normal tree. The runtime still prints MiKTeX's update-check notice, which is an installation-maintenance notice rather than a TeX compilation error.

Probe processes that stalled in automatic maintenance or format generation were terminated through their own command sessions. The direct initialization and final minimal compilation both completed normally. Existing user MATLAB and TeX processes were not targeted.
