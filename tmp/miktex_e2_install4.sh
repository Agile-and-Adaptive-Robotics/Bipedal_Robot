#!/usr/bin/env bash
# MiKTeX take 4 — install-only from the completed local repo (5.9 GB, done).
# Root rule learned: --user-config/--user-data must NOT nest inside
# --user-install. Keep install root on D:, leave config/data at user defaults.
set -x
cd /d/temp/miktex_setup_e2/standalone || exit 1

rm -rf /d/MiKTeX

./miktexsetup_standalone.exe --shared=no --package-set=complete \
  --local-package-repository="D:\\temp\\miktex_repo" \
  --user-install="D:\\MiKTeX" --modify-path=yes \
  install
echo "INSTALL-EXIT=$?"

ls -la "/d/MiKTeX/miktex/bin/x64/pdflatex.exe" "/d/MiKTeX/miktex/bin/x64/bibtex.exe" || { echo "PDFLATEX-MISSING"; exit 3; }

"/d/MiKTeX/miktex/bin/x64/initexmf.exe" --set-config-value=[MPM]AutoInstall=1 || echo "AUTOINSTALL-SET-FAILED"
echo "MIKTEX-SETUP-DONE"
