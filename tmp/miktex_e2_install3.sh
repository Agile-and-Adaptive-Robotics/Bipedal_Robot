#!/usr/bin/env bash
# MiKTeX take 3 — headless miktexsetup_standalone (the CI-grade path; the GUI
# wizard hung on its own twice). Complete package set per Ben's "install the
# whole thing on D:". User-scope to D:\MiKTeX, PATH modified by the installer.
set -x
cd /d/temp/miktex_setup_e2/standalone || exit 1

# Pick a repository URL from the official list; fall back to a CTAN mirror path.
REMOTE=$(./miktexsetup_standalone.exe --list-repositories 2>/dev/null | grep -m1 "^https")
[ -z "$REMOTE" ] && REMOTE="https://mirror.ctan.org/systems/win32/miktex/pub/"
echo "REMOTE=$REMOTE"

mkdir -p /d/temp/miktex_repo

./miktexsetup_standalone.exe --package-set=complete \
  --local-package-repository="D:\\temp\\miktex_repo" \
  --remote-package-repository="$REMOTE" \
  download
echo "DOWNLOAD-EXIT=$?"
[ -d /d/temp/miktex_repo ] || { echo "REPO-MISSING"; exit 2; }

# Clean slate: remove the partial tree from the failed wizard runs (created by
# this session's earlier attempts only), then install fresh from the local repo.
rm -rf /d/MiKTeX

./miktexsetup_standalone.exe --shared=no --package-set=complete \
  --local-package-repository="D:\\temp\\miktex_repo" \
  --user-install="D:\\MiKTeX" --user-config="D:\\MiKTeX\\Config" \
  --user-data="D:\\MiKTeX\\Data" --modify-path=yes \
  install
echo "INSTALL-EXIT=$?"

ls -la "/d/MiKTeX/miktex/bin/x64/pdflatex.exe" "/d/MiKTeX/miktex/bin/x64/bibtex.exe" || { echo "PDFLATEX-MISSING"; exit 3; }

"/d/MiKTeX/miktex/bin/x64/initexmf.exe" --set-config-value=[MPM]AutoInstall=1 || echo "AUTOINSTALL-SET-FAILED"
echo "MIKTEX-SETUP-DONE"
