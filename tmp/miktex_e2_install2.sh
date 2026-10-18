#!/usr/bin/env bash
# MiKTeX user-scope install to D:\MiKTeX — take 2. The unattended wizard returns
# before its child finishes copying; the previous task exited early and the
# child was killed mid-copy (289 files, no engines). This version keeps the
# task alive and polls until the engines land, THEN runs the post-install steps.
set -x
cd /d/temp/miktex_setup_e2 || exit 1

./basic-miktex-x64.exe --unattended \
  --user-install="D:\\MiKTeX" --user-config="D:\\MiKTeX\\Config" \
  --user-data="D:\\MiKTeX\\Data"
echo "WIZARD-EXIT=$?"

# Keep this task alive while the wizard's child copies; do NOT exit early.
for i in $(seq 1 60); do
  if [ -f "/d/MiKTeX/miktex/bin/x64/pdflatex.exe" ] && [ -f "/d/MiKTeX/miktex/bin/x64/bibtex.exe" ]; then
    echo "ENGINES-PRESENT after ~$((i*15))s of polling"
    break
  fi
  sleep 15
done

ls -la "/d/MiKTeX/miktex/bin/x64/pdflatex.exe" "/d/MiKTeX/miktex/bin/x64/bibtex.exe" || { echo "PDFLATEX-STILL-MISSING"; exit 3; }

MIK="/d/MiKTeX/miktex/bin/x64"
"$MIK/initexmf.exe" --set-config-value=[MPM]AutoInstall=1 || echo "AUTOINSTALL-SET-FAILED"
"$MIK/miktex.exe" packages update-package-database || echo "DB-UPDATE-FAILED"
"$MIK/miktex.exe" packages install-all >/dev/null 2>&1 && echo "INSTALL-ALL-OK" || echo "INSTALL-ALL-UNSUPPORTED (AutoInstall=1 fallback active)"

powershell.exe -NoProfile -Command "$p=[Environment]::GetEnvironmentVariable('Path','User'); if ($p -notlike '*D:\MiKTeX\miktex\bin\x64*') { [Environment]::SetEnvironmentVariable('Path', $p + ';D:\MiKTeX\miktex\bin\x64', 'User'); Write-Output 'PATH-ADDED' } else { Write-Output 'PATH-PRESENT' }"
echo "MIKTEX-SETUP-DONE"
