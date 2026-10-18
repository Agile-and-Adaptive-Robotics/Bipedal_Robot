#!/usr/bin/env bash
# MiKTeX user-scope install to D:\MiKTeX on easteregg2 (Ben's go, 2026-10-17).
# Downloads the basic installer from CTAN, installs unattended user-scope,
# updates the package database, then attempts a complete package set with
# AutoInstall=1 as the fallback so compiles never stall on missing packages.
set -x
mkdir -p /d/temp/miktex_setup_e2
cd /d/temp/miktex_setup_e2 || exit 1

if [ ! -f basic-miktex-x64.exe ]; then
  curl -L --fail -o basic-miktex-x64.exe \
    "https://mirror.ctan.org/systems/win32/miktex/setup/windows-x64/basic-miktex-25.12-x64.exe" || { echo "DOWNLOAD-FAILED"; exit 2; }
fi
ls -la basic-miktex-x64.exe

./basic-miktex-x64.exe --unattended \
  --user-install="D:\\MiKTeX" --user-config="D:\\MiKTeX\\Config" \
  --user-data="D:\\MiKTeX\\Data"
INST=$?
echo "INSTALLER-EXIT=$INST"
ls -la "/d/MiKTeX/miktex/bin/x64/pdflatex.exe" || { echo "PDFLATEX-MISSING"; exit 3; }

MIK="/d/MiKTeX/miktex/bin/x64"
"$MIK/miktex.exe" packages update-package-database || echo "DB-UPDATE-FAILED"
"$MIK/initexmf.exe" --set-config-value=[MPM]AutoInstall=1 || echo "AUTOINSTALL-SET-FAILED"
"$MIK/miktex.exe" packages install-all >/dev/null 2>&1 && echo "INSTALL-ALL-OK" || echo "INSTALL-ALL-UNSUPPORTED (AutoInstall=1 fallback active)"

powershell.exe -NoProfile -Command "$p=[Environment]::GetEnvironmentVariable('Path','User'); if ($p -notlike '*D:\MiKTeX\miktex\bin\x64*') { [Environment]::SetEnvironmentVariable('Path', $p + ';D:\MiKTeX\miktex\bin\x64', 'User'); Write-Output 'PATH-ADDED' } else { Write-Output 'PATH-PRESENT' }"
echo "MIKTEX-SETUP-DONE"
