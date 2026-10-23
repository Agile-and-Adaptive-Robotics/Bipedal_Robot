@echo off
setlocal
set SRC=D:\Github\Bipedal_Robot\tmp\overleaf_live_20261023
set OUT=D:\temp\diss_u14u21_20261023
cd /d "%SRC%"
"D:\MiKTeX\miktex\bin\x64\pdflatex.exe" -interaction=nonstopmode -output-directory="%OUT%" main.tex > "%OUT%\pass6.log" 2>&1
"D:\MiKTeX\miktex\bin\x64\pdflatex.exe" -interaction=nonstopmode -output-directory="%OUT%" main.tex > "%OUT%\pass7.log" 2>&1
echo REBUILD DONE
