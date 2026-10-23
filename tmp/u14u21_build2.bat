@echo off
setlocal
set SRC=D:\Github\Bipedal_Robot\tmp\overleaf_live_20261023
set OUT=D:\temp\diss_u14u21_20261023
set BIBINPUTS=%SRC%
set BSTINPUTS=%SRC%
cd /d "%OUT%"
"D:\MiKTeX\miktex\bin\x64\bibtex.exe" main > "%OUT%\bibtex2.log" 2>&1
cd /d "%SRC%"
"D:\MiKTeX\miktex\bin\x64\pdflatex.exe" -interaction=nonstopmode -output-directory="%OUT%" main.tex > "%OUT%\pass4.log" 2>&1
"D:\MiKTeX\miktex\bin\x64\pdflatex.exe" -interaction=nonstopmode -output-directory="%OUT%" main.tex > "%OUT%\pass5.log" 2>&1
echo BIBTEX RERUN DONE
