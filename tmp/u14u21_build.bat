@echo off
setlocal
set SRC=D:\Github\Bipedal_Robot\tmp\overleaf_live_20261023
set OUT=D:\temp\diss_u14u21_20261023
if exist "%OUT%" rmdir /s /q "%OUT%"
mkdir "%OUT%"
cd /d "%SRC%"
set BIBINPUTS=%SRC%
set BSTINPUTS=%SRC%
"D:\MiKTeX\miktex\bin\x64\pdflatex.exe" -interaction=nonstopmode -output-directory="%OUT%" main.tex > "%OUT%\pass1.log" 2>&1
"D:\MiKTeX\miktex\bin\x64\bibtex.exe" "%OUT%\main" > "%OUT%\bibtex.log" 2>&1
"D:\MiKTeX\miktex\bin\x64\pdflatex.exe" -interaction=nonstopmode -output-directory="%OUT%" main.tex > "%OUT%\pass2.log" 2>&1
"D:\MiKTeX\miktex\bin\x64\pdflatex.exe" -interaction=nonstopmode -output-directory="%OUT%" main.tex > "%OUT%\pass3.log" 2>&1
echo BUILD SCRIPT DONE
