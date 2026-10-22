@echo off
rem Build one dissertation tree with MiKTeX (pdflatex, bibtex, pdflatex x2)
rem usage: build_one.bat <source-dir> <build-dir>
setlocal
set SRCDIR=%~1
set BUILD=%~2
set PDFLATEX=D:\MiKTeX\miktex\bin\x64\pdflatex.exe
set BIBTEX=D:\MiKTeX\miktex\bin\x64\bibtex.exe
set SRCBIB=%SRCDIR%

if not exist "%BUILD%" mkdir "%BUILD%"
cd /d "%SRCDIR%"

echo === PASS 1 ===
"%PDFLATEX%" -interaction=nonstopmode -halt-on-error -output-directory="%BUILD%" main.tex
if errorlevel 1 goto :fail

echo === PASS 2 (populate aux) ===
"%PDFLATEX%" -interaction=nonstopmode -halt-on-error -output-directory="%BUILD%" main.tex
if errorlevel 1 goto :fail

echo === BIBTEX ===
set BIBINPUTS=%SRCBIB%
set BSTINPUTS=%SRCBIB%
pushd "%BUILD%"
"%BIBTEX%" main
if errorlevel 1 goto :fail
popd

echo === PASS 3 ===
"%PDFLATEX%" -interaction=nonstopmode -halt-on-error -output-directory="%BUILD%" main.tex
if errorlevel 1 goto :fail

echo === PASS 4 ===
"%PDFLATEX%" -interaction=nonstopmode -halt-on-error -output-directory="%BUILD%" main.tex
if errorlevel 1 goto :fail

echo === BUILD OK ===
exit /b 0

:fail
echo === BUILD FAILED ===
exit /b 1
