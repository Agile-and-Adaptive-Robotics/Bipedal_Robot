@echo off
rem Lock-safe dissertation build v2: pdflatex+bibtex into D:\temp\diss_build, then copy back to build_r8.
taskkill /IM acrotray.exe /F >nul 2>&1
set PF=D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal
set OUT=D:\temp\diss_build
set TEX=D:\MiKTeX\miktex\bin\x64
rmdir /s /q "%OUT%" >nul 2>&1
mkdir "%OUT%" >nul 2>&1
mkdir "%OUT%\chapters" >nul 2>&1
cd /d "%PF%"
"%TEX%\pdflatex.exe" -interaction=nonstopmode -halt-on-error -output-directory="%OUT%" main.tex > "%OUT%\p1.log" 2>&1
if errorlevel 1 goto fail
cd /d "%OUT%"
set BIBINPUTS=%PF%
set BSTINPUTS=%PF%
"%TEX%\bibtex.exe" main > "%OUT%\bib.log" 2>&1
cd /d "%PF%"
"%TEX%\pdflatex.exe" -interaction=nonstopmode -output-directory="%OUT%" main.tex > "%OUT%\p2.log" 2>&1
"%TEX%\pdflatex.exe" -interaction=nonstopmode -output-directory="%OUT%" main.tex > "%OUT%\p3.log" 2>&1
findstr /c:"^!" "%OUT%\p3.log" > "%OUT%\errors.txt"
findstr /c:"undefined" "%OUT%\p3.log" > "%OUT%\undef.txt"
echo --- errors --- & type "%OUT%\errors.txt"
echo --- undefined --- & type "%OUT%\undef.txt"
for /f "tokens=2 delims=:" %%a in ('findstr /c:"Output written" "%OUT%\p3.log"') do set PAGES=%%a
echo PAGES:%PAGES%
copy /y "%OUT%\main.pdf" "%PF%\build_r8\main.pdf" >nul
copy /y "%OUT%\main.aux" "%OUT%\main.bbl" "%OUT%\main.blg" "%OUT%\main.log" "%OUT%\main.toc" "%OUT%\main.lof" "%OUT%\main.lot" "%OUT%\main.out" "%PF%\build_r8\" >nul 2>&1
xcopy /e /i /y "%OUT%\chapters" "%PF%\build_r8\chapters" >nul 2>&1
echo BUILD_OK
goto :eof
:fail
echo BUILD_FAIL - see %OUT%\p1.log
