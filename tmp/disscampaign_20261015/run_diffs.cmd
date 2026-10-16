@echo off
set GIT=C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.7\resources\app\git\cmd\git.exe
set OL=D:\Github\Bipedal_Robot\tmp\disscampaign_20261015\overleaf_pull
set PF=D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal
set OUT=D:\Github\Bipedal_Robot\tmp\disscampaign_20261015\diffs
mkdir %OUT% 2>nul

for %%F in (main.tex) do "%GIT%" diff --no-index --unified=3 -- "%PF%\%%F" "%OL}\%%F" > "%OUT%\%%F.diff" 2>&1
"%GIT%" diff --no-index --unified=3 -- "%PF%\main.tex" "%OL%\main.tex" > "%OUT%\main.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\10-introduction.tex" "%OL%\chapters\10-introduction.tex" > "%OUT%\10-introduction.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\15-background.tex" "%OL%\chapters\15-background.tex" > "%OUT%\15-background.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\20-methods.tex" "%OL%\chapters\20-methods.tex" > "%OUT%\20-methods.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\30-results.tex" "%OL%\chapters\30-results.tex" > "%OUT%\30-results.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\40-discussion.tex" "%OL%\chapters\40-discussion.tex" > "%OUT%\40-discussion.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\50-futurework.tex" "%OL%\chapters\50-futurework.tex" > "%OUT%\50-futurework.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\60-conclusion.tex" "%OL%\chapters\60-conclusion.tex" > "%OUT%\60-conclusion.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\90-bib.tex" "%OL%\chapters\90-bib.tex" > "%OUT%\90-bib.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\92-AppendixA.tex" "%OL%\chapters\92-AppendixA.tex" > "%OUT%\92-AppendixA.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\93-AppendixB.tex" "%OL%\chapters\93-AppendixB.tex" > "%OUT%\93-AppendixB.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\94-AppendixC.tex" "%OL%\chapters\94-AppendixC.tex" > "%OUT%\94-AppendixC.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\01-copyright.tex" "%OL%\chapters\01-copyright.tex" > "%OUT%\01-copyright.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\02-abstract.tex" "%OL%\chapters\02-abstract.tex" > "%OUT%\02-abstract.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\04-dedication.tex" "%OL%\chapters\04-dedication.tex" > "%OUT%\04-dedication.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\chapters\06-ack.tex" "%OL%\chapters\06-ack.tex" > "%OUT%\06-ack.tex.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\BolenFrontiers22.bib" "%OL%\BolenFrontiers22.bib" > "%OUT%\BolenFrontiers22.bib.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\thesis.bib" "%OL%\thesis.bib" > "%OUT%\thesis.bib.diff"
"%GIT%" diff --no-index --unified=3 -- "%PF%\psutex.cls" "%OL%\psutex.cls" > "%OUT%\psutex.cls.diff"
echo DIFF_BATCH_DONE
for %%F in ("%OUT%\*.diff") do @echo %%~nxF : %%~zF bytes
