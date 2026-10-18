$candidates = @(
  'C:\Program Files\MiKTeX\miktex\bin\x64\pdflatex.exe',
  'C:\Users\Ben Bolen\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe',
  'C:\Users\Ben Bolen\AppData\Local\MiKTeX\miktex\bin\x64\pdflatex.exe',
  'D:\MiKTeX\miktex\bin\x64\pdflatex.exe',
  'D:\Program Files\MiKTeX\miktex\bin\x64\pdflatex.exe',
  'C:\texlive\2026\bin\windows\pdflatex.exe',
  'D:\texlive\2026\bin\windows\pdflatex.exe'
)
foreach ($x in $candidates) { if (Test-Path $x) { Write-Output "HIT $x" } }
$cmd = Get-Command pdflatex -ErrorAction SilentlyContinue
if ($cmd) { Write-Output ("ONPATH " + $cmd.Source) }
$userPath = [Environment]::GetEnvironmentVariable('Path','User') + ';' + [Environment]::GetEnvironmentVariable('Path','Machine')
foreach ($seg in ($userPath -split ';' | Where-Object { $_ -match 'miktex|texlive' })) { Write-Output ("PATHSEG " + $seg) }
Write-Output 'SEARCH-DONE'
