$ErrorActionPreference = 'Stop'
$dissertationRoot = 'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation'
$sourceRoot = Join-Path $dissertationRoot 'ProofFinal'
$buildRoot = Join-Path $dissertationRoot 'Archive\build_session1_sns_20261031_v5'
if (Test-Path -LiteralPath $buildRoot) { throw 'Fresh build directory already exists; use a new revision name.' }
New-Item -ItemType Directory -Path (Join-Path $buildRoot 'chapters') -Force | Out-Null
$env:BIBINPUTS = "$sourceRoot;"
$env:BSTINPUTS = "$sourceRoot;"
$texExe = 'D:\MiKTeX\miktex\bin\x64\pdflatex.exe'
$bibExe = 'D:\MiKTeX\miktex\bin\x64\bibtex.exe'
Push-Location $sourceRoot
try {
    & $texExe '-interaction=nonstopmode' '-halt-on-error' "-output-directory=$buildRoot" 'main.tex' *> (Join-Path $buildRoot 'pass1.console.log')
    if ($LASTEXITCODE -ne 0) { Get-Content (Join-Path $buildRoot 'pass1.console.log') -Tail 28; throw 'First LaTeX pass failed' }
    Push-Location $buildRoot
    try {
        & $bibExe 'main' *> 'bibtex.console.log'
        if ($LASTEXITCODE -ne 0) { Get-Content 'bibtex.console.log' -Tail 28; throw 'BibTeX failed' }
    } finally { Pop-Location }
    foreach ($pass in 2..4) {
        & $texExe '-interaction=nonstopmode' '-halt-on-error' "-output-directory=$buildRoot" 'main.tex' *> (Join-Path $buildRoot "pass$pass.console.log")
        if ($LASTEXITCODE -ne 0) { Get-Content (Join-Path $buildRoot "pass$pass.console.log") -Tail 28; throw "LaTeX pass $pass failed" }
    }
    Get-Content (Join-Path $buildRoot 'main.log') -Tail 5
    Write-Output "Build: $buildRoot\main.pdf"
} finally { Pop-Location }
