$ErrorActionPreference='Stop'
$root='D:\GitHub\Bipedal_Robot'
$log=Join-Path $root 'tmp/torque_review_20261027'
$tex='D:\MiKTeX\miktex\bin\x64\pdflatex.exe'
$paperSource=Join-Path $root 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal'
$paper=Join-Path $root 'Documentation/Reports and Papers/Knee_Torque_Test/build_review_20261027'
$diss=Join-Path $root 'Documentation/Reports and Papers/Dissertation/ProofFinal'
New-Item -ItemType Directory -Path $paper -Force | Out-Null
Get-ChildItem -LiteralPath $paperSource -File | Where-Object { $_.Extension -in @('.tex','.cls','.sty','.bib','.bbl','.jpg','.png','.bst','.eps') -or $_.Name -like '*eps-converted-to.pdf' } | Copy-Item -Destination $paper
Set-Location -LiteralPath $paper
foreach($entry in @('frontiers_SupplementaryMaterial','Main','frontiers_SupplementaryMaterial','Main')) {
    $count=if(Test-Path "$entry.aux"){'pass2'}else{'pass1'}
    & $tex --disable-installer -interaction=nonstopmode -halt-on-error "$entry.tex" *> (Join-Path $log "$entry-$count.txt")
    if($LASTEXITCODE -ne 0){throw "Compile failed: $entry (see $count)"}
    Write-Output "Built $entry ($count)"
}
Set-Location -LiteralPath $diss
$build='build_r12_torque_review_20261027'
New-Item -ItemType Directory -Path "$build/chapters" -Force | Out-Null
Copy-Item -LiteralPath 'build_r11_20261027/main.bbl' -Destination "$build/main.bbl"
foreach($pass in 1..3) {
    & $tex --disable-installer -interaction=nonstopmode -halt-on-error "-output-directory=$build" main.tex *> (Join-Path $log "dissertation-pass$pass.txt")
    if($LASTEXITCODE -ne 0){throw "Dissertation compile failed ($pass)"}
    Write-Output "Built dissertation ($pass)"
    if($pass -eq 1) {
        $env:BIBINPUTS="$diss;"
        $env:BSTINPUTS="$diss;"
        Push-Location -LiteralPath $build
        try {
            & 'D:\MiKTeX\miktex\bin\x64\bibtex.exe' main *> (Join-Path $log 'dissertation-bibtex.txt')
            if($LASTEXITCODE -ne 0){throw 'Dissertation bibliography check failed'}
        } finally { Pop-Location }
    }
}
