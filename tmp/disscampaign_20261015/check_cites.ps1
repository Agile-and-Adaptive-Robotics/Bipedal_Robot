$pf = 'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal'
$ol = 'D:\Github\Bipedal_Robot\tmp\disscampaign_20261015\overleaf_pull'

Write-Output '--- steele_experimental_2023 citations in PF tex ---'
Get-ChildItem -Recurse -File $pf -Filter '*.tex' | Select-String -Pattern 'steele_experimental_2023' | ForEach-Object { $_.Filename + ':' + $_.LineNumber }

Write-Output '--- steele_experimental_2023 citations in OL tex ---'
Get-ChildItem -Recurse -File $ol -Filter '*.tex' | Select-String -Pattern 'steele_experimental_2023' | ForEach-Object { $_.Filename + ':' + $_.LineNumber }

Write-Output '--- all steele cites in PF (context) ---'
Get-ChildItem -Recurse -File $pf -Filter '*.tex' | Select-String -Pattern 'steele_' | ForEach-Object { $_.Filename + ':' + $_.LineNumber + ': ' + $_.Line.Substring([Math]::Max(0,$_.Line.IndexOf('steele_')-40), [Math]::Min(120, $_.Line.Length - [Math]::Max(0,$_.Line.IndexOf('steele_')-40))) }

Write-Output '--- graphicspath lines in PF tex files ---'
Get-ChildItem -Recurse -File $pf -Filter '*.tex' | Select-String -Pattern '\\graphicspath' | ForEach-Object { $_.Filename + ':' + $_.LineNumber + ': ' + $_.Line.Trim() }
