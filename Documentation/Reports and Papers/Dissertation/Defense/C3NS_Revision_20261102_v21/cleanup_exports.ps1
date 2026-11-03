$ErrorActionPreference='Stop'
$taskRoot=(Resolve-Path -LiteralPath $PSScriptRoot).Path
$defenseRoot=(Resolve-Path -LiteralPath (Join-Path $taskRoot '..')).Path
$check=Get-Content -LiteralPath (Join-Path $taskRoot 'build\delivery_checks.json') -Raw | ConvertFrom-Json
if(!$check.passed -or $check.review_pdf_pages -ne 145){throw 'Final delivery verification is incomplete'}
if((Get-FileHash -LiteralPath $check.deck -Algorithm SHA256).Hash.ToLower() -ne $check.sha256){throw 'Final hash mismatch'}
$relative=@(
 'C3NS_Visual_20261029_v11\output\Bolen_PhD_Defense_C3NS_20261029_v11.pdf',
 'C3NS_Revision_20261030_v14\output\Bolen_PhD_Defense_RoseGreen_20261030_v14.pdf',
 'C3NS_Revision_20261030_v15\output\Bolen_PhD_Defense_RoseGreen_20261030_v15.pdf',
 'C3NS_Revision_20261031_v17\output\Bolen_PhD_Defense_RoseGreen_20261031_v17.pdf',
 'C3NS_Revision_20261031_v18\output\Bolen_PhD_Defense_RoseGreen_20261031_v18.pdf',
 'C3NS_Revision_20261102_v19\output\Bolen_PhD_Defense_RoseGreen_20261102_v19.pdf',
 'C3NS_Revision_20261102_v19\output\Bolen_PhD_Defense_RoseGreen_20261102_v19.pptx',
 'C3NS_Revision_20261102_v21\build\artifact_draft.pptx',
 'C3NS_Revision_20261102_v21\build\candidate.pptx',
 'C3NS_Revision_20261102_v21\build\pre_conclusion_cleanup.pptx',
 'C3NS_Revision_20261102_v21\build\validated_without_receipt.pptx',
 'C3NS_Revision_20261102_v21\build\route_gif_native_playback.mp4'
)
$targets=@()
foreach($rel in $relative){
 $candidate=Join-Path $defenseRoot $rel
 if(Test-Path -LiteralPath $candidate){
  $resolved=(Resolve-Path -LiteralPath $candidate).Path
  if(!$resolved.StartsWith($defenseRoot+'\',[StringComparison]::OrdinalIgnoreCase)){throw 'Cleanup target outside defense workspace'}
  if((Get-Item -LiteralPath $resolved).PSIsContainer){throw 'Cleanup expects files only'}
  $targets+=$resolved
 }
}
$report=@()
foreach($target in $targets){
 $entry=[pscustomobject]@{path=$target;bytes=(Get-Item -LiteralPath $target).Length;sha256=(Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLower();reason='Superseded generated export or verified regenerable intermediate'}
 Remove-Item -LiteralPath $target
 $report+=$entry
}
$report | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (Join-Path $taskRoot 'build\cleanup.json') -Encoding utf8
[pscustomobject]@{files_removed=$report.Count;bytes_removed=($report | Measure-Object bytes -Sum).Sum;retained_powerpoint_revisions=@(20,21)} | ConvertTo-Json
