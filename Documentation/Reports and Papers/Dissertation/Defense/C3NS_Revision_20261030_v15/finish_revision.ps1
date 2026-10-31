$ErrorActionPreference='Stop'
$parent=Split-Path $PSScriptRoot
$baselineRoot=Join-Path $parent 'C3NS_Revision_20261030_v14'
$final=Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261030_v15.pptx'
$preservation=Get-Content -LiteralPath (Join-Path $PSScriptRoot 'build\preservation.json') -Raw|ConvertFrom-Json
$visual=Get-Content -LiteralPath (Join-Path $PSScriptRoot 'build\visual_review.json') -Raw|ConvertFrom-Json
$playback=Get-Content -LiteralPath (Join-Path $PSScriptRoot 'build\native_video_playback.json') -Raw|ConvertFrom-Json
if(!$preservation.passed -or $visual.unchanged_pixel_matches -ne 81 -or !$playback.native_playback_advances -or !$playback.loop_restarted){throw 'Incomplete verification'}
$visual|Add-Member -NotePropertyName full_size_inspected_changed_slides -NotePropertyValue @(1,19) -Force
$visual|Add-Member -NotePropertyName pdf_pages -NotePropertyValue 83 -Force
$visual|Add-Member -NotePropertyName powerpoint_sha256 -NotePropertyValue (Get-FileHash -LiteralPath $final).Hash -Force
$visual|ConvertTo-Json -Depth 7|Set-Content -LiteralPath (Join-Path $PSScriptRoot 'build\visual_review.json') -Encoding utf8
$meta=Get-Content -LiteralPath (Join-Path $baselineRoot 'current_deck.json') -Raw|ConvertFrom-Json
$meta.revision=15
$meta.current_deck=$final
$meta.pdf=Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261030_v15.pdf'
$meta.sha256=(Get-FileHash -LiteralPath $final).Hash.ToLowerInvariant()
$meta.baseline=Join-Path $baselineRoot 'output\Bolen_PhD_Defense_RoseGreen_20261030_v14.pptx'
$meta.baseline_sha256=(Get-FileHash -LiteralPath $meta.baseline).Hash.ToLowerInvariant()
$meta.status='review_copy_pending_released_video_and_flexor_photo'
$meta.retained_final_revisions=@(14,15)
$meta.offline_app=Join-Path $PSScriptRoot 'output\sadb_app.html'
$meta.validation='build/validation_delivery2.json'
$meta|Add-Member -NotePropertyName original_title_honeycomb_restored -NotePropertyValue $true -Force
$meta|Add-Member -NotePropertyName honeycomb_source -NotePropertyValue 'Original C3NS template artwork, ppt/media/image3.png, preserved byte for byte; native Rose green recoloring.' -Force
$meta|Add-Member -NotePropertyName archived_walker_video_embedded -NotePropertyValue $true -Force
$meta|Add-Member -NotePropertyName video_playback_verified -NotePropertyValue $true -Force
$meta|Add-Member -NotePropertyName video_duration_seconds -NotePropertyValue 5.4 -Force
$meta|Add-Member -NotePropertyName released_walker_video_pending -NotePropertyValue $true -Force
$pointer=Join-Path $parent 'current_deck.json'
$prior=Get-Content -LiteralPath $pointer -Raw|ConvertFrom-Json
if($prior.revision -notin @(14,15)){throw 'A different current revision appeared; reconcile before replacing pointer'}
$json=$meta|ConvertTo-Json -Depth 8
$json|Set-Content -LiteralPath (Join-Path $PSScriptRoot 'current_deck.json') -Encoding utf8
$json|Set-Content -LiteralPath $pointer -Encoding utf8
$cleanup=Get-Content -LiteralPath (Join-Path $baselineRoot 'cleanup_versions.ps1') -Raw
$retained=@'
$retained=@(
    (Join-Path $defenseRoot 'C3NS_Revision_20261030_v14\output\Bolen_PhD_Defense_RoseGreen_20261030_v14.pptx'),
    (Join-Path $defenseRoot 'C3NS_Revision_20261030_v15\output\Bolen_PhD_Defense_RoseGreen_20261030_v15.pptx')
)
'@
$cleanup=[regex]::Replace($cleanup,'(?s)\$retained=@\(.*?\r?\n\)',[System.Text.RegularExpressions.MatchEvaluator]{param($m)$retained})
$cleanup=$cleanup.Replace("'C3NS_Revision_20261030_v14'","'C3NS_Revision_20261030_v14','C3NS_Revision_20261030_v15'")
Set-Content -LiteralPath (Join-Path $PSScriptRoot 'cleanup_versions.ps1') -Value $cleanup -Encoding utf8
[pscustomobject]@{revision=15;sha256=$meta.sha256;bytes=(Get-Item -LiteralPath $final).Length;core_slides=70;question_slides=13;planned_minutes=35;native_video_verified=$true}|ConvertTo-Json
