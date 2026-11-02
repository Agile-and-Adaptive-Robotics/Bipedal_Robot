$ErrorActionPreference='Stop'
$app=New-Object -ComObject PowerPoint.Application
$initial=$app.Presentations.Count
$deck=$null
$sourceOrder = @(Get-Content -LiteralPath (Join-Path $PSScriptRoot 'build/slide_order.json') -Raw | ConvertFrom-Json)
$gifSlides = @(46,47,49,67 | ForEach-Object { [array]::IndexOf($sourceOrder,$_)+1 })
try {
  $deck=$app.Presentations.Open((Join-Path $PSScriptRoot 'output/Bolen_PhD_Defense_RoseGreen_20261102_v21.pptx'),-1,0,0)
  # Keep original native slide objects and their animations in this unsaved test.
  foreach($n in ($deck.Slides.Count..1)) {if($n -notin $gifSlides){$deck.Slides.Item($n).Delete()}}
  foreach($s in $deck.Slides){$s.SlideShowTransition.Hidden=0}
  $path=Join-Path $PSScriptRoot 'build/route_gif_native_playback.mp4'
  $deck.CreateVideo($path,$false,12,720,20,80)
  $start=Get-Date
  do {
    Start-Sleep -Seconds 2
    $status=$deck.CreateVideoStatus
    if(((Get-Date)-$start).TotalSeconds -gt 300){throw 'Native animation export exceeded five minutes'}
  } while($status -in @(1,2))
  if($status -ne 3){throw "Native animation export status $status"}
  [pscustomobject]@{source_slides=$gifSlides;baseline_slides=@(46,47,49,67);seconds_per_slide=12;status=$status;output=$path;bytes=(Get-Item -LiteralPath $path).Length} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $PSScriptRoot 'build/native_gif_export.json')
  Get-Content -LiteralPath (Join-Path $PSScriptRoot 'build/native_gif_export.json')
} finally {
  if($deck){$deck.Saved=-1;$deck.Close();[void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck)}
  if($initial -eq 0 -and $app.Presentations.Count -eq 0){$app.Quit()}
  [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}
