$ErrorActionPreference='Stop'
$app=New-Object -ComObject PowerPoint.Application
$initialCount=$app.Presentations.Count
$deck=$null;$show=$null
try {
 $deck=$app.Presentations.Open((Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261031_v17.pptx'),-1,0,0)
 $settings=$deck.SlideShowSettings;$settings.ShowType=2;$settings.RangeType=2;$settings.StartingSlide=52;$settings.EndingSlide=52;$settings.AdvanceMode=1;$settings.ShowPresenterView=0
 $show=$settings.Run();Start-Sleep -Milliseconds 500
 $samples=@()
 for($i=0;$i -lt 9;$i++){
  $samples+=,[pscustomobject]@{t=$i;treadmill=[double]$show.View.Player(4).CurrentPosition;free_drop=[double]$show.View.Player(14).CurrentPosition}
  Start-Sleep -Milliseconds 1000
 }
 $report=@()
 foreach($id in @(4,14)){
  $key=if($id -eq 4){'treadmill'}else{'free_drop'};$values=@($samples|ForEach-Object{$_.$key});$advances=$false;$restarted=$false
  for($j=1;$j -lt $values.Count;$j++){if($values[$j]-gt $values[$j-1]){$advances=$true};if($values[$j]-lt $values[$j-1]-100){$restarted=$true}}
  $shape=$deck.Slides.Item(52).Shapes|Where-Object{$_.Id -eq $id}
  $report+=,[pscustomobject]@{shape_id=$id;length_ms=$shape.MediaFormat.Length;positions_ms=$values;advances=$advances;loop_restarted=$restarted;loop_until_stopped=$shape.AnimationSettings.PlaySettings.LoopUntilStopped;play_on_entry=$shape.AnimationSettings.PlaySettings.PlayOnEntry;volume=$shape.MediaFormat.Volume}
 }
 $report|ConvertTo-Json -Depth 5|Set-Content (Join-Path $PSScriptRoot 'build\native_video_playback.json') -Encoding utf8
 $report|ConvertTo-Json -Depth 5
 if(@($report|Where-Object{!$_.advances -or !$_.loop_restarted}).Count){throw 'A recording failed native playback or loop check'}
}finally{
 if($show){try{$show.View.Exit()}catch{}}
 if($deck){try{$deck.Saved=-1;$deck.Close()}catch{};[void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck)}
 if($initialCount -eq 0 -and $app.Presentations.Count -eq 0){$app.Quit()}
 [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}
