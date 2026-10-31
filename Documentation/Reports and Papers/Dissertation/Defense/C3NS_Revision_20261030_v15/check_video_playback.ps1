$ErrorActionPreference='Stop'
$app=New-Object -ComObject PowerPoint.Application
$initialCount=$app.Presentations.Count
$deck=$null;$show=$null
try {
    $deck=$app.Presentations.Open((Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261030_v15.pptx'),-1,0,0)
    $settings=$deck.SlideShowSettings
    $settings.ShowType=2
    $settings.RangeType=2
    $settings.StartingSlide=19
    $settings.EndingSlide=19
    $settings.AdvanceMode=1
    $settings.ShowPresenterView=0
    $show=$settings.Run()
    Start-Sleep -Milliseconds 1500
    $player=$show.View.Player(4)
    $first=[double]$player.CurrentPosition
    $firstState=[int]$player.State
    Start-Sleep -Milliseconds 1400
    $second=[double]$player.CurrentPosition
    Start-Sleep -Milliseconds 4000
    $third=[double]$player.CurrentPosition
    $shape=$deck.Slides.Item(19).Shapes|Where-Object{$_.Id -eq 4}
    $report=[ordered]@{slide=19;shape_id=4;positions_ms=@($first,$second,$third);initial_player_state=$firstState;length_ms=$shape.MediaFormat.Length;loop_until_stopped=$shape.AnimationSettings.PlaySettings.LoopUntilStopped;play_on_entry=$shape.AnimationSettings.PlaySettings.PlayOnEntry;volume=$shape.MediaFormat.Volume;native_playback_advances=($second -gt $first);loop_restarted=($third -lt $second)}
    $report|ConvertTo-Json -Depth 4|Set-Content -LiteralPath (Join-Path $PSScriptRoot 'build\native_video_playback.json') -Encoding utf8
    $report|ConvertTo-Json
    if(!$report.native_playback_advances){throw 'Native video did not advance automatically'}
} catch {
    Write-Output ('Playback check error: '+$_.Exception.Message+' at '+$_.ScriptStackTrace)
    throw
} finally {
    if($show){try{$show.View.Exit()}catch{}}
    if($deck){try{$deck.Saved=-1;$deck.Close()}catch{};[void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck)}
    try{if($initialCount -eq 0 -and $app.Presentations.Count -eq 0){$app.Quit()}}catch{}
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}
