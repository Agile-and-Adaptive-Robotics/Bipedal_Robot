$ErrorActionPreference='Stop'
$app=New-Object -ComObject PowerPoint.Application
$initialCount=$app.Presentations.Count
$deck=$null
try {
    $deck=$app.Presentations.Add(0)
    $deck.PageSetup.SlideWidth=720
    $deck.PageSetup.SlideHeight=405
    $slide=$deck.Slides.Add(1,12)
    $height=210
    $width=$height*650/688
    $video=$slide.Shapes.AddMediaObject2((Join-Path $PSScriptRoot 'build\w2l_archived_motion.mp4'),0,-1,66,108,$width,$height)
    $video.Name='Original W2L AnimatLab recording'
    $video.AlternativeText='Archived AnimatLab recording of the W2L model moving its legs.'
    $video.MediaFormat.Volume=0
    $video.MediaFormat.SetDisplayPictureFromFile((Join-Path $PSScriptRoot 'build\w2l_poster.png'))
    $video.AnimationSettings.PlaySettings.PlayOnEntry=-1
    $video.AnimationSettings.PlaySettings.LoopUntilStopped=-1
    $video.AnimationSettings.PlaySettings.StopAfterSlides=1
    $slide.TimeLine.MainSequence.Item(1).Timing.TriggerType=2
    $slide.TimeLine.MainSequence.Item(1).Timing.TriggerDelayTime=0
    $video.MediaFormat.Volume=0
    $deck.SaveAs((Join-Path $PSScriptRoot 'build\video_component.pptx'),24)
    Write-Output 'Embedded archived recording; muted, automatic playback, looping.'
} finally {
    if($deck){$deck.Close();[void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck)}
    if($initialCount -eq 0 -and $app.Presentations.Count -eq 0){$app.Quit()}
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}
