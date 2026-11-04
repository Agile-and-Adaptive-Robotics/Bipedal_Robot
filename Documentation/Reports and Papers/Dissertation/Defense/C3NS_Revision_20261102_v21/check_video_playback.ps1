$ErrorActionPreference = 'Stop'
$deckPath = Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261102_v21.pptx'
$reportPath = Join-Path $PSScriptRoot 'build\native_video_playback.json'
$originalDeckHash = (Get-FileHash -LiteralPath $deckPath -Algorithm SHA256).Hash
$sourceOrder = @(Get-Content -LiteralPath (Join-Path $PSScriptRoot 'build\slide_order.json') -Raw | ConvertFrom-Json)

function Invoke-MediaSession {
    param([ValidateSet('walker','sadb')][string]$Mode)
    $app = New-Object -ComObject PowerPoint.Application
    $initialCount = $app.Presentations.Count
    $deck = $null
    $show = $null
    $slide = $null
    $originalHidden = $null
    try {
        # A no-window presentation may close when its slideshow exits. Each test
        # opens its own read-only presentation and completes before Exit.
        $deck = $app.Presentations.Open($deckPath, -1, 0, 0)
        $sourceNumber = if ($Mode -eq 'walker') { 56 } else { 127 }
        $number = [array]::IndexOf($sourceOrder,$sourceNumber)+1
        $ids = if ($Mode -eq 'walker') { @(4,14) } else { @(15) }
        $slide = $deck.Slides.Item($number)
        $originalHidden = $slide.SlideShowTransition.Hidden
        if ($Mode -eq 'sadb') { $slide.SlideShowTransition.Hidden = 0 }
        $shapes = @{}
        foreach ($id in $ids) {
            $shapes[$id] = $slide.Shapes | Where-Object { $_.Id -eq $id }
            if (!$shapes[$id]) { throw "Movie shape $id is absent from slide $number" }
        }
        $settings = $deck.SlideShowSettings
        $settings.ShowType = 2
        $settings.RangeType = 2
        $settings.StartingSlide = $number
        $settings.EndingSlide = $number
        $settings.AdvanceMode = 1
        $settings.ShowPresenterView = 0
        $show = $settings.Run()
        Start-Sleep -Milliseconds 500
        # Walker movies must autoplay from their saved settings. The authentic
        # offline walkthrough is allowed to be click-to-play.
        if ($Mode -eq 'sadb') { $show.View.Player(15).Play() }
        $positions = @{}
        foreach ($id in $ids) { $positions[$id] = @() }
        $sampleCount = if ($Mode -eq 'walker') { 9 } else { 5 }
        for ($i = 0; $i -lt $sampleCount; $i++) {
            foreach ($id in $ids) {
                $positions[$id] += [double]$show.View.Player($id).CurrentPosition
            }
            if ($i -lt $sampleCount - 1) { Start-Sleep -Milliseconds 1000 }
        }
        $rows = @()
        foreach ($id in $ids) {
            $values = $positions[$id]
            $advances = $false
            $restarted = $false
            for ($j = 1; $j -lt $values.Count; $j++) {
                if ($values[$j] -gt $values[$j-1]) { $advances = $true }
                if ($values[$j] -lt $values[$j-1] - 100) { $restarted = $true }
            }
            $name = switch ($id) { 4 {'treadmill'} 14 {'free_drop'} 15 {'sadb_offline_walkthrough'} }
            $shape = $shapes[$id]
            $rows += ,[pscustomobject]@{
                slide_number=$number
                recording=$name
                shape_id=$id
                length_ms=$shape.MediaFormat.Length
                positions_ms=$values
                advances=$advances
                loop_required=($Mode -eq 'walker')
                loop_restarted=$restarted
                loop_until_stopped=$shape.AnimationSettings.PlaySettings.LoopUntilStopped
                play_on_entry=$shape.AnimationSettings.PlaySettings.PlayOnEntry
                volume=$shape.MediaFormat.Volume
                originally_hidden=$originalHidden
                presentations_before_session=$initialCount
            }
        }
        # Persist walker evidence before Exit. Append SADb after its own session.
        $previous = if ($Mode -eq 'sadb') { @(Get-Content -LiteralPath $reportPath -Raw | ConvertFrom-Json) } else { @() }
        @($previous + $rows) | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $reportPath -Encoding utf8
    } finally {
        if ($slide -and $null -ne $originalHidden) {
            try { $slide.SlideShowTransition.Hidden = $originalHidden } catch {}
        }
        if ($deck) { try { $deck.Saved = -1 } catch {} }
        if ($show) { try { $show.View.Exit() } catch {} }
        if ($deck) {
            try { $deck.Close() } catch {}
            try { [void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck) } catch {}
        }
        # Never close an unrelated presentation or quit an app holding one.
        try { if ($initialCount -eq 0 -and $app.Presentations.Count -eq 0) { $app.Quit() } } catch {}
        try { [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app) } catch {}
    }
}

Invoke-MediaSession -Mode walker
Invoke-MediaSession -Mode sadb
$report = @(Get-Content -LiteralPath $reportPath -Raw | ConvertFrom-Json)
$report | ConvertTo-Json -Depth 5
if ((Get-FileHash -LiteralPath $deckPath -Algorithm SHA256).Hash -ne $originalDeckHash) {
    throw 'The validation session unexpectedly changed the saved deck'
}
if ($report.Count -ne 3 -or @($report | Where-Object { !$_.advances -or ($_.loop_required -and !$_.loop_restarted) }).Count) {
    throw 'A recording failed its native playback or required loop check'
}

