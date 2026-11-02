$ErrorActionPreference = 'Stop'
$baseline = Join-Path $PSScriptRoot '..\C3NS_Revision_20261031_v17\output\Bolen_PhD_Defense_RoseGreen_20261031_v17.pptx'
$sourcePath = Join-Path $PSScriptRoot 'build\Bolen_PhD_Defense_RoseGreen_20261031_v18_artifact_draft.pptx'
$candidate = Join-Path $PSScriptRoot 'build\candidate.pptx'
Copy-Item -LiteralPath $baseline -Destination $candidate -Force

$app = New-Object -ComObject PowerPoint.Application
$initialCount = $app.Presentations.Count
$dest = $null
try {
    $dest = $app.Presentations.Open($candidate, -1, 0, 0)

    # Replace only the visually rebuilt slides. The baseline package remains the
    # owner of its native charts, workbooks, hidden-slide state, GIFs, and movies.
    $replace = @(69,65,64,60,59,58,57,56,55,54,53,49,11,8,6,4,3)
    foreach ($number in $replace) {
        [void]$dest.Slides.InsertFromFile($sourcePath, $number, $number, $number)
        $dest.Slides.Item($number).Delete()
        Write-Output ("Replaced slide {0}" -f $number)
    }

    # Slide 44 keeps the validated native PowerPoint chart and embedded workbook,
    # but uses the committee-submitted dissertation values.
    $slide = $dest.Slides.Item(44)
    $title = $slide.Shapes | Where-Object { $_.Top -lt 40 -and $_.Width -gt 500 -and $_.HasTextFrame -eq -1 } | Select-Object -First 1
    $title.TextFrame.TextRange.Text = 'Route redesign raises mean torque coverage from 37.35% to 65.81%'
    $title.TextFrame.TextRange.Font.Name = 'Arial'
    $title.TextFrame.TextRange.Font.Size = 30
    $title.TextFrame.TextRange.Font.Bold = -1
    $title.TextFrame.TextRange.Font.Color.RGB = 2379031

    $chartShape = $slide.Shapes | Where-Object { $_.HasChart -eq -1 } | Select-Object -First 1
    $chart = $chartShape.Chart
    $series1 = $chart.SeriesCollection(1)
    $series2 = $chart.SeriesCollection(2)
    $series1.Values = [double[]]@(0.6581, 0.9931)
    $series2.Values = [double[]]@(0.3735, 0.9300)
    $series1.Name = 'Redesign'
    $series2.Name = 'Baseline'
    [void]$series1.ApplyDataLabels()
    [void]$series2.ApplyDataLabels()
    $series1.DataLabels().Position = 3
    $series2.DataLabels().Position = 3
    $series1.DataLabels().NumberFormat = '0.00%'
    $series2.DataLabels().NumberFormat = '0.00%'
    $series1.DataLabels().Font.Color = 16777215
    $series2.DataLabels().Font.Color = 0
    $series1.DataLabels().Font.Bold = -1
    $series2.DataLabels().Font.Bold = -1

    try {
        $notes = $slide.NotesPage.Shapes.Placeholders.Item(2).TextFrame.TextRange
        $notes.Text = @'
The submitted dissertation corrects the whole-leg comparison: mean torque coverage rises from 37.35% to 65.81%, while mean active-force coverage rises from 93.00% to 99.31%. Emphasize that route redesign expands the joint-angle conditions in which the robot model meets or exceeds the human benchmark.

Source: Submitted dissertation, Chapter 4, Table 4.12, and Appendix C.
'@
    } catch {
        Write-Warning ('Could not update slide 44 notes: ' + $_.Exception.Message)
    }

    $dest.Save()
    Write-Output ([pscustomobject]@{
        path = $candidate
        slides = $dest.Slides.Count
        hidden = @($dest.Slides | Where-Object { $_.SlideShowTransition.Hidden -eq -1 }).Count
        charts = @($dest.Slides | ForEach-Object { $_.Shapes | Where-Object { $_.HasChart -eq -1 } }).Count
    } | ConvertTo-Json)
} finally {
    if ($null -ne $dest) { $dest.Close(); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($dest) }
    if ($initialCount -eq 0 -and $app.Presentations.Count -eq 0) { $app.Quit() }
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}
