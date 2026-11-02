$ErrorActionPreference = 'Stop'
$inputPath = Join-Path $PSScriptRoot 'build\candidate.pptx'
$outputPath = Join-Path $PSScriptRoot 'build\candidate_powerpoint.pptx'
$app = New-Object -ComObject PowerPoint.Application
$initialCount = $app.Presentations.Count
$deck = $null
try {
    $deck = $app.Presentations.Open($inputPath, -1, 0, 0)
    foreach ($slide in $deck.Slides) {
        foreach ($shape in $slide.Shapes) {
            if ($shape.HasChart -eq -1) {
                try {
                    $shape.Chart.ChartData.Activate()
                    $workbook = $shape.Chart.ChartData.Workbook
                    $workbook.Close($true)
                    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($workbook)
                } catch {
                    Write-Warning ("Chart data activation failed on slide {0}: {1}" -f $slide.SlideIndex, $_.Exception.Message)
                }
            }
        }
    }
    $deck.SaveAs($outputPath, 24)
    Write-Output $outputPath
} finally {
    if ($null -ne $deck) { $deck.Close(); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck) }
    if ($initialCount -eq 0 -and $app.Presentations.Count -eq 0) { $app.Quit() }
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}
