$ErrorActionPreference = 'Stop'
$app = New-Object -ComObject PowerPoint.Application
$decks = @()
try {
    $paths = @(
        (Join-Path $PSScriptRoot '..\C3NS_Revision_20261031_v17\output\Bolen_PhD_Defense_RoseGreen_20261031_v17.pptx'),
        (Join-Path $PSScriptRoot 'build\Bolen_PhD_Defense_RoseGreen_20261031_v18_artifact_draft.pptx')
    )
    foreach ($path in $paths) {
        $deck = $app.Presentations.Open($path, -1, 0, 0)
        $decks += $deck
        Write-Output "DECK $path"
        foreach ($shape in $deck.Slides.Item(44).Shapes) {
            $text = ''
            if ($shape.HasTextFrame -eq -1 -and $shape.TextFrame.HasText -eq -1) {
                $text = ($shape.TextFrame.TextRange.Text -replace "`r|`n", ' ')
            }
            [pscustomobject]@{Id=$shape.Id;Name=$shape.Name;Type=$shape.Type;HasChart=$shape.HasChart;Left=[math]::Round($shape.Left,1);Top=[math]::Round($shape.Top,1);Width=[math]::Round($shape.Width,1);Height=[math]::Round($shape.Height,1);Text=$text}
        }
        $deck.Close()
    }
} finally {
    foreach ($deck in $decks) { try { [void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck) } catch {} }
    if ($app.Presentations.Count -eq 0) { $app.Quit() }
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}
