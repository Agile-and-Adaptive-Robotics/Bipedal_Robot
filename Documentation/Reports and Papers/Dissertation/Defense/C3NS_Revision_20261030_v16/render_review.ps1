param(
    [string]$DeckPath = (Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261030_v16.pptx'),
    [int[]]$Selection = @(),
    [switch]$MakePdf
)
$ErrorActionPreference='Stop'
$renderPath=Join-Path $PSScriptRoot 'build\render'
New-Item -ItemType Directory -Force -Path $renderPath | Out-Null
$app=New-Object -ComObject PowerPoint.Application
$initialPresentationCount=$app.Presentations.Count
$deck=$null
$report=[System.Collections.Generic.List[object]]::new()
try {
    $deck=$app.Presentations.Open($DeckPath,-1,0,0)
    $slideCount=$deck.Slides.Count
    if(!$Selection.Count){$Selection=1..$slideCount}
    foreach($number in $Selection){
        $slide=$deck.Slides.Item($number)
        $slide.Export((Join-Path $renderPath ('slide_{0:D2}.png' -f $number)),'PNG',1600,900)
        foreach($shape in $slide.Shapes){
            if($shape.HasTextFrame -eq -1 -and $shape.TextFrame.HasText -eq -1){
                $range=$shape.TextFrame2.TextRange
                $report.Add([pscustomobject]@{
                    slide=$number;name=$shape.Name;text=$shape.TextFrame.TextRange.Text;
                    x=[double]$shape.Left;y=[double]$shape.Top;w=[double]$shape.Width;h=[double]$shape.Height;
                    bound_x=[double]$range.BoundLeft;bound_y=[double]$range.BoundTop;
                    bound_w=[double]$range.BoundWidth;bound_h=[double]$range.BoundHeight
                })
            }
        }
        if($number % 10 -eq 0){Write-Output ('Rendered '+$number+' / '+$slideCount)}
    }
    $geometryName=if($Selection.Count -eq $slideCount){'native_text_geometry.json'}else{'native_text_geometry_subset.json'}
    $report | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $PSScriptRoot ('build\'+$geometryName)) -Encoding utf8
    if($MakePdf){
        $pdfPath=Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261030_v16.pdf'
        $hiddenSlides=@($deck.Slides|Where-Object{$_.SlideShowTransition.Hidden -eq -1}|ForEach-Object{$_.SlideIndex})
        foreach($hiddenNumber in $hiddenSlides){$deck.Slides.Item($hiddenNumber).SlideShowTransition.Hidden=0}
        $deck.PrintOptions.PrintHiddenSlides=-1
        $printRange=$deck.PrintOptions.Ranges.Add(1,$slideCount)
        $deck.ExportAsFixedFormat($pdfPath,2,1,0,1,1,-1,$printRange,1,'',-1,-1,-1,-1,0)
        foreach($hiddenNumber in $hiddenSlides){$deck.Slides.Item($hiddenNumber).SlideShowTransition.Hidden=-1}
        $deck.Saved=-1
        Write-Output ('PDF: '+$pdfPath)
    }
    [pscustomobject]@{rendered=$Selection.Count;slides=$slideCount;hidden=(@($deck.Slides|Where-Object{$_.SlideShowTransition.Hidden -eq -1})).Count;render_path=$renderPath} | ConvertTo-Json
}finally{
    if($null -ne $deck){$deck.Close();[void][Runtime.InteropServices.Marshal]::ReleaseComObject($deck)}
    if($initialPresentationCount -eq 0 -and $app.Presentations.Count -eq 0){$app.Quit()}
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}


