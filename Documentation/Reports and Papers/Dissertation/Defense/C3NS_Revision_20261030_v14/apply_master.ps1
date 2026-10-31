param([string]$SourceDeck = 'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\C3NS_Revision_20261030_v13\output\Bolen_PhD_Defense_C3NS_20261030_v13_review.pptx')
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.IO.Compression.FileSystem
$workRoot = $PSScriptRoot
$draftPath = Join-Path $workRoot 'build\candidate.pptx'
$drawingNs = 'http://schemas.openxmlformats.org/drawingml/2006/main'
$presentationNs = 'http://schemas.openxmlformats.org/presentationml/2006/main'
$relationshipNs = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
$deepGreen = '174D24'
$leafGreen = '58D51B'
$huntMinimal = 'A. J. Hunt, INRC2 Meeting, Intel Labs (Feb 15, 2019).'
$huntFull = 'Source: Hunt, A. J. (2019, February 15). Real-time control of dynamic walking in a dog-like robot. Presented at the 2nd Intel Neuromorphic Research Community (INRC2) Meeting.'
$emuPerPoint = 12700
$changes = [System.Collections.Generic.List[object]]::new()
function Get-XmlDocument([string]$text) {
    $doc = [xml]::new()
    $doc.PreserveWhitespace = $true
    $doc.LoadXml($text)
    return ,$doc
}
function Set-Transform($node, [double]$x, [double]$y, [double]$width, [double]$height) {
    $xf = $node.SelectSingleNode('./*[local-name()="spPr"]/*[local-name()="xfrm"]')
    if (!$xf) { throw "Missing transform on $($node.OuterXml.Substring(0,80))" }
    $off = $xf.SelectSingleNode('./*[local-name()="off"]')
    $ext = $xf.SelectSingleNode('./*[local-name()="ext"]')
    $off.SetAttribute('x', [string][long][math]::Round($x * $emuPerPoint))
    $off.SetAttribute('y', [string][long][math]::Round($y * $emuPerPoint))
    $ext.SetAttribute('cx', [string][long][math]::Round($width * $emuPerPoint))
    $ext.SetAttribute('cy', [string][long][math]::Round($height * $emuPerPoint))
}
function Set-TextColor($node, [string]$color, [int]$size = 0) {
    foreach ($textRun in $node.SelectNodes('.//*[local-name()="r" or local-name()="fld"]')) {
        if (!$textRun.SelectSingleNode('./*[local-name()="rPr"]')) {
            $properties=$node.OwnerDocument.CreateElement('a','rPr',$drawingNs)
            [void]$textRun.PrependChild($properties)
        }
    }
    foreach ($run in $node.SelectNodes('.//*[local-name()="rPr" or local-name()="defRPr" or local-name()="endParaRPr"]')) {
        foreach ($oldFill in @($run.SelectNodes('./*[local-name()="solidFill" or local-name()="gradFill" or local-name()="noFill"]'))) { [void]$run.RemoveChild($oldFill) }
        $fill = $node.OwnerDocument.CreateElement('a','solidFill',$drawingNs)
        $rgb = $node.OwnerDocument.CreateElement('a','srgbClr',$drawingNs)
        $rgb.SetAttribute('val',$color)
        [void]$fill.AppendChild($rgb)
        $firstAfterLine=$run.FirstChild
        if($firstAfterLine -and $firstAfterLine.LocalName -eq 'ln'){$firstAfterLine=$firstAfterLine.NextSibling}
        if($firstAfterLine){[void]$run.InsertBefore($fill,$firstAfterLine)}else{[void]$run.AppendChild($fill)}
        if ($size -gt 0) { $run.SetAttribute('sz',[string]$size) }
    }
}
function Add-Fragment($doc, $parent, [string]$fragment) {
    $container = $doc.CreateDocumentFragment()
    $container.InnerXml = $fragment
    [void]$parent.AppendChild($container)
}
function New-BandXml([int]$id, [string]$name, [double]$y, [double]$height, [string]$color) {
    $yEmu = [long][math]::Round($y * $emuPerPoint)
    $heightEmu = [long][math]::Round($height * $emuPerPoint)
    return "<p:sp xmlns:p='$presentationNs' xmlns:a='$drawingNs'><p:nvSpPr><p:cNvPr id='$id' name='$name'/><p:cNvSpPr/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x='0' y='$yEmu'/><a:ext cx='9144000' cy='$heightEmu'/></a:xfrm><a:prstGeom prst='rect'><a:avLst/></a:prstGeom><a:solidFill><a:srgbClr val='$color'/></a:solidFill><a:ln><a:noFill/></a:ln></p:spPr></p:sp>"
}
function New-CroppedFooterXml([int]$id,[string]$name,[double]$x,[double]$y,[double]$height,[int]$cropLeft,[int]$cropRight,[int]$cropTop,[int]$cropBottom) {
    $width = $height * (2048.0 * (100000-$cropLeft-$cropRight)) / (1152.0 * (100000-$cropTop-$cropBottom))
    $xEmu=[long][math]::Round($x*$emuPerPoint); $yEmu=[long][math]::Round($y*$emuPerPoint)
    $wEmu=[long][math]::Round($width*$emuPerPoint); $hEmu=[long][math]::Round($height*$emuPerPoint)
    return "<p:pic xmlns:p='$presentationNs' xmlns:a='$drawingNs' xmlns:r='$relationshipNs'><p:nvPicPr><p:cNvPr id='$id' name='$name'/><p:cNvPicPr><a:picLocks noChangeAspect='1'/></p:cNvPicPr><p:nvPr/></p:nvPicPr><p:blipFill><a:blip r:embed='rId13'><a:biLevel thresh='90000'/><a:duotone><a:srgbClr val='$deepGreen'/><a:srgbClr val='FFFFFF'/></a:duotone></a:blip><a:srcRect l='$cropLeft' r='$cropRight' t='$cropTop' b='$cropBottom'/><a:stretch><a:fillRect/></a:stretch></p:blipFill><p:spPr><a:xfrm><a:off x='$xEmu' y='$yEmu'/><a:ext cx='$wEmu' cy='$hEmu'/></a:xfrm><a:prstGeom prst='rect'><a:avLst/></a:prstGeom><a:ln><a:noFill/></a:ln></p:spPr></p:pic>"
}

if (Test-Path -LiteralPath $draftPath) { throw "Draft already exists: $draftPath" }
$source = [IO.Compression.ZipFile]::OpenRead($SourceDeck)
$output = [IO.Compression.ZipFile]::Open($draftPath,[IO.Compression.ZipArchiveMode]::Create)
try {
    foreach ($entry in $source.Entries) {
        $name = $entry.FullName
        $modified = $false
        $data = $null
        if ($name -match '^ppt/(slides/slide\d+|slideMasters/slideMaster\d+|slideLayouts/slideLayout\d+)\.xml$' -or $name -match '^ppt/notesSlides/notesSlide(5|53)\.xml$') {
            $reader = [IO.StreamReader]::new($entry.Open())
            try { $originalText = $reader.ReadToEnd() } finally { $reader.Dispose() }
            $doc = Get-XmlDocument $originalText
            if($name -match '^ppt/notesSlides/notesSlide(5|53)\.xml$'){
                $notesBody=$doc.SelectSingleNode('//*[local-name()="sp"][.//*[local-name()="ph" and @type="body"]]/*[local-name()="txBody"]')
                if(!$notesBody){throw "Notes body missing in $name"}
                Add-Fragment $doc $notesBody "<a:p xmlns:a='$drawingNs'><a:r><a:t>$huntFull</a:t></a:r></a:p>"
                $modified=$true
                $changes.Add(@{part=$name;operation='Hunt_meeting_lecture_source_added'})
            }
            foreach($attr in @($doc.SelectNodes('//@*[namespace-uri()="http://schemas.openxmlformats.org/officeDocument/2006/relationships" and string(.)=""]'))){
                # PowerPoint writes an empty relationship on local media actions.
                # The action remains; an empty relationship ID has no target.
                [void]$attr.OwnerElement.RemoveAttributeNode($attr)
                $modified=$true
                $changes.Add(@{part=$name;operation='remove_empty_relationship_attribute'})
            }
            $paletteChanges = 0
            foreach ($color in $doc.SelectNodes('//*[local-name()="spTree"]//*[local-name()="srgbClr"]')) {
                # Chart XML and all raster/vector evidence assets remain byte-identical.
                if ($color.GetAttribute('val').ToUpperInvariant() -in @('B57C50','006D77','158158')) {
                    $color.SetAttribute('val',$deepGreen); $paletteChanges++
                } elseif ($color.GetAttribute('val').ToUpperInvariant() -eq '92CF48') {
                    $color.SetAttribute('val',$leafGreen); $paletteChanges++
                }
            }
            if ($paletteChanges) { $modified=$true; $changes.Add(@{part=$name;operation='theme_colors';count=$paletteChanges}) }
            if ($name -eq 'ppt/slideMasters/slideMaster2.xml') {
                $bg = $doc.SelectSingleNode('//*[local-name()="cSld"]/*[local-name()="bg"]')
                $bg.InnerXml = "<p:bgPr xmlns:p='$presentationNs' xmlns:a='$drawingNs'><a:solidFill><a:srgbClr val='FFFFFF'/></a:solidFill><a:effectLst/></p:bgPr>"
                $tree = $doc.SelectSingleNode('//*[local-name()="spTree"]')
                foreach ($shape in @($tree.SelectNodes('./*[local-name()="sp"]'))) {
                    if ($shape.SelectSingleNode('.//*[local-name()="ph" and @type="sldNum"]')) { [void]$tree.RemoveChild($shape) }
                }
                $logo = $tree.SelectSingleNode('./*[local-name()="pic"]')
                $extent = $logo.SelectSingleNode('.//*[local-name()="ext"]')
                $logoHeight = 38.0 * [double]$extent.GetAttribute('cy') / [double]$extent.GetAttribute('cx')
                Set-Transform $logo 675 6 38 $logoHeight
                Add-Fragment $doc $tree (New-BandXml 900 'Existing footer green band' 389 16 $deepGreen)
                Add-Fragment $doc $tree (New-BandXml 901 'Rose green footer accent' 388 1 $leafGreen)
                Add-Fragment $doc $tree (New-CroppedFooterXml 902 'Original C3NS NeuroNex footer wordmark, smaller' 6 390 14 600 74000 92200 300)
                Add-Fragment $doc $tree (New-CroppedFooterXml 903 'Original footer neural icons, smaller' 644 390 14 82700 4100 92000 0)
                $modified=$true
                $changes.Add(@{part=$name;operation='existing_branding_scaled';logo_width_pt=38;footer_height_pt=16;original_footer_height_pt=32.7})
            }
            if ($name -eq 'ppt/slideLayouts/slideLayout12.xml') {
                $doc.DocumentElement.SetAttribute('showMasterSp','1')
                $bg=$doc.SelectSingleNode('//*[local-name()="cSld"]/*[local-name()="bg"]')
                $bg.InnerXml="<p:bgPr xmlns:p='$presentationNs' xmlns:a='$drawingNs'><a:solidFill><a:srgbClr val='FFFFFF'/></a:solidFill><a:effectLst/></p:bgPr>"
                $modified=$true
                $changes.Add(@{part=$name;operation='cover_uses_reduced_existing_branding'})
            }
            if ($name -match '^ppt/slides/slide(\d+)\.xml$') {
                $slideNumber=[int]$Matches[1]
                if($slideNumber -in @(5,53)){
                    $creditUpdated=$false
                    foreach($textNode in $doc.SelectNodes('//*[local-name()="t"]')){
                        if($slideNumber -eq 5 -and $textNode.InnerText.Contains('Hunt, INRC (2019)')){
                            $textNode.InnerText=$textNode.InnerText.Replace('Hunt, INRC (2019)',$huntMinimal)
                            $creditUpdated=$true
                        }elseif($slideNumber -eq 53 -and $textNode.InnerText -eq 'Table 6.1 · contribution 6'){
                            $textNode.InnerText=$huntMinimal+' · Dissertation Table 6.1'
                            $creditUpdated=$true
                        }
                    }
                    if(!$creditUpdated){throw "Expected Hunt credit location missing on slide $slideNumber"}
                    $modified=$true
                    $changes.Add(@{part=$name;operation='Hunt_meeting_lecture_credit_corrected'})
                }
                if ($slideNumber -eq 1) {
                    $doc.DocumentElement.SetAttribute('showMasterSp','1')
                    $tree=$doc.SelectSingleNode('//*[local-name()="spTree"]')
                    $title=$tree.SelectSingleNode('./*[local-name()="sp"][1]')
                    Set-Transform $title 45 94 620 160
                    Set-TextColor $title $deepGreen 3600
                    foreach($para in $title.SelectNodes('.//*[local-name()="pPr"]')){$para.SetAttribute('algn','l')}
                    $bodyPr=$title.SelectSingleNode('.//*[local-name()="bodyPr"]');$bodyPr.SetAttribute('anchor','t')
                    $psu=$tree.SelectSingleNode('./*[local-name()="pic"][1]')
                    Set-Transform $psu 45 326 160 (160*1122068.0/3681170.0)
                    $subtitle=$tree.SelectSingleNode('./*[local-name()="sp"][2]')
                    Set-Transform $subtitle 45 272 620 42
                    Set-TextColor $subtitle '202C33' 1600
                    foreach($para in $subtitle.SelectNodes('.//*[local-name()="pPr"]')){$para.SetAttribute('algn','l')}
                    $modified=$true
                    $changes.Add(@{part=$name;operation='cover_repositioned_existing_text_and_PSU_logo'})
                } else {
                    foreach($shape in $doc.SelectNodes('//*[local-name()="spTree"]/*[local-name()="sp"]')) {
                        $text=($shape.SelectNodes('.//*[local-name()="t"]')|ForEach-Object{$_.InnerText}) -join ''
                        $off=$shape.SelectSingleNode('./*[local-name()="spPr"]/*[local-name()="xfrm"]/*[local-name()="off"]')
                        if($text -match '^\d+$' -and $off -and [long]$off.GetAttribute('y') -gt 4500000) {
                            Set-Transform $shape 697 392 16 10
                            Set-TextColor $shape 'FFFFFF' 800
                            $modified=$true
                        }
                    }
                }
                if($slideNumber -in @(36,66)){
                    # Trim only stray, partially captured panel captions. All circuit
                    # elements and quantitative marks keep their original scale and position.
                    $pics=$doc.SelectNodes('//*[local-name()="spTree"]/*[local-name()="pic"]')
                    $pic=if($slideNumber -eq 36){$pics.Item(0)}else{$pics.Item($pics.Count-1)}
                    $fill=$pic.SelectSingleNode('./*[local-name()="blipFill"]')
                    $crop=$doc.CreateElement('a','srcRect',$drawingNs)
                    $crop.SetAttribute('l','0');$crop.SetAttribute('r','0')
                    $trimTop=if($slideNumber -eq 36){6500}else{0}
                    $trimBottom=if($slideNumber -eq 66){5000}else{0}
                    $crop.SetAttribute('t',[string]$trimTop);$crop.SetAttribute('b',[string]$trimBottom)
                    [void]$fill.InsertAfter($crop,$fill.SelectSingleNode('./*[local-name()="blip"]'))
                    $off=$pic.SelectSingleNode('./*[local-name()="spPr"]/*[local-name()="xfrm"]/*[local-name()="off"]')
                    $ext=$pic.SelectSingleNode('./*[local-name()="spPr"]/*[local-name()="xfrm"]/*[local-name()="ext"]')
                    $height=[long]$ext.GetAttribute('cy')
                    $off.SetAttribute('y',[string]([long]$off.GetAttribute('y')+[long][math]::Round($height*$trimTop/100000.0)))
                    $ext.SetAttribute('cy',[string][long][math]::Round($height*(100000-$trimTop-$trimBottom)/100000.0))
                    $modified=$true
                    $changes.Add(@{part=$name;operation='trim_incomplete_raster_panel_caption';crop_top=$trimTop;crop_bottom=$trimBottom;quantitative_marks_unchanged=$true})
                }
            }
            if ($modified) { $data=[Text.UTF8Encoding]::new($false).GetBytes($doc.OuterXml) }
        }
        $target=$output.CreateEntry($name,[IO.Compression.CompressionLevel]::Optimal)
        $target.LastWriteTime=$entry.LastWriteTime
        $targetStream=$target.Open()
        try {
            if ($null -ne $data) { $targetStream.Write($data,0,$data.Length) }
            else { $sourceStream=$entry.Open();try{$sourceStream.CopyTo($targetStream)}finally{$sourceStream.Dispose()} }
        } finally { $targetStream.Dispose() }
    }
} finally { $output.Dispose();$source.Dispose() }
$changes | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath (Join-Path $workRoot 'build\master_changes.json') -Encoding utf8
[pscustomobject]@{draft=$draftPath;source_sha256=(Get-FileHash -LiteralPath $SourceDeck -Algorithm SHA256).Hash;changes=$changes.Count} | ConvertTo-Json
