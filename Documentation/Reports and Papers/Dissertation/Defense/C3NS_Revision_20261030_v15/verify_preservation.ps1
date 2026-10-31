$ErrorActionPreference='Stop'
Add-Type -AssemblyName System.IO.Compression.FileSystem
$parent=Split-Path $PSScriptRoot
$baseline=Join-Path $parent 'C3NS_Revision_20261030_v14\output\Bolen_PhD_Defense_RoseGreen_20261030_v14.pptx'
$final=Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261030_v15.pptx'
$allowed=@('[Content_Types].xml','ppt/slides/slide1.xml','ppt/slides/_rels/slide1.xml.rels','ppt/slides/slide19.xml','ppt/slides/_rels/slide19.xml.rels','ppt/notesSlides/notesSlide19.xml')
function Hash-Entry($entry){$s=$entry.Open();try{[Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($s))}finally{$s.Dispose()}}
function Xml-Entry($entry){$r=[IO.StreamReader]::new($entry.Open());try{[xml]$x=$r.ReadToEnd();return ,$x}finally{$r.Dispose()}}
$a=[IO.Compression.ZipFile]::OpenRead($baseline)
$b=[IO.Compression.ZipFile]::OpenRead($final)
$failures=[System.Collections.Generic.List[string]]::new()
$changed=[System.Collections.Generic.List[string]]::new()
$added=[System.Collections.Generic.List[string]]::new()
$identity=0
try {
    foreach($entry in $a.Entries){$other=$b.GetEntry($entry.FullName);if(!$other){$failures.Add('Missing: '+$entry.FullName);continue};if((Hash-Entry $entry) -cne (Hash-Entry $other)){$changed.Add($entry.FullName);if($entry.FullName -notin $allowed){$failures.Add('Unexpected change: '+$entry.FullName)}}else{$identity++}}
    foreach($entry in $b.Entries){if(!$a.GetEntry($entry.FullName)){$added.Add($entry.FullName);if($entry.FullName -notlike 'ppt/media/w2l20210217_*'){$failures.Add('Unexpected addition: '+$entry.FullName)}}}
    foreach($number in 1..83){$old=Xml-Entry ($a.GetEntry("ppt/slides/slide$number.xml"));$new=Xml-Entry ($b.GetEntry("ppt/slides/slide$number.xml"));if($old.DocumentElement.GetAttribute('show') -cne $new.DocumentElement.GetAttribute('show')){$failures.Add('Visibility changed: '+$number)}}
    $cover=Xml-Entry ($b.GetEntry('ppt/slides/slide1.xml'))
    if(!$cover.SelectSingleNode('//*[local-name()="pic"][.//*[local-name()="cNvPr" and @name="Original C3NS honeycomb artwork"]]')){$failures.Add('Original honeycomb missing')}
    $slide=Xml-Entry ($b.GetEntry('ppt/slides/slide19.xml'))
    if(!$slide.SelectSingleNode('//*[local-name()="videoFile"]') -or !$slide.SelectSingleNode('//*[local-name()="timing"]')){$failures.Add('Native video or timing missing')}
    $oldNotes=Xml-Entry ($a.GetEntry('ppt/notesSlides/notesSlide19.xml'));$newNotes=Xml-Entry ($b.GetEntry('ppt/notesSlides/notesSlide19.xml'))
    $ta=($oldNotes.SelectNodes('//*[local-name()="t"]')|ForEach-Object{$_.InnerText}) -join ''
    $tb=($newNotes.SelectNodes('//*[local-name()="t"]')|ForEach-Object{$_.InnerText}) -join ''
    if(!$tb.StartsWith($ta,[StringComparison]::Ordinal)){$failures.Add('Original speaker notes changed')}
    $counts=[ordered]@{original_media=@($a.Entries|Where-Object{$_.FullName -like 'ppt/media/*'}).Count;charts=@($a.Entries|Where-Object{$_.FullName -like 'ppt/charts/*'}).Count;workbooks=@($a.Entries|Where-Object{$_.FullName -like 'ppt/embeddings/*'}).Count;themes=@($a.Entries|Where-Object{$_.FullName -like 'ppt/theme/*'}).Count;masters=@($a.Entries|Where-Object{$_.FullName -like 'ppt/slideMasters/*'}).Count}
    $result=[ordered]@{passed=$failures.Count -eq 0;unchanged_package_parts=$identity;unchanged_original_assets=$counts;changed_parts=@($changed);added_parts=@($added);slides=83;visibility_preserved=83;chart_owners=@(8,14,58,59);failures=@($failures);sha256=(Get-FileHash -LiteralPath $final).Hash}
    $result|ConvertTo-Json -Depth 7|Set-Content -LiteralPath (Join-Path $PSScriptRoot 'build\preservation.json') -Encoding utf8
    $result|ConvertTo-Json -Depth 5
    if($failures.Count){throw 'Preservation failed'}
} finally {$a.Dispose();$b.Dispose()}
