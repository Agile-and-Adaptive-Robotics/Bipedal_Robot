$ErrorActionPreference='Stop'
Add-Type -AssemblyName System.IO.Compression.FileSystem
$baselinePath=Join-Path (Split-Path $PSScriptRoot -Parent) 'C3NS_Revision_20261030_v13\output\Bolen_PhD_Defense_C3NS_20261030_v13_review.pptx'
$finalPath=Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261030_v14.pptx'
function EntryHash($entry){
    $stream=$entry.Open()
    try{return [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($stream))}finally{$stream.Dispose()}
}
function EntryXml($entry){
    $reader=[IO.StreamReader]::new($entry.Open())
    try{[xml]$doc=$reader.ReadToEnd();return ,$doc}finally{$reader.Dispose()}
}
$before=[IO.Compression.ZipFile]::OpenRead($baselinePath)
$after=[IO.Compression.ZipFile]::OpenRead($finalPath)
$checked=@{media=0;charts=0;workbooks=0;notes=0;themes=0;slide_text=0;slide_visibility=0;plot_placement=0}
$huntMinimal='A. J. Hunt, INRC2 Meeting, Intel Labs (Feb 15, 2019).'
$huntFull='Source: Hunt, A. J. (2019, February 15). Real-time control of dynamic walking in a dog-like robot. Presented at the 2nd Intel Neuromorphic Research Community (INRC2) Meeting.'
$failures=[System.Collections.Generic.List[object]]::new()
try{
    foreach($entry in $before.Entries){
        $category=switch -Regex ($entry.FullName){
            '^ppt/media/'{'media';break}
            '^ppt/charts/'{'charts';break}
            '^ppt/embeddings/'{'workbooks';break}
            '^ppt/notesSlides/'{'notes';break}
            '^ppt/theme/'{'themes';break}
            default{$null}
        }
        if($category){
            $other=$after.GetEntry($entry.FullName)
            if($entry.FullName -match '^ppt/notesSlides/notesSlide(5|53)\.xml$'){
                $oldNotes=EntryXml $entry;$newNotes=EntryXml $other
                $oldText=($oldNotes.SelectNodes('//*[local-name()="t"]')|ForEach-Object{$_.InnerText}) -join ''
                $newText=($newNotes.SelectNodes('//*[local-name()="t"]')|ForEach-Object{$_.InnerText}) -join ''
                if($newText.Replace($huntFull,'') -cne $oldText -or !$newText.Contains($huntFull)){$failures.Add(@{part=$entry.FullName;check='bounded_Hunt_source_addition'})}
            }elseif(!$other -or (EntryHash $entry) -ne (EntryHash $other)){$failures.Add(@{part=$entry.FullName;check='byte_identity'})}
            $checked[$category]++
        }
        if($entry.FullName -match '^ppt/slides/slide(\d+)\.xml$'){
            $a=EntryXml $entry;$b=EntryXml ($after.GetEntry($entry.FullName))
            $ta=@($a.SelectNodes('//*[local-name()="t"]')|ForEach-Object{$_.InnerText})
            $tb=@($b.SelectNodes('//*[local-name()="t"]')|ForEach-Object{$_.InnerText})
            $beforeText=$ta -join [char]31;$afterText=$tb -join [char]31
            if($entry.FullName -eq 'ppt/slides/slide5.xml'){$beforeText=$beforeText.Replace('Hunt, INRC (2019)',$huntMinimal)}
            if($entry.FullName -eq 'ppt/slides/slide53.xml'){$beforeText=$beforeText.Replace('Table 6.1 · contribution 6',$huntMinimal+' · Dissertation Table 6.1')}
            if($beforeText -cne $afterText){$failures.Add(@{part=$entry.FullName;check='slide_text'})}
            $checked.slide_text++
            if($a.DocumentElement.GetAttribute('show') -ne $b.DocumentElement.GetAttribute('show')){$failures.Add(@{part=$entry.FullName;check='visibility'})}
            $checked.slide_visibility++
            $ga=@($a.SelectNodes('//*[local-name()="graphicFrame"]/*[local-name()="xfrm"]')|ForEach-Object{$_.OuterXml})
            $gb=@($b.SelectNodes('//*[local-name()="graphicFrame"]/*[local-name()="xfrm"]')|ForEach-Object{$_.OuterXml})
            if(($ga -join '') -cne ($gb -join '')){$failures.Add(@{part=$entry.FullName;check='native_plot_placement'})}
            $checked.plot_placement++
        }
    }
}finally{$before.Dispose();$after.Dispose()}
$result=@{passed=$failures.Count -eq 0;checked=$checked;authorized_text_changes=@('Hunt meeting-lecture credit and full notes source on slides 5 and 53');failures=@($failures);baseline_sha256=(Get-FileHash -LiteralPath $baselinePath -Algorithm SHA256).Hash;final_sha256=(Get-FileHash -LiteralPath $finalPath -Algorithm SHA256).Hash}
$result|ConvertTo-Json -Depth 8|Set-Content -LiteralPath (Join-Path $PSScriptRoot 'build\preservation.json') -Encoding utf8
$result|ConvertTo-Json -Depth 5
if($failures.Count){throw 'Preservation check failed'}
