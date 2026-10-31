$ErrorActionPreference='Stop'
Add-Type -AssemblyName System.IO.Compression.FileSystem
$parent=Split-Path $PSScriptRoot
New-Item -ItemType Directory -Force -Path (Join-Path $PSScriptRoot 'build'),(Join-Path $PSScriptRoot 'output')|Out-Null
$source=Join-Path $parent 'C3NS_Revision_20261030_v15\output\Bolen_PhD_Defense_RoseGreen_20261030_v15.pptx'
$candidate=Join-Path $PSScriptRoot 'build\candidate.pptx'
if(Test-Path -LiteralPath $candidate){Remove-Item -LiteralPath $candidate}
Copy-Item -LiteralPath $source -Destination $candidate
$order=@(1..22)+@(29..70)+@(72..83)+@(71)+@(23..28)
if($order.Count -ne 83 -or @($order|Sort-Object -Unique).Count -ne 83){throw 'Invalid order'}
$z=[IO.Compression.ZipFile]::Open($candidate,[IO.Compression.ZipArchiveMode]::Update)
function ReadXml($name){$rd=[IO.StreamReader]::new($z.GetEntry($name).Open());try{$d=[xml]::new();$d.PreserveWhitespace=$true;$d.LoadXml($rd.ReadToEnd());return ,$d}finally{$rd.Dispose()}}
function WriteXml($name,$doc){$z.GetEntry($name).Delete();$e=$z.CreateEntry($name);$w=[IO.StreamWriter]::new($e.Open(),[Text.UTF8Encoding]::new($false));try{$w.Write($doc.OuterXml)}finally{$w.Dispose()}}
try {
 $presentation=ReadXml 'ppt/presentation.xml';$lst=$presentation.SelectSingleNode('//*[local-name()="sldIdLst"]')
 $old=@($lst.ChildNodes)
 foreach($n in $old){[void]$lst.RemoveChild($n)}
 $map=@()
 for($i=0;$i -lt $order.Count;$i++){
  $sourceNumber=$order[$i];[void]$lst.AppendChild($old[$sourceNumber-1])
  $slide=ReadXml "ppt/slides/slide$sourceNumber.xml"
  $slide.DocumentElement.SetAttribute('show',$(if($i -ge 64){'0'}else{'1'}))
  foreach($shape in $slide.SelectNodes('//*[local-name()="spTree"]/*[local-name()="sp"]')){
   $off=$shape.SelectSingleNode('.//*[local-name()="xfrm"]/*[local-name()="off"]');$texts=@($shape.SelectNodes('.//*[local-name()="t"]'))
   if($off -and [long]$off.y -gt 4900000 -and $texts.Count -eq 1 -and $texts[0].InnerText -match '^\d+$'){$texts[0].InnerText=[string]($i+1)}
  }
  WriteXml "ppt/slides/slide$sourceNumber.xml" $slide
  $map+=[pscustomobject]@{position=$i+1;source_slide=$sourceNumber;hidden=$i -ge 64;sadb=$sourceNumber -in (@(23..28)+@(71))}
 }
 foreach($sections in @($presentation.SelectNodes('//*[local-name()="sectionLst"]'))){[void]$sections.ParentNode.RemoveChild($sections)}
 WriteXml 'ppt/presentation.xml' $presentation
 $map|ConvertTo-Json|Set-Content -LiteralPath (Join-Path $PSScriptRoot 'build\slide_order.json') -Encoding utf8
}finally{$z.Dispose()}
Copy-Item -LiteralPath (Join-Path $parent 'C3NS_Revision_20261030_v15\output\sadb_app.html') -Destination (Join-Path $PSScriptRoot 'output\sadb_app.html')
