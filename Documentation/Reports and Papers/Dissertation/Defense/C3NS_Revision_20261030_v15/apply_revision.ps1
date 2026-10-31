$ErrorActionPreference='Stop'
Add-Type -AssemblyName System.IO.Compression.FileSystem
$parent=Split-Path $PSScriptRoot
$source=Join-Path $parent 'C3NS_Revision_20261030_v14\output\Bolen_PhD_Defense_RoseGreen_20261030_v14.pptx'
$candidate=Join-Path $PSScriptRoot 'build\candidate.pptx'
if(Test-Path -LiteralPath $candidate){Remove-Item -LiteralPath $candidate}
Copy-Item -LiteralPath $source -Destination $candidate
$z=[IO.Compression.ZipFile]::Open($candidate,[IO.Compression.ZipArchiveMode]::Update)
$component=[IO.Compression.ZipFile]::OpenRead((Join-Path $PSScriptRoot 'build\video_component.pptx'))
$p='http://schemas.openxmlformats.org/presentationml/2006/main'
$a='http://schemas.openxmlformats.org/drawingml/2006/main'
$r='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
function Read-Xml($zip,$name){$reader=[IO.StreamReader]::new($zip.GetEntry($name).Open());try{$doc=[xml]::new();$doc.PreserveWhitespace=$true;$doc.LoadXml($reader.ReadToEnd());return ,$doc}finally{$reader.Dispose()}}
function Write-Xml($name,$doc){$entry=$z.GetEntry($name);if($entry){$entry.Delete()};$entry=$z.CreateEntry($name,[IO.Compression.CompressionLevel]::Optimal);$writer=[IO.StreamWriter]::new($entry.Open(),[Text.UTF8Encoding]::new($false));try{$writer.Write($doc.OuterXml)}finally{$writer.Dispose()}}
function Fragment($doc,$xml){$f=$doc.CreateDocumentFragment();$f.InnerXml=$xml;return $f}
function Set-Box($shape,$x,$y,$w,$h){$t=$shape.SelectSingleNode('.//*[local-name()="xfrm"]');$off=$t.SelectSingleNode('*[local-name()="off"]');$ext=$t.SelectSingleNode('*[local-name()="ext"]');$off.x=[string][long]($x*12700);$off.y=[string][long]($y*12700);$ext.cx=[string][long]($w*12700);$ext.cy=[string][long]($h*12700)}
function White-Text($shape){foreach($pr in $shape.SelectNodes('.//*[local-name()="rPr" or local-name()="defRPr" or local-name()="endParaRPr"]')){$fill=$pr.SelectSingleNode('*[local-name()="solidFill"]');if(!$fill){$fill=$pr.OwnerDocument.CreateElement('a','solidFill',$a);$before=$pr.SelectSingleNode('*[local-name()="latin" or local-name()="ea" or local-name()="cs" or local-name()="sym" or local-name()="hlinkClick" or local-name()="extLst"]');if($before){[void]$pr.InsertBefore($fill,$before)}else{[void]$pr.AppendChild($fill)}};$fill.RemoveAll();$color=$pr.OwnerDocument.CreateElement('a','srgbClr',$a);$color.SetAttribute('val','FFFFFF');[void]$fill.AppendChild($color)}}
function Add-Rel($doc,$id,$type,$target){$e=$doc.CreateElement('Relationship',$doc.DocumentElement.NamespaceURI);$e.SetAttribute('Id',$id);$e.SetAttribute('Type',$type);$e.SetAttribute('Target',$target);[void]$doc.DocumentElement.AppendChild($e)}
try {
    $cover=Read-Xml $z 'ppt/slides/slide1.xml'
    $cover.DocumentElement.SetAttribute('showMasterSp','0')
    $tree=$cover.SelectSingleNode('//*[local-name()="spTree"]')
    $bg=@"
<p:pic xmlns:p="$p" xmlns:a="$a" xmlns:r="$r"><p:nvPicPr><p:cNvPr id="800" name="Original C3NS honeycomb artwork" descr="Original C3NS honeycomb pattern recolored in Rose green."/><p:cNvPicPr/><p:nvPr/></p:nvPicPr><p:blipFill><a:blip r:embed="rIdCoverArt"><a:biLevel thresh="90000"/><a:duotone><a:srgbClr val="174D24"/><a:srgbClr val="58D51B"/></a:duotone></a:blip><a:stretch><a:fillRect/></a:stretch></p:blipFill><p:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="9144000" cy="5143500"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>
"@
    [void]$tree.InsertBefore((Fragment $cover $bg),$tree.SelectSingleNode('*[local-name()="sp" or local-name()="pic"]'))
    $mask=@"
<p:sp xmlns:p="$p" xmlns:a="$a"><p:nvSpPr><p:cNvPr id="801" name="Cover branding reduction mask"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr><p:spPr><a:xfrm><a:off x="5588000" y="4064000"/><a:ext cx="3556000" cy="889000"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:solidFill><a:srgbClr val="174D24"/></a:solidFill><a:ln><a:noFill/></a:ln></p:spPr></p:sp>
"@
    [void]$tree.AppendChild((Fragment $cover $mask))
    $title=$tree.SelectSingleNode('*[local-name()="sp"][.//*[local-name()="cNvPr" and @id="117"]]')
    Set-Box $title 285 108 390 155
    White-Text $title
    $author=$tree.SelectSingleNode('*[local-name()="sp"][.//*[local-name()="cNvPr" and @id="2"]]')
    Set-Box $author 285 280 390 42
    White-Text $author
    $psu=$tree.SelectSingleNode('*[local-name()="pic"][.//*[local-name()="cNvPr" and @id="118"]]')
    $oldExt=$psu.SelectSingleNode('.//*[local-name()="xfrm"]/*[local-name()="ext"]')
    Set-Box $psu 285 337 140 (140*[double]$oldExt.cy/[double]$oldExt.cx)
    $seal=@"
<p:pic xmlns:p="$p" xmlns:a="$a" xmlns:r="$r"><p:nvPicPr><p:cNvPr id="802" name="Small original white C3NS seal"/><p:cNvPicPr/><p:nvPr/></p:nvPicPr><p:blipFill><a:blip r:embed="rIdCoverArt"><a:biLevel thresh="90000"/><a:duotone><a:srgbClr val="174D24"/><a:srgbClr val="FFFFFF"/></a:duotone></a:blip><a:srcRect l="61500" r="27000" t="79500" b="500"/><a:stretch><a:fillRect/></a:stretch></p:blipFill><p:spPr><a:xfrm><a:off x="8572500" y="76200"/><a:ext cx="482600" cy="472108"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>
"@
    [void]$tree.AppendChild((Fragment $cover $seal))
    $master=Read-Xml $z 'ppt/slideMasters/slideMaster2.xml'
    foreach($id in 900..903){$shape=$master.SelectSingleNode('//*[local-name()="spTree"]/*[.//*[local-name()="cNvPr" and @id="'+$id+'"]]');$copy=$cover.ImportNode($shape,$true);foreach($attr in $copy.SelectNodes('.//@*[local-name()="embed"]')){$attr.Value='rIdCoverFooter'};[void]$tree.AppendChild($copy)}
    $coverRels=Read-Xml $z 'ppt/slides/_rels/slide1.xml.rels'
    Add-Rel $coverRels 'rIdCoverArt' ($r+'/image') '../media/image3.png'
    Add-Rel $coverRels 'rIdCoverFooter' ($r+'/image') '../media/image1.png'
    Write-Xml 'ppt/slides/slide1.xml' $cover
    Write-Xml 'ppt/slides/_rels/slide1.xml.rels' $coverRels

    $slide=Read-Xml $z 'ppt/slides/slide19.xml'
    if($slide.SelectSingleNode('/*/*[local-name()="timing"]')){throw 'Unexpected existing slide19 animation; merge required.'}
    $rels=Read-Xml $z 'ppt/slides/_rels/slide19.xml.rels'
    $compSlide=Read-Xml $component 'ppt/slides/slide1.xml'
    $compRels=Read-Xml $component 'ppt/slides/_rels/slide1.xml.rels'
    $media=$compSlide.SelectSingleNode('//*[local-name()="pic"][.//*[local-name()="videoFile"]]')
    if(!$media){throw 'Native component does not contain a video.'}
    $oldId=$media.SelectSingleNode('.//*[local-name()="cNvPr"]').id
    $idMap=@{}
    foreach($rel in $compRels.DocumentElement.ChildNodes){if($rel.Target -match '^../media/'){$newId='rIdW2L'+$rel.Id;$idMap[$rel.Id]=$newId;$file=[IO.Path]::GetFileName($rel.Target);Add-Rel $rels $newId $rel.Type ('../media/w2l20210217_'+$file)}}
    foreach($entry in $component.Entries|Where-Object{$_.FullName -like 'ppt/media/*' -and $_.Length -gt 0}){$new=$z.CreateEntry('ppt/media/w2l20210217_'+[IO.Path]::GetFileName($entry.FullName),[IO.Compression.CompressionLevel]::Optimal);$from=$entry.Open();$to=$new.Open();try{$from.CopyTo($to)}finally{$to.Dispose();$from.Dispose()}}
    $copy=$slide.ImportNode($media,$true)
    $copy.SelectSingleNode('.//*[local-name()="cNvPr"]').id='4'
    foreach($attr in @($copy.SelectNodes('.//@*'))){if($attr.NamespaceURI -eq $r){if($idMap.ContainsKey($attr.Value)){$attr.Value=$idMap[$attr.Value]}elseif([string]::IsNullOrEmpty($attr.Value)){[void]$attr.OwnerElement.RemoveAttributeNode($attr)}}}
    $old=$slide.SelectSingleNode('//*[local-name()="spTree"]/*[local-name()="pic"][.//*[local-name()="cNvPr" and @id="4"]]')
    [void]$old.ParentNode.ReplaceChild($copy,$old)
    $timing=$compSlide.SelectSingleNode('/*/*[local-name()="timing"]')
    if(!$timing){throw 'Automatic video playback timing missing.'}
    $copyTiming=$slide.ImportNode($timing,$true)
    foreach($attr in $copyTiming.SelectNodes('.//@spid')){if($attr.Value -eq $oldId){$attr.Value='4'}}
    $before=$slide.DocumentElement.SelectSingleNode('*[local-name()="extLst"]')
    if($before){[void]$slide.DocumentElement.InsertBefore($copyTiming,$before)}else{[void]$slide.DocumentElement.AppendChild($copyTiming)}
    foreach($pair in @(@(5,'Original recording · 2021'),@(6,'Archived AnimatLab motion'),@(15,'Released recording pending'))){$shape=$slide.SelectSingleNode('//*[local-name()="spTree"]/*[.//*[local-name()="cNvPr" and @id="'+$pair[0]+'"]]');$text=$shape.SelectSingleNode('.//*[local-name()="t"]');$text.InnerText=$pair[1]}
    Write-Xml 'ppt/slides/slide19.xml' $slide
    Write-Xml 'ppt/slides/_rels/slide19.xml.rels' $rels
    $types=Read-Xml $z '[Content_Types].xml'
    $compTypes=Read-Xml $component '[Content_Types].xml'
    foreach($type in $compTypes.DocumentElement.ChildNodes|Where-Object{$_.LocalName -eq 'Default'}){if(!$types.DocumentElement.SelectSingleNode('*[local-name()="Default" and @Extension="'+$type.Extension+'"]')){[void]$types.DocumentElement.AppendChild($types.ImportNode($type,$true))}}
    Write-Xml '[Content_Types].xml' $types
    $notes=Read-Xml $z 'ppt/notesSlides/notesSlide19.xml'
    $body=$notes.SelectSingleNode('//*[local-name()="sp"][.//*[local-name()="ph" and @type="body"]]/*[local-name()="txBody"]')
    $para=$notes.CreateElement('a','p',$a);$run=$notes.CreateElement('a','r',$a);$text=$notes.CreateElement('a','t',$a)
    $text.InnerText='Media: Ben supplied the original 2021-02-17-20-58-24.mp4 recording. The full 5.4-second sequence is retained, cropped around the model and platform, muted, and looped. This is the archived-motion panel; a native released/falling recording remains pending.'
    [void]$run.AppendChild($text);[void]$para.AppendChild($run);[void]$body.AppendChild($para)
    Write-Xml 'ppt/notesSlides/notesSlide19.xml' $notes
} finally {$component.Dispose();$z.Dispose()}
Copy-Item -LiteralPath (Join-Path $parent 'C3NS_Revision_20261030_v14\output\sadb_app.html') -Destination (Join-Path $PSScriptRoot 'output\sadb_app.html')
Write-Output $candidate
