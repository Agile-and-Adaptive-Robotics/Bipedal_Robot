param([switch]$Apply)
$ErrorActionPreference='Stop'
$defenseRoot=(Resolve-Path -LiteralPath (Split-Path $PSScriptRoot -Parent)).Path
$retained=@(
    (Join-Path $defenseRoot 'C3NS_Revision_20261030_v14\output\Bolen_PhD_Defense_RoseGreen_20261030_v14.pptx'),
    (Join-Path $defenseRoot 'C3NS_Revision_20261030_v15\output\Bolen_PhD_Defense_RoseGreen_20261030_v15.pptx')
)
foreach($keptPath in $retained){if(!(Test-Path -LiteralPath $keptPath)){throw "Retained revision is missing: $keptPath"}}
$generatedFolders=@(
    'AE_Revision_20261023','AE_Rebuild_20261028_v9','C3NS_Rebuild_20261028_v10',
    'C3NS_Visual_20261029_v11','C3NS_Revision_20261029_v12',
    'C3NS_Revision_20261030_v13','C3NS_Revision_20261030_v14','C3NS_Revision_20261030_v15'
)
$targets=[System.Collections.Generic.List[object]]::new()
$outdatedHashes=[System.Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
foreach($folderName in $generatedFolders){
    $folder=Join-Path $defenseRoot $folderName
    if(!(Test-Path -LiteralPath $folder)){continue}
    foreach($file in Get-ChildItem -LiteralPath $folder -Recurse -Force -File -Filter '*.pptx'){
        if($file.Name.StartsWith('~$') -or $file.FullName -in $retained -or $file.FullName -match '\\reference\\'){continue}
        $resolved=(Resolve-Path -LiteralPath $file.FullName).Path
        if(!$resolved.StartsWith($defenseRoot+'\',[StringComparison]::OrdinalIgnoreCase)){throw "File outside defense directory: $resolved"}
        $hash=(Get-FileHash -LiteralPath $resolved -Algorithm SHA256).Hash
        [void]$outdatedHashes.Add($hash)
        $targets.Add([pscustomobject]@{path=$resolved;bytes=$file.Length;sha256=$hash})
    }
}
# Remove only proven copies of superseded defense decks from shared locations.
foreach($duplicatePath in @(
    (Join-Path $defenseRoot 'Bolen_Dissertation_Defense.pptx'),
    (Join-Path $defenseRoot 'slideshow\AH001_20261028\.build\combined_candidate.pptx')
)){
    if(Test-Path -LiteralPath $duplicatePath){
        $hash=(Get-FileHash -LiteralPath $duplicatePath -Algorithm SHA256).Hash
        if($outdatedHashes.Contains($hash)){
            $file=Get-Item -LiteralPath $duplicatePath
            $targets.Add([pscustomobject]@{path=$file.FullName;bytes=$file.Length;sha256=$hash})
        }
    }
}
$report=[ordered]@{applied=$false;retained=$retained;file_count=$targets.Count;bytes=($targets|Measure-Object -Property bytes -Sum).Sum;files=@($targets)}
$report|ConvertTo-Json -Depth 6|Set-Content -LiteralPath (Join-Path $PSScriptRoot 'build\cleanup_manifest.json') -Encoding utf8
if($Apply){
    foreach($item in $targets){
        $resolved=(Resolve-Path -LiteralPath $item.path).Path
        if(!$resolved.StartsWith($defenseRoot+'\',[StringComparison]::OrdinalIgnoreCase) -or $resolved -in $retained){throw 'Unsafe cleanup target'}
        if((Get-FileHash -LiteralPath $resolved -Algorithm SHA256).Hash -ne $item.sha256){throw "File changed since inventory: $resolved"}
        Remove-Item -LiteralPath $resolved
    }
    $report.applied=$true
    $report|ConvertTo-Json -Depth 6|Set-Content -LiteralPath (Join-Path $PSScriptRoot 'build\cleanup_manifest.json') -Encoding utf8
}
[pscustomobject]@{applied=[bool]$Apply;file_count=$targets.Count;GiB=[math]::Round($report.bytes/1GB,3);retained=$retained}|ConvertTo-Json -Depth 3

