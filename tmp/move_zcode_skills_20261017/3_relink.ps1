# Rebuild .agents\skills junctions from old D: ZCode_Skills path to the new C: location (2026-10-17, easteregg2)
$root = 'C:\Users\Ben Bolen\.agents\skills'
$oldPrefix = 'D:\GitHub\ZCode_Skills'
$newRoot = 'C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills'

Get-ChildItem -LiteralPath $root -Force | ForEach-Object {
  $tstr = ($_.Target -join ',')
  if ($_.LinkType -eq 'Junction' -and $tstr -and $tstr.StartsWith($oldPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
    $full = $_.FullName
    $rel = $tstr.Substring($oldPrefix.Length).TrimStart('\')
    $newT = if ($rel) { Join-Path $newRoot $rel } else { $newRoot }
    & cmd.exe /c rmdir "$full"
    if (Test-Path -LiteralPath $full) { Write-Output "FAILED-REMOVE $full"; return }
    & cmd.exe /c mklink /J "$full" "$newT" | Out-Null
    if (Test-Path -LiteralPath (Join-Path $full 'SKILL.md')) { Write-Output "RELINKED $($_.Name) -> $newT" }
    else { Write-Output "JUNCTION-CREATED-NO-SKILLMD $($_.Name) -> $newT" }
  } else {
    Write-Output "SKIP $($_.Name) (LinkType=$($_.LinkType) Target=$tstr)"
  }
}
