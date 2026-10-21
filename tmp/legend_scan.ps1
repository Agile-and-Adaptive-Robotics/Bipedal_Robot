Add-Type -AssemblyName System.Drawing
$src = "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\10-introduction\components\look_roadmap.png"
$bmp = New-Object System.Drawing.Bitmap($src)
Write-Output "size: $($bmp.Width) x $($bmp.Height)"

# Legend strip approx rows 2..40. Find dark (<120) pixels per column, split by row:
# line rows = rows where the left swatch (cols 10..55) is dark.
$lineRows = @()
for ($y = 0; $y -lt 45; $y++) {
  $cnt = 0
  for ($x = 10; $x -le 55; $x++) {
    $c = $bmp.GetPixel($x, $y)
    if ((($c.R + $c.G + $c.B) / 3) -lt 120) { $cnt++ }
  }
  if ($cnt -ge 10) { $lineRows += $y }
}
Write-Output ("legend swatch line rows: " + ($lineRows -join ","))

# For columns 50..420, find dark pixels NOT in line rows (text pixels) -> text extent
$textCols = @{}
for ($x = 50; $x -lt 430; $x++) {
  $cnt = 0
  for ($y = 2; $y -lt 42; $y++) {
    if ($lineRows -contains $y) { continue }
    $c = $bmp.GetPixel($x, $y)
    if ((($c.R + $c.G + $c.B) / 3) -lt 120) { $cnt++ }
  }
  if ($cnt -ge 2) { $textCols[$x] = $cnt }
}
$keys = $textCols.Keys | Sort-Object
if ($keys.Count -gt 0) {
  Write-Output ("text columns span: " + $keys[0] + " .. " + $keys[$keys.Count-1])
  # find the gap between entry-1 text and entry-2 text
  $prev = -10; $gaps = @()
  foreach ($k in $keys) {
    if ($k - $prev -gt 8) { $gaps += ("gap before col " + $k + " (prev " + $prev + ")") }
    $prev = $k
  }
  Write-Output ("gaps: " + ($gaps -join " | "))
}

# Dashed swatch 2 location: find dark runs in line rows for cols 260..420
$dashCols = @()
foreach ($y in $lineRows) {
  for ($x = 260; $x -lt 430; $x++) {
    $c = $bmp.GetPixel($x, $y)
    if ((($c.R + $c.G + $c.B) / 3) -lt 120) { $dashCols += $x }
  }
}
if ($dashCols.Count -gt 0) {
  $dk = $dashCols | Sort-Object -Unique
  Write-Output ("dash-2 columns: " + $dk[0] + " .. " + $dk[$dk.Count-1])
}
$bmp.Dispose()
