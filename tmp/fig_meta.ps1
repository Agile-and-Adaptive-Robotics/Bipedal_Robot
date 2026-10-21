$files = @(
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\10-introduction\components\look_roadmap.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns\components\look_pipeline.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns\components\look_curriculum.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\70-overall\components\look_workflow.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\70-overall\components\look_sensorium.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\70-overall\components\look_cutaway.png"
)
foreach ($f in $files) {
  $bytes = [System.IO.File]::ReadAllBytes($f)
  # PNG: 8-byte signature, then chunks. Find pHYs (pHYs) and tEXt chunks.
  $pos = 8
  $found = @()
  while ($pos -lt $bytes.Length - 4) {
    $len = ($bytes[$pos] -shl 24) -bor ($bytes[$pos+1] -shl 16) -bor ($bytes[$pos+2] -shl 8) -bor $bytes[$pos+3]
    $type = [System.Text.Encoding]::ASCII.GetString($bytes[($pos+4)..($pos+7)])
    if ($type -eq "pHYs") {
      $ppuX = ($bytes[$pos+8] -shl 24) -bor ($bytes[$pos+9] -shl 16) -bor ($bytes[$pos+10] -shl 8) -bor $bytes[$pos+11]
      $ppuY = ($bytes[$pos+12] -shl 24) -bor ($bytes[$pos+13] -shl 16) -bor ($bytes[$pos+14] -shl 8) -bor $bytes[$pos+15]
      $unit = $bytes[$pos+16]
      $dpi = [math]::Round($ppuX * 0.0254, 1)
      $found += "pHYs: ${ppuX}x${ppuY} ppunit unit=$unit -> ${dpi} dpi"
    }
    if ($type -eq "tEXt") {
      $chunk = [System.Text.Encoding]::ASCII.GetString($bytes[($pos+8)..($pos+7+$len)])
      if ($chunk.Length -gt 80) { $chunk = $chunk.Substring(0,80) }
      $found += "tEXt: $chunk"
    }
    if ($type -eq "IEND") { break }
    $pos = $pos + 8 + $len + 4
    if ($len -lt 0) { break }
  }
  $name = Split-Path (Split-Path $f -Parent) -Leaf
  $leaf = Split-Path $f -Leaf
  if ($found.Count -eq 0) { $found = @("(no pHYs/tEXt chunks)") }
  Write-Output ("{0}/{1}: {2}" -f $name, $leaf, ($found -join " | "))
}
