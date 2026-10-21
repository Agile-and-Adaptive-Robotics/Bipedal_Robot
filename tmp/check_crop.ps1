Add-Type -AssemblyName System.Drawing
foreach ($f in @("crop_pipe_loop","crop_pipe_band","crop_pipe_left","crop_wf_no","crop_wf_route")) {
  $p = "D:\Github\Bipedal_Robot\tmp\$f.png"
  $img = [System.Drawing.Image]::FromFile($p)
  $bmp = New-Object System.Drawing.Bitmap($img)
  $tot = 0; $n = 0; $dark = 0
  for ($x = 0; $x -lt $bmp.Width; $x += 7) {
    for ($y = 0; $y -lt $bmp.Height; $y += 7) {
      $c = $bmp.GetPixel($x, $y)
      $lum = ($c.R + $c.G + $c.B) / 3
      $tot += $lum; $n++
      if ($lum -lt 128) { $dark++ }
    }
  }
  $mean = [math]::Round($tot / $n, 1)
  $pct = [math]::Round(100.0 * $dark / $n, 1)
  Write-Output "$f : $($bmp.Width)x$($bmp.Height) meanLum=$mean darkPx%=$pct fmt=$($img.PixelFormat)"
  $bmp.Dispose(); $img.Dispose()
}
