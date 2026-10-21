Add-Type -AssemblyName System.Drawing
foreach ($n in @("crop_pipe_loop","crop_wf_no","crop_pipe_band","crop_pipe_left","crop_wf_route")) {
  $src = "D:\Github\Bipedal_Robot\tmp\$n.png"
  $img = [System.Drawing.Image]::FromFile($src)
  $dst = "D:\Github\Bipedal_Robot\tmp\$n.jpg"
  $jpegCodec = [System.Drawing.Imaging.ImageCodecInfo]::GetImageEncoders() | Where-Object { $_.MimeType -eq "image/jpeg" }
  $ep = New-Object System.Drawing.Imaging.EncoderParameters(1)
  $ep.Param[0] = New-Object System.Drawing.Imaging.EncoderParameter([System.Drawing.Imaging.Encoder]::Quality, [long]95)
  $img.Save($dst, $jpegCodec, $ep)
  $img.Dispose()
  Write-Output "saved $dst"
}
