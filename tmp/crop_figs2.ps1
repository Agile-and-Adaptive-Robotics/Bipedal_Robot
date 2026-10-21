Add-Type -AssemblyName System.Drawing

function CropAndSave($src, $x, $y, $w, $h, $dst, [int]$z) {
  $img = [System.Drawing.Image]::FromFile($src)
  $srcRect = New-Object System.Drawing.Rectangle([int]$x, [int]$y, [int]$w, [int]$h)
  $dw = [int]($w * $z); $dh = [int]($h * $z)
  $bmp = New-Object System.Drawing.Bitmap($dw, $dh)
  $g = [System.Drawing.Graphics]::FromImage($bmp)
  $g.Clear([System.Drawing.Color]::White)
  $g.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
  $g.PixelOffsetMode = [System.Drawing.Drawing2D.PixelOffsetMode]::Half
  $dstRect = New-Object System.Drawing.Rectangle(0, 0, $dw, $dh)
  $g.DrawImage($img, $dstRect, $srcRect, [System.Drawing.GraphicsUnit]::Pixel)
  $g.Dispose()
  $bmp.Save($dst, [System.Drawing.Imaging.ImageFormat]::Png)
  $bmp.Dispose()
  $img.Dispose()
  Write-Output "saved $dst ($dw x $dh)"
}

$pipe = "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns\components\look_pipeline.png"
$wf   = "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\70-overall\components\look_workflow.png"

CropAndSave $pipe 430 150 480 160 "D:\Github\Bipedal_Robot\tmp\crop_pipe_loop.png" 3
CropAndSave $pipe 700 150 311 160 "D:\Github\Bipedal_Robot\tmp\crop_pipe_band.png" 3
CropAndSave $pipe 230 150 260 160 "D:\Github\Bipedal_Robot\tmp\crop_pipe_left.png" 3
CropAndSave $wf 520 230 430 227 "D:\Github\Bipedal_Robot\tmp\crop_wf_no.png" 3
CropAndSave $wf 300 230 350 200 "D:\Github\Bipedal_Robot\tmp\crop_wf_route.png" 3
