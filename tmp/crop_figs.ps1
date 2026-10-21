Add-Type -AssemblyName System.Drawing

function CropAndSave($src, $x, $y, $w, $h, $dst, $scale) {
  $img = [System.Drawing.Image]::FromFile($src)
  $rect = New-Object System.Drawing.Rectangle($x, $y, $w, $h)
  $bmp = New-Object System.Drawing.Bitmap($w, $h)
  $g = [System.Drawing.Graphics]::FromImage($bmp)
  $g.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
  $g.PixelOffsetMode = [System.Drawing.Drawing2D.PixelOffsetMode]::Half
  $g.DrawImage($img, (New-Object System.Drawing.Rectangle(0,0,$w*$scale,$h*$scale)), $rect, [System.Drawing.GraphicsUnit]::Pixel)
  $g.Dispose()
  $bmp.Save($dst, [System.Drawing.Imaging.ImageFormat]::Png)
  $bmp.Dispose()
  $img.Dispose()
  Write-Output "saved $dst ($($w*$scale) x $($h*$scale))"
}

$pipe = "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns\components\look_pipeline.png"
$wf   = "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\70-overall\components\look_workflow.png"

# A1: pipeline - afferents label zone + connector vertical (x~568px) crossing
CropAndSave $pipe 430 150 480 160 "D:\Github\Bipedal_Robot\tmp\crop_pipe_loop.png" 3
# A2: pipeline - right band edge overflow zone
CropAndSave $pipe 700 150 311 160 "D:\Github\Bipedal_Robot\tmp\crop_pipe_band.png" 3
# A3: pipeline - per-muscle label / left band edge
CropAndSave $pipe 230 150 260 160 "D:\Github\Bipedal_Robot\tmp\crop_pipe_left.png" 3
# B: workflow - no-branch elbow + verify box top
CropAndSave $wf 520 230 430 227 "D:\Github\Bipedal_Robot\tmp\crop_wf_no.png" 3
# B2: workflow - route box south + horizontal segment arrival
CropAndSave $wf 300 230 350 200 "D:\Github\Bipedal_Robot\tmp\crop_wf_route.png" 3
