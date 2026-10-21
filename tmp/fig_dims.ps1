$files = @(
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\10-introduction\components\look_roadmap.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns\components\look_pipeline.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns\components\look_curriculum.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\70-overall\components\look_workflow.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\70-overall\components\look_sensorium.png",
  "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\70-overall\components\look_cutaway.png"
)
Add-Type -AssemblyName System.Drawing
foreach ($f in $files) {
  $img = [System.Drawing.Image]::FromFile($f)
  $w = $img.Width; $h = $img.Height
  $img.Dispose()
  $dpiAtPrint = [math]::Round($w / (16.5 / 2.54), 0)
  $name = Split-Path $f -Leaf
  $dir = Split-Path (Split-Path $f -Parent) -Leaf
  Write-Output ("{0}/{1}: {2} x {3} px  -> {4} dpi at 16.5 cm print width; aspect {5}" -f $dir, $name, $w, $h, $dpiAtPrint, [math]::Round($w/$h,2))
}
