$sim = 'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe'
$repo = 'D:\GitHub\Bipedal_Robot'
$out = "$repo\Code\MuJoCo_SNS\spinal\campaigns\20260930\animatlab"
New-Item -ItemType Directory -Force -Path $out | Out-Null
$models = @(
  "$repo\Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\results\Walker_2_Layer_CPG_BilateralRG_Ground_Standalone.asim",
  "$repo\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG_Standalone_modern.asim",
  "$repo\Neuromechanical_Models\Biped_2xCPG_wSubs\Biped_2xCPG_wSubs_Standalone.asim"
)
foreach ($m in $models) {
  if (-not (Test-Path $m)) { Write-Host "MISSING: $m"; continue }
  $name = [IO.Path]::GetFileNameWithoutExtension($m)
  $wd = "$out\$name"
  New-Item -ItemType Directory -Force -Path $wd | Out-Null
  Copy-Item $m "$wd\" -Force
  $p = Start-Process -FilePath $sim -ArgumentList "`"$wd\$name.asim`"" -WorkingDirectory $wd -PassThru -WindowStyle Hidden
  if ($p.WaitForExit(180000)) { Write-Host "$name exit=$($p.ExitCode)" } else { $p.Kill(); Write-Host "$name TIMEOUT killed" }
}
Get-ChildItem $out -Recurse -Filter *.txt | ForEach-Object { Write-Host ($_.FullName + ' ' + $_.Length + 'B') }
