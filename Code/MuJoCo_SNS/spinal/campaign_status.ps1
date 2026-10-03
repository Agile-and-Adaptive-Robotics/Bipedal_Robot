Set-Location 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
Write-Host ("py procs: " + (Get-Process python -ErrorAction SilentlyContinue | Measure-Object).Count + "  matlab: " + (Get-Process MATLAB -ErrorAction SilentlyContinue | Measure-Object).Count)
Write-Host "--- robust ---"
Get-ChildItem robust_results_*.jsonl -ErrorAction SilentlyContinue | ForEach-Object {
  Write-Host ($_.Name + ": " + (Get-Content $_.FullName | Measure-Object -Line).Lines + " cells")
}
Write-Host "--- retune ---"
foreach ($f in 'pruned_air','pruned_w1','pruned_w2','pruned_w3') {
  if (Test-Path "$f.log") {
    $tail = (Get-Content "$f.log" -Tail 2) -join ' | '
    Write-Host ("$f : $tail")
  } else { Write-Host "$f : no log" }
}
Write-Host "--- validation ---"
if (Test-Path gait_validate.log) { (Get-Content gait_validate.log | Select-String 'val:|wrote') | ForEach-Object { Write-Host $_.Line } }
Write-Host "--- li chain ---"
if (Test-Path w2l_mujoco\li_chain.log) { Get-Content w2l_mujoco\li_chain.log -Tail 10 | ForEach-Object { Write-Host $_ } }
Write-Host "--- simscape ---"
$p = 'D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\logs\easteregg2_runs_20260930\runs.log'
if (Test-Path $p) { Get-Content $p -Tail 10 | ForEach-Object { Write-Host $_ } } else { Write-Host "no log" }
