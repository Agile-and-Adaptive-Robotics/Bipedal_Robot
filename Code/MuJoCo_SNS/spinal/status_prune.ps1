Set-Location 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
Write-Host ("procs: " + (Get-Process python -ErrorAction SilentlyContinue | Measure-Object).Count)
foreach ($f in 's3k','w2lvar','syn6') {
  $n = 0
  if (Test-Path "prune_results_$f.jsonl") { $n = (Get-Content "prune_results_$f.jsonl" | Measure-Object -Line).Lines }
  Write-Host ("{0}: {1} cells" -f $f, $n)
}
