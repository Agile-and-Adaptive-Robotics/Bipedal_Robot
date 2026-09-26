param([int]$ExitCode = -9999)
# Writes the curr_syn6_s4 done-marker AFTER the python process has exited.
# First line = exit code; then the last 3 lines of the run log.
$log  = "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs\curr_syn6_s4.log"
$done = "D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs\curr_syn6_s4.done"
$tail = @()
try {
    $tail = @(Get-Content -Path $log -Tail 3 -Encoding UTF8 -ErrorAction Stop)
} catch {
    $tail = @("<log unreadable: $($_.Exception.Message)>")
}
$lines = @([string]$ExitCode) + @($tail)
# .NET WriteAllLines -> UTF-8 WITHOUT BOM, so line 1 parses as a bare int.
[System.IO.File]::WriteAllLines($done, $lines)
