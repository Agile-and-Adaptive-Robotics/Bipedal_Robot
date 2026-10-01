$ErrorActionPreference = 'Stop'
$wd  = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\campaigns\20260930\animatlab\Li_rerun'
$exe = 'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe'
$asim = Join-Path $wd 'walk new new tester_Standalone.asim'
$out = Join-Path $wd 'simulator_stdout.log'
$err = Join-Path $wd 'simulator_stderr.log'
$t0 = Get-Date
$p = Start-Process -FilePath $exe -ArgumentList "`"$asim`"" -WorkingDirectory $wd `
     -RedirectStandardOutput $out -RedirectStandardError $err -PassThru
$exited = $p.WaitForExit(180000)
$dt = [math]::Round(((Get-Date) - $t0).TotalSeconds, 1)
if ($exited) {
    "RESULT exited code=$($p.ExitCode) after $dt s"
} else {
    "RESULT timeout after $dt s - killing"
    Stop-Process -Id $p.Id -Force
    Start-Sleep -Seconds 2
    "RESULT killed"
}
"--- files in workdir ---"
Get-ChildItem $wd | ForEach-Object { '{0}`t{1}' -f $_.Name, $_.Length }
"--- stdout ---"
if (Test-Path $out) { Get-Content $out | Select-Object -First 40 }
"--- stderr ---"
if (Test-Path $err) { Get-Content $err | Select-Object -First 40 }
