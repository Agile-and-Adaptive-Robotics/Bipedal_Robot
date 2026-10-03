$wd = 'D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape'
$log = "$wd\logs\easteregg2_runs_20260930"
New-Item -ItemType Directory -Force -Path $log | Out-Null
$cmd = "cmd /c cd /d $wd && matlab -batch `"try, run('sns_units_test_2n.m'); catch e, disp(getReport(e)); end, try, sns_verify_from_json; catch e, disp(getReport(e)); end, try, open_system('KneeReflexDemo'); sim('KneeReflexDemo'); disp('KneeReflexDemo OK'); catch e, disp(getReport(e)); end, try, sim('BPACPGLegDemo'); disp('BPACPGLegDemo OK'); catch e, disp(getReport(e)); end`" >> $log\runs.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
"simscape launched"
