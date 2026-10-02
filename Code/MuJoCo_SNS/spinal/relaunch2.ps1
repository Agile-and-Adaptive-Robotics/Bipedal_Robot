$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
$py = 'D:\Anaconda\envs\myo\python.exe'
# retune (per-study DBs now)
foreach ($j in @('air','walk 1','walk 2','walk 3')) {
  $n = $j.Split(' ')[0] + (@{'walk'='w'}[$j.Split(' ')[0]]) + $j.Split(' ')[-1]
  $n = if ($j -eq 'air') {'pruned_air'} else {'pruned_w' + $j.Split(' ')[1]}
  $cmd = "cmd /c cd /d $wd && $py _curriculum_pruned.py $j >> $n.log 2>&1"
  Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
  Start-Sleep -Milliseconds 500
}
# syn6 robust shards (bug fixed)
foreach ($s in 0,1) {
  $cmd = "cmd /c cd /d $wd && $py prune_robust.py syn6 $s 2 >> robust_syn6_$s.log 2>&1"
  Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
  Start-Sleep -Milliseconds 400
}
# li chain with GROUP redirection
$wdm = "$wd\w2l_mujoco"
$aproj = 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj'
$cmd = "cmd /c (cd /d $wdm && set W2L_APROJ=$aproj&& $py make_w2l_mjcf.py && $py fix_joint_axes.py && $py validate_body.py && $py test_li_stepping.py && $py test_w2l_air.py) >> $wdm\li_chain.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wdm } | Out-Null
# simscape with R2025a files + demos path
$wds = 'D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape'
$ml = "cmd /c cd /d $wds && matlab -batch `"addpath(pwd); addpath(fullfile(pwd,'demos')); addpath(fullfile(pwd,'figures')); try, run('sns_units_test_2n.m'); catch e, disp(getReport(e)); end, try, sim('KneeReflexDemo_R2025a'); disp('KneeReflexDemo_R2025a OK'); catch e, disp(getReport(e)); end`" >> $wds\logs\easteregg2_runs_20260930\runs2.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $ml; CurrentDirectory = $wds } | Out-Null
"relaunched: retune 4, syn6 robust 2, li chain 1, simscape 1"
