$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal'
$py = 'D:\Anaconda\envs\myo\python.exe'
foreach ($j in @(@('air','pruned_air'), @('walk 1','pruned_w1'), @('walk 2','pruned_w2'), @('walk 3','pruned_w3'))) {
  $cmd = "cmd /c cd /d $wd && $py _curriculum_pruned.py $($j[0]) >> $($j[1]).log 2>&1"
  Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
  Start-Sleep -Milliseconds 500
}
$cmd = "cmd /c cd /d $wd && $py gait_validate_all.py >> gait_validate.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
$wdm = "$wd\w2l_mujoco"
$aproj = 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj'
$cmd = "cmd /c (cd /d $wdm && set W2L_APROJ=$aproj&& $py make_w2l_mjcf.py && $py fix_joint_axes.py && $py validate_body.py && $py test_li_stepping.py && $py test_w2l_air.py) >> $wdm\li_chain.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wdm } | Out-Null
"relaunched: retune 4, validation 1, li chain 1"
