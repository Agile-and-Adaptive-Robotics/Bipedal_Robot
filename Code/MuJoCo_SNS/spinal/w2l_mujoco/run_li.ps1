$wd = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco'
$py = 'D:\Anaconda\envs\myo\python.exe'
$env:W2L_APROJ = 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj'
$cmd = "cmd /c cd /d $wd && set W2L_APROJ=D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj&& $py make_w2l_mjcf.py && $py fix_joint_axes.py && $py validate_body.py && $py test_li_stepping.py && $py test_w2l_air.py >> li_chain.log 2>&1"
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $cmd; CurrentDirectory = $wd } | Out-Null
"li chain launched"
