param(
    [double]$Duration = 40,
    [double]$ActivationFloor = 0.05,
    [double]$RgTauH = 0.35,
    [ValidateSet('exact_fast','reference')][string]$Backend = 'exact_fast'
)
$ErrorActionPreference = 'Stop'
if ((hostname) -ine 'easteregg2') { throw 'This preset is for easteregg2.' }
$robotDirectory = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso'
$output = Join-Path $robotDirectory ('runs\support65_' + (Get-Date -Format 'yyyyMMdd_HHmmss'))
$runnerArguments = @(
    '--mode','step','--duration',"$Duration",'--stand-time','5',
    '--support','0.8','--support-cap-final','0.65',
    '--support-wean-start','5','--support-wean-duration','10',
    '--horizontal-stiffness','200','--lateral-shift','0.03',
    '--lateral-stiffness','1000','--lateral-damping','120',
    '--knee-amplitude','0.6','--knee-amplitude-right','0.85',
    '--hip-flexion','0.35','--hip-extension','-0.12',
    '--pattern','event','--unload-swing','--load-compensation','1',
    '--activation-floor',"$ActivationFloor",'--rg-tau-h',"$RgTauH",
    '--walk-drive','3','--backend',$Backend,'--output',$output
)
& 'D:\Anaconda\envs\myo\python.exe' (Join-Path $robotDirectory 'run_neural_robot.py') @runnerArguments
if ($LASTEXITCODE -ne 0) { throw 'Simulation process failed; inspect its output.' }
& 'D:\Anaconda\envs\myo\python.exe' (Join-Path $robotDirectory 'assess_supported_step.py') ($output + '.npz')
