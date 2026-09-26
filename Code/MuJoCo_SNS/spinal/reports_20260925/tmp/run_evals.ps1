# Goal-1 evaluation batch: never-run SCONE tutorials, Hyfydy engine, fallback OpenSim.
# Writes: logs\<tag>.log (sconecmd -l 2 output) + logs\motion_<tag>.sto (motion, via -r)
$ErrorActionPreference = 'Continue'
$SCONE = 'C:\Program Files\SCONE\bin\sconecmd.exe'
$TUT   = 'C:\Users\Ben Bolen\Documents\SCONE\Tutorials3'
$LOGS  = 'D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs'

$evals = @(
    @{ tag='tut1_intro';              hyf='Tutorial 1 - Introduction - Hyfydy.scone';              os='Tutorial 1 - Introduction - OpenSim.scone' },
    @{ tag='tut2a_highjump';          hyf='Tutorial 2a - High Jump - Hyfydy.scone';                os='Tutorial 2a - High Jump - OpenSim.scone' },
    @{ tag='tut2b_jump_poly';         hyf='Tutorial 2b - High Jump Polynomial - Hyfydy.scone';     os='Tutorial 2b - High Jump Polynomial - OpenSim.scone' },
    @{ tag='tut2c_straight_pose';     hyf='Tutorial 2c - Straight Pose Jump - Hyfydy.scone';       os='Tutorial 2c - Straight Pose Jump - OpenSim.scone' },
    @{ tag='tut3b_motor_noise_bal';   hyf='Tutorial 3b - Motor Noise Balance - Hyfydy.scone';      os='Tutorial 3b - Motor Noise Balance - OpenSim.scone' },
    @{ tag='tut4b_fast_gait';         hyf='Tutorial 4b - Fast Gait - Hyfydy.scone';                os='Tutorial 4b - Fast Gait - OpenSim.scone' },
    @{ tag='tut4c_perturbed_gait';    hyf='Tutorial 4c - Perturbed Gait - Hyfydy.scone';           os='Tutorial 4c - Perturbed Gait - OpenSim.scone' },
    @{ tag='tut4d_slippery_slope';    hyf='Tutorial 4d - Slippery Slope - Hyfydy.scone';           os='Tutorial 4d - Slippery Slope - OpenSim.scone' },
    @{ tag='tut5a_pf_weakness';       hyf='Tutorial 5a - Plantarflexor Weakness - Hyfydy.scone';   os='Tutorial 5a - Plantarflexor Weakness - OpenSim.scone' },
    @{ tag='tut5b_short_hamstrings';  hyf='Tutorial 5b - Short Hamstrings - Hyfydy.scone';         os='Tutorial 5b - Short Hamstrings - OpenSim.scone' },
    @{ tag='tut5c_hyper_reflexia';    hyf='Tutorial 5c - Hyper-reflexia - Hyfydy.scone';           os='Tutorial 5c - Hyper-reflexia - OpenSim.scone' },
    @{ tag='tut6a_script_bodyheight'; hyf='Tutorial 6a - Script - Body Height - Hyfydy.scone';     os='Tutorial 6a - Script - Body Height - OpenSim.scone' }
)

foreach ($e in $evals) {
    $tag = $e.tag
    $scen = Join-Path $TUT $e.hyf
    $log  = Join-Path $LOGS ($tag + '.log')
    $mot  = Join-Path $LOGS ('motion_' + $tag)
    Write-Host ("##### EVAL {0}  [{1}]" -f $tag, (Get-Date -Format 'HH:mm:ss'))
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    & $SCONE -e $scen -l 2 -r $mot 2>&1 | Tee-Object -FilePath $log
    $code = $LASTEXITCODE
    $sw.Stop()
    Write-Host ("##### {0} exit={1} wall={2:N1}s" -f $tag, $code, $sw.Elapsed.TotalSeconds)
    if ($code -ne 0) {
        Write-Host ("##### FALLBACK to OpenSim variant for {0}" -f $tag)
        $scen2 = Join-Path $TUT $e.os
        $log2  = Join-Path $LOGS ($tag + '_OS_FALLBACK.log')
        $mot2  = Join-Path $LOGS ('motion_' + $tag + '_OS')
        & $SCONE -e $scen2 -l 2 -r $mot2 2>&1 | Tee-Object -FilePath $log2
        Write-Host ("##### {0}_OS exit={1}" -f $tag, $LASTEXITCODE)
    }
}
Write-Host '##### ALL EVALS DONE'
