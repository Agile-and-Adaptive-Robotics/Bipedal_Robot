# Goal-1 re-evaluations of best pars + optional 6b/6c/6d evals.
$ErrorActionPreference = 'Continue'
$SCONE = 'C:\Program Files\SCONE\bin\sconecmd.exe'
$RES   = 'C:\Users\Ben Bolen\Documents\SCONE\results'
$TUT   = 'C:\Users\Ben Bolen\Documents\SCONE\Tutorials3'
$LOGS  = 'D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\logs'
$GDIR  = Join-Path $env:TEMP 'scone_opt_gait'
$BDIR  = Join-Path $env:TEMP 'scone_opt_balance'

function Reeval-InPlace($dir, $par, $tag, $motion) {
    Push-Location $dir
    Write-Host ("##### REEVAL {0}" -f $tag)
    if ($motion) {
        & $SCONE -e $par -l 2 -r (Join-Path $LOGS ('motion_' + $tag)) 2>&1 |
            Tee-Object -FilePath (Join-Path $LOGS ($tag + '.log'))
    } else {
        & $SCONE -e $par -l 2 2>&1 |
            Tee-Object -FilePath (Join-Path $LOGS ($tag + '.log'))
    }
    Write-Host ("##### {0} exit={1}" -f $tag, $LASTEXITCODE)
    Pop-Location
}

# --- R1/R2: gait results dir (config has max_duration=4) ---
Reeval-InPlace (Join-Path $RES '260925.130052.H0918v3.RS2.S10WA3K1G14.D4.R42') '0044_0.873_0.858.par' 'reeval_gait4a_best_d4' $true
Reeval-InPlace (Join-Path $RES '260925.130052.H0918v3.RS2.S10WA3K1G14.D4.R42') '0000_23.911_0.901.par' 'reeval_gait4a_init_d4' $false
# --- R3: balance 40-gen best at its D12 ---
Reeval-InPlace (Join-Path $RES '260925.130211.H0918v3.R36.BW.D12.R42') '0037_94.802_80.445.par' 'reeval_bal40_best_d12' $true
# --- R4: balance 300-gen best at its D12 ---
Reeval-InPlace (Join-Path $RES '260925.130428.H0918v3.R36.BW.D12.R42') '0293_37.931_2.207.par' 'reeval_bal300_best_d12' $true
# --- R5: high-jump best at its D2 ---
Reeval-InPlace (Join-Path $RES '260925.130348.H0918v3.FC2.Jump.D2.R42') '0019_82.742_93.606.par' 'reeval_jump2a_best' $true

# --- R6/R7: gait best + init at FULL scenario default duration 20 (scratch config.scone) ---
Copy-Item (Join-Path $RES '260925.130052.H0918v3.RS2.S10WA3K1G14.D4.R42\0044_0.873_0.858.par') (Join-Path $GDIR 'best44.par') -Force
Copy-Item (Join-Path $RES '260925.130052.H0918v3.RS2.S10WA3K1G14.D4.R42\0000_23.911_0.901.par') (Join-Path $GDIR 'init00.par') -Force
Copy-Item (Join-Path $GDIR 'Tutorial 4a - Gait - Hyfydy.scone') (Join-Path $GDIR 'config.scone') -Force
Reeval-InPlace $GDIR 'best44.par' 'reeval_gait4a_best_d20' $true
Reeval-InPlace $GDIR 'init00.par' 'reeval_gait4a_init_d20' $true

# --- R8/R9/R10: balance bests + init at FULL scenario default duration 30 ---
Copy-Item (Join-Path $RES '260925.130211.H0918v3.R36.BW.D12.R42\0037_94.802_80.445.par') (Join-Path $BDIR 'best37.par') -Force
Copy-Item (Join-Path $RES '260925.130428.H0918v3.R36.BW.D12.R42\0293_37.931_2.207.par') (Join-Path $BDIR 'best293.par') -Force
Copy-Item (Join-Path $RES '260925.130428.H0918v3.R36.BW.D12.R42\0000_99.012_96.562.par') (Join-Path $BDIR 'init00.par') -Force
Copy-Item (Join-Path $BDIR 'Tutorial 3a - Balance - Hyfydy.scone') (Join-Path $BDIR 'config.scone') -Force
Reeval-InPlace $BDIR 'best37.par' 'reeval_bal40_best_d30' $true
Reeval-InPlace $BDIR 'best293.par' 'reeval_bal300_best_d30' $true
Reeval-InPlace $BDIR 'init00.par' 'reeval_bal_init_d30' $true

# --- R11-R13: optional Lua tutorials 6b/6c/6d (Hyfydy, fallback OpenSim) ---
$lua = @(
    @{ tag='tut6b_gyro_balance'; hyf='Tutorial 6b - Script - Gyro Balance Gait - Hyfydy.scone'; os='Tutorial 6b - Script - Gyro Balance Gait - OpenSim.scone' },
    @{ tag='tut6c_reflex_mod';   hyf='Tutorial 6c - Script - Reflex Modulation - Hyfydy.scone'; os=$null },
    @{ tag='tut6d_neural_delays';hyf='Tutorial 6d - Script - Neural Delays - Hyfydy.scone';     os=$null }
)
foreach ($e in $lua) {
    $tag  = $e.tag
    $scen = Join-Path $TUT $e.hyf
    $log  = Join-Path $LOGS ($tag + '.log')
    $mot  = Join-Path $LOGS ('motion_' + $tag)
    Write-Host ("##### EVAL {0}" -f $tag)
    & $SCONE -e $scen -l 2 -r $mot 2>&1 | Tee-Object -FilePath $log
    Write-Host ("##### {0} exit={1}" -f $tag, $LASTEXITCODE)
    if ($LASTEXITCODE -ne 0 -and $e.os) {
        Write-Host ("##### FALLBACK to OpenSim variant for {0}" -f $tag)
        & $SCONE -e (Join-Path $TUT $e.os) -l 2 -r (Join-Path $LOGS ('motion_' + $tag + '_OS')) 2>&1 |
            Tee-Object -FilePath (Join-Path $LOGS ($tag + '_OS_FALLBACK.log'))
        Write-Host ("##### {0}_OS exit={1}" -f $tag, $LASTEXITCODE)
    }
}
Write-Host '##### ALL REEVALS DONE'
