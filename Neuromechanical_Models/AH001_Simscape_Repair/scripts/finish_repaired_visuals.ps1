$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..')).Path
Set-Location -LiteralPath $repoRoot
$auditRoot = Join-Path $repoRoot 'Code\Matlab\Mesh_Optimization\Functional_Routing_20261027\campaign_20261103\audit_20261105'
$buildRoot = Join-Path $repoRoot 'Neuromechanical_Models\AH001_Simscape_Repair\packaging\build'
$packagingRoot = Split-Path -Parent $buildRoot
$pythonExe = 'C:\Users\Ben\.anaconda3\envs\myoconv\python.exe'
$blenderExe = Join-Path $buildRoot 'blender_runtime\blender.exe'
$deadline = (Get-Date).AddHours(3)
while (-not (Test-Path -LiteralPath (Join-Path $auditRoot 'reoptimization\repaired_routes.json'))) {
    if ((Get-Date) -gt $deadline) { throw 'Timed out waiting for the MATLAB repaired-route export. Inspect repair_optimization4.log.' }
    Start-Sleep -Seconds 5
}
& $pythonExe (Join-Path $PSScriptRoot 'export_latest_routes_opensim.py') --repaired *> (Join-Path $buildRoot 'repaired_opensim_export.log')
if ($LASTEXITCODE -ne 0) { throw 'OpenSim export/verification failed. See repaired_opensim_export.log.' }
$env:BLENDER_USER_CONFIG = Join-Path $buildRoot 'blender_config'
& $blenderExe -b -t 4 -P (Join-Path $PSScriptRoot 'build_latest_opensim_blender.py') -- --repaired *> (Join-Path $buildRoot 'repaired_blender_build.log')
if ($LASTEXITCODE -ne 0) { throw 'Blender build failed. See repaired_blender_build.log.' }
& $blenderExe -b (Join-Path $auditRoot 'Repaired_Routes_OpenSim_Matched.blend') -t 4 -P (Join-Path $PSScriptRoot 'verify_repaired_blender.py') *> (Join-Path $buildRoot 'repaired_blender_verification.log')
if ($LASTEXITCODE -ne 0) { throw 'Saved Blender geometry verification failed. See repaired_blender_verification.log.' }
& $pythonExe (Join-Path $PSScriptRoot 'write_routing_audit.py') *> (Join-Path $buildRoot 'repaired_report.log')
if ($LASTEXITCODE -ne 0) { throw 'Audit report failed. See repaired_report.log.' }
Copy-Item -LiteralPath (Join-Path $auditRoot 'Repaired_Routes_OpenSim_Matched.blend') -Destination (Join-Path $packagingRoot 'Repaired_Robot_20261105.blend')
Copy-Item -LiteralPath (Join-Path $auditRoot 'Repaired_Routes_OpenSim_Matched.png') -Destination (Join-Path $packagingRoot 'Repaired_Robot_20261105.png')
Copy-Item -LiteralPath (Join-Path $auditRoot 'knee_comparison.mot') -Destination (Join-Path $packagingRoot 'knee_comparison_20261105.mot')
foreach ($branch in @('compliant','rigid')) {
    Copy-Item -LiteralPath (Join-Path $auditRoot "gait2392_robotbody_repaired_${branch}_blue.osim") -Destination (Join-Path $packagingRoot "gait2392_robotbody_repaired_${branch}_blue.osim")
}
Write-Output "READY: $(Join-Path $packagingRoot 'Repaired_Robot_20261105.blend')"
