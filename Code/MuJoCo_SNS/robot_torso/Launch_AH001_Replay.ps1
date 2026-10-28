param([switch]$RunNew)
$ErrorActionPreference = 'Stop'
if ([Environment]::MachineName -ine 'easteregg2') { throw 'This launcher is configured for easteregg2.' }
$robotDirectory = $PSScriptRoot
$taskPython = 'D:\Anaconda\envs\myo\python.exe'
$sourceRun = Join-Path $robotDirectory 'runs\abd_defense_exact'
if ($RunNew) {
    $settings = (Get-Content -LiteralPath ($sourceRun + '.json') -Raw | ConvertFrom-Json).run_options
    $sourceRun = Join-Path $robotDirectory ('runs\interactive_abductor_' + (Get-Date -Format 'yyyyMMdd_HHmmss'))
    $settings.output = $sourceRun
    $runnerArguments = [Collections.Generic.List[string]]::new()
    foreach ($property in $settings.PSObject.Properties) {
        if ($null -eq $property.Value) { continue }
        $flag = '--' + $property.Name.Replace('_','-')
        if ($property.Value -is [bool]) {
            if ($property.Value) { $runnerArguments.Add($flag) }
        } else {
            $runnerArguments.Add($flag)
            $runnerArguments.Add([Convert]::ToString($property.Value,[Globalization.CultureInfo]::InvariantCulture))
        }
    }
    & $taskPython (Join-Path $robotDirectory 'run_neural_robot.py') @runnerArguments
    if ($LASTEXITCODE -ne 0) { throw 'The simulation failed; inspect its output.' }
}
Write-Host 'Opening the recorded AH001 model in MuJoCo. Close the viewer to finish.'
Write-Host 'Supported hybrid control with reference tensions; this is not a Simscape model.'
& $taskPython (Join-Path $robotDirectory 'replay_neural_robot.py') --run ($sourceRun + '.npz') --view
if ($LASTEXITCODE -ne 0) { throw 'The replay viewer failed; inspect its output.' }
