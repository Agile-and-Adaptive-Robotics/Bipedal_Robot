# export_run.ps1 - headless AnimatLab standalone-sim export via AnimatGUI API (32-bit)
param(
  [string]$Aproj,
  [string]$OutAsim,
  [int]$WaitLoadSec = 40
)
$ErrorActionPreference = 'Continue'
$bin = 'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin'
Set-Location $bin
[Environment]::CurrentDirectory = $bin   # PowerShell does NOT sync this with Set-Location
Write-Output ("cwd: " + (Get-Location) + " ; env: " + [Environment]::CurrentDirectory)

# probe bin for dependent assemblies (string-only handler, reentrancy-guarded)
$script:inResolve = $false
$handler = [ResolveEventHandler] {
  param($s, $e)
  if ($script:inResolve) { return $null }
  $script:inResolve = $true
  try {
    $simple = $e.Name.Split(',')[0]
    $p = Join-Path 'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin' ($simple + '.dll')
    if (Test-Path $p) { return [Reflection.Assembly]::LoadFrom($p) }
    return $null
  } finally { $script:inResolve = $false }
}
[AppDomain]::CurrentDomain.add_AssemblyResolve($handler)

foreach ($d in @('AnimatGUI.dll','AnimatGUICtrls.dll','ManagedAnimatTools.dll','ManagedAnimatInterfaces.dll',
                 'LicensedAnimatGUI.dll','IntegrateFireGUI.dll','FiringRateGUI.dll','AnimatServer.dll')) {
  try { [void][Reflection.Assembly]::LoadFrom((Join-Path $bin $d)); Write-Output ("loaded " + $d) }
  catch { Write-Output ("LOAD-ERR " + $d + ": " + $_.Exception.Message) }
}

try {
  Write-Output "instantiating AnimatApplication..."
  $app = New-Object AnimatGUI.Forms.AnimatApplication
  Write-Output ("app OK: " + $app.GetType().FullName)
  Write-Output "Initialize()..."
  $app.Initialize()
  Write-Output "Initialize OK"
} catch {
  Write-Output ("APP-INIT-FAILED: " + $_.Exception.ToString())
  exit 2
}

try {
  Write-Output ("LoadProject: " + $Aproj)
  $app.LoadProject($Aproj)
  Write-Output ("ProjectIsOpen: " + $app.ProjectIsOpen)
} catch {
  Write-Output ("LOAD-FAILED: " + $_.Exception.ToString())
  exit 3
}

Start-Sleep -Seconds 5

try {
  Write-Output "ExportStandAloneSim..."
  $app.ExportStandAloneSim($Aproj)
  Write-Output "export call returned"
} catch {
  Write-Output ("EXPORT-FAILED: " + $_.Exception.ToString())
}

# report what new .asim files exist in likely spots
foreach ($dir in @((Split-Path $Aproj), (Get-Location).Path)) {
  Get-ChildItem (Join-Path $dir '*.asim') -ErrorAction SilentlyContinue |
    ForEach-Object { Write-Output ("ASIM: " + $_.FullName + " " + $_.Length + " " + $_.LastWriteTime) }
}
try { $app.CloseQuiet() } catch {}
Write-Output "DONE"
