# uia_dump.ps1 - launch AnimatLab2 on a project, dump top windows + full menu tree.
param([string]$Proj, [int]$WaitSec = 45)
$exe = "D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatLab2.exe"
$p = Start-Process -FilePath $exe -ArgumentList "`"$Proj`"" -PassThru
Start-Sleep -Seconds $WaitSec

Add-Type -AssemblyName UIAutomationClient
Add-Type -AssemblyName UIAutomationTypes

$root = [System.Windows.Automation.AutomationElement]::RootElement
$cond = New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ProcessIdProperty, $p.Id)
$wins = $root.FindAll([System.Windows.Automation.TreeScope]::Children, $cond)
foreach ($w in $wins) {
  Write-Output ("TOP: " + $w.Current.ClassName + " | " + $w.Current.Name + " | ctrl=" + $w.Current.ControlType.ProgrammaticName)
}

function DumpMenu($el, $depth) {
  $kids = $el.FindAll([System.Windows.Automation.TreeScope]::Children, [System.Windows.Automation.Condition]::TrueCondition)
  foreach ($k in $kids) {
    $ct = $k.Current.ControlType.ProgrammaticName
    $nm = $k.Current.Name
    if ($nm) { Write-Output (("  " * $depth) + $ct + " : " + $nm) }
    if ($ct -like "*MenuItem*") { DumpMenu $k ($depth + 1) }
  }
}

$main = $wins | Where-Object { $_.Current.ControlType.ProgrammaticName -eq "ControlType.Window" } | Select-Object -First 1
if ($main) {
  Write-Output "=== MENU TREE of main window ==="
  DumpMenu $main 0
}
Write-Output ("PID=" + $p.Id)
