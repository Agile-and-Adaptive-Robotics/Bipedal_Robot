# uia_tree.ps1 - full descendant dump of the AnimatLab main window
param([int]$ProcId, [string]$Out = "tree.txt")
Add-Type -AssemblyName UIAutomationClient
Add-Type -AssemblyName UIAutomationTypes

$root = [System.Windows.Automation.AutomationElement]::RootElement
$cond = New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ProcessIdProperty, $ProcId)
$wins = $root.FindAll([System.Windows.Automation.TreeScope]::Children, $cond)
$main = $wins | Where-Object { $_.Current.ControlType.ProgrammaticName -eq "ControlType.Window" } | Select-Object -First 1
if (-not $main) { Write-Output "NO MAIN WINDOW"; exit 1 }

$sb = New-Object System.Text.StringBuilder
function Walk($el, $depth) {
  if ($depth -gt 10) { return }
  $kids = $el.FindAll([System.Windows.Automation.TreeScope]::Children, [System.Windows.Automation.Condition]::TrueCondition)
  foreach ($k in $kids) {
    $c = $k.Current
    $line = ("." * $depth) + $c.ControlType.ProgrammaticName + " | " + $c.ClassName + " | '" + $c.Name + "'"
    [void]$sb.AppendLine($line)
    Walk $k ($depth + 1)
  }
}
Walk $main 0
$sb.ToString() | Out-File -FilePath $Out -Encoding utf8
Write-Output ("lines: " + ($sb.ToString() -split "`n").Count)
