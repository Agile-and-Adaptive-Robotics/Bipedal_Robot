# uia_menus.ps1 - attach to an AnimatLab2 PID and dump every MenuBar/MenuItem found.
param([int]$ProcId)
Add-Type -AssemblyName UIAutomationClient
Add-Type -AssemblyName UIAutomationTypes

$root = [System.Windows.Automation.AutomationElement]::RootElement
$cond = New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ProcessIdProperty, $ProcId)
$wins = $root.FindAll([System.Windows.Automation.TreeScope]::Children, $cond)
foreach ($w in $wins) { Write-Output ("TOP: " + $w.Current.ClassName + " | " + $w.Current.Name) }

$main = $wins | Where-Object { $_.Current.ControlType.ProgrammaticName -eq "ControlType.Window" } | Select-Object -First 1
if (-not $main) { Write-Output "NO MAIN WINDOW"; exit 1 }

# find menu items anywhere (descendants), bounded walk
$menuCond = New-Object System.Windows.Automation.PropertyCondition([System.Windows.Automation.AutomationElement]::ControlTypeProperty, [System.Windows.Automation.ControlType]::MenuItem)
$items = $main.FindAll([System.Windows.Automation.TreeScope]::Descendants, $menuCond)
Write-Output ("MenuItems found: " + $items.Count)

# Build parent chain map by walking tree once, depth-limited
function Walk($el, $path, $depth) {
  if ($depth -gt 8) { return }
  $kids = $el.FindAll([System.Windows.Automation.TreeScope]::Children, [System.Windows.Automation.Condition]::TrueCondition)
  foreach ($k in $kids) {
    $ct = $k.Current.ControlType.ProgrammaticName
    if ($ct -eq "ControlType.MenuItem" -or $ct -eq "ControlType.MenuBar" -or $ct -eq "ControlType.Menu") {
      Write-Output (($path) + " > " + $k.Current.Name + " [" + $ct + "]")
      Walk $k ($path + " > " + $k.Current.Name) ($depth + 1)
    } else {
      Walk $k $path ($depth + 1)
    }
  }
}
Walk $main "ROOT" 0
