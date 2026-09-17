# gui_check.ps1 — launch AnimatLab2.exe on a project, detect error dialogs, kill.
# Usage: powershell -NoProfile -ExecutionPolicy Bypass -File gui_check.ps1 <aproj> [waitSeconds]
param([string]$Proj, [int]$WaitSec = 30)
$exe = "D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatLab2.exe"
$p = Start-Process -FilePath $exe -ArgumentList "`"$Proj`"" -PassThru
Start-Sleep -Seconds $WaitSec

Add-Type @"
using System;
using System.Text;
using System.Collections.Generic;
using System.Runtime.InteropServices;
public class WinEnum {
  [DllImport("user32.dll")] static extern bool EnumWindows(EnumWindowsProc cb, IntPtr lp);
  [DllImport("user32.dll")] static extern int GetClassName(IntPtr h, StringBuilder s, int n);
  [DllImport("user32.dll")] static extern int GetWindowText(IntPtr h, StringBuilder s, int n);
  [DllImport("user32.dll")] static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
  [DllImport("user32.dll")] static extern bool IsWindowVisible(IntPtr h);
  delegate bool EnumWindowsProc(IntPtr h, IntPtr lp);
  public static List<string> GetWindows(uint targetPid) {
    var list = new List<string>();
    EnumWindows((h, lp) => {
      uint pid; GetWindowThreadProcessId(h, out pid);
      if (pid == targetPid && IsWindowVisible(h)) {
        var cn = new StringBuilder(256); GetClassName(h, cn, 256);
        var ti = new StringBuilder(512); GetWindowText(h, ti, 512);
        list.Add(cn.ToString() + " | " + ti.ToString());
      }
      return true;
    }, IntPtr.Zero);
    return list;
  }
}
"@

$wins = [WinEnum]::GetWindows($p.Id)
$dialogs = @($wins | Where-Object { $_ -match '^#32770' })
$errors  = @($dialogs  | Where-Object { $_ -match '(?i)error|exception|fail' })
Write-Output ("PROCESS_ALIVE=" + (-not $p.HasExited))
Write-Output ("WINDOWS=" + $wins.Count)
foreach ($w in $wins) { Write-Output ("  WIN: " + $w) }
if ($dialogs.Count -gt 0) { Write-Output ("DIALOGS=" + $dialogs.Count) } else { Write-Output "DIALOGS=0" }
if ($errors.Count -gt 0)  { Write-Output ("ERROR_DIALOGS=" + $errors.Count) } else { Write-Output "ERROR_DIALOGS=0" }

# let a second wave appear (staggered load errors), then re-check
if ($dialogs.Count -gt 0) {
  Start-Sleep -Seconds 4
  $wins2 = [WinEnum]::GetWindows($p.Id)
  $d2 = @($wins2 | Where-Object { $_ -match '^#32770' });
  Write-Output ("DIALOGS_AFTER_4S=" + $d2.Count)
  foreach ($w in $d2) { Write-Output ("  DLG: " + $w) }
}
Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
Start-Sleep -Seconds 1
Get-Process -Name AnimatLab2 -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
Write-Output "KILLED"
