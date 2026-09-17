# gui_err_text.ps1 — launch AnimatLab2 on a project and dump the text of any
# Error/Exception window's child controls. Kills the app afterwards.
param([string]$Proj, [int]$WaitSec = 30)
$exe = "D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatLab2.exe"
$p = Start-Process -FilePath $exe -ArgumentList "`"$Proj`"" -PassThru
Start-Sleep -Seconds $WaitSec

Add-Type @"
using System;
using System.Text;
using System.Collections.Generic;
using System.Runtime.InteropServices;
public class WinDump {
  delegate bool EnumWindowsProc(IntPtr h, IntPtr lp);
  [DllImport("user32.dll")] static extern bool EnumWindows(EnumWindowsProc cb, IntPtr lp);
  [DllImport("user32.dll")] static extern bool EnumChildWindows(IntPtr h, EnumWindowsProc cb, IntPtr lp);
  [DllImport("user32.dll")] static extern int GetClassName(IntPtr h, StringBuilder s, int n);
  [DllImport("user32.dll")] static extern int GetWindowText(IntPtr h, StringBuilder s, int n);
  [DllImport("user32.dll")] static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
  [DllImport("user32.dll")] static extern bool IsWindowVisible(IntPtr h);
  public static List<IntPtr> TopWindows(uint targetPid) {
    var list = new List<IntPtr>();
    EnumWindows((h, lp) => { uint pid; GetWindowThreadProcessId(h, out pid);
      if (pid == targetPid && IsWindowVisible(h)) list.Add(h); return true; }, IntPtr.Zero);
    return list;
  }
  public static List<string> ChildrenText(IntPtr top) {
    var list = new List<string>();
    EnumChildWindows(top, (h, lp) => {
      var cn = new StringBuilder(256); GetClassName(h, cn, 256);
      var ti = new StringBuilder(1024); GetWindowText(h, ti, 1024);
      if (ti.Length > 0) list.Add(cn.ToString() + ": " + ti.ToString());
      return true; }, IntPtr.Zero);
    return list;
  }
  public static string ClassAndTitle(IntPtr h) {
    var cn = new StringBuilder(256); GetClassName(h, cn, 256);
    var ti = new StringBuilder(512); GetWindowText(h, ti, 512);
    return cn.ToString() + " | " + ti.ToString();
  }
}
"@

foreach ($h in [WinDump]::TopWindows($p.Id)) {
  $ct = [WinDump]::ClassAndTitle($h)
  Write-Output ("TOP: " + $ct)
  if ($ct -match '(?i)error|exception|warn') {
    foreach ($c in [WinDump]::ChildrenText($h)) { Write-Output ("   CHILD " + $c) }
  }
}
Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
Get-Process -Name AnimatLab2 -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
Write-Output "KILLED"
