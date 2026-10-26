# watch_dialogs.ps1 - poll a PID's top-level windows; log any dialog-looking ones
param([int]$ProcId, [int]$Seconds = 120, [string]$Out = "dialogs.txt")
Add-Type @"
using System;
using System.Text;
using System.Collections.Generic;
using System.Runtime.InteropServices;
public class WinDump2 {
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
$lines = New-Object System.Collections.Generic.List[string]
$seen = @{}
$sw = [Diagnostics.Stopwatch]::StartNew()
Add-Type @"
using System;
using System.Text;
using System.Collections.Generic;
using System.Runtime.InteropServices;
public class WinClick {
  [DllImport("user32.dll")] public static extern IntPtr SendMessage(IntPtr h, uint msg, IntPtr wp, IntPtr lp);
  [DllImport("user32.dll")] public static extern bool EnumChildWindows(IntPtr h, EnumChildWindowsProc cb, IntPtr lp);
  public delegate bool EnumChildWindowsProc(IntPtr h, IntPtr lp);
  [DllImport("user32.dll")] public static extern int GetClassName(IntPtr h, StringBuilder s, int n);
  [DllImport("user32.dll")] public static extern int GetWindowText(IntPtr h, StringBuilder s, int n);
  public static List<IntPtr> ChildButtons(IntPtr top) {
    var list = new List<IntPtr>();
    EnumChildWindows(top, (h, lp) => {
      var cn = new StringBuilder(256); GetClassName(h, cn, 256);
      if (cn.ToString() == "Button") list.Add(h);
      return true; }, IntPtr.Zero);
    return list;
  }
}
"@
while ($sw.Elapsed.TotalSeconds -lt $Seconds) {
  foreach ($h in [WinDump2]::TopWindows($ProcId)) {
    $ct = [WinDump2]::ClassAndTitle($h)
    if (-not $seen.ContainsKey($ct)) {
      $seen[$ct] = $true
      $lines.Add("TOP: " + $ct)
      foreach ($c in [WinDump2]::ChildrenText($h)) { $lines.Add("   CHILD " + $c) }
      # auto-dismiss Error dialogs (BM_CLICK = 0x00F5 on their OK button)
      if ($ct -match '\| Error$') {
        foreach ($b in [WinClick]::ChildButtons($h)) {
          $txt = (New-Object System.Text.StringBuilder(256))
          [void][WinClick]::SendMessage($b, 0x000D, [IntPtr]256, $txt)  # WM_GETTEXT
          if ($txt.ToString() -eq 'OK') {
            [void][WinClick]::SendMessage($b, 0x00F5, [IntPtr]::Zero, [IntPtr]::Zero)
            $lines.Add("   CLICKED OK on Error dialog")
          }
        }
      }
    }
  }
  Start-Sleep -Milliseconds 700
}
$lines | Out-File -FilePath $Out -Encoding utf8
Write-Output ("watched " + $ProcId + " for " + $Seconds + "s; windows: " + $lines.Count)
