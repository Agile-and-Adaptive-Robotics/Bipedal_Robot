# il_calls_exe.ps1 - ordered IL call list for AnimatLab2.exe methods (32-bit)
param([string]$TypeName, [string]$MethodName)
$bin = 'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin'
$a = [Reflection.Assembly]::LoadFrom((Join-Path $bin 'AnimatLab2.exe'))
$t = $a.GetType($TypeName)
$mod = $a.ManifestModule
$flags = [Reflection.BindingFlags]'Public,NonPublic,Static,Instance,DeclaredOnly'
foreach ($m in ($t.GetMethods($flags) | Where-Object { $_.Name -eq $MethodName })) {
  Write-Output ("=== " + $m.ToString())
  $body = $m.GetMethodBody()
  if (-not $body) { Write-Output "  (no body)"; continue }
  $il = $body.GetILAsByteArray()
  for ($i = 0; $i -lt $il.Length; $i++) {
    if (($il[$i] -eq 0x28 -or $il[$i] -eq 0x6F) -and ($i + 4) -lt $il.Length) {
      $tok = [BitConverter]::ToInt32($il, $i + 1)
      try {
        $mi = $mod.ResolveMethod($tok)
        Write-Output ("  -> " + $mi.DeclaringType.FullName + '::' + $mi.Name)
      } catch {}
    }
    if ($il[$i] -eq 0x72 -and ($i + 4) -lt $il.Length) {
      $tok = [BitConverter]::ToInt32($il, $i + 1)
      try { Write-Output ("  str: " + $mod.ResolveString($tok)) } catch {}
    }
  }
}
