# il_calls.ps1 - list method calls inside a method body (32-bit PS, AnimatGUI.dll)
param([string]$TypeName, [string]$MethodName)
$bin = 'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin'
$a = [Reflection.Assembly]::LoadFrom((Join-Path $bin 'AnimatGUI.dll'))
$t = $a.GetType($TypeName)
$flags = [Reflection.BindingFlags]'Public,NonPublic,Instance,Static,DeclaredOnly'
$all = $t.GetMethods($flags) | Where-Object { $_.Name -eq $MethodName }
$mod = $a.ManifestModule
foreach ($m in $all) {
  Write-Output ("=== " + $m.ToString())
  try {
    $body = $m.GetMethodBody()
    if (-not $body) { Write-Output "  (no body)"; continue }
    $il = $body.GetILAsByteArray()
    $seen = @{}
    for ($i = 0; $i -lt $il.Length; $i++) {
      $op = $il[$i]
      $isCall = ($op -eq 0x28 -or $op -eq 0x6F)   # call, callvirt
      if ($isCall -and ($i + 4) -lt $il.Length) {
        $tok = [BitConverter]::ToInt32($il, $i + 1)
        try {
          $mi = $mod.ResolveMethod($tok)
          $key = $mi.DeclaringType.FullName + '::' + $mi.Name
          Write-Output ("  -> " + $key)
        } catch {}
      }
    }
    Write-Output ("  strings:")
    foreach ($s in $body.LocalVariables) { }
    # dump ldstr contents: rescan for 0x72 tokens via ResolveString
    $seenS = @{}
    for ($i = 0; $i -lt $il.Length; $i++) {
      if ($il[$i] -eq 0x72 -and ($i + 4) -lt $il.Length) {
        $tok = [BitConverter]::ToInt32($il, $i + 1)
        try {
          $s = $mod.ResolveString($tok)
          if (-not $seenS.ContainsKey($s)) { $seenS[$s] = $true; Write-Output ("  str: " + $s) }
        } catch {}
      }
    }
  } catch { Write-Output ("  ERR " + $_.Exception.Message) }
}
