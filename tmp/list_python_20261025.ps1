Get-CimInstance Win32_Process | Where-Object { $_.Name -match 'python' } |
  ForEach-Object {
    $c = $_.CommandLine
    if ($c -and $c.Length -gt 130) { $c = $c.Substring(0, 130) }
    "{0}  {1}" -f $_.ProcessId, $c
  }
