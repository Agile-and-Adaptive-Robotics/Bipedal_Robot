$ErrorActionPreference = 'Stop'
$reportPath = 'D:\GitHub\Bipedal_Robot\CHATGPT_REPORT.md'
$evidenceDir = $PSScriptRoot
$blockPath = Join-Path $evidenceDir 'SESSION_1_UPDATE_20261027.md'
$startMarker = '<!-- SESSION_1_SIMULATION_RECOVERY_20261027_START -->'
$endMarker = '<!-- SESSION_1_SIMULATION_RECOVERY_20261027_END -->'
$utf8 = [Text.UTF8Encoding]::new($false, $true)
$replacement = [IO.File]::ReadAllText($blockPath, $utf8).TrimEnd("`r", "`n")
if (-not $replacement.StartsWith($startMarker) -or -not $replacement.EndsWith($endMarker)) {
    throw 'Replacement must contain exactly the bounded Session 1 block.'
}
$stream = $null
try {
    # Read the latest shared document under a brief exclusive-write lock.
    # Preserve every byte outside our own marker-bounded section.
    $stream = [IO.FileStream]::new($reportPath, [IO.FileMode]::Open, [IO.FileAccess]::ReadWrite, [IO.FileShare]::Read)
    $before = [byte[]]::new($stream.Length)
    $stream.ReadExactly($before)
    $current = $utf8.GetString($before)
    $start = $current.IndexOf($startMarker, [StringComparison]::Ordinal)
    $end = $current.IndexOf($endMarker, [StringComparison]::Ordinal)
    if ($start -lt 0 -or $end -lt $start -or $current.IndexOf($startMarker, $start + 1, [StringComparison]::Ordinal) -ge 0) {
        throw 'Expected one existing Session 1 block; no report changes made.'
    }
    $end += $endMarker.Length
    $prefix = $current.Substring(0, $start)
    $suffix = $current.Substring($end)
    $updated = $prefix + $replacement + $suffix
    $bytes = $utf8.GetBytes($updated)
    $prefixBytes = $utf8.GetBytes($prefix)
    $suffixBytes = $utf8.GetBytes($suffix)
    $originalPrefix = [byte[]]::new($prefixBytes.Length)
    $originalSuffix = [byte[]]::new($suffixBytes.Length)
    [Buffer]::BlockCopy($before, 0, $originalPrefix, 0, $originalPrefix.Length)
    [Buffer]::BlockCopy($before, $before.Length - $suffixBytes.Length, $originalSuffix, 0, $originalSuffix.Length)
    if ([Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($originalPrefix)) -cne
        [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($prefixBytes)) -or
        [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($originalSuffix)) -cne
        [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($suffixBytes))) {
        throw 'UTF-8 round-trip changed bytes outside Session 1; no report changes made.'
    }
    $stream.Position = 0
    $stream.Write($bytes)
    $stream.SetLength($bytes.Length)
    $stream.Flush($true)
    $stream.Position = 0
    $verified = [byte[]]::new($stream.Length)
    $stream.ReadExactly($verified)
    $verifiedText = $utf8.GetString($verified)
    if ($verifiedText -cne $updated) { throw 'Report verification failed.' }
    $record = [ordered]@{
        machine = [Environment]::MachineName
        session = 'Session 1'
        update_scope = 'Only existing Session 1 marker-bounded block'
        outside_block_bytes_preserved = $true
        preserved_prefix_bytes = $prefixBytes.Length
        preserved_suffix_bytes = $suffixBytes.Length
        prefix_sha256 = [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($prefixBytes))
        suffix_sha256 = [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($suffixBytes))
        report_bytes = $verified.Length
    }
    [IO.File]::WriteAllText((Join-Path $evidenceDir 'session1_report_update_verification.json'), ($record | ConvertTo-Json), $utf8)
    $record | ConvertTo-Json
} finally {
    if ($null -ne $stream) { $stream.Dispose() }
}
