$ErrorActionPreference = 'Stop'
$reportPath = 'D:\GitHub\Bipedal_Robot\CHATGPT_REPORT.md'
$evidenceDir = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261027_robot_control'
$blockPath = Join-Path $evidenceDir 'SESSION_1_UPDATE_20261027.md'
$marker = '<!-- SESSION_1_SIMULATION_RECOVERY_20261027_START -->'
$utf8 = [Text.UTF8Encoding]::new($false)
$block = [IO.File]::ReadAllText($blockPath, $utf8)
$stream = $null
try {
    # Briefly exclude simultaneous writes. Read, append, and verify on the same
    # handle; never rewrite Session 2's existing content from a stale snapshot.
    $stream = [IO.FileStream]::new($reportPath, [IO.FileMode]::Open, [IO.FileAccess]::ReadWrite, [IO.FileShare]::Read)
    $before = [byte[]]::new($stream.Length)
    [void]$stream.ReadExactly($before)
    if ($utf8.GetString($before).Contains($marker)) { throw 'Session 1 block already exists; use a bounded edit to that block.' }
    $backup = Join-Path $evidenceDir 'CHATGPT_REPORT_before_session1_append.md'
    if (-not [IO.File]::Exists($backup)) { [IO.File]::WriteAllBytes($backup, $before) }
    $addition = $utf8.GetBytes("`n`n" + $block + "`n")
    [void]$stream.Seek(0, [IO.SeekOrigin]::End)
    $stream.Write($addition, 0, $addition.Length)
    $stream.Flush($true)
    [void]$stream.Seek(0, [IO.SeekOrigin]::Begin)
    $prefix = [byte[]]::new($before.Length)
    [void]$stream.ReadExactly($prefix)
    $beforeHash = [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($before))
    $prefixHash = [Convert]::ToHexString([Security.Cryptography.SHA256]::HashData($prefix))
    if ($prefixHash -ne $beforeHash) { throw 'Existing report bytes changed unexpectedly; do not restore over another session.' }
    $verification = [ordered]@{ machine='easteregg2'; session='Session 1'; report=$reportPath; before_bytes=$before.Length; appended_bytes=$addition.Length; prior_sha256=$beforeHash; prior_bytes_preserved=$true; update_scope='append only' }
    [IO.File]::WriteAllText((Join-Path $evidenceDir 'session1_report_append_verification.json'), ($verification | ConvertTo-Json), $utf8)
    $verification | ConvertTo-Json
} finally {
    if ($null -ne $stream) { $stream.Dispose() }
}
