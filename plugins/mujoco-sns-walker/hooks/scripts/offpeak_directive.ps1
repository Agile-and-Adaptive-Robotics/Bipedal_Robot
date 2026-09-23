# UserPromptSubmit hook (mujoco-sns-walker plugin).
# Reads the prompt payload from stdin. When the user asks for unattended /
# off-peak / keep-running execution, injects a standing directive so the
# agent keeps iterating autonomously instead of pausing for a go-ahead.
# No phrase match -> silent exit 0 (no context injected).
# NOTE: keep this file PURE ASCII - Windows PowerShell 5.1 misparses
# BOM-less UTF-8 CJK literals. Chinese match phrases are built from
# [char] codepoints at runtime instead.

$ErrorActionPreference = 'SilentlyContinue'
$raw = [Console]::In.ReadToEnd()
if (-not $raw) { exit 0 }

try { $obj = $raw | ConvertFrom-Json } catch { exit 0 }

$parts = @()
foreach ($f in @('prompt', 'message', 'text')) {
  $v = $obj.$f
  if ($v) { $parts += [string]$v }
}
$text = $parts -join ' '
if (-not $text) { $text = $raw }

# Chinese phrases (built from codepoints so the file stays ASCII):
#   xianShi = "off-peak/idle-time", diFeng = "low-peak", guoYe = "overnight",
#   jiXuPao = "keep running", ziZhu = "autonomous"
$cjk = @(
  ([string][char]0x95F2 + [char]0x65F6),
  ([string][char]0x4F4E + [char]0x5CF0),
  ([string][char]0x8FC7 + [char]0x591C),
  ([string][char]0x7EE7 + [char]0x7EED + [char]0x8DD1),
  ([string][char]0x81EA + [char]0x4E3B)
) -join '|'

$pattern = '(?i)keep running|off[- ]?peak|over.?night|free compute|idle[- ]?time|unattended|keep going|while i.?m (away|asleep|out)|don.t (wait|stop)|autonomous|continue (running|iterating)|' + $cjk
if ($text -notmatch $pattern) { exit 0 }

$directive = @'
[STANDING DIRECTIVE - OFF-PEAK / UNATTENDED RUN]
Ben asked for continuous unattended execution. Until he explicitly says otherwise in a later message:
- Do NOT pause for go-ahead between rounds. When one tuning/eval round completes, analyze the results, decide the next round yourself, and launch it.
- Launch long runs as background tasks; poll their logs, never babysit turn-by-turn.
- Iterate the full loop autonomously: probe results -> diagnose -> adjust -> relaunch. Keep a running summary of what was tried and what won.
- If genuinely blocked on a decision ONLY Ben can make (connectome/wiring changes, destructive operations, dissertation .tex edits), finish every non-blocked item first, leave the question for the final report, and keep working around it.
- When idle-time/off-peak execution is what he asked for, use the OffPeakCreate tool instead of idling.
- End with a consolidated report: best result + numbers, full lineage of rounds, what remains open.
'@

[PSCustomObject]@{
  hookSpecificOutput = [PSCustomObject]@{
    hookEventName     = 'UserPromptSubmit'
    additionalContext = $directive
  }
} | ConvertTo-Json -Depth 4 -Compress
exit 0
