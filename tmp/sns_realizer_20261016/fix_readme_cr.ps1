$p = 'D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\README_SNS_Simscape.md'
$c = [IO.File]::ReadAllText($p)
$bad = "spinal" + [char]13 + "eports_2026106"
$bad2 = "spinal" + [char]13 + "eports_20261016"
$good = "spinal\reports_20261016"
$n = 0
if ($c.Contains($bad2)) { $c = $c.Replace($bad2, $good); $n++ }
if ($c.Contains($bad))  { $c = $c.Replace($bad,  $good); $n++ }
# also handle the already-CR-stripped misspelling without backslash
$bad3 = "spinaleports_20261016"
if ($c.Contains($bad3)) { $c = $c.Replace($bad3, $good); $n++ }
[IO.File]::WriteAllText($p, $c, (New-Object System.Text.UTF8Encoding($false)))
Write-Output ("FIXED occurrences=" + $n)
