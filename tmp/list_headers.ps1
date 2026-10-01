$root = 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\campaigns\20260930\animatlab'
Get-ChildItem -Recurse $root -Filter *.txt | ForEach-Object {
    $h = Get-Content -TotalCount 1 $_.FullName
    $rel = $_.FullName.Substring($root.Length + 1)
    '{0} | {1} bytes | {2}' -f $rel, $_.Length, $h
}
