$ErrorActionPreference = "Stop"
$out = "D:\Github\Bipedal_Robot\SAD_audit"

function Norm-Loose([string]$s) {
  if ([string]::IsNullOrWhiteSpace($s)) { return "" }
  $t = ([regex]::Replace($s, '\s+', ' ')).Trim().ToLower()
  $t = [regex]::Replace($t, '[^a-z0-9]+', ' ')
  return ([regex]::Replace($t, '\s+', ' ')).Trim()
}
function Read-CsvUtf8([string]$path) {
  return (Get-Content -LiteralPath $path -Encoding UTF8 | ConvertFrom-Csv)
}

# ---- Biology side: drop the 10 rows whose misses were adjudicated as matches ----
$adjKeys = @('HBAEC9DM','E2G75RYG','VKC5THEJ','HJUK3CTI','PUVMGVWH','HQZAC8XN','FIN3YLAZ','BGU8UZJE','JA7MEXI6','YZSP4ZPD')
$bioMiss = Read-CsvUtf8 "$out\Biology_not_in_airtable.csv"
$bioCorr = @($bioMiss | Where-Object { $adjKeys -notcontains $_.key })
$bioMiss | Where-Object { $adjKeys -contains $_.key } | ForEach-Object { Write-Output ("  dropped from Biology gap: " + $_.key + " :: " + $_.title.Substring(0, [Math]::Min(80, $_.title.Length))) }
$bioCorr | Export-Csv "$out\Biology_not_in_airtable_corrected.csv" -NoTypeInformation -Encoding utf8
Write-Output ("Biology missing: raw " + @($bioMiss).Count + " -> corrected " + @($bioCorr).Count)

# ---- Personal side: drop rows matching the 7 adjudicated papers by loose title ----
$adjTitles = @(
  'A role for hip position in initiating the swing-to-stance transition in walking cats.',
  'Leg Coordination Mechanisms in the Stick Insect Applied to Hexapod Robot Locomotion',
  'Speed dependency in α-motoneuron activity and locomotor modules in human locomotion: Indirect evidence for phylogenetically conserved spinal circuits',
  'Five basic muscle activation patterns account for muscle activity during human locomotion: Basic muscle activation patterns',
  'Contribution of hind limb flexor muscle afferents to the timing of phase transitions in the cat step cycle',
  'Bio-inspired controller achieving forward speed modulation with a 3D bipedal walker',
  'Relative Contribution of Proprioceptive and Vestibular Sensory Systems to Locomotion: Opportunities for Discovery in the Age of Molecular Science.'
)
$adjLoose = @{}
foreach ($t in $adjTitles) { $adjLoose[(Norm-Loose $t)] = 1 }
$perMiss = Read-CsvUtf8 "$out\personalSAD_not_in_airtable.csv"
$perCorr = New-Object System.Collections.ArrayList
foreach ($p in $perMiss) {
  $kl = Norm-Loose $p.title
  if ($kl -and $adjLoose.ContainsKey($kl)) {
    Write-Output ("  dropped from personal gap: " + $p.key + " :: " + $p.title.Substring(0, [Math]::Min(80, $p.title.Length)))
  } else { [void]$perCorr.Add($p) }
}
$perCorr | Export-Csv "$out\personalSAD_not_in_airtable_corrected.csv" -NoTypeInformation -Encoding utf8
Write-Output ("Personal missing: raw " + @($perMiss).Count + " -> corrected " + @($perCorr).Count)
Write-Output "DONE"
