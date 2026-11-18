param(
    [Parameter(Mandatory)][string]$Deck,
    [Parameter(Mandatory)][string]$Pdf
)
$ErrorActionPreference = 'Stop'
$sourcePath = (Resolve-Path -LiteralPath $Deck).Path
$outputPath = [IO.Path]::GetFullPath($Pdf)
$originalHash = (Get-FileHash -LiteralPath $sourcePath).Hash
if ($originalHash -ne '9AF4C36A3EF2C70FC19D8801A687BBE03CCFC2A35094D404FBD097960F624370') {
    throw 'The input must be the verified EB475WS4 thumb-drive defense deck.'
}
if (Test-Path -LiteralPath $outputPath) { throw 'Choose a new output PDF.' }
$powerPoint = New-Object -ComObject PowerPoint.Application
$presentation = $null
try {
    $presentation = $powerPoint.Presentations.Open($sourcePath, -1, 0, 0)
    $presentation.SaveAs($outputPath, 32)
} finally {
    if ($null -ne $presentation) { $presentation.Close() }
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($powerPoint)
}
if ((Get-FileHash -LiteralPath $sourcePath).Hash -ne $originalHash) { throw 'Source hash changed.' }
