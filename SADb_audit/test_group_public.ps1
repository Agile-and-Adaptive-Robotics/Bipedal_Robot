# Test: is an AARL group PDF anonymously fetchable from api.zotero.org?
$ErrorActionPreference = 'Continue'
$row = Import-Csv 'D:\Github\Bipedal_Robot\SADb_audit\pdf_inventory.csv' |
    Where-Object { $_.library -eq 'groups/735051' -and $_.file_exists -eq 'True' } |
    Select-Object -First 1
if (-not $row) { Write-Output 'no group row found'; exit }
Write-Output ("att_key: " + $row.att_key + "  doi: " + $row.doi + "  surname: " + $row.surname)
curl.exe -s -L -m 30 "https://api.zotero.org/groups/735051/items/$($row.att_key)/file" -o zt_pub_test.bin -w "http=%{http_code}`n"
$fs = [System.IO.File]::OpenRead('D:\Github\Bipedal_Robot\SADb_audit\zt_pub_test.bin')
$head = New-Object byte[] 8
$null = $fs.Read($head, 0, 8)
$fs.Close()
$txt = [System.Text.Encoding]::ASCII.GetString($head)
Write-Output ("first bytes: " + $txt)
