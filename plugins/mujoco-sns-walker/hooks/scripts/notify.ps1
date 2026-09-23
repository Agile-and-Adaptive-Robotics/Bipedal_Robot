# Windows notification with fallback chain: toast -> msg.exe -> beep.
# Used by the PreCompact hook (high-context warning) and reusable from
# any other hook or script. Always beeps so the user notices something.
param(
  [string]$Title = 'ZCode',
  [string]$Message = 'notification'
)
$ErrorActionPreference = 'SilentlyContinue'
$ok = $false
try {
  $null = [Windows.UI.Notifications.ToastNotificationManager, Windows.UI.Notifications, ContentType = WindowsRuntime]
  $null = [Windows.Data.Xml.Dom.XmlDocument, Windows.Data.Xml.Dom.XmlDocument, ContentType = WindowsRuntime]
  $xml = [Windows.UI.Notifications.ToastNotificationManager]::GetTemplateContent([Windows.UI.Notifications.ToastTemplateType]::ToastText02)
  $texts = $xml.GetElementsByTagName('text')
  $null = $texts.Item(0).AppendChild($xml.CreateTextNode($Title))
  $null = $texts.Item(1).AppendChild($xml.CreateTextNode($Message))
  $toast = [Windows.UI.Notifications.ToastNotification]::new($xml)
  # PowerShell's own AppUserModelID gives the toast an identity on Win10/11.
  $appId = '{1AC14E77-02E7-4E5D-B744-2EB1AE5198B7}\WindowsPowerShell\v1.0\powershell.exe'
  [Windows.UI.Notifications.ToastNotificationManager]::CreateToastNotifier($appId).Show($toast)
  $ok = $true
} catch { $ok = $false }
if (-not $ok) {
  & msg.exe * "$Title : $Message" 2>$null
  if ($LASTEXITCODE -eq 0) { $ok = $true }
}
[Console]::Beep(1000, 250)
exit 0
