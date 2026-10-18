try {
    Start-Process -FilePath 'D:\Anaconda\python.exe' `
        -ArgumentList '"D:\GitHub\Bipedal_Robot\tmp\sw_cam\export_tibia008.py"' `
        -Verb RunAs -Wait
    Write-Output 'ELEVATED RUN FINISHED'
} catch {
    Write-Output ('UAC/run failed: ' + $_.Exception.Message)
}
