try {
    Start-Process -FilePath 'cmd.exe' `
        -ArgumentList '/c', '"D:\Anaconda\python.exe" "D:\GitHub\Bipedal_Robot\tmp\sw_cam\export_tibia008.py" > "D:\GitHub\Bipedal_Robot\tmp\sw_cam\elevated_run.log" 2>&1' `
        -Verb RunAs -Wait
    Write-Output 'ELEVATED RUN FINISHED'
} catch {
    Write-Output ('UAC/run failed: ' + $_.Exception.Message)
}
