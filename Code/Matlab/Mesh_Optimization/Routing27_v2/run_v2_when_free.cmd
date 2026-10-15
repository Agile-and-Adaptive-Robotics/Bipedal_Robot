@echo off
rem Watcher: wait for the v1 campaign's MATLAB to exit, then launch v2.
rem v1 (Routing27) is never touched; v2 starts only when the box is free.
setlocal
set LOG=D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\Results\watcher.log
if not exist "D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\Results" mkdir "D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\Results"
echo [%date% %time%] watcher armed, waiting for v1 MATLAB to exit >> "%LOG%"
:waitloop
tasklist /FI "IMAGENAME eq MATLAB.exe" 2>nul | find /i "MATLAB.exe" >nul
if %ERRORLEVEL%==0 (
    ping -n 120 127.0.0.1 >nul
    goto waitloop
)
echo [%date% %time%] v1 MATLAB gone; validating v2 harness first >> "%LOG%"
call "D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\run_routing27_v2.cmd" validate
echo [%date% %time%] validate exit %ERRORLEVEL%; checking feasibility gate >> "%LOG%"
findstr /r /c:",1,[0-9-]*,[0-9.]*,[0-9]*,[0-9.]*,[0-9.]*,[0-9]*,.*,.*,1,.*0$" "D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\Results\routing27_summary_validate.csv" >nul 2>&1
if not %ERRORLEVEL%==0 (
    findstr /i "FEASIBLE" "D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\Results\routing27_validate.log" >nul 2>&1
)
echo [%date% %time%] launching v2 full >> "%LOG%"
start "" /b cmd /c "D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\run_routing27_v2.cmd full" >nul 2>&1
echo [%date% %time%] v2 full launched >> "%LOG%"
endlocal
