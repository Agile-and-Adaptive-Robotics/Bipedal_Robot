@echo off
setlocal
set "BLENDER_EXE=%USERPROFILE%\Desktop\blender-5.3\blender.exe"
if not exist "%BLENDER_EXE%" (
  echo Blender was not found at "%BLENDER_EXE%".
  echo Update BLENDER_EXE in this file to your Blender executable.
  pause
  exit /b 1
)
if not exist "%~dp0build" mkdir "%~dp0build"
set "BLENDER_USER_CONFIG=%~dp0build\blender_config"
"%BLENDER_EXE%" -b --factory-startup -t 4 -P "%~dp0..\scripts\build_packaging_blender.py" > "%~dp0build\blender_build.log" 2>&1
if errorlevel 1 (
  type "%~dp0build\blender_build.log"
  echo Blender could not complete the scene. Check whether Blender itself opens.
  pause
  exit /b 1
)
echo Saved output\AH001_Packaging.blend and output\Blender_Packaging.png.
pause
