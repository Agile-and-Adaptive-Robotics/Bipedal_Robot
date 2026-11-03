@echo off
setlocal
set "EDITOR_PYTHON=D:\Users\%USERNAME%\ConnectomeEditorRuntime\Scripts\pythonw.exe"
if defined CONNECTOME_EDITOR_PYTHON set "EDITOR_PYTHON=%CONNECTOME_EDITOR_PYTHON%"
if not exist "%EDITOR_PYTHON%" (
  echo Connectome Editor runtime was not found.
  echo Set CONNECTOME_EDITOR_PYTHON to a Python environment with requirements.txt installed.
  pause
  exit /b 1
)
start "Connectome Editor" "%EDITOR_PYTHON%" "%~dp0app.py" %*
