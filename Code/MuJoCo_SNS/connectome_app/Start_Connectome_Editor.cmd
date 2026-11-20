@echo off
setlocal
if defined CONNECTOME_EDITOR_PYTHON (
    "%CONNECTOME_EDITOR_PYTHON%" "%~dp0app.py" %*
) else (
    python "%~dp0app.py" %*
)
exit /b %errorlevel%
