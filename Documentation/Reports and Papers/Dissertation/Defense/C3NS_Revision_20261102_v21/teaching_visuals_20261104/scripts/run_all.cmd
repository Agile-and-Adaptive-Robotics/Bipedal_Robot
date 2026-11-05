@echo off
rem run_all.cmd — rebuild every teaching visual end to end (this machine).
rem Paths are the laptop (DESKTOP-5Q16KE9) defaults; edit if they move.
setlocal
set PY=C:\Users\Ben\.anaconda3\python.exe
set BLENDER=C:\Users\Ben\.zcode\cli\plugins\data\blender@zcode-plugins-official\engine\blender-5.2.1-windows-x64\blender.exe
set HERE=%~dp0

pushd "%HERE%.."

echo == geometry prep ==
"%PY%" -I scripts\prep_geometry.py || goto :eof

echo == visual 3 ==
"%BLENDER%" --background --factory-startup --python scripts\build_visual3_frames.py || goto :eof
echo == visual 2 ==
"%BLENDER%" --background --factory-startup --python scripts\build_visual2_pulley.py || goto :eof
echo == visual 1 ==
"%BLENDER%" --background --factory-startup --python scripts\build_visual1_cam.py || goto :eof

echo == GIF assembly ==
"%PY%" -I scripts\make_gifs.py
popd
echo DONE
