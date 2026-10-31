param([string]$SourceVideo='C:\Users\Ben Bolen\Downloads\2021-02-17-20-58-24.mp4')
$ErrorActionPreference='Stop'
$ffmpeg='D:\Anaconda\envs\myo\Library\bin\ffmpeg.exe'
New-Item -ItemType Directory -Force -Path (Join-Path $PSScriptRoot 'build')|Out-Null
& $ffmpeg -hide_banner -loglevel error -y -i $SourceVideo -vf 'crop=650:688:350:200' -an -c:v libx264 -crf 18 -preset fast -pix_fmt yuv420p -movflags +faststart (Join-Path $PSScriptRoot 'build\w2l_archived_motion.mp4')
if($LASTEXITCODE){throw 'Video preparation failed'}
& $ffmpeg -hide_banner -loglevel error -y -ss 0.4 -i (Join-Path $PSScriptRoot 'build\w2l_archived_motion.mp4') -frames:v 1 (Join-Path $PSScriptRoot 'build\w2l_poster.png')
if($LASTEXITCODE){throw 'Poster extraction failed'}
