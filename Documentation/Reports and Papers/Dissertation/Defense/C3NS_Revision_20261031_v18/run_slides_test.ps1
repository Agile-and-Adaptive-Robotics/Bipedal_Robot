$ErrorActionPreference = 'Stop'
$testTemp = Join-Path $PSScriptRoot 'build\slides_test_temp'
New-Item -ItemType Directory -Force -Path $testTemp | Out-Null
$env:TEMP = $testTemp
$env:TMP = $testTemp
$env:RUNTIME_NODE = 'C:\Users\Ben Bolen\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe'
$env:RUNTIME_NODE_MODULES = 'C:\Users\Ben Bolen\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\node_modules'
& 'C:\Users\Ben Bolen\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  'C:\Users\Ben Bolen\.codex\plugins\cache\openai-primary-runtime\presentations\26.1017.22136\skills\presentations\container_tools\slides_test.py' `
  (Join-Path $PSScriptRoot 'output\Bolen_PhD_Defense_RoseGreen_20261031_v18.pptx')
