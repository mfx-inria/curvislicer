@echo off
setlocal enabledelayedexpansion

set gurobi=0


set volumic=0
set nozzle=0.4
set layer=0.3
set filament=1.75
set ironing=0
set model=""


set arg=none

for %%A in (%*) do call :Loop %%A
goto :EndLoop

:Loop
  if "%arg%" EQU "none" (
    set arg=%1
  ) else (
    set %arg%=%1
    set arg=none
  )
  goto :End

:EndLoop
if "%arg%" EQU "none" (
  echo Error in arguments
  exit
)

set model=%arg%
set model=%model:\=/%

echo Generate tetmesh "from %model%.stl" ...
call toTetmesh.bat %model%
echo Done!

echo Optimize...

if "%gurobi%" EQU "1" (
  .\bin\curvislice_grb.exe ./ %model%.msh -l %layer%
) else (
  .\bin\curvislice_osqp.exe ./ %model%.msh -l %layer%
)
echo Done!

echo Prepare lua for IceSL
call luaGenerator.bat %model% %volumic% %nozzle% %layer% %filament% %ironing%

if not exist %appdata%\IceSL\icesl-printers\fff\curvi (
  echo Create 'curvi' printer profile for IceSL
  xcopy /S /I /Q /Y .\resources\curvi "%AppData%\IceSL\icesl-printers\fff\curvi"
)

.\tools\icesl\bin\icesl-slicer.exe settings.lua --service

.\bin\uncurve.exe -l %layer% --gcode %model%


:End
