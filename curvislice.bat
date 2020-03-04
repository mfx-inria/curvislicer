@echo off
setlocal enabledelayedexpansion

set volumic=0
set nozzle=0.4
set layer=0.3
set filament=1.75
set gurobi=0


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

echo Generate tetmesh "from %model%.stl" ...
call toTetmesh.bat %model%
echo Done!

set optimize="optimize_osqp.bat"
if "%gurobi%" EQU "1" (
  set optimize="optimize_grb.bat"
)

echo Optimize...
call %optimize% %model% -l %layer%
echo Done!

:End