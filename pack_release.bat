@echo off

set zip_name=curvislice

set bin_folder="bin"
set icesl_folder="tools/icesl"
set tetwild_folder="tools/tetwild-windows"
set profile_folder="resources/curvi"
set models_folder="models"

set curvi_script="curvislice.bat"
set lua_script="luaGenerator.bat"
set grb_script="optimize_grb.bat"
set osqp_script="optimize_osqp.bat"
set tetmesh_script="toTetmesh.bat"

set _7zip="C:\Program Files\7-Zip\7z.exe"

echo Zipping ...

%_7zip% a -tzip "%zip_name%.zip" %bin_folder% %icesl_folder% %tetwild_folder% %profile_folder% %models_folder% %curvi_script% %lua_script% %grb_script% %osqp_script% %tetmesh_script%

echo Done !
echo Release Zip is %zip_name%.zip