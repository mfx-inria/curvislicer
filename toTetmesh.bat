@echo off

if not exist %1.msh ( 
.\tools\tetwild-windows\TetWild.exe %1.stl -e 0.025 --targeted-num-v 1000 --is-laplacian
.\tools\tetwild-windows\TetWild.exe --save-mid-result 2 %1__sf.obj -e 0.025 --targeted-num-v 1000 --is-laplacian

move %1__sf__mid2.000000.msh %1.msh
)

REM Clean files
for /f "delims=" %%a in ('dir /b /s /a-d %1*  ^| findstr /i /r ".*\.obj"') do (
  del "%%a"
)

for /f "delims=" %%a in ('dir /b /s /a-d %1*  ^| findstr /i /r ".*\.csv"') do (
  del "%%a"
)

for /f "delims=" %%a in ('dir /b /s /a-d %1*  ^| findstr /i /r ".*\.csv"') do (
  del "%%a"
)

for /f "delims=" %%a in ('dir /b /s /a-d %1*  ^| findstr /i /r ".*_\.msh"') do (
  del "%%a"
)