@echo off
set a=%1
set a=%a:/=\%

if not exist %a%.msh (
  .\tools\tetwild-windows\TetWild.exe --save-mid-result 2 %a%.stl -e 0.001 --targeted-num-v 1000 --is-laplacian
  move "%a%__mid2.000000.msh" "%a%.msh"
)

REM Clean files
for /f "delims=" %%a in ('dir /b /s /a-d %a%*  ^| findstr /i /r ".*\.obj"') do (
  del "%%a"
)

for /f "delims=" %%a in ('dir /b /s /a-d %a%*  ^| findstr /i /r ".*\.csv"') do (
  del "%%a"
)

for /f "delims=" %%a in ('dir /b /s /a-d %a%*  ^| findstr /i /r ".*\.csv"') do (
  del "%%a"
)

for /f "delims=" %%a in ('dir /b /s /a-d %a%*  ^| findstr /i /r ".*_\.msh"') do (
  del "%%a"
)
