@echo off
set a=%1
set a=%a:/=\%

if not exist %a%.msh (
  .\tools\tetwild-windows\TetWild.exe --save-mid-result 2 %a%.stl -e 0.025 --targeted-num-v 1000 --is-laplacian
  move "%a%__mid2.000000.msh" "%a%.msh"
)
