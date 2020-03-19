REM call luaGenerator.bat %model% %volumic% %nozzle% %layer% %filament% %ironing%

@echo off
@echo set_setting_value('printer', 'curvi')> settings.lua

@echo set_setting_value('filament_diameter_mm_0', %5)>> settings.lua
@echo set_setting_value('nozzle_diameter_mm_0', %3)>> settings.lua
@echo set_setting_value('z_layer_height_mm', %4)>> settings.lua
@echo set_setting_value('gcode_volumic', %2)>> settings.lua

@echo emit(load('%1/after.stl'))>> settings.lua

@echo set_service('FilamentSlicer')>> settings.lua
@echo run_service('%1.gcode')>> settings.lua
