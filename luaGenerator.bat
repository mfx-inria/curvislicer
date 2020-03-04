@echo off
@echo set_setting_value('printer', curvi)> settings.lua

@echo set_setting_value('filament_diameter_mm_0', 3.85)>> settings.lua
@echo set_setting_value('nozzle_diameter_mm_0', 2.0)>> settings.lua
@echo set_setting_value('z_layer_height_mm', 1.0)>> settings.lua
@echo set_setting_value('gcode_volumic', false)>> settings.lua

@echo emit(difference(ccube(40),sphere(25)))>> settings.lua

@echo set_service('FilamentSlicer')>> settings.lua
@echo run_service(Path..'scratch/tmp.gcode')>> settings.lua