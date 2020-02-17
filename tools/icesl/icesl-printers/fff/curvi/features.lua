version = 2

bed_size_x_mm = 220
bed_size_y_mm = 220
bed_size_z_mm = 240
nozzle_diameter_mm_0 = 0.4

extruder_count = 1

z_offset   = 0.0

priming_mm_per_sec = 40


z_layer_height_mm = 0.3
z_layer_height_mm_min = 0.01
z_layer_height_mm_max = 10.0

print_speed_mm_per_sec_min = 5
print_speed_mm_per_sec_max = 80

bed_temp_degree_c = 50
bed_temp_degree_c_min = 0
bed_temp_degree_c_max = 120

perimeter_print_speed_mm_per_sec_min = 5
perimeter_print_speed_mm_per_sec_max = 80

first_layer_print_speed_mm_per_sec = 10
first_layer_print_speed_mm_per_sec_min = 1
first_layer_print_speed_mm_per_sec_max = 80

for i=0,63,1 do
  _G['filament_diameter_mm_'..i] = 1.75
  _G['filament_priming_mm_'..i] = 4.0
  _G['extruder_temp_degree_c_' ..i] = 210
  _G['extruder_temp_degree_c_'..i..'_min'] = 150
  _G['extruder_temp_degree_c_'..i..'_max'] = 270
  _G['extruder_mix_count_'..i] = 1
end

add_checkbox_setting('gcode_ultimaker2','ultimaker2')
add_checkbox_setting('gcode_delta','delta')
