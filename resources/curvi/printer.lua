version = 2

function comment(text)
  output('; ' .. text)
end

extruder_e = 0
extruder_e_restart = 0

function prep_extruder(extruder)
end

function header()
  output('o X ' .. gcode_to_model_x .. ' Y ' .. gcode_to_model_y .. ' Z ' .. gcode_to_model_z)
  output('t ' .. z_layer_height_mm)
end

function footer()
end

function layer_start(zheight)
  output(';(<layer ' .. layer_id .. '>)')
  if layer_id == 0 then
    output('G0 F600 Z' .. ff(zheight))
  else
    output('G0 F100 Z' .. ff(zheight))
  end
  -- output('G1 Z' .. f(zheight))
end

function layer_stop()
  comment('</layer>')
end

function retract(extruder,e)
  -- speed = priming_mm_per_sec * 60;
  -- letter = 'E'
  -- output('G1 F' .. speed .. ' ' .. letter .. f(e - len - extruder_e_restart))
  output('R')
  return e
end

function prime(extruder,e)
  -- speed = priming_mm_per_sec * 60;
  -- letter = 'E'
  -- output('G1 F' .. speed .. ' ' .. letter .. f(e + len - extruder_e_restart))
  output('P')
  return e
end

current_extruder = 0
current_frate = 0

function select_extruder(extruder)
end

function swap_extruder(from,to,x,y,z)
end

function move_xyz(x,y,z)
  output('G1 X' .. f(x) .. ' Y' .. f(y) .. ' Z' .. f(z+z_offset))
end

function move_xyze(x,y,z,e)
  to_mm_cube = 1.0
  if gcode_ultimaker2 then
    r = filament_diameter_mm[extruders[0]] / 2
    to_mm_cube = 3.14159 * r * r
  end
  if traveling == 1 then
    traveling = 0 -- start path
    if      path_is_perimeter then output(';perimeter')
    elseif  path_is_shell     then output(';shell')
    elseif  path_is_infill    then output(';infill')
    elseif  path_is_raft      then output(';raft')
    elseif  path_is_brim      then output(';brim')
    elseif  path_is_shield    then output(';shield')
    elseif  path_is_support   then output(';support')
    elseif  path_is_tower     then output(';tower')
    elseif  path_is_ironing   then output(';ironing')
    end
  end

  extruder_e = e
  letter = 'E'
  instr = 'G'
  if path_is_ironing then
    instr = 'I'
  end
  if z == current_z then
    output(instr .. '1 F' .. f(current_frate) .. ' X' .. f(x) .. ' Y' .. f(y) .. ' ' .. letter .. ff((e-extruder_e_restart)*to_mm_cube))
  else
    output(instr .. '1 F' .. f(current_frate) .. ' X' .. f(x) .. ' Y' .. f(y) .. ' Z' .. ff(z) .. ' ' .. letter .. ff((e-extruder_e_restart)*to_mm_cube))
    current_z = z
  end
end

function move_e(e)
  to_mm_cube = 1.0
  if gcode_ultimaker2 then
    r = filament_diameter_mm[extruders[0]] / 2
    to_mm_cube = 3.14159 * r * r
  end
  extruder_e = e
  letter = 'E'
  output('G1 ' .. letter .. f((e - extruder_e_restart)*to_mm_cube))
end

function set_feedrate(feedrate)
  current_frate = feedrate
end

function extruder_start()
end

function extruder_stop()
end

function progress(percent)
end

function set_extruder_temperature(extruder,temperature)
end

function set_fan_speed(speed)
end
