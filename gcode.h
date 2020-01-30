// SL 2018-07-03
#pragma once

#include <algorithm>
#include <LibSL.h>
#include <LibSL/CppHelpers/BasicParser.h>

// start interpreting the gcode
void gcode_start(const char *gcode);

v3f  gcode_offset();

// advances to the next position
// return false if none exists (end of gcode)
bool gcode_advance();

// returns the next position to reach (x,y,z,e)
v4f  gcode_next_pos();

// returns the speed in mm/sec
float gcode_speed();

// return current line in gcode stream
int gcode_line();

// restart from scratch
void gcode_reset();

// returns true in a reading error occured
bool gcode_error();

// returns true if has to prime
bool gcode_do_prime();

// returns true if has to retract
bool gcode_do_retract();

// returns true if ironing is active
bool gcode_is_ironing();

// returns the layer thickness
float gcode_layer_thickness();
