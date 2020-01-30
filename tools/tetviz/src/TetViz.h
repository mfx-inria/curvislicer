#pragma once

#include <LibSL/LibSL.h>
#include <LibSL/LibSL_gl4.h>

#include "imgui\imgui.h"
#define IMGUI_DEFINE_MATH_OPERATORS true
#include "imgui\imgui_internal.h"

class TetViz{
private:
  TetViz & operator= (const TetViz&) {};
  TetViz(const TetViz&) {};

  static TetViz *s_instance;
  TetViz();
  ~TetViz() {}

public:
  static void TetViz::launch();

  static TetViz *Instance() {
    if (!s_instance)
      s_instance = new TetViz;
    return s_instance;
  }


protected:
  const float ZN = 1.0f;
  const float ZF = 2000.0f;

public:
  /**
  *  Called at each frame
  */
  static void mainRender();

  /**
  *  Called on resize
  *  @param width the new width
  *  @param height the new height
  */
  static void mainOnResize(uint _width, uint _height);


  /**
  *  Called when a key is pressed
  *  @param k the typed character
  */
  static void mainKeyPressed(uchar _k);
  static void mainScanCodePressed(uint _sc);
  static void mainScanCodeUnpressed(uint _sc);
  static void mainMouseMoved(uint _x, uint _y);
  static void mainMousePressed(uint _x, uint _y, uint _button, uint _flags);





  bool draw();

private:
  ImVec2 m_size;
  Trackball m_Trackball;
};