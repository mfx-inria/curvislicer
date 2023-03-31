#!/bin/bash

mkdir BUILD
cd BUILD

/mingw64/bin/cmake -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles" ..
mingw32-make install

cd ..
