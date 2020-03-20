gurobi=0

volumic=0
nozzle=0.4
layer=0.3
filament=1.75
ironing=0
model=


arg="none"

for elem in $*
do
  if [ $arg = "none" ]
  then
    arg=$elem
  else
    eval "$arg=$elem"
    arg="none"
  fi
done
  
if [ $arg = "none" ]
then
  echo Error in arguments
  exit
fi

model=${arg%.*}

echo Generate tetmesh "from $model.stl" ...
./toTetmesh.sh $model
echo Done!

echo Optimize...

if [ $gurobi = "1" ]
then
  ./bin/curvislice_grb $model.msh -l $layer
else
  ./bin/curvislice_osqp $model.msh -l $layer
fi
echo Done!

echo Prepare lua for IceSL
./luaGenerator.sh $model $volumic $nozzle $layer $filament $ironing

if [ -e "~/.icesl/icesl-printers/fff/curvi" ]
then
  echo "'curvi' printer profile already exist"
else
  echo Create 'curvi' printer profile for IceSL
  cp -r "./resources/curvi" "~/icesl/icesl-printers/fff/curvi"
fi

./tools/icesl/bin/icesl-slicer settings.lua --service

./bin/uncurve -l $layer --gcode $model

clear

echo '  ______                                  __            __  __                     '
echo ' /      \                                /  |          /  |/  |                    '
echo '/$$$$$$  | __    __   ______   __     __ $$/   _______ $$ |$$/   _______   ______  '
echo '$$ |  $$/ /  |  /  | /      \ /  \   /  |/  | /       |$$ |/  | /       | /      \ '
echo '$$ |      $$ |  $$ |/$$$$$$  |$$  \ /$$/ $$ |/$$$$$$$/ $$ |$$ |/$$$$$$$/ /$$$$$$  |'
echo '$$ |   __ $$ |  $$ |$$ |  $$/  $$  /$$/  $$ |$$      \ $$ |$$ |$$ |      $$    $$ |'
echo '$$ \__/  |$$ \__$$ |$$ |        $$ $$/   $$ | $$$$$$  |$$ |$$ |$$ \_____ $$$$$$$$/ '
echo '$$    $$/ $$    $$/ $$ |         $$$/    $$ |/     $$/ $$ |$$ |$$       |$$       |'
echo ' $$$$$$/   $$$$$$/  $$/           $/     $$/ $$$$$$$/  $$/ $$/  $$$$$$$/  $$$$$$$/ '

echo '==================================================================================='
echo '==>'
echo "     Gcode generated at: $model.gcode"
echo '==>'
echo '===================================================================================