# CurviSlicer

CurviSlicer is research project about achieving curved printing on standard, off-the-shelf, 3-axis FDM printers. 
Our goal is to improve surface quality and accuracy by curving the layers not only at the top, but also throughout the part. This avoids leaving porosity inside, and allows to accurately position the top curved surfaces. 

The image below compares adaptive slicing with flat layer (top) to the same number of curved layers using CurviSlicer (bottom). The adaptive slicer concentrates the thin slices around the car hood (as it should) but then it has to use thick slices everywhere else. The curved slices nicely follow the car outlines.

![](https://github.com/mfx-inria/curvislicer/blob/master/resources/car.png "Comparing flat and curved layers.")

The method was developed by an international team of researchers: Jimmy Etienne, Nicolas Ray, Daniele Panozzo, Samuel Hornus, Charlie C.L. Wang, Jonàs Martínez, Sara Mcmains, Marc Alexa, Brian Wyvill and Sylvain Lefebvre ; you can find the academic paper here: https://hal.archives-ouvertes.fr/hal-02120033/document
The work finds its origin in a brainstorm session during the 2018 Computational Geometry workshop at the Bellairs Research Institute, co-organized by Sue Whitesides and Sylvain Lazard.

The implementation was done by Jimmy Etienne and Sylvain Lefebvre, with guidance from colleagues. Adrien Bedel helped greatly to modify the code to support OSQP.
Please don't expect high quality, production ready code, this is a research prototype. The code depends on many other great projects such as [TetWild](https://github.com/Yixin-Hu/TetWild) and [OSQP](https://github.com/oxfordcontrol/osqp).

## Important note

The original implementation in the paper uses the [Gurobi](https://www.gurobi.com/) commercial solver. This initial implementation is in the SIGGRAPH 2019 branch, please us it for reproducibility of the paper results (speed and quality). The master branch is modified to use OSQP and while it works great, there are differences and limitations compared to the Gurobi version. Of course we'll keep improving it!

# How to use

This repository is meant to be built from source, and includes binaries (Windows) of some required external tools. It is primarily developed under Windows with Visual Studio Community 2019. There is not reason this would not work under Linux, but we did not have time yet to make the scripts, external binaires, etc. Contributions are welcome!

We will also provide a binary release package, so check the files there.

## Prerequisites:

You need to have Visual C++ Studio and CMake installed.

### IceSL

Install the latest version of [IceSL](https://icesl.loria.fr/download/) (the latest version adds a small feature to work with curvislice).

Once installed, copy the folder "curvi" (in the /resources folder) into IceSL printer profiles folder ; on Windows this is **%appdata%/IceSL/icesl-printers/fff/**

(under Linux it should be in **~/.icesl/icesl-printers/fff/**)

## Download, build and run

### Download
**git clone --recurse-submodules https://github.com/mfx-inria/curvislicer.git**

This will download other repositories as:
	SolverWrapper (wrapper API around Gurobi and OSQP).
	OSQP
	LibSL-small

### Build

By default, only OSQP solver will be used. If you want to use Gurobi instead, you'll have to enable the cmake flag "BUILD_WITH_GRB" and choose the "GRB_VERSION".

Then, you need to build the **INSTALL** project, it will generate the executables and put them in the **bin** folder

### Run

On Windows, from a command line run:

curvislice.bat <volumic=0> <nozzle=0.4> <layer=0.3> <filament=1.75> <ironing=0> [stl_filename]

It will automagically generate your gcode files.

For example, a great starting point is to simply run

curvislice.bat models/wing.stl

The gcode is then found in models/wing.gcode

Linux: you have to run every step by hand for now. A setting example for IceSL is provided (settings.lua).

### Printing

The produced gcode is standard Marlin style for 1.75 mm filament and 0.4 mm nozzle. In our experience it works best on delta-style printers, as the Z axis is comparably efficient to the X,Y axes. On other types of printers some adaptation of flow is required ; our tool **uncurve** has some command line parameters for this purpose, but these are mostly experimental.

# Caution, this software generates complex curved trajectories that may result in collisions between the printer carriage and the print. This could damage your printer.

*We are expecting a certain clearance around the nozzle, so make sure there is space around -- basically a 45 degree cone going up from the nozzle tip on at least 5 centimeters, but larger parts may require more clearance. The angle to optimize for can be controlled from the command line.*

![](https://github.com/mfx-inria/curvislicer/blob/master/resources/nozzle-clearance.jpg "Typical space required around the nozzle.")

### License

[Affero GPL 3.0](https://www.gnu.org/licenses/agpl-3.0.en.html)
