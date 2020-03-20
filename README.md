# CurviSlicer

CurviSlicer is research project about achieving curved printing on standard, off-the-shelf, 3-axis FDM printers. 
Our goal is to improve surface quality and accuracy by curving the layers not only at the top, but also throughout the part. This reduces internal porosity and fragilities, and allows to accurately position the top curved surfaces. 

The image below compares adaptive slicing with flat layers (top) to the same number of layers using CurviSlicer (bottom). The adaptive slicer concentrates the thin slices around the car hood (as it should) but then it has to use thick slices everywhere else. Instead, CurviSlicer outputs curved slices that nicely follow the car outlines. These are roughly the same print time.

![](https://github.com/mfx-inria/curvislicer/blob/master/resources/car.png "Comparing flat and curved layers.")

The method was developed by an international team of researchers: Jimmy Etienne, Nicolas Ray, Daniele Panozzo, Samuel Hornus, Charlie C.L. Wang, Jonàs Martínez, Sara Mcmains, Marc Alexa, Brian Wyvill and Sylvain Lefebvre ; you can [find the academic paper here](https://hal.archives-ouvertes.fr/hal-02120033/document).
The work finds its origin in a brainstorm session during the 2018 Computational Geometry workshop at the Bellairs Research Institute, co-organized by Sue Whitesides and Sylvain Lazard.

The implementation was done by Jimmy Etienne and Sylvain Lefebvre, with guidance from colleagues. Adrien Bedel helped greatly to modify the code to support OSQP.
Please don't expect high quality, production ready code, this is a research prototype. The code depends on many other great projects such as [TetWild](https://github.com/Yixin-Hu/TetWild) and [OSQP](https://github.com/oxfordcontrol/osqp).

## Important note

The original implementation in the paper uses the [Gurobi](https://www.gurobi.com/) commercial solver. This initial implementation is in the SIGGRAPH 2019 branch. **Please use it for reproducibility of the paper results (speed and quality)**. The master branch is modified to use OSQP and while it works great, there are differences and limitations compared to the Gurobi version. We'll keep improving it!

# How to use (Windows)

This repository is meant to be built from source, and includes Windows binaries of some required external tools. Sources are meant to be compiled with Visual Studio C++ 2019. 

We will also provide a binary release package, so check the available files there.

## Prerequisites:

You need to have Visual Studio C++ and CMake latest installed.

### IceSL

Install the latest version of [IceSL](https://icesl.loria.fr/download/) (adds a small feature to work with curvislice).

Once installed, copy the folder "curvi" (in the /resources folder) into IceSL printer profiles folder ; on Windows this is **%appdata%/IceSL/icesl-printers/fff/**

## Download

```git clone --recurse-submodules https://github.com/mfx-inria/curvislicer.git```

This will automatically download other repositories:
	SolverWrapper (wrapper API around Gurobi and OSQP),
	OSQP,
	LibSL-small.

## Build

Then, you need to build the **INSTALL** project, it will generate the executables and put them in the **bin** folder.

By default, the OSQP solver version will be built. If you want to use Gurobi instead, you'll have to enable the CMake flag "BUILD_WITH_GRB" and choose the "GRB_VERSION" (and quite obviously you need to have Gurobi installed with a license).

## Run

From a Windows command line run:

```curvislice.bat <volumic=0> <nozzle=0.4> <layer=0.3> <filament=1.75> <ironing=0> [stl_filename]```

It will automagically generate your gcode files.

For example, a great starting point is to simply run

```curvislice.bat models/wing.stl```

The gcode is then found in models/wing.gcode

# How to use (Linux)

There is not reason this would not work under Linux, but we did not have time yet to make the scripts and the build system for all dependencies. Contributions are welcome!

You can follow the Windows procedure, but will have to manually compile dependencies (TetWild) and create shell scripts from the Windows batch files.

# Printing

The produced GCode is standard Marlin style for 1.75 mm filament and 0.4 mm nozzle. 
Note that it has **no header and no footer**. These you will have to add manually to fit your printer. Also please make sur the produced GCode properly fits your bed as we use an 'average' print bed configuration.

In our experience the GCode prints best on delta-style printers, as the Z axis is comparably efficient to the X,Y axes. On other types of printers some adaptation of flow is required ; our tool **uncurve** has some command line parameters for this purpose, but these are mostly experimental.

# Caution, this software generates complex curved trajectories that may result in collisions between the printer carriage and the print. This could damage your printer.

*We are expecting a certain clearance around the nozzle, so make sure there is space around -- basically a 45 degree cone going up from the nozzle tip on at least 5 centimeters, but larger parts may require more clearance. The angle to optimize for can be controlled from the command line.*

![](https://github.com/mfx-inria/curvislicer/blob/master/resources/nozzle-clearance.jpg "Typical space required around the nozzle.")

# Slicing parameters

We slice with IceSL using default parameters that may not be best for your models.
You can change these parameters by opening IceSL, selecting the *curvi* printer, 
changing parameters and slicing any object (slicing will save the settings for next time).

We encourage you to play with ironing and the type of top covers (curved covers or zigzag covers).

# Integrating in another slicer

CurviSlicer was developed using IceSL, however it can easily be used with different slicers.
The optimizer generates a model that has to be sliced 'flat'. The model is called *after.stl*
and can be found in the sub-directory having the name of your model and created during processing.

However the slicer has to output a special GCode format, see *printer.lua* in the *curvi* printer profile directory (in /resources). Our tool *uncurve* also needs to know how the 3D mesh spatially relates to the trajectories produced by the slicer, as well as the layer thickness.
This requires two special lines at the top, here is an example:
```o X-10.3 Y-5.7 Z0.0
t 0.2```
Here, it means that given a trajectory point, we have to add the offset (-10.3,-5.7,0) to
locate this same point in the input 3D model space. The line starting with *t* gives the slicing layer height.

You are welcome to use CurviSlicer within the scope of the license (see below). Please cite our paper in your publications and the credits of your software!

```@article{curvislicer,
author = {Etienne, Jimmy and Ray, Nicolas and Panozzo, Daniele and Hornus, Samuel and Wang, Charlie C. L. and Mart\'{\i}nez, Jon\`{a}s and McMains, Sara and Alexa, Marc and Wyvill, Brian and Lefebvre, Sylvain},
title = {CurviSlicer: Slightly Curved Slicing for 3-Axis Printers},
year = {2019},
volume = {38},
number = {4},
journal = {ACM Transactions on Graphics},
articleno = {Article 81},
numpages = {11},
}```

### License

[Affero GPL 3.0](https://www.gnu.org/licenses/agpl-3.0.en.html)
