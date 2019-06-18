# Fitting module readme

Yorck Ramachers (Warwick)
Last updated June 18, 2019

The Fitter module is a SuperNEMO reconstruction module. It attempts to
fit clustered tracker hits and fill the TTD data bank in Falaise in the current
event data model. Unfortunately, that is not possible since the TTD data 
bank does not offer to hold neither required data on fit characteristics 
nor fit errors.

The fitter_module code is therefore left empty for now - this is not a working Falaise module yet.

## Files:

- fitter_module.cpp
- fitter_module.h
- fitter_library.cpp
- fitter_library.h
- CMakeLists.txt
- fit.conf.in


## Description

Add to an flreconstruct pipeline to cluster tracker hits for reconstruction. To build it, do

``` console
$ ls
CMakeLists.txt  fit.conf fitter_module.cpp  fitter_module.h  README.md  testing

$ mkdir build
$ cd build
$ cmake -DCMAKE_PREFIX_PATH=$(brew --prefix) ..
...
$ make
...
... If you are developing the module, you can test it by calling
$ make test
...
... or obtain more detail on the tests and launch in the build directory
$ ctest -V
```

Note: if you get a QT5 error, you may need to specify the QT5 path when you run the cmake line, as given by `brew --prefix qt5-base`. For example, you can run:
``` console
$ cmake -DCMAKE_PREFIX_PATH="$(brew --prefix qt5-base);$(brew --prefix)" ..
``` 

The build will create the `libFitter.so` shared library. Assuming that you have an `input.brio` file that contains a `TCD` bank (i.e. clustered data) from the reconstruction, this can be run after editing the configuration file to 
point at the library location. No module configuration parameter are required.:

``` console
...
[name="flreconstruct.plugins" type="flreconstruct::section"]
plugins : string[1] = "Fitter"
Fitter.directory : string = "/my/path/to/build/directory"

# Define the modules in the pipeline:
[name="pipeline" type="fitter_module"]
...
$ flreconstruct -i /path/to/input.brio -p fit.conf -o /tmp/fitted_data.brio
```

## Implementation

This reconstruction module attempts to use the bare required minimum to run in the Falaise pipeline and use the existing data model. 

Input is the clustered data bank, default name 'TCD'. Output is put into the default 'TTD' data bank, the tracker trajectory data bank. Additional required  data to write can not be put into TTD as it does not provide storage for 
fit errors or fit characteristics. Therefore, a different persistence data model has to be used.

This module requires only ROOT as external library and C++11 standard libraries for the methods it uses.

This module uses the Catch test framework for Falaise modules. All unit tests address the algorithm only, not the Falaise module structure nor the embedding of this module into flreconstruct.

## Fitting process:

Each cluster of Geiger Hits is fit with three expected models: a line model in 3D, a broken line mode in 3D and a helix model in 3D. The former requires four fit parameters, the latter five. The broken line model also features 4 
fit parameters but delivers additionally all relative angles between data points, i.e. a perfect line would have zero angles, showing no deviation from a perfect line model. The fitter object, SNFitter, has methods to calculate 
initial best guesses for the fit parameters of each model. 

Like for every high-dimensional fit problem, good initial values are crucial. Tests on idealized data shows that for fitting to rings (cylinders in 3D), each model can show multiple shallow minima in the parameter space and fits 
might not necessarily find the best (deepest) minimum even with fairly decent initial values. This turns fitting potentially into a hard problem, requiring subsequent evaluation of solutions.

The chosen minimizer is TMinuit2 from the ROOT framework as a tried and tested method. It offers an interface to using a functor (a fit function object as opposed to a function) for fitting. This allows for flexible handling of required data for the fitter. There are 
correspondingly fit objects for the line and the helix, respectively. All they are required to do is to calculate the total distance between all data objects (the GeigerRing) and model. The minimizer will then attempt to minimize 
that (weighted) distance. This is therefore a textbook least squares minimization problem. 

Consequently, best fit values for all parameters as well as their errors require storage. Likewise, the goodness-of-fit value characterising the process and any potential messages from the process require storage. The fit results 
need to be evaluated subsequently and these numbers are hence important. Likewise any fit extrapolation of a given model can only be done properly with fit errors known. Vertex determination needs to result in an extended area 
rather than an intersection point. This can only take place with model uncertainties known to the vertex extrapolator.

While the above models merely require Euclidian geometry to obtain the distances of cylinders to models, an extension to a broken line model (V. Blobel, NIM A566 (2006) 14) is more suitable for broken tracks due to multiple scattering. 

The broken line algorithm relies on being able to define separated planes in propagation direction. Data points (with propagation vector) must fix non-overlapping planes such that a comparison to neighbouring direction vectors 
becomes possible. That appears to be impossible for ring-data in the tracker where the suitable tangent points depend on the momentary model solution and non-overlapping planes can change order at each iteration. 

The rings therefore have to be replaced first with a collection of tangent points as data to fit a broken line. This fitting capability hence introduces a large section of code dedicated purely to tangent point calculations and 
path finding since not all tangent points are viable for a potential charged particle path. Every possible path should then be fit separately and sorted as a candidate charge path according to fit quality. Sorting only to shortest 
paths can yield surprising kinks since the origin of the potential path points are after all rings.

## Algorithm 

The distance calculations for the line in 3D follow a ROOT example from L. Moneta with the distance line to point as D= | (xp-x0) cross  u |, where xp is the data point and x0 an arbitrary point on the line, while u is the unit 
vector defining the line direction. Geometrically, this distance is the magnitude of the vector perpendicular to the line through the data point. Therefore, if the point xp is the centre of a circle (the ring) the distance to the 
circle circumference is merely the distance to the centre minus the ring radius.

This geometry calculation does not change much when swapping the line with a helix. Since tracker helices are nothing else than circles winding around a cylinder surface with axis fixed by the magnetic field direction, the real 
challenge is to calculate the shortest total distance between a model circle and a set of circles (the data). Again, the vector connecting centres defines all the important points. It passes through the data circle circumference 
and the helix (being the circumference of a circle too). The distance between those two points then contributes to the total sum of distances.

The broken line fit algorithm is documented in (V. Blobel, NIM A566 (2006) 14) and maps well to SuperNEMO data iff the data to fit consists of suitable path points made from tangent points on the collection of rings. The natural 
propagation direction for the algorithm to walk from point to point is the x-coordinate axis in both trackers. The line model for the broken line allows for vector to vector angle deviations at every step, each adding to the Chi 
square sum hence straight lines are favoured but not enforced. The collection of individual angles for a best fit line is stored and can be analysed once the best fit has been obtained. A Kink finder algorithm is included which 
searches for a significant angle above angle noise and signals, if found, a scatter event along the track. The longer the track, the more data there is, the more sensitive the kink finder gets. Should a kink be located, the track 
solution is split into at least two tracks, the nearest track stub to the foil and the calo. Each is then the best estimate for extrapolation to the foil and calo, respectively.

## Utilities

Initial values for the line model are swiftly obtained from a 2D line fit to Geiger ring centres using the TLinearFitter in ROOT. The problem with that approach is the multiple degeneracy of tracker trajectories. The idealized  
event data samples were fit successfully only if at least four different initial value sets of slopes and intercepts is offered to the full 3D fitter. These are obtained from allowing an arbitrary but fixed slope deviation to add 
(subtract) from the best fit slope and additionally shift the intercept solution by the nearest ring radius, again to both sides. At least one of the four fits converged in all test cases.

The inital values for the helix are surprisingly less critical given that the minimization challenge should be harder with five parameter. It appears that the helix parameter space simply offers far fewer localized minima for he 
fitter to get stuck in and convergence appears to be easier, regardless of the starting values. 

Analytic calculations assuming data points (here the Geiger ring centres) on a ring allow to calculate the helix radius and centre position in 2D. A linear fit to the z-coordinates gives the helix pitch and the mean value between 
the extreme z-coordinate positions gives the z-centre coordinate. These values (radius, pitch, xc, yc, zc) are one possible representation of the five helix parameters. They are often far from the true values but give the minimizer 
a good starting point. For the degenerate (purely horizontal track) case, they are useless. 

A backup function for inital values appears to work just as fine. It suggests that the y and z centre locations (y: how the helix curves, up or down, z: pitch direction - up or down) are more important than staying close to true values.

The broken line algorithm uses a lot of code to calculate tangent points on all rings in relation to each other and subsequently to find suitable paths from one end of the track to the other. While the tangent point calculations 
follow well documented Euclidean geometry formulae, the path finding required a few utility objects such as an Interval object that respects intervals on a circle rather than purely on a straight line. Another significant utility 
object is a weighted graph object that hosts the classic Dijkstra shortest path algorithm between a node source and a node target. This weighted graph is quite different to the unweighted graph used in the cluster library. Weights 
are defined as doubles (here Euclidean distances) while nodes remain integer identifiers with the same argument as explained in the cluster library.

## Data model

LineFit structure: 
- 4 best fit parameter (double), slopes and intercepts
- 4 error values (double)
- Chi2 (double)
- Probability (double)
- Fit status (integer)

HelixFit structure:
- 5 best fit parameter (double), radius, pitch, centre (x,y,z)
- 5 error values (double)
- Chi2 (double)
- Probability (double)
- Fit status (integer)

BrokenLineFit structure:
- 2 LineFit structures with the second potentially empty
- breakpoint collection (vector<int>); potentially empty
- angle collection (vector<double>)
- path point collection (vector<PathPoint>); see below
- path length (double)
- Chi2 (double)
- Probability (double)
- Fit status (integer)

PathPoint stucture:
- 3 tangent point coordinates (double)
- 3 errors in Cartesian coordinates (double)
- 1 ID pair (pair<int, size_t>)

