# Fitting module readme

Yorck Ramachers (Warwick)
Last updated April 17, 2019

The Fitter module is a SuperNEMO reconstruction module. It attempts to
fit clustered tracker hits and fill the TTD data bank in Falaise in the current
event data model.


## Files:

- Fitter_module.cpp
- Fitter_module.h
- CMakeLists.txt
- fit.conf


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

The build will create the `libFitter.so` shared library. Assuming that you have an `input.brio` file that contains
a `TCD` bank from the reconstruction, this can be run after editing
the configuration file to point at the library location.:

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

This reconstruction module attempts to use the bare required minimum to run in
the Falaise pipeline and use the existing data model. None of the pre-defined
'base' classes are used other than the dpp::base_module which renders this
code into a Falaise module.

Input is the calibrated data bank, default name 'CD'. Output is put into the
default 'TCD' data bank, the tracker cluster data bank. 

This module requires only ROOT as external library and C++11 standard 
libraries for the methods it uses.

This module uses the Catch test framework for Falaise
modules. All unit tests address the algorithm only, not the Falaise module
structure nor the embedding of this module into flreconstruct.

## Fitting process:


## Algorithm 


## Utilities

