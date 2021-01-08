# Swig-SoFiA 2

This version is used to create a dynamic shared library \_sofia.so\ that can be called from Python. It allows data sources to be memory arrays, rather than reading in a FITS file.

Ensure you have the dependency WCSLIB mentioned below, as well as SWIG installed. You will also need numpy and in particular the numpy swig extensions:

    pip3 install numpy

This will install numpy in your ~/.local/lib/python3/site-packages directory. For me, using pip did not install the swig extensions successfully, so I downloaded numpy-1.19.4, found the numpy.i and pyfragments.swg files (under tools/swig) and copied them across to my build directory.

Now compile as follows:

1. "make clean"
2. "make"

Read the comments in the example wrapper file "python_spark.py" to see how to run it.

To compile the stand-alone version, ignore the above, just run "./compile.sh -fopenmp"

# SoFiA 2

This is version 2 of the HI Source Finding Application (SoFiA). SoFiA 2 is a reimplementation of the original SoFiA pipeline in the C programming language. It is intended for use in HI data analysis pipelines and will be developed and maintained in parallel to SoFiA 1.x. While SoFiA 2 is still under active development, a **stable release** is already available at this point in time and can be used in production mode. In addition, users will be able to continue using **SoFiA 1.x** (https://github.com/SoFiA-Admin/SoFiA) for processing their data.

## Improvements in SoFiA 2

* Being written in C and making extensive use of multi-threading, SoFiA 2 is much faster than SoFiA 1.x.
* SoFiA 2 requires significantly less memory than SoFiA 1.x (down from > 5 × cube size to ~ 2.3 × cube size).
* SoFiA 2 currently has only a single external dependency (WCSLIB) and should therefore compile and run on any machine with a Linux or Unix operating system and the GCC compiler and WCSLIB installed.

## Installation

Please ensure that the **GNU C compiler** (GCC) and **WCSLIB** (https://www.atnf.csiro.au/people/mcalabre/WCS/) are installed on your machine, before downloading and extracting the SoFiA 2 **source code** into a directory of your choice. You may want to first check if WCSLIB is either already installed or available from your operating system’s software repository (`wcslib-dev` package on Ubuntu/Debian; `wcslib-devel` on Fedora/CentOS/Red Hat; `wcslib` on MacOS/Homebrew) before attempting to install it manually.

Once WCSLIB is installed, simply execute the `compile.sh` script to compile SoFiA 2 using the GCC compiler:

`./compile.sh -fopenmp`

Note that the `-fopenmp` parameter is optional and will enable **multi-threading** using OpenMP. If your compiler does not support OpenMP, this parameter can simply be omitted to install a single-threaded version of SoFiA 2. Please ensure that you read and follow the **instructions** printed at the end of the compilation process to finalise the installation. If a compiler error related to WCSLIB shows up, please ensure that WCSLIB is installed in a standard location where it can be found by the GCC. Then run the compilation script again to see if the error message has disappeared.

As an alternative, we provide a `Makefile` for those who prefer to use `make` to install SoFiA 2. The `Makefile` itself contains a few examples of how to invoke it with different compilers (with or without OpenMP support). Note that if you prefer to use `make`, you may still want to create a symbolic link or alias to the `sofia` executable file in the end to make SoFiA 2 easily accessible across your system.

## Documentation

An overview of all control parameters as well as a PDF copy of the SoFiA 2 User Manual can be found on the wiki at https://github.com/SoFiA-Admin/SoFiA-2/wiki. The wiki also contains a small test data cube and parameter file that users can download to verify their SoFiA 2 installation.

## Feedback

Should you decide to run SoFiA 2 on your own data cubes, we would welcome any feedback on how well SoFiA 2 works for you and whether any improvements could be made. If you have a GitHub account, you can directly create a new issue (https://github.com/SoFiA-Admin/SoFiA-2/issues/new) on GitHub for questions, feature requests or bug reports. Alternatively, please feel free to directly contact the project leader, Tobias Westmeier, via e-mail at `tobias.westmeier (at) uwa.edu.au` to provide feedback on your experience with SoFiA 2. Note that the main purpose of SoFiA 2 is to facilitate the processing of HI data from SKA precursor surveys, and we are unlikely to be able to accommodate requests for additional features beyond this primary scope.

## Copyright and licence

Copyright (C) 2020 Tobias Westmeier

SoFiA 2 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License  along with this program. If not, see http://www.gnu.org/licenses/.
