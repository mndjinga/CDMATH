CDMATH
======

CDMATH is a geometrical and numerical toolbox designed for numerical analysts who work on the discretisation of partial differential equations on general shapes and meshes and who would rather focus on high-level scripting. The library originates from [CDMATH](http://cdmath.jimdo.com), a collaborative workgroup with the same name. It is based on the [MEDcoupling](http://docs.salome-platform.org/latest/dev/MEDCoupling/index.html) library of the [SALOME](http://www.salome-platform.org/) project for the handling of meshes and fields, and on the library [PETSC](https://www.mcs.anl.gov/petsc/) for the handling of matrices and linear solvers. The library currently developed for linux distributions and is maintained on Ubuntu 16.04 LTS, as well as on Fedora 24, 25 and 26.


Download CDMATH sources to compile
----------------------------------

Create your source directory. For instance:
* `mkdir ~/workspace/cdmath`
* `cd ~/workspace/cdmath`

Download from GitHub. For instance:
* `wget https://github.com/mndjinga/CDMATH/archive/master.zip`

Then unzip the file to a directory cdmath-master
* `unzip master.zip`


Set environment for the compilation of CDMATH
---------------------------------------------
Dependencies. The following packages list is sufficient on Ubuntu 14.04, Ubuntu 16.04 :

 - `cmake` (mandatory)
 - `g++` or another C++ compiler (mandatory)
 - `libhdf5-dev` (mandatory)
 - `python-dev`, `python-numpy` and `swig`, if you want to generate Python executables and libraries of CDMATH (highly recommended). Use the compilation option `-DCDMATH_WITH_PYTHON=ON`.
 - `python-matplotlib` and `paraview` for postprocessing tools such as plotting curves (matplotlib) or generating 3D views (paraview). Use the compilation option `-DCDMATH_WITH_POSTPRO=ON` (recommended).
 - `petsc` if you want to solve large spase linear system. Typically required for implicit methods (recommended).
 - `jupyter`, in order to generate nice reports from test case simulations
 - `doxygen`, `graphviz` and `mscgen`, if you want to generate a nice code source documentation in `~/workspace/cdmath/cdmath_install/doc/` (recommended). Use the compilation option `-DCDMATH_WITH_DOCUMENTATION=ON`.
 - `libcppunit-dev`, if you want to generate unit tests. Use the compilation option `-DCDMATH_WITH_TESTS=ON` (optional).
 - `libopenmpi-dev`, in particular if you need to use the compilation option `-DMEDFILE_USE_MPI=ON` (optional).
 - `rpm`, if you want to generate RPM installation packages. Use the compilation option `-DCDMATH_WITH_PACKAGE=ON` (optional).

Directories. Create the suggested build and installation folders:
* `cd ~/workspace/cdmath`
* `mkdir cdmath_build`
* `mkdir cdmath_install`
* `cd cdmath_build`


Compile and install CDMATH
--------------------------
Generate makefiles for a minimum version:
* `cmake ../cdmath-master/ -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release`

Or generate makefiles for an all-options version:
* `cmake ../cdmath-master -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -G"Eclipse CDT4 - Unix Makefiles" -D_ECLIPSE_VERSION=4.3 -DMEDFILE_USE_MPI=ON -DCDMATH_WITH_PETSC=ON -DCDMATH_WITH_PYTHON=ON  -DCDMATH_WITH_POSTPRO=ON -DCDMATH_WITH_TESTS=ON -DCDMATH_WITH_DOCUMENTATION=ON -DCDMATH_WITH_PACKAGE=ON`

Compile and install:
* `make`
* `make install`

Run unit and example tests:
* make check

Run validation tests:
* make validation

Notes for compilation options:
* Eclipse: The Cmake options `-G"Eclipse CDT4 - Unix Makefiles" -D_ECLIPSE_VERSION=4.3` create project files if you want to develop CDMATH with Eclipse Kepler or higher.
* HDF5: On some systems (not Ubuntu 14.04 nor Ubuntu 16.04), you may have to use the compilation option `-DHDF5_ROOT_DIR=/path/to/hdf5/library` too.
* MPI: On some systems (not Ubuntu 14.04, nor Ubuntu 16.04), you may have to use the compilation option `-DMPI_ROOT_DIR=/path/to/mpi/library` too. You may also have to set the environment variable `export MPI_ROOT_DIR=/path/to/mpi/library`. Moreover, on some systems (not Ubuntu 14.04, nor Ubuntu 16.04), the compilation option `-DMEDFILE_USE_MPI=ON` may be mandatory and be set to `ON`.
* PETSc: If the library Petsc is already installed in your system (packages libpetsc-dev for ubuntu and petsc-devel for fedora 25 and 26), you may save time and disk space by using the installed library instead of installing a new one. In order to do so use the compilation options `-DPETSC_DIR=/path/to/petsc/installation/petsc -DPETSC_ARCH=arch-linux2-c-opt`. If you prefer to compile PETSc yourself from the sources you may follow the instructions given in [the official documentation](http://www.mcs.anl.gov/petsc/documentation/installation.html).


Use CDMATH
----------
To use CDMATH with your C++ code `main.cxx`:
 * C++ libraries: `export LD_LIBRARY_PATH=~/workspace/cdmath/cdmath_install/lib`
 * To know how to include the right libraries for compilation, see the makefiles of the examples. They include the list `-linterpkernel -lmedC -lmedloader -lmedcoupling -lbase -lmesh -llinearsolver`.

To use CDMATH with your Python code `main.py`, you can load the CDMATH environment in your terminal using the command
 * source `~/workspace/cdmath/cdmath_install/env_CDMATH.sh`

The CDMATH environment variables consist in :
 * C++ libraries: `export LD_LIBRARY_PATH=~/workspace/cdmath/cdmath_install/lib`
 * Python libraries: `export PYTHONPATH=~/workspace/cdmath/cdmath_install/lib/cdmath:~/workspace/cdmath/cdmath_install/bin/cdmath`

Create Linux installation packages for CDMATH
---------------------------------------------
After popular request, here is how you can create packages for Ubuntu 14.04 and Ubuntu 16.04 and Red Hat-based Linux distributions:

1. Download CDMATH as explained hereabove.
2. Set the environment as explained hereabove (in particular, make sure you have `rpm` installed).
3. Generate a makefile with `cmake -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -DCDMATH_WITH_PACKAGE=ON ../cdmath-master/` and eventually other options (documentation, tests, swig, etc).
4. Compile with `make package`.

You will then find a Debian package in the build directory; you may install it on Ubuntu 14.04 or Ubuntu 16.04. You will also find an RPM package, which you may install on Red Hat-based distributions. This way, the packages you generate may include all the compilation options you want.

Unfortunately, the Debian package may be said to be of “bad quality” for Debian standards as far as ownership is concerned. This is true and due to limitations in CMake/CPack. The package should still install nonetheless.
