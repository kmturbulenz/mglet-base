MGLET-base flow solver
======================

MGLET-base is an incompressible Navier-Stokes flow solver written in Fortran,
provided by [KM Turbulenz GmbH](https://km-turbulenz.de/).

Relationship to commercial MGLET
--------------------------------

The present MGLET-base is a subset of the commercial product MGLET developed
and sold by [KM Turbulenz GmbH](https://km-turbulenz.de/) since 2005. This
open source code base was established in 2023 to provide a framework for
contributions and interactions with academic user communities. The present
MGLET-base contains the following features:

- Incompressible flow solver
- Ghost cell immersed boundary method for DNS and low-Re wall-resolved LES
of stationary geometries

The commercial MGLET contains the following additional features and
advantages:

- Unique cut-cell immersed boundary method for high-Re and wall-modeled
simulations, with the possibility of doing moving geometries
- Aeroacoustic solver capable of predicting flow noise, with a huge variety of
options on boundary conditions, loudspeakers, damping materials, acoustic
feedback and various acoustic evaluation tools
- Structural solver for acoustic-fluid-structure-interaction
- Advanced and efficient grid-generator with many refinement options such as
refinement boxes, surfaces and volumes, giving large flexibility in grid shapes
and setups. With this tool, grids are generated in seconds.
- A comprehensive set of integrated run-time sampling and postprocessing
options, including surface loads, wall-shear stresses, Fourier transforms
of 3-D volumetric fields, sound sources and more
- Postprocessing options with integration and output to Paraview, generation of
volumetric- and surface-fields, automatic image generation, spectra generation
and plotting
- A comprehensive report tool to generate automatic reports from a simulation
- Workflow management solutions that streamline the configuration, execution and
post-processing of a simulation
- A finished distribution package ready to use on any HPC system including
comprehensive documentation
- User training and support

If you need a turnkey solution for flow simulations, KM Turbulenz GmbH
provides MGLET as a ready-to-use commercial flow and aeroacoustic
solver including everything listed above and more. Please
[contact us](https://km-turbulenz.de/about/) for a discussion about the
possibilities MGLET will give you.


Purpose of MGLET-base
---------------------

The present code is a new and modern re-implementation of the "traditional"
and well-known MGLET numerics that date back to the 1980's. The code is
re-implemented Fortran 2018 making use of modern programming techniques to
improve flexibility and expandability.

The purpose of the MGLET-base is to provide high-quality building blocks for
academic work in flow simulations. The code is built to be easily expandable
for specific projects and needs. Implementing new physical models for any
transport phenomena should be easy. Although the present code compiles and runs
certain basic flow cases, it is not intended as the final tool for any
specific task or need. Instead, it serves as a starting point for further
developments of specific features.


### Brief history

- 1980's: Initial version developed by Heinrich Werner within the scope of his
dissertation "Grobstruktursimulation der turbulenten Strömung über eine
querliegende Rippe in einem Plattenkanal bei hoher Reynoldszahl". Performing
high-Reynolds number large-eddy simulations of turbulent channel flows with
obstacles on Cray X-MP and Y-MP vector computers.

- 1990's: Michael Manhart introduced message passing with the MPI library and
the multi-grid approach, including local grid refinement. Since this time the
code was named **M**ulti-**G**rid-**L**arge-**E**ddy-**T**urbulence.

- 2000's: Frederic Tremblay, Nikolaus Peller and Johannes Kreuzinger: Ghost-cell
immersed boundary implementation.

- 2005: [KM Turbulenz GmbH](https://km-turbulenz.de/) is founded and acquired
necessary rights to the MGLET code from TU Munich.

- 2005 -> present: KM Turbulenz develops MGLET into a commercial state-of-the
art flow and aeroacoustics solver

- 2023: MGLET-base open source core released


Numerical properties
-------------------

MGLET uses staggered, Cartesian grids organized in a overlapping, hierarchical
multi-level pattern allowing local grid refinement.

The staggered arrangement of the flow variables allows for energy
conserving spatial schemes. Coupled with explicit Runge-Kutta time integration
methods it is tailored for DNS and LES simulations of flows at a wide range
of Reynolds numbers.


Performance and scaling
-----------------------

MGLET's handling of many independent Cartesian grids is very efficient, and
simulations have been conducted with over 250 000 grids of size 24x24x24 - that
is nearly 3.5 billion grid cells. It has also been proven to run on tens of
thousands of MPI ranks in big supercomputers throughout Europe.


Building MGLET-base
-------------------

The most important build requirement is a set with C, C++ and Fortran
compiler with Fortran 2008 + TS 29113 support. Then you will need an MPI
library and a HDF5 library with MPI support.

The code is currently tested with GNU compilers version >= 11.2 and Intel
compilers from the oneAPI toolkits version >= 2022.2.1. Older compilers
than these will typically not work.

Since version 7.2 build 7213, the
[NAG Fortran compiler](https://nag.com/fortran-compiler/) also compiles a
working binary that runs all testcases successfully. This is only tested on
x64 Linux systems.

### Dependencies

* An MPI implementation of your choice that provides the `MPI_f08` Fortran
bindings
* HDF5: https://github.com/HDFGroup/hdf5
* CMake: https://cmake.org/download/

The following dependencies are fetched and built automatically by CMake:

* Nlohman JSON: https://github.com/nlohmann/json
* Exprtk: https://github.com/ArashPartow/exprtk

### Build instructions

MGLET-base make use of
[CMake presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)
to store a default set of build settings for the most common environments.
There are currently four pre-defined presets:

* `gnu-debug`: GNU compilers `gcc`, `g++` and `gfortran` with common debug flags
* `gnu-release`: GNU compilers `gcc`, `g++` and `gfortran` with flags
for release (performance-optimized) builds
* `intel-debug`: Intel compilers `icx`, `icpx` and `ifx` with common debug flags
* `intel-release`: Intel compilers `icx`, `icpx` and `ifx` with flags for
release (performance-optimized) builds
* `nag-debug`: NAG Fortran compiler `nagfor` with GNU C and C++ compilers `gcc`
and `g++` with debug flags. This configuration is very effective of discovering
errors that the Intel and GNU compilers cannot detect.
* `nag-release`: NAG Fortran compiler `nagfor` with GNU C and C++ compilers
`gcc` and `g++` with flags for release (performance-optimized) builds.

In order to build MGLET there is a few simple steps to follow:

1. Check out the source code

2. Create a separate `build` directory in the source code root and enter that
directory: `mkdir build && cd build`

3. Run CMake: `cmake --preset=gnu-release ..` (replace the preset with your
desired one)

4. Compile: `make`

5. Run tests: `ctest --output-on-failure --test-dir tests`


License and trademarks
----------------------

MGLET is a registered trademark of Kreuzinger und Manhart Turbulenz GmbH.

The present code is licensed under the Apache License, Version 2.0
(the "License"). A copy of the license is contained in this repository or at:

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
