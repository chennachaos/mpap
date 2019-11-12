Multi-Physics Analysis Program (MPAP)

Fluid-Structure Interaction analysis using the CutFEM approach on hierarchical b-spline grids.

The finite element formulation is published in 3 papers in [CMAME](https://www.sciencedirect.com/journal/computer-methods-in-applied-mechanics-and-engineering) journal.
[Paper 1](https://www.sciencedirect.com/science/article/pii/S0045782516304893)
[Paper 2](https://www.sciencedirect.com/science/article/pii/S0045782516313706)
[Paper 3](https://www.sciencedirect.com/science/article/pii/S0045782518301026)


Animations of some of the simulations performed using this code are available on my [YouTube](https://www.youtube.com/playlist?list=PL9IBrbGcgPbK2frCdXdOvCH4YwHRhiVts) channel.


* Programming languages: **C++** and **Fortran**
* C++ standard: **C++11** (or above)
* Required third-party libraries:
  1.) [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  2.) [PETSc](https://www.mcs.anl.gov/petsc/)
  3.) [VTK](https://vtk.org/)
  4.) [CGAL](https://www.cgal.org/)


## Compilation and building
* Create **build** and **bin** directories.
* Modify the CMakeLists.txt file accordingly.
  * Change the paths to the compilers, Eigen, CGAL, PETSc and VTK libraries.
  * Change the path to *bin*.
* Run the following commands (from the repository folder)
  * `mkdir build`
  * `cd build`
  * `cmake ..`
  * `make install`


## Execution
./mpap  \<pathtodirectory\>  \<inpfilename\>

