Multi-Physics Analysis Program (MPAP)

Fluid-Structure Interaction analysis using the CutFEM approach on hierarchical b-spline grids.

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

