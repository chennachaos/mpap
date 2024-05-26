Multi-Physics Analysis Program (MPAP) for Fluid-Structure Interaction (FSI) using the CutFEM approach on hierarchical b-spline grids.

The finite element formulation is published in 3 papers in [CMAME](https://www.sciencedirect.com/journal/computer-methods-in-applied-mechanics-and-engineering) journal.
[Paper 1](https://www.sciencedirect.com/science/article/pii/S0045782516304893)
[Paper 2](https://www.sciencedirect.com/science/article/pii/S0045782516313706)
[Paper 3](https://www.sciencedirect.com/science/article/pii/S0045782518301026)


Animations of some of the simulations performed using this code are available on my [YouTube](https://www.youtube.com/playlist?list=PL9IBrbGcgPbK2frCdXdOvCH4YwHRhiVts) channel.


* Programming languages: **C++** and **Fortran**
* C++ standard: **C++14** (or above)
* Required third-party libraries:
  1. CMake
  2. Blas
  3. Lapack
  4. Boost
  5. MPI (OpenMPI, MPICH or Intel MPI. Your choice!)
  6. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  7. [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
  8. [PETSc](https://www.mcs.anl.gov/petsc/)
  9. [VTK](https://vtk.org/)
  10. [CGAL](https://www.cgal.org/)
  11. Metis and ParMetis
  12. SuperLU


## Compilation and building
1. Clone the repository or download the zip file and extract its contents.
2. Go to the directory of the repository in a terminal.
3. Create **build** and **bin** directories.
    * `mkdir build`
    * `mkdir bin`
4. Modify the CMake file accordingly.
    * Copy `CMakeLists-chennalaptop.txt` to a new file for your machine, say `CMakeLists-local.txt`.
    * Change the paths to the compilers, Eigen, CGAL, PETSc and VTK libraries.
    * Change the path in the `install` function.
    * Create a symbolic link to the local CMake file.
      * `ln -sf CMakeLists-local.txt CMakeLists.txt`
5. Enter the `build` directory.
    * `cd build`
6. Configure using the CMake file
    * `cmake ..`
7. Compile, build and install the executable `mpap`. This step will also copy the exe to the `bin` folder.
    * `make install`

## Execution
* Simulations are usually run from the `bin` folder.
* Copy `petsc_options.dat` file from the `project/sampleinputs` folder to the `bin` folder.
* For serial run: `./mpap  <pathtodirectory>  <inpfilename>` (White space between each entry)
* For parallel run: `mpirun -n <nprocs>./mpap  <pathtodirectory>  <inpfilename>`
* Example:
  * `./mpap  ../project/sampleinputs  IsquareFixed`
  * `mpirun -n 4 ./mpap  ../project/sampleinputs  IsquareFixed`

## Understanding and using output(s)
For each simulation the output contains two items:

1. Files for visualisation in ParaView. Generated in the `bin` folder.
    * One **.pvtu** file for each time step. Contains velocity, pressure and vorticity for the flow field. This file is an amalgamation of **.vtu** files for individual subdomains. You need to open the **.pvtu** files in ParaView for visualising the whole fluid domain.
    * One **.vtu** file for each solid, if the solids are allowed to move. Only done for flexible solids.

2. A file with the letter prefixed `T` to the name of the input file containing the data for forces and displacements. (This will be created in the same directory as the input file.)


## Citing MPAP
If you use MPAP in your research, please cite it as:

```
@misc{mpap,
  author = {Chennakesava Kadapa},
  title = {{MPAP}: A Fluid-Structure Interaction simulation framework based on CutFEM},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/chennachaos/mpap}},
}
```

