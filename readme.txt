Cosmo++ is a numerical library for cosmology. The functionality is described in detail in the paper arXiv:1312.4961. Any published material that uses Cosmo++ should cite that paper. The library implements many features that have been developed by other authors. It is the responsibility of the user to cite all of the relevant papers. Please see the paper arXiv:1312.4961 for some references. It is also the responsibility of the user to cite any papers that are required to be cited for using the publicly available libraries needed for Cosmo++ (see the COMPILATION section below).

|---------------|
| DOCUMENTATION |
|---------------|

The library is fully documented using doxygen. Please run “doxygen” in the main directory to generate the documentation in html format.

|-------------|
| COMPILATION |
|-------------|

Cosmo++ can be compiled by itself without any other libraries. However, some of the functionality will be missing. To get the full functionality of Cosmo++ please install the following libraries before compiling Cosmo++. The user can include only some of the libraries. Only those modules of Cosmo++ will be compiled for which all of the required libraries have been included.

1. CFITSIO: Please download from http://heasarc.gsfc.nasa.gov/fitsio/ then compile and install.
2. HEALPIX: Please download from http://healpix.jpl.nasa.gov then compile and install. The C and C++ libraries need to be included in the installation.
3. LAPACK++: This is the C++ version of Lapack. Please download from http://math.nist.gov/lapack++/ then compile and install.
4. CLASS: Please download from http://class-code.net then compile and install.
5. MINUIT: Please download from http://seal.web.cern.ch/seal/work-packages/mathlibs/minuit/ then compile and install.
6. MULTINEST: Please download from http://ccpforge.cse.rl.ac.uk/gf/project/multinest/ then compile and install.
7. PLANCK likelihood code and data files: Please download from http://pla.esac.esa.int/pla/aio/planckProducts.html then compile and install the likelihood code.
8. WMAP9 likelihood code and data files: Please download from http://lambda.gsfc.nasa.gov/product/map/current/likelihood_get.cfm then compile.

You need to create a file make.inc and give the references to the directories of the libraries described above. Some (or all) of the entries can be left empty. An example file is included, make_template.inc. If the HEALPIX library is included then you must also include CFITSIO. If the Planck likelihood is included, you need to specify the location of a directory containing all of the Planck likelihood files you are planning to use. Please define this in PLANCKDATADIR in make.inc before compiling the library. If the WMAP9 likelihood is included (this should be included using the WMAP9DIR flag) then you need to also specify the compiler used for it (currently supported gfortran or ifort) using the WMAP9COMPILER flag, and the library linking flags using WMAP9LIBFLAGS. The library flags should include whatever is needed for linking to the WMAP9 likelihood code (e.g. lapack, blas), as well as the WMAP9 library file (-lwmap9) which needs to be made when compiling the WMAP9 likelihood code.

In order to turn on OpenMP you need to specify the OMP_FLAG in make.inc (this is usually -fopenmp). If this is left empty then the library will be compiled without OpenMP.

In order to compile Cosmo++ with MPI you need to specify the MPI_COMP flag in make.inc (this should be the name of the MPI C++ compiler, e.g. mpic++). If this is left empty then the library will be compiled without MPI.

No changes to the Makefile should be required. Just run “make” in the main directory.

The gcc compiler is needed.

|---------|
| TESTING |
|---------|

The library has multiple unit tests that can be run after compiling. To get a list of tests call "./bin/test list" after compiling. The tests are divided into two categories - fast and slow. You can run all of the fast tests by calling "./bin/test fast" and all of the slow tests by calling "./bin/test slow". Each individual test can be run by calling "./bin/test test_name" where test_name is the name of the particular tests. All of the names can be obtained using "./bin/test list".

Calling "make test" will run all of the fast tests. It is recommended to do this after compiling the library to make sure everything is ok. All of the fast tests take only about one minute.

The slow tests are much slower and it is recommended to run them individually. If you are using some functionality of the library that has a slow test it is recommended to run that test before using the library. Some of the tests may take a few days. Using many MPI nodes may accelerate the slow tests.

|-------|
| USAGE |
|-------|

The header files are in the include directory, the compiled library will be in the lib directory. The examples will be compiled into executables in the bin directory. You need to point to the include directory from your own project (i.e. using the -I flag) and then link your code to the libcosmopp.a static library. You also need to point to the include directories for the libraries given above and link to their libraries as well.

|--------------------------|
| COMMENTS AND SUGGESTIONS |
|--------------------------|

Please e-mail any comments, suggestions, or bugs to Grigor Aslanyan, g.aslanyan@auckland.ac.nz.

I will appreciate any suggestions for any tools that you would like to see in Cosmo++. If you want to include your own code into Cosmo++, please let me know.
