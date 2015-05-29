Cosmo++ is a numerical library for cosmology. The functionality is described in detail in the paper arXiv:1312.4961. Any published material that uses Cosmo++ should cite that paper. The library implements many features that have been developed by other authors. It is the responsibility of the user to cite all of the relevant papers. Please see the paper arXiv:1312.4961 for some references. It is also the responsibility of the user to cite any papers that are required to be cited for using the publicly available libraries needed for Cosmo++ (see the COMPILATION section below).

Some of the newer functionality of Cosmo++, such as the k-d tree, the learn-as-you-go emulator, and the accelerated Planck likelihood code are described in the paper arXiv:xxxx.xxxx. If you are using this functionality please make sure to cite that paper as well.

|---------------|
| DOCUMENTATION |
|---------------|

The library is fully documented using doxygen. Please run “doxygen” in the main directory to generate the documentation in html format.

|-------------|
| COMPILATION |
|-------------|

The latest version of Cosmo++ uses cmake for compilation. cmake is available on most unix/linux systems and on macos. If your system does not have cmake installed you can download it from http://www.cmake.org.

Cosmo++ can be compiled by following these steps:
1. Create a cmake_settings.txt file by copying from cmake_settings_template.txt. Edit the file to specify some compile settings (see below).
2. Create a directory where it will be built (for example, you can create a directory named 'build' inside the Cosmo++ main directory) and change into that directory.
3. Run 'cmake <cosmo++_dir>' where <cosmo++_dir> is the main directory of Cosmo++ (in the example above you would run 'cmake ..').
4. Run 'make'. You can also comile on multiple threads by running 'make -jN' with N replaced by the number of threads you want to use. For example, 'make -j4' will use 4 threads.
5. Run 'make test' to test your build (optional).
6. Run 'make install'

The comments inside cmake_settings_template.txt explain the different settings. We give a more detailed explanation below.

You can edit cmake_settings.txt to specify the following settings:
1. Installation directory: You need to specify the installation directory. You can do it by either specifying the variable CMAKE_INSTALL_PREFIX in cmake_settings.txt or by passing -DCMAKE_INSTALL_PREFIX=<DIR> to cmake (in step 3).
2. Build configuration: This is release by default, i.e. with optimisations turned on. If you would like to change the configuration to debug, set the variable CMAKE_BUILD_TYPE to DEBUG by uncommenting the corresponding line in cmake_settings.txt. Alternatively, you can pass -DCMAKE_BUILD_TYPE=DEBUG to cmake (step 3).

Cosmo++ can be compiled by itself without any other libraries. However, some of the functionality will be missing. To get the full functionality of Cosmo++ please install the following libraries before compiling Cosmo++. The user can include only some of the libraries. Only those modules of Cosmo++ will be compiled for which all of the required libraries have been included.

1. LAPACK: There are different versions of the library, and some systems may already include it. It can be downloaded from http://www.netlib.org/lapack/. If you would like to include lapack functionality in Cosmo++, set the variable LAPACK_LIB_FLAGS in cmake_settings.txt. These are the flags for linking the lapack library. A few different examples are given in the comments in cmake_settings_template.txt.
2. CFITSIO: Please download from http://heasarc.gsfc.nasa.gov/fitsio/ then compile and install. Then set the variable CFITSIO_DIR to the directory where cfitsio is installed.
3. HEALPIX: Please download from http://healpix.jpl.nasa.gov then compile and install. The C and C++ libraries need to be included in the installation. Set the variable HEALPIX_DIR to point to the directory where Healpix is installed, and HEALPIXPP_DIR to the directory where the Healpix C++ libraries are installed.
4. CLASS: Please download from http://class-code.net then compile and install. Set the variable CLASS_DIR to point to the directory where Class is installed.
5. MULTINEST: Please download from http://ccpforge.cse.rl.ac.uk/gf/project/multinest/ then compile and install. Set the variable MULTINEST_DIR to point to the directory where MultiNest is installed.
6. POLYCHORD: Please download from https://ccpforge.cse.rl.ac.uk/gf/project/polychord/ then compile and install. Set the variable POLYCHORD_DIR to point to the directory where PolyChord is installed. NOTE: The PolyChord support in Cosmo++ is in the experimental stage and is not well tested yet.
7. PLANCK likelihood code and data files: Please download from http://pla.esac.esa.int/pla/aio/planckProducts.html then compile and install the likelihood code. Set the variable PLANCK_DIR to point to the directory where Planck likelihood code is installed, PLANCK_DATA_DIR to point to the directory where the Planck likelihood files are, and the PLANCK_COMPILE_FLAGS variable to the specific compilation flags used in compiling the likelihood code.
8. WMAP9 likelihood code and data files: Please download from http://lambda.gsfc.nasa.gov/product/map/current/likelihood_get.cfm then compile. Set the variable WMAP9_DIR to point to the directory where it's installed.
9. MINUIT: Please download from http://seal.web.cern.ch/seal/work-packages/mathlibs/minuit/ then compile and install. Set the variable MINUIT_DIR to the directory where Minuit is installed.

You can specify any combination of the libraries, or just specify none. However, note that CFITSIO must be specified if HEALPIX and/or WMAP9 are specified. Also, LAPACK must be specified if WMAP9 is specified.

Cmake will choose the default C/C++ and Fortran compilers of your system. If you would like to choose different compilers you can do so through setting some environment variables. Namely, the variable CC should specify the C compiler, the variable CXX should specify the C++ compiler, and FC should be the Fortran compiler. For example, running 'export CXX=g++' will set your C++ compiler to g++. You need to specify these environment variables BEFORE running cmake (step 3).

No Fortran compiler is needed unless you're compiling Cosmo++ with PolyChord.

Cmake will automatically check if MPI libraries are installed on your system and will compile Cosmo++ with MPI if possible.

If your C++ compiler supports openmp, then Cosmo++ will be compiled with openmp.

To avoid linking problems with other libraries, it is recommended to compile the other libraries with the same compilers and compiling options as Cosmo++. For example, if MultNest is compiled with MPI and Cosmo++ is compiled without MPI then the build will fail. Also, if Healpix is compiled with openmp but the C++ compiler used for Cosmo++ does not support openmp then the build will fail.

|---------|
| TESTING |
|---------|

The library has multiple unit tests that can be run after compiling.

To get a list of tests call "./bin/cosmo_test list" after compiling. The tests are divided into two categories - fast and slow. You can run all of the fast tests by calling "./bin/cosmo_test fast" and all of the slow tests by calling "./bin/cosmo_test slow". Each individual test can be run by calling "./bin/cosmo_test test_name" where test_name is the name of the particular tests. All of the names can be obtained using "./bin/cosmo_test list".

Calling "make test" will run all of the fast tests. It is recommended to do this after compiling the library to make sure everything is ok. All of the fast tests take only a few minutes.

If a test fails and you would like to see a more detailed output, you can set the environment variable "export CTEST_OUTPUT_ON_FAILURE=1" before running "make test". This will make sure that the output is displayed for the tests that failed. You can also re-run the individual tests that failed as described above.

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
