Cosmo.cpp is a numerical library for cosmology. The functionality is described in detail in the paper TBD. Any published material that uses cosmo.cpp should cite that paper. The library implements many features that have been developed by other authors. It is the responsibility of the user to cite all of the relevant papers. Please see the paper TBD for some references. It is also the responsibility of the user to cite any papers that are required to be cited for using the publicly available libraries needed for cosmo.cpp (see the COMPILATION section below).

|---------------|
| DOCUMENTATION |
|---------------|

The library is fully documented using doxygen. Please run “doxygen” in the main directory to generate the documentation in html format.

|-------------|
| COMPILATION |
|-------------|

Several publicly available libraries are required for the successful compilation of cosmo.cpp.

1. CFITSIO: Please download from http://heasarc.gsfc.nasa.gov/fitsio/ then compile and install.
2. HEALPIX: Please download from http://healpix.jpl.nasa.gov then compile and install. The C and C++ libraries need to be included in the installation.
3. BOOST: Please download from http://www.boost.org. No compilation or installation is required.
4. LAPACK++: This is the C++ version of Lapack. Please download from http://math.nist.gov/lapack++/ then compile and install.
5. CLASS: Please download from http://class-code.net then compile and install.

Optional libraries (you may choose to not include them in which case some of the functionality of cosmo.cpp will miss).

6. MINUIT: Please download from http://seal.web.cern.ch/seal/work-packages/mathlibs/minuit/ then compile and install.
7. MULTINEST: Please download from http://ccpforge.cse.rl.ac.uk/gf/project/multinest/ then compile and install.
8. PLANCK likelihood code and data files: Please download from http://pla.esac.esa.int/pla/aio/planckProducts.html then compile and install the likelihood code.

You need to create a file make.inc and give the references to the directories of the libraries described above. The optional libraries may be left empty. An example file is included, make_template.inc. No changes to the Makefile should be required. Just run “make” in the main directory. The gcc compiler is needed.

|-------|
| USAGE |
|-------|

The header files are in the include directory, the compiled library will be in the lib directory. The examples will be compiled into executables in the bin directory. You need to point to the include directory from your own project (i.e. using the -I flag) and then link your code to the libcosmocpp.a static library. You also need to point to the include directories for the libraries given above and link to their libraries as well.

|--------------------------|
| COMMENTS AND SUGGESTIONS |
|--------------------------|

Please e-mail any comments, suggestions, or bugs to Grigor Aslanyan, g.aslanyan@auckland.ac.nz.

I will appreciate any suggestions for any tools that you would like to see in cosmo.cpp. If you want to include your own code into cosmo.cpp, please let me know.
