OMP_FLAG = -fopenmp
MPI_COMP = mpic++

COMPILER_FLAGS = -g -O4 -ffast-math -fpic -std=c++11

CFITSIO = /usr/local
HEALPIX = /Users/aslanyan/programs/Healpix_3.11
HEALPIXPP = $(HEALPIX)/src/cxx/osx
LAPACKPPLIBDIR = /usr/local/lib
LAPACKPPINCDIR = /Users/aslanyan/C++/library/lapackpp-2.5.4/include
CLASSDIR = /Users/aslanyan/programs/class_v2.3.2

MINUIT = /usr/local
MULTINEST = /Users/aslanyan/programs/MultiNest_v3.2_CMake/multinest
PLANCKDIR = /Users/aslanyan/programs/planck/COM_Code_Likelihood-v1.0_R1.10.ext/plc-1.0
PLANCKDATADIR = /Users/aslanyan/programs/planck
PLANCKCOMPILEFLAGS = -D LAPACK_CLIK
WMAP9DIR = /Users/aslanyan/programs/wmap_likelihood_v5
WMAP9COMPILER = gfortran # supported gfortran or ifort
WMAP9LIBFLAGS = -lwmap9 -lgfortran -framework Accelerate

# MY_COMPILE_FLAGS below defines some compiling options.
# include -D CREATE_LOG to create log file. CREATE_LOG1 will create a log file with more output and CREATE_LOG2 will have even more detailed output.
# include -D VERBOSE to compile in verbose mode. VERBOSE1 will have more output on the screen and VERBOSE2 will have even more output.
# include -D CHECKS_ON to turn on the checks (slows down, use it only for debugging)
# include -D SLOW_CHECKS_ON to turn on all the checks, including slow ones (may extremely slow down, use it only for debugging)
MY_COMPILE_FLAGS = -D VERBOSE1 -D CREATE_LOG -D CHECKS_ON
