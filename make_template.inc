CFITSIO = /usr/local
HEALPIX = /Users/aslanyan/programs/Healpix_3.11
HEALPIXPP = $(HEALPIX)/src/cxx/osx
BOOSTDIR = /Users/aslanyan/C++/library/boost_1_51_0
LAPACKPPLIBDIR = /usr/local/lib
LAPACKPPINCDIR = /Users/aslanyan/C++/library/lapackpp-2.5.4/include
CLASSDIR = /Users/aslanyan/programs/class_v2.0.2

MINUIT =
MULTINEST =
PLANCKDIR =
PLANCKDATADIR =

# MY_COMPILE_FLAGS below defines some compiling options.
# include -D CREATE_LOG to create log file
# include -D VERBOSE to compile in verbose mode
# include -D CHECKS_ON to turn on the checks (slows down, use it only for debugging)
# include -D SLOW_CHECKS_ON to turn on all the checks, including slow ones (extremely slows down, use it only for debugging)
MY_COMPILE_FLAGS = -D VERBOSE -D CHECKS_ON
