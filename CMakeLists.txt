cmake_minimum_required (VERSION 2.8.10)

project ("Cosmo++" CXX)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE RELEASE)
endif (NOT CMAKE_BUILD_TYPE)

include(cmake_settings.txt)

# compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb")
if(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
	set(CMAKE_CXX_FLAGS_RELEASE "-O4 -ffast-math")
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
	set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
else(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)

#check for openmp
find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	add_definitions(-DCOSMO_OMP)
endif(OPENMP_FOUND)

#enable fortran if polychord is specified
if(POLYCHORD_DIR)
	enable_language (Fortran)
endif(POLYCHORD_DIR)

#check for MPI
find_package(MPI)
if(MPI_FOUND)
	add_definitions(-DCOSMO_MPI)
	#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
	#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
	include_directories(SYSTEM ${MPI_INCLUDE_PATH})
	if(NOT MPI_CXX_LIBRARIES)
		message(FATAL_ERROR "MPI found but the MPI C++ libraries not found!")
	endif(NOT MPI_CXX_LIBRARIES)
endif(MPI_FOUND)

#lapack
if(LAPACK_LIB_FLAGS)
	add_definitions(-DCOSMO_LAPACK)
endif(LAPACK_LIB_FLAGS)

#cfitsio
if(CFITSIO_DIR)
	include_directories(${CFITSIO_DIR}/include)
	find_library(CFITSIOLIB cfitsio ${CFITSIO_DIR}/lib)
	if(NOT CFITSIOLIB)
		message(FATAL_ERROR "cfitsio library not found!")
	endif(NOT CFITSIOLIB)
endif(CFITSIO_DIR)

#healpix
if(HEALPIX_DIR)
	if(NOT CFITSIO_DIR)
		message(FATAL_ERROR "CFITSIO_DIR must be specified if HEALPIX_DIR is specified")
	endif(NOT CFITSIO_DIR)
	if(NOT HEALPIXPP_DIR)
		message(FATAL_ERROR "HEALPIXPP_DIR must be specified if HEALPIX_DIR is specified")
	endif(NOT HEALPIXPP_DIR)

	include_directories(${HEALPIX_DIR}/include ${HEALPIXPP_DIR}/include)
	add_definitions(-DCOSMO_HEALPIX -DHEALPIX_DATA_DIR=${HEALPIX_DIR}/data)

	find_library(HEALPIXCXXLIB healpix_cxx ${HEALPIXPP_DIR}/lib)
	if(NOT HEALPIXCXXLIB)
		message(FATAL_ERROR "healpix_cxx library not found!")
	endif(NOT HEALPIXCXXLIB)

	find_library(CXXSUPPORTLIB cxxsupport ${HEALPIXPP_DIR}/lib)
	if(NOT CXXSUPPORTLIB)
		message(FATAL_ERROR "cxxsupport library not found!")
	endif(NOT CXXSUPPORTLIB)

	find_library(SHARPLIB sharp ${HEALPIXPP_DIR}/lib)
	if(NOT SHARPLIB)
		message(FATAL_ERROR "sharp library not found!")
	endif(NOT SHARPLIB)

	find_library(FFTPACKLIB fftpack ${HEALPIXPP_DIR}/lib)
	if(NOT FFTPACKLIB)
		message(FATAL_ERROR "fftpack library not found!")
	endif(NOT FFTPACKLIB)

	find_library(CHEALPIXLIB chealpix ${HEALPIX_DIR}/lib)
	if(NOT CHEALPIXLIB)
		message(FATAL_ERROR "chealpix library not found!")
	endif(NOT CHEALPIXLIB)

	find_library(CUTILSLIB c_utils ${HEALPIXPP_DIR}/lib)
	if(NOT CUTILSLIB)
		message(FATAL_ERROR "c_utils library not found!")
	endif(NOT CUTILSLIB)
endif(HEALPIX_DIR)

#class
if(CLASS_DIR)
	include_directories(${CLASS_DIR}/include)
	add_definitions(-DCOSMO_CLASS)

	find_library(CLASSLIB class ${CLASS_DIR})
	if(NOT CLASSLIB)
		message(FATAL_ERROR "class library not found!")
	endif(NOT CLASSLIB)
endif(CLASS_DIR)

#multinest
if(MULTINEST_DIR)
	include_directories(${MULTINEST_DIR}/include)
	add_definitions(-DCOSMO_MULTINEST)

	if(MPI_FOUND)
		find_library(MULTINESTLIB multinest_mpi ${MULTINEST_DIR}/lib)
	else(MPI_FOUND)
		find_library(MULTINESTLIB multinest ${MULTINEST_DIR}/lib)
	endif(MPI_FOUND)

	if(NOT MULTINESTLIB)
		message(FATAL_ERROR "multinest library not found!")
	endif(NOT MULTINESTLIB)
endif(MULTINEST_DIR)

#polychord
if(POLYCHORD_DIR)
	include_directories(${POLYCHORD_DIR}/src)
	add_definitions(-DCOSMO_POLYCHORD)
	set(CMAKE_Fortran_FLAGS "-cpp")
	get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
	if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
		add_definitions(-DPOLYCHORD_GFORT)
		set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
		set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -ggdb")
	elseif(Fortran_COMPILER_NAME MATCHES "ifort.*")
		add_definitions(-DPOLYCHORD_IFORT)
		set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
		set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -ggdb")
	else(Fortran_COMPILER_NAME MATCHES "gfortran.*")
		message(FATAL_ERROR "The fortran compiler supported for Polychord can only be gfortran or ifort. Please use one of these compilers by setting the environment variable FC or turn off Polychord")
	endif(Fortran_COMPILER_NAME MATCHES "gfortran.*")
	find_library(POLYCHORDLIB chord ${POLYCHORD_DIR}/src)
	if(NOT POLYCHORDLIB)
		message(FATAL_ERROR "chord library not found!")
	endif(NOT POLYCHORDLIB)
	if(MPI_FOUND)
		add_definitions(${MPI_Fortran_COMPILE_FLAGS})
		include_directories(${MPI_Fortran_INCLUDE_PATH})
		set(CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} ${MPI_LINK_FLAGS}")
	endif(MPI_FOUND)
endif(POLYCHORD_DIR)

#planck
if(PLANCK_DIR)
	if(NOT PLANCK_DATA_DIR)
		message(FATAL_ERROR "PLANCK_DATA_DIR must be specified if PLANCK_DIR is specified")
	endif(NOT PLANCK_DATA_DIR)

	include_directories(${PLANCK_DIR}/include)
	add_definitions(-DCOSMO_PLANCK -DPLANCK_DATA_DIR=${PLANCK_DATA_DIR})
	if(PLANCK_COMPILE_FLAGS)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PLANCK_COMPILE_FLAGS}")
	endif(PLANCK_COMPILE_FLAGS)
	find_library(PLANCKLIB clik ${PLANCK_DIR}/lib)
	if(NOT PLANCKLIB)
		message(FATAL_ERROR "clik library not found!")
	endif(NOT PLANCKLIB)
endif(PLANCK_DIR)

#wmap9
if(WMAP9_DIR)
	if(NOT WMAP9_COMPILER)
		message(FATAL_ERROR "WMAP9_COMPILER must be specified if WMAP9_DIR is specified")
	endif(NOT WMAP9_COMPILER)
	if(NOT LAPACK_LIB_FLAGS)
		message(FATAL_ERROR "LAPACK_LIB_FLAGS must be specified if WMAP9_DIR is specified")
	endif(NOT LAPACK_LIB_FLAGS)
	if(NOT CFITSIO_DIR)
		message(FATAL_ERROR "CFITSIO_DIR must be specified if WMAP9_DIR is specified")
	endif(NOT CFITSIO_DIR)

	add_definitions(-DCOSMO_WMAP9)
	if(WMAP9_COMPILER MATCHES "gfortran.*")
		add_definitions(-DWMAP9_GFORT)
	elseif(WMAP9_COMPILER MATCHES "ifort.*")
		add_definitions(-DWMAP9_IFORT)
	else(WMAP9_COMPILER MATCHES "gfortran.*")
		message(FATAL_ERROR "WMAP9_COMPILER can only be gfortran or ifort, others not supported at the time")
	endif(WMAP9_COMPILER MATCHES "gfortran.*")

	find_library(WMAP9LIB wmap9 ${WMAP9_DIR})
	if(NOT WMAP9LIB)
		message(FATAL_ERROR "wmap9 library not found!")
	endif(NOT WMAP9LIB)
endif(WMAP9_DIR)

#minuit
if(MINUIT_DIR)
	include_directories(${MINUIT_DIR}/include)
	add_definitions(-DCOSMO_MINUIT)

	find_library(MINUITLIB Minuit2 ${MINUIT_DIR}/lib)
	if(NOT MINUITLIB)
		message(FATAL_ERROR "Minuit2 library not found!")
	endif(NOT MINUITLIB)
endif(MINUIT_DIR)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/test_files)
file(COPY ${PROJECT_SOURCE_DIR}/slow_test_files DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR})

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(source)
add_subdirectory(examples)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
#install(DIRECTORY ${PROJECT_BINARY_DIR}/bin DESTINATION bin)
#install(DIRECTORY ${PROJECT_BINARY_DIR}/lib DESTINATION lib)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/data DESTINATION ${CMAKE_INSTALL_PREFIX})

