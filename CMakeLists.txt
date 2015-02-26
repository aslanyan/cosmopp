cmake_minimum_required (VERSION 2.8)

project ("Cosmo++" CXX)

# compiler options
set(CMAKE_BUILD_TYPE RELEASE)

set(CMAKE_CXX_FLAGS "-std=c++11 -fpic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb")
if(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
	set(CMAKE_CXX_FLAGS_RELEASE "-O4 -ffast-math")
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
	set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
else(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
	set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)

# specific options
SET(CMAKE_INSTALL_PREFIX /Users/aslanyan/code/install_cosmo_pp)

# cosmo++ specific compiler flags
# include -D CREATE_LOG to create log file. CREATE_LOG1 will create a log file with more output and CREATE_LOG2 will have even more detailed output.
# include -D VERBOSE to compile in verbose mode. VERBOSE1 will have more output on the screen and VERBOSE2 will have even more output.
# include -D CHECKS_ON to turn on the checks (slows down, use it only for debugging)
# include -D SLOW_CHECKS_ON to turn on all the checks, including slow ones (may extremely slow down, use it only for debugging)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D VERBOSE1 -D CREATE_LOG -D CHECKS_ON")

#add openmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -D COSMO_OMP")

#check for MPI
find_package(MPI)
if(MPI_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D COSMO_MPI")
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    MESSAGE("${MPI_C_LIBRARIES}")
endif(MPI_FOUND)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/test_files)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/slow_test_files)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(source)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin/ DESTINATION bin)
install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/ DESTINATION lib)

