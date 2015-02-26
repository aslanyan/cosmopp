cmake_minimum_required (VERSION 2.8)

set(LIB_FILES macros.cpp test_framework.cpp whole_matrix.cpp scale_factor.cpp markov_chain.cpp)

set(TEST_FILES test_unit_conversions.cpp test_int_operations.cpp test_integral.cpp test_conjugate_gradient.cpp test_polynomial.cpp test_legendre.cpp test_spherical_harmonics.cpp test_matrix.cpp test_wigner_3j.cpp test_table_function.cpp test_cubic_spline.cpp test_three_rotation.cpp)

add_library(cosmopp STATIC ${LIB_FILES})

add_executable(test test.cpp ${TEST_FILES})
target_link_libraries(test cosmopp)
if(MPI_FOUND)
	MESSAGE("${MPI_C_LIBRARIES}")
	target_link_libraries(test ${MPI_C_LIBRARIES})
endif(MPI_FOUND)

