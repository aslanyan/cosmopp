include make.inc

CC=gcc

ifdef OMP_FLAG
OMP_COMPILE_FLAGS = -D COSMO_OMP
endif

ifdef MPI_COMP
CC = $(MPI_COMP)
MPI_COMPILE_FLAGS = -D COSMO_MPI
endif

HEALPIX_AND_LAPACKPP_OBJ = obj/simulate.o obj/likelihood.o obj/master.o
HEALPIX_AND_LAPACKPP_TARGET = bin/generate_white_noise

HEALPIX_AND_LAPACK_AND_CLASS_TEST_OBJ = obj/test_cmb_gibbs.o obj/test_like_high.o

PLANCK_AND_CLASS_OBJ = obj/planck_like.o
PLANCK_AND_CLASS_TEST_OBJ = obj/test_mcmc_planck.o
PLANCK_AND_CLASS_TARGET = bin/test_planck

WMAP_AND_CLASS_OBJ = obj/wmap9_like.o
WMAP_AND_CLASS_TEST_OBJ = obj/test_wmap9_like.o

PLANCK_AND_MULTINEST_AND_CLASS_TEST_OBJ = obj/test_multinest_planck.o
PLANCK_AND_MULTINEST_AND_CLASS_TARGET = bin/mn_scan_planck

ifdef HEALPIX
ifndef CFITSIO
$(error Need to specify CFITSIO if HEALPIX is specified)
endif
ifndef HEALPIXPP
$(error Need to specify HEALPIXPP if HEALPIX is specified)
endif
HEALPIX_INCLUDE_FLAGS = -I $(CFITSIO)/include -I $(HEALPIX)/include -I $(HEALPIXPP)/include
HEALPIX_COMPILE_FLAGS = -D COSMO_HEALPIX -D HEALPIX_DATA_DIR=$(HEALPIX)/data
HEALPIX_LIB_FLAGS1 = -L $(CFITSIO)/lib -L $(HEALPIX)/lib -L $(HEALPIXPP)/lib
HEALPIX_LIB_FLAGS2 = -lchealpix -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio
HEALPIX_OBJ = obj/utils.o obj/c_matrix.o obj/c_matrix_generator.o obj/mode_directions.o obj/cmb_gibbs.o obj/mask_apodizer.o
HEALPIX_TARGET = bin/apodize_mask
else
HEALPIX_INCLUDE_FLAGS =
HEALPIX_COMPILE_FLAGS =
HEALPIX_LIB_FLAGS1 =
HEALPIX_LIB_FLAGS2 =
HEALPIX_OBJ =
HEALPIX_TARGET =
HEALPIX_AND_LAPACKPP_OBJ =
HEALPIX_AND_LAPACKPP_TARGET =
HEALPIX_AND_LAPACK_AND_CLASS_TEST_OBJ =
endif

ifdef LAPACKPPINCDIR
ifndef LAPACKPPLIBDIR
$(error Need to specify LAPACKPPLIBDIR if LAPACKPPINCDIR is specified)
endif
LAPACKPP_INCLUDE_FLAGS = -I $(LAPACKPPINCDIR)
LAPACKPP_COMPILE_FLAGS = -D COSMO_LAPACKPP
LAPACKPP_LIB_FLAGS1 = -L $(LAPACKPPLIBDIR)
LAPACKPP_LIB_FLAGS2 = -llapackpp
else
LAPACKPP_INCLUDE_FLAGS =
LAPACKPP_COMPILE_FLAGS =
LAPACKPP_LIB_FLAGS1 =
LAPACKPP_LIB_FLAGS2 =
HEALPIX_AND_LAPACKPP_OBJ =
HEALPIX_AND_LAPACKPP_TARGET =
HEALPIX_AND_LAPACK_AND_CLASS_TEST_OBJ =
endif

ifdef CLASSDIR
CLASS_INCLUDE_FLAGS = -I $(CLASSDIR)/include
CLASS_COMPILE_FLAGS = -D COSMO_CLASS
CLASS_LIB_FLAGS1 = -L $(CLASSDIR)
CLASS_LIB_FLAGS2 = -lclass
CLASS_OBJ = obj/cmb.o
CLASS_TEST_OBJ = obj/test_cmb.o
else
CLASS_INCLUDE_FLAGS =
CLASS_COMPILE_FLAGS =
CLASS_LIB_FLAGS1 =
CLASS_LIB_FLAGS2 =
CLASS_OBJ =
CLASS_TEST_OBJ =
PLANCK_AND_CLASS_OBJ =
PLANCK_AND_CLASS_TEST_OBJ =
PLANCK_AND_CLASS_TARGET =
PLANCK_AND_MULTINEST_AND_CLASS_TEST_OBJ =
PLANCK_AND_MULTINEST_AND_CLASS_TARGET =
WMAP_AND_CLASS_OBJ =
HEALPIX_AND_LAPACK_AND_CLASS_TEST_OBJ =
WMAP_AND_CLASS_TEST_OBJ =
endif

ifdef MINUIT
MINUIT_INCLUDE_FLAGS = -I $(MINUIT)/include
MINUIT_COMPILE_FLAGS = -D COSMO_MINUIT
MINUIT_LIB_FLAGS1 = -L $(MINUIT)/lib
MINUIT_LIB_FLAGS2 = -lMinuit2
MINUIT_TEST_OBJ = obj/test_fit.o
else
MINUIT_INCLUDE_FLAGS =
MINUIT_COMPILE_FLAGS =
MINUIT_LIB_FLAGS1 =
MINUIT_LIB_FLAGS2 =
MINUIT_TEST_OBJ =
endif

ifdef MULTINEST
MULTINEST_INCLUDE_FLAGS = -I $(MULTINEST)/includes
MULTINEST_COMPILE_FLAGS = -D COSMO_MULTINEST
MULTINEST_LIB_FLAGS1 = -L $(MULTINEST)/lib
MULTINEST_LIB_FLAGS2 = -lmultinest_mpi
MULTINEST_OBJ = obj/mn_scanner.o
MULTINEST_TEST_OBJ = obj/test_multinest.o
else
MULTINEST_INCLUDE_FLAGS =
MULTINEST_COMPILE_FLAGS =
MULTINEST_LIB_FLAGS1 =
MULTINEST_LIB_FLAGS2 =
MULTINEST_OBJ =
MULTINEST_TEST_OBJ =
PLANCK_AND_MULTINEST_AND_CLASS_TEST_OBJ =
PLANCK_AND_MULTINEST_AND_CLASS_TARGET =
endif

ifdef PLANCKDIR
ifndef PLANCKDATADIR
$(error Need to specify PLANCKDATADIR if PLANCKDIR is specified) 
endif
PLANCK_COMPILE_FLAGS = -D COSMO_PLANCK -D LAPACK_CLIK -D PLANCK_DATA_DIR=$(PLANCKDATADIR)
PLANCK_INCLUDE_FLAGS = -I $(PLANCKDIR)/include
PLANCK_LIB_FLAGS1 = -L $(PLANCKDIR)/lib
PLANCK_LIB_FLAGS2 = -lclik
PLANCK_OBJ =
PLANCK_TEST_OBJ =
else
PLANCK_COMPILE_FLAGS =
PLANCK_INCLUDE_FLAGS =
PLANCK_LIB_FLAGS1 =
PLANCK_LIB_FLAGS2 =
PLANCK_OBJ =
PLANCK_TEST_OBJ =
PLANCK_AND_MULTINEST_AND_CLASS_TEST_OBJ =
PLANCK_AND_MULTINEST_AND_CLASS_TARGET =
PLANCK_AND_CLASS_OBJ =
PLANCK_AND_CLASS_TEST_OBJ =
PLANCK_AND_CLASS_TARGET =
endif

ifdef WMAP9DIR
ifndef WMAP9COMPILER
$(error Need to specify WMAP9COMPILER if WMAP9DIR is specified)
endif
ifndef WMAP9LIBFLAGS
$(error Need to specify WMAP9LIBFLAGS if WMAP9DIR is specified)
endif
ifeq ($(WMAP9COMPILER),gfortran)
WMAP_COMPILE_FLAGS = -D COSMO_WMAP9 -D WMAP9_GFORT
else
WMAP_COMPILE_FLAGS = -D COSMO_WMAP9 -D WMAP9_IFORT
endif
WMAP_INCLUDE_FLAGS =
WMAP_LIB_FLAGS1 = -L $(WMAP9DIR)
WMAP_LIB_FLAGS2 = $(WMAP9LIBFLAGS)
WMAP_OBJ =
else
WMAP_COMPILE_FLAGS =
WMAP_INCLUDE_FLAGS =
WMAP_LIB_FLAGS1 =
WMAP_LIB_FLAGS2 =
WMAP_OBJ =
WMAP_AND_CLASS_OBJ =
WMAP_AND_CLASS_TEST_OBJ =
endif

LFLAGS1 = $(HEALPIX_LIB_FLAGS1) $(LAPACKPP_LIB_FLAGS1) $(CLASS_LIB_FLAGS1) $(MINUIT_LIB_FLAGS1) $(MULTINEST_LIB_FLAGS1) $(PLANCK_LIB_FLAGS1) $(WMAP_LIB_FLAGS1) -fopenmp

LFLAGS2 = -lstdc++ $(HEALPIX_LIB_FLAGS2) $(LAPACKPP_LIB_FLAGS2) $(CLASS_LIB_FLAGS2) $(MINUIT_LIB_FLAGS2) $(MULTINEST_LIB_FLAGS2) $(PLANCK_LIB_FLAGS2) $(WMAP_LIB_FLAGS2)

MYFLAGS = $(MY_COMPILE_FLAGS) $(HEALPIX_COMPILE_FLAGS) $(LAPACKPP_COMPILE_FLAGS) $(CLASS_COMPILE_FLAGS) $(MINUIT_COMPILE_FLAGS) $(MULTINEST_COMPILE_FLAGS) $(PLANCK_COMPILE_FLAGS) $(WMAP_COMPILE_FLAGS) $(MPI_COMPILE_FLAGS) $(OMP_COMPILE_FLAGS)

INCLUDE_FLAGS = -I include $(HEALPIX_INCLUDE_FLAGS) $(LAPACKPP_INCLUDE_FLAGS) $(CLASS_INCLUDE_FLAGS) $(MINUIT_INCLUDE_FLAGS) $(MULTINEST_INCLUDE_FLAGS) $(PLANCK_INCLUDE_FLAGS) $(WMAP_INCLUDE_FLAGS)

CFLAGS = $(MYFLAGS) -c -g -O2 -fpic -Werror -std=c++11 $(INCLUDE_FLAGS) $(OMP_FLAG)
CFLAGS_NO11 = $(MYFLAGS) -c -g -O2 -Werror -fpic $(INCLUDE_FLAGS) $(OMP_FLAG)

#Header file tree
MACROS_HPP = include/macros.hpp
EXCEPTION_HANDLER_HPP = include/exception_handler.hpp
PROGRESS_METER_HPP = include/progress_meter.hpp $(MACROS_HPP)
NUMERICS_HPP = include/numerics.hpp
TEST_FRAMEWORK_HPP = include/test_framework.hpp

ANGULAR_COORDINATES_HPP = include/angular_coordinates.hpp
COMPLEX_TYPES_HPP = include/complex_types.hpp
MATH_CONSTANTS_HPP = include/math_constants.hpp
PHYS_CONSTANTS_HPP = include/phys_constants.hpp $(MATH_CONSTANTS_HPP)
UNIT_CONVERSIONS_HPP = include/unit_conversions.hpp $(PHYS_CONSTANTS_HPP)
FUNCTION_HPP = include/function.hpp
HISTOGRAM_HPP = include/histogram.hpp
CHI_SQUARED_HPP = include/chi_squared.hpp $(MARCOR_HPP) $(FUNCTION_HPP)
RANDOM_HPP = include/random.hpp
INT_OPERATION_HPP = include/int_operations.hpp
LIKELIHOOD_FUNCTION_HPP = include/likelihood_function.hpp
INTEGRAL_HPP = include/integral.hpp $(MACROS_HPP) $(FUNCTION_HPP)
PARAMETRIC_FUNCTION_HPP = include/parametric_function.hpp $(MACROS_HPP) $(FUNCTION_HPP)
FIT_HPP = include/fit.hpp $(MACROS_HPP) $(PARAMETRIC_FUNCTION_HPP)
POLYNOMIAL_HPP = include/polynomial.hpp $(PARAMETRIC_FUNCTION_HPP)
LEGENDRE_HPP = include/legendre.hpp $(POLYNOMIAL_HPP)
TABLE_FUNCTION_HPP = include/table_function.hpp $(MACROS_HPP) $(FUNCTION_HPP)
CUBIC_SPLINE_HPP = include/cubic_spline.hpp $(MACROS_HPP) $(FUNCTION_HPP)
GAUSS_SMOOTH_HPP = include/gauss_smooth.hpp $(MACROS_HPP) $(FUNCTION_HPP)
THREE_VECTOR_HPP = include/three_vector.hpp $(ANGULAR_COORDINATES_HPP) $(MATH_CONSTANTS_HPP)
THREE_ROTATION_HPP = include/three_rotation.hpp $(MACROS_HPP) $(THREE_VECTOR_HPP) $(NUMERICS_HPP)
MCMC_HPP = include/mcmc.hpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(LIKELIHOOD_FUNCTION_HPP) $(RANDOM_HPP)
MN_SCANNER_HPP = include/mn_scanner.hpp $(LIKELIHOOD_FUNCTION_HPP) $(TABLE_FUNCTION_HPP)
WIGNER_3J_HPP = include/wigner_3j.hpp $(MACROS_HPP)
SPIN_SPHERICAL_HARMONICS_HPP = include/spin_spherical_harmonics.hpp $(MACROS_HPP) $(ANGULAR_COORDINATES_HPP) $(COMPLEX_TYPES_HPP) $(MATH_CONSTANTS_HPP)
CONJUGATE_GRADIENT_HPP = include/conjugate_gradient.hpp $(MACROS_HPP)
MARKOV_CHAIN_HPP = include/markov_chain.hpp $(FUNCTION_HPP) $(TABLE_FUNCTION_HPP)

C_MATRIX_HPP = include/c_matrix.hpp
WHOLE_MATRIX_HPP = include/whole_matrix.hpp
C_MATRIX_GENERATOR_HPP = include/c_matrix_generator.hpp $(C_MATRIX_HPP) $(WHOLE_MATRIX_HPP)
LIKELIHOOD_HPP = include/likelihood.hpp $(COMPLEX_TYPES_HPP) $(C_MATRIX_HPP) $(WHOLE_MATRIX_HPP)
MASTER_HPP = include/master.hpp $(MACROS_HPP)
MODE_DIRECTIONS_HPP = include/mode_directions.hpp
SIMULATE_HPP = include/simulate.hpp $(WHOLE_MATRIX_HPP)
UTILS_HPP = include/utils.hpp
POWER_SPECTRUM_HPP = include/power_spectrum.hpp $(FUNCTION_HPP) $(TABLE_FUNCTION_HPP) $(CUBIC_SPLINE_HPP)
COSMOLOGICAL_PARAMS_HPP = include/cosmological_params.hpp $(MACROS_HPP) $(PHYS_CONSTANTS_HPP) $(UNIT_CONVERSIONS_HPP) $(POWER_SPECTRUM_HPP)
CMB_HPP = include/cmb.hpp $(TABLE_FUNCTION_HPP) $(COSMOLOGICAL_PARAMS_HPP)
PLANCK_LIKE_HPP = include/planck_like.hpp $(LIKELIHOOD_FUNCTION_HPP) $(CMB_HPP)
WMAP9_LIKE_HPP = include/wmap9_like.hpp $(MACROS_HPP) $(LIKELIHOOD_FUNCTION_HPP) $(CMB_HPP)
SCALE_FACTOR_HPP = include/scale_factor.hpp $(MACROS_HPP) $(TABLE_FUNCTION_HPP) $(COSMOLOGICAL_PARAMS_HPP)
CMB_GIBBS_HPP = include/cmb_gibbs.hpp $(RANDOM_HPP)
MASK_APODIZER_HPP = include/mask_apodizer.hpp
INFLATION_HPP = include/inflation.hpp $(FUNCTION_HPP) $(TABLE_FUNCTION_HPP)

TEST_UNIT_CONVERSIONS_HPP = include/test_unit_conversions.hpp $(TEST_FRAMEWORK_HPP)
TEST_INT_OPERATIONS_HPP = include/test_int_operations.hpp $(TEST_FRAMEWORK_HPP)
TEST_INTEGRAL_HPP = include/test_integral.hpp $(TEST_FRAMEWORK_HPP)
TEST_CONJUGATE_GRADIENT_HPP = include/test_conjugate_gradient.hpp $(TEST_FRAMEWORK_HPP)
TEST_POLYNOMIAL_HPP = include/test_polynomial.hpp $(TEST_FRAMEWORK_HPP)
TEST_LEGENDRE_HPP = include/test_legendre.hpp $(TEST_FRAMEWORK_HPP)
TEST_MCMC_HPP = include/test_mcmc.hpp $(TEST_FRAMEWORK_HPP)
TEST_MULTINEST_HPP = include/test_multinest.hpp $(TEST_FRAMEWORK_HPP)
TEST_MCMC_PLANCK_HPP = include/test_mcmc_planck.hpp $(TEST_FRAMEWORK_HPP)
TEST_MULTINEST_PLANCK_HPP = include/test_multinest_planck.hpp $(TEST_FRAMEWORK_HPP)
TEST_CMB_HPP = include/test_cmb.hpp $(TEST_FRAMEWORK_HPP)
TEST_CMB_GIBBS_HPP = include/test_cmb_gibbs.hpp $(TEST_FRAMEWORK_HPP)
TEST_FIT = include/test_fit.hpp $(TEST_FRAMEWORK_HPP)
TEST_WMAP9_LIKE = include/test_wmap9_like.hpp $(TEST_FRAMEWORK_HPP)
TEST_LIKE_HIGH = include/test_like_high.hpp $(TEST_FRAMEWORK_HPP)

all: lib/libcosmopp.a bin/test $(HEALPIX_AND_LAPACKPP_TARGET) $(HEALPIX_TARGET) $(PLANCK_AND_CLASS_TARGET) $(PLANCK_AND_MULTINEST_AND_CLASS_TARGET)

GENERAL_OBJ = obj/macros.o obj/test_framework.o obj/mcmc.o obj/whole_matrix.o obj/scale_factor.o obj/markov_chain.o obj/inflation.o

OBJ_LIBRARY = $(GENERAL_OBJ) $(HEALPIX_OBJ) $(HEALPIX_AND_LAPACKPP_OBJ) $(CLASS_OBJ) $(MULTINEST_OBJ) $(PLANCK_OBJ) $(PLANCK_AND_CLASS_OBJ) $(WMAP_OBJ) $(WMAP_AND_CLASS_OBJ)
lib/libcosmopp.a: $(OBJ_LIBRARY)
	ar rcs $@ $(OBJ_LIBRARY)

GENERAL_TEST_OBJ = obj/test.o obj/test_unit_conversions.o obj/test_int_operations.o obj/test_integral.o obj/test_conjugate_gradient.o obj/test_polynomial.o obj/test_legendre.o obj/test_mcmc.o

OBJ_TEST = $(OBJ_LIBRARY) $(GENERAL_TEST_OBJ) $(MULTINEST_TEST_OBJ) $(PLANCK_TEST_OBJ) $(PLANCK_AND_MULTINEST_AND_CLASS_TEST_OBJ) $(CLASS_TEST_OBJ) $(PLANCK_AND_CLASS_TEST_OBJ) $(HEALPIX_AND_LAPACK_AND_CLASS_TEST_OBJ) $(MINUIT_TEST_OBJ) $(WMAP_AND_CLASS_TEST_OBJ)
bin/test: $(OBJ_TEST)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_TEST) $(LFLAGS2)

ifdef HEALPIX_AND_LAPACKPP_TARGET
OBJ_GENERATE_WHITE_NOISE = obj/generate_white_noise.o obj/macros.o obj/simulate.o obj/whole_matrix.o
bin/generate_white_noise: $(OBJ_GENERATE_WHITE_NOISE)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_GENERATE_WHITE_NOISE) $(LFLAGS2)
endif

ifdef HEALPIX_TARGET
OBJ_APODIZE_MASK = obj/apodize_mask.o obj/macros.o obj/mask_apodizer.o
bin/apodize_mask: $(OBJ_APODIZE_MASK)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_APODIZE_MASK) $(LFLAGS2)
endif

ifdef PLANCK_AND_CLASS_TARGET
OBJ_TEST_PLANCK = obj/test_planck.o obj/macros.o obj/planck_like.o obj/cmb.o
bin/test_planck: $(OBJ_TEST_PLANCK)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_TEST_PLANCK) $(LFLAGS2)
endif

ifdef PLANCK_AND_MULTINEST_AND_CLASS_TARGET
OBJ_MN_SCAN_PLANCK = obj/mn_scan_planck.o obj/macros.o obj/planck_like.o obj/cmb.o obj/mn_scanner.o
bin/mn_scan_planck: $(OBJ_MN_SCAN_PLANCK)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_MN_SCAN_PLANCK) $(LFLAGS2)
endif

obj/macros.o: source/macros.cpp $(MACROS_HPP)
	$(CC) $(CFLAGS) source/macros.cpp -o $@

obj/c_matrix.o: source/c_matrix.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(UTILS_HPP) $(C_MATRIX_HPP)
	$(CC) $(CFLAGS) source/c_matrix.cpp -o $@

obj/c_matrix_generator.o: source/c_matrix_generator.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(ANGULAR_COORDINATES_HPP) $(THREE_ROTATION_HPP) $(PROGRESS_METER_HPP) $(C_MATRIX_GENERATOR_HPP) $(UTILS_HPP) $(RANDOM_HPP) $(LEGENDRE_HPP)
	$(CC) $(CFLAGS) source/c_matrix_generator.cpp -o $@

obj/cmb.o: source/cmb.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(CMB_HPP)
	$(CC) $(CFLAGS) source/cmb.cpp -o $@

obj/likelihood.o: source/likelihood.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(PROGRESS_METER_HPP) $(THREE_ROTATION_HPP) $(LIKELIHOOD_HPP) $(UTILS_HPP)
	$(CC) $(CFLAGS) source/likelihood.cpp -o $@

obj/master.o: source/master.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(PROGRESS_METER_HPP) $(WIGNER_3J_HPP) $(MASTER_HPP)
	$(CC) $(CFLAGS) source/master.cpp -o $@

ifdef PLANCK_AND_MULTINEST_AND_CLASS_TARGET
obj/mn_scan_planck.o: source/mn_scan_planck.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(LIKELIHOOD_FUNCTION_HPP) $(PLANCK_LIKE_HPP) $(MN_SCANNER_HPP)
	$(CC) $(CFLAGS) source/mn_scan_planck.cpp -o $@
endif

ifdef MULTINEST
obj/mn_scanner.o: source/mn_scanner.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(MN_SCANNER_HPP)
	$(CC) $(CFLAGS) source/mn_scanner.cpp -o $@
endif

obj/mode_directions.o: source/mode_directions.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(ANGULAR_COORDINATES_HPP) $(THREE_ROTATION_HPP) $(PROGRESS_METER_HPP) $(MATH_CONSTANTS_HPP) $(NUMERICS_HPP) $(MODE_DIRECTIONS_HPP)
	$(CC) $(CFLAGS) source/mode_directions.cpp -o $@

ifdef PLANCKDIR
obj/planck_like.o: source/planck_like.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(PLANCK_LIKE_HPP)
	$(CC) $(CFLAGS_NO11) source/planck_like.cpp -o $@
endif

ifdef WMAP9DIR
obj/wmap9_like.o: source/wmap9_like.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(WMAP9_LIKE_HPP)
	$(CC) $(CFLAGS) source/wmap9_like.cpp -o $@
endif

obj/scale_factor.o: source/scale_factor.cpp $(UNIT_CONVERSIONS_HPP) $(TABLE_FUNCTION_HPP) $(SCALE_FACTOR_HPP)
	$(CC) $(CFLAGS) source/scale_factor.cpp -o $@

obj/simulate.o: source/simulate.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(WHOLE_MATRIX_HPP) $(RANDOM_HPP) $(SIMULATE_HPP)
	$(CC) $(CFLAGS) source/simulate.cpp -o $@

obj/test.o: source/test.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(TEST_FRAMEWORK_HPP) $(TEST_UNIT_CONVERSIONS_HPP) $(TEST_INT_OPERATIONS_HPP) $(TEST_INTEGRAL_HPP) $(TEST_CONJUGATE_GRADIENT_HPP) $(TEST_POLYNOMIAL_HPP) $(TEST_LEGENDRE_HPP) $(TEST_MCMC_HPP) $(TEST_MULTINEST_HPP) $(TEST_MCMC_PLANCK_HPP) $(TEST_CMB_HPP) $(TEST_CMB_GIBBS_HPP) $(TEST_FIT_HPP) $(TEST_WMAP9_LIKE_HPP) $(TEST_LIKE_HIGH_HPP)
	$(CC) $(CFLAGS) source/test.cpp -o $@

obj/test_framework.o: source/test_framework.cpp $(MACROS_HPP) $(NUMERICS_HPP) $(TEST_FRAMEWORK_HPP)
	$(CC) $(CFLAGS) source/test_framework.cpp -o $@

ifdef PLANCK_AND_CLASS_TARGET
obj/test_planck.o: source/test_planck.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(PLANCK_LIKE_HPP)
	$(CC) $(CFLAGS) source/test_planck.cpp -o $@
endif

obj/utils.o: source/utils.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(UTILS_HPP)
	$(CC) $(CFLAGS) source/utils.cpp -o $@

obj/whole_matrix.o: source/whole_matrix.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(THREE_ROTATION_HPP) $(WHOLE_MATRIX_HPP)
	$(CC) $(CFLAGS) source/whole_matrix.cpp -o $@

obj/cmb_gibbs.o: source/cmb_gibbs.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(PROGRESS_METER_HPP) $(UTILS_HPP) $(CONJUGATE_GRADIENT_HPP) $(NUMERICS_HPP) $(CMB_GIBBS_HPP)
	$(CC) $(CFLAGS) source/cmb_gibbs.cpp -o $@

ifdef HEALPIX_AND_LAPACKPP_TARGET
obj/generate_white_noise.o: source/generate_white_noise.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(SIMULATE_HPP)
	$(CC) $(CFLAGS) source/generate_white_noise.cpp -o $@
endif

obj/mask_apodizer.o: source/mask_apodizer.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(PROGRESS_METER_HPP) $(THREE_VECTOR_HPP) $(MATH_CONSTANTS_HPP) $(MASK_APODIZER_HPP)
	$(CC) $(CFLAGS) source/mask_apodizer.cpp -o $@

ifdef HEALPIX_TARGET
obj/apodize_mask.o: source/apodize_mask.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MASK_APODIZER_HPP)
	$(CC) $(CFLAGS) source/apodize_mask.cpp -o $@
endif

obj/mcmc.o: source/mcmc.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MCMC_HPP)
	$(CC) $(CFLAGS) source/mcmc.cpp -o $@

obj/markov_chain.o: source/markov_chain.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(CUBIC_SPLINE_HPP) $(GAUSS_SMOOTH_HPP) $(PROGRESS_METER_HPP) $(MARKOV_CHAIN_HPP)
	$(CC) $(CFLAGS) source/markov_chain.cpp -o $@

obj/inflation.o: source/inflation.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(PROGRESS_METER_HPP) $(UNIT_CONVERSIONS_HPP) $(INFLATION_HPP)
	$(CC) $(CFLAGS) source/inflation.cpp -o $@

obj/test_unit_conversions.o: source/test_unit_conversions.cpp $(MACROS_HPP) $(TEST_UNIT_CONVERSIONS_HPP) $(UNIT_CONVERSIONS_HPP)
	$(CC) $(CFLAGS) source/test_unit_conversions.cpp -o $@

obj/test_int_operations.o: source/test_int_operations.cpp $(MACROS_HPP) $(TEST_INT_OPERATIONS_HPP) $(INT_OPERATIONS_HPP)
	$(CC) $(CFLAGS) source/test_int_operations.cpp -o $@

obj/test_integral.o: source/test_integral.cpp $(TEST_INTEGRAL_HPP) $(INTEGRAL_HPP)
	$(CC) $(CFLAGS) source/test_integral.cpp -o $@

obj/test_conjugate_gradient.o: source/test_conjugate_gradient.cpp $(TEST_CONJUGATE_GRADIENT_HPP) $(CONJUGATE_GRADIENT_HPP)
	$(CC) $(CFLAGS) source/test_conjugate_gradient.cpp -o $@

obj/test_polynomial.o: source/test_polynomial.cpp $(TEST_POLYNOMIAL_HPP) $(POLYNOMIAL_HPP)
	$(CC) $(CFLAGS) source/test_polynomial.cpp -o $@

obj/test_legendre.o: source/test_legendre.cpp $(TEST_LEGENDRE_HPP) $(LEGENDRE_HPP)
	$(CC) $(CFLAGS) source/test_legendre.cpp -o $@

obj/test_mcmc.o: source/test_mcmc.cpp $(TEST_MCMC_HPP) $(MCMC_HPP) $(MARKOV_CHAIN_HPP) $(NUMERICS_HPP)
	$(CC) $(CFLAGS) source/test_mcmc.cpp -o $@

obj/test_multinest.o: source/test_multinest.cpp $(TEST_MULTINEST_HPP) $(MN_SCANNER_HPP) $(MARKOV_CHAIN_HPP) $(NUMERICS_HPP)
	$(CC) $(CFLAGS) source/test_multinest.cpp -o $@

obj/test_mcmc_planck.o: source/test_mcmc_planck.cpp $(TEST_MCMC_PLANCK_HPP) $(MCMC_HPP) $(PLANCK_LIKE_HPP) $(MARKOV_CHAIN_HPP) $(NUMERICS_HPP)
	$(CC) $(CFLAGS) source/test_mcmc_planck.cpp -o $@

obj/test_multinest_planck.o: source/test_multinest_planck.cpp $(TEST_MULTINEST_PLANCK_HPP) $(MN_SCANNER_HPP) $(PLANCK_LIKE_HPP) $(MARKOV_CHAIN_HPP) $(NUMERICS_HPP)
	$(CC) $(CFLAGS) source/test_multinest_planck.cpp -o $@

obj/test_cmb.o: source/test_cmb.cpp $(CMB_HPP) $(TEST_CMB_HPP)
	$(CC) $(CFLAGS) source/test_cmb.cpp -o $@

obj/test_cmb_gibbs.o: source/test_cmb_gibbs.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(CMB_GIBBS_HPP) $(CMB_HPP) $(SIMULATE_HPP) $(UTILS_HPP) $(RANDOM_HPP) $(MARKOV_CHAIN_HPP) $(NUMERICS_HPP) $(PROGRESS_METER_HPP) $(TEST_CMB_GIBBS_HPP)
	$(CC) $(CFLAGS) source/test_cmb_gibbs.cpp -o $@

obj/test_fit.o: source/test_fit.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(FIT_HPP) $(POLYNOMIAL_HPP) $(TEST_FIT_HPP)
	$(CC) $(CFLAGS) source/test_fit.cpp -o $@

obj/test_wmap9_like.o: source/test_wmap9_like.cpp $(WMAP9_LIKE_HPP) $(NUMERICS_HPP) $(TEST_WMAP9_LIKE_HPP)
	$(CC) $(CFLAGS) source/test_wmap9_like.cpp -o $@

obj/test_like_high.o: source/test_like_high.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(UTILS_HPP) $(HISTOGRAM_HPP) $(LIKELIHOOD_HPP) $(CMB_HPP) $(MASTER_HPP) $(SIMULATE_HPP) $(PROGRESS_METER_HPP) $(CHI_SQUARED_HPP) $(TEST_LIKE_HIGH_HPP)
	$(CC) $(CFLAGS) source/test_like_high.cpp -o $@

clean:
	rm obj/* bin/* lib/*

test: all
	./bin/test fast
