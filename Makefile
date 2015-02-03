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

HEALPIX_AND_LAPACK_AND_CLASS_TEST_OBJ = obj/test_cmb_gibbs.o obj/test_like_high.o obj/test_like_low.o

CLASS_AND_LAPACK_OBJ = obj/matter_likelihood.o
CLASS_AND_LAPACK_TEST_OBJ = obj/test_matter_likelihood.o

PLANCK_AND_CLASS_OBJ = obj/planck_like.o
PLANCK_AND_CLASS_TEST_OBJ = obj/test_planck_like.o
PLANCK_AND_CLASS_TARGET = bin/example_planck

PLANCK_AND_CLASS_AND_LAPACK_TEST_OBJ = obj/test_mcmc_planck.o

WMAP_AND_CLASS_OBJ = obj/wmap9_like.o
WMAP_AND_CLASS_TEST_OBJ = obj/test_wmap9_like.o

PLANCK_AND_MULTINEST_AND_CLASS_TEST_OBJ = obj/test_multinest_planck.o
PLANCK_AND_MULTINEST_AND_CLASS_TARGET = bin/example_mn_planck

ANN_AND_LAPACK_OBJ = obj/fast_approximator.o obj/fast_approximator_error.o obj/learn_as_you_go.o
ANN_AND_LAPACK_TEST_OBJ = obj/test_fast_approximator.o obj/test_fast_approximator_error.o

MINUIT_AND_LAPACKPP_OBJ = obj/gaussian_process.o
MINUIT_AND_LAPACKPP_TEST_OBJ = obj/test_gaussian_process.o


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
HEALPIX_TEST_OBJ = obj/test_mask_apodizer.o
HEALPIX_TARGET = bin/apodize_mask
else
HEALPIX_INCLUDE_FLAGS =
HEALPIX_COMPILE_FLAGS =
HEALPIX_LIB_FLAGS1 =
HEALPIX_LIB_FLAGS2 =
HEALPIX_OBJ =
HEALPIX_TEST_OBJ =
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
LAPACKPP_OBJ = obj/mcmc.o
LAPACKPP_TEST_OBJ = obj/test_mcmc.o
LAPACKPP_TARGET = bin/example_metropolis_hastings
else
LAPACKPP_INCLUDE_FLAGS =
LAPACKPP_COMPILE_FLAGS =
LAPACKPP_LIB_FLAGS1 =
LAPACKPP_LIB_FLAGS2 =
LAPACKPP_OBJ =
CLASS_AND_LAPACK_OBJ =
CLASS_AND_LAPACK_TEST_OBJ =
HEALPIX_AND_LAPACKPP_OBJ =
HEALPIX_AND_LAPACKPP_TARGET =
HEALPIX_AND_LAPACK_AND_CLASS_TEST_OBJ =
PLANCK_AND_CLASS_AND_LAPACK_TEST_OBJ =
ANN_AND_LAPACK_OBJ =
ANN_AND_LAPACK_TEST_OBJ =
MINUIT_AND_LAPACKPP_OBJ =
MINUIT_AND_LAPACKPP_TEST_OBJ =
endif

ifdef ANNDIR
ANN_INCLUDE_FLAGS = -I $(ANNDIR)/include
ANN_COMPILE_FLAGS = -D COSMO_ANN
ANN_LIB_FLAGS1 = -L $(ANNDIR)/lib
ANN_LIB_FLAGS2 = -lANN
ANN_OBJ = obj/k_nearest_neighbors.o
ANN_TEST_OBJ = obj/test_k_nearest_neighbors.o
else
ANN_INCLUDE_FLAGS =
ANN_LIB_FLAGS1 =
ANN_LIB_FLAGS2 =
ANN_OBJ =
ANN_TEST_OBJ =
ANN_AND_LAPACK_OBJ =
ANN_AND_LAPACK_TEST_OBJ =
endif

ifdef CLASSDIR
CLASS_INCLUDE_FLAGS = -I $(CLASSDIR)/include
CLASS_COMPILE_FLAGS = -D COSMO_CLASS
CLASS_LIB_FLAGS1 = -L $(CLASSDIR)
CLASS_LIB_FLAGS2 = -lclass
CLASS_OBJ = obj/cmb.o
CLASS_TEST_OBJ = obj/test_cmb.o
CLASS_TARGET = bin/example_cl
else
CLASS_INCLUDE_FLAGS =
CLASS_COMPILE_FLAGS =
CLASS_LIB_FLAGS1 =
CLASS_LIB_FLAGS2 =
CLASS_OBJ =
CLASS_TEST_OBJ =
CLASS_TARGET =
CLASS_AND_LAPACK_OBJ =
CLASS_AND_LAPACK_TEST_OBJ =
PLANCK_AND_CLASS_OBJ =
PLANCK_AND_CLASS_TEST_OBJ =
PLANCK_AND_CLASS_TARGET =
PLANCK_AND_CLASS_AND_LAPACK_TEST_OBJ =
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
MINUIT_AND_LAPACKPP_OBJ =
MINUIT_AND_LAPACKPP_TEST_OBJ =
endif

ifdef MULTINEST
MULTINEST_INCLUDE_FLAGS = -I $(MULTINEST)/include
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
PLANCK_COMPILE_FLAGS = -D COSMO_PLANCK -D PLANCK_DATA_DIR=$(PLANCKDATADIR) $(PLANCKCOMPILEFLAGS)
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
PLANCK_AND_CLASS_AND_LAPACK_TEST_OBJ =
endif

ifdef WMAP9DIR
ifndef WMAP9COMPILER
$(error Need to specify WMAP9COMPILER if WMAP9DIR is specified)
endif
ifndef WMAP9LIBFLAGS
$(error Need to specify WMAP9LIBFLAGS if WMAP9DIR is specified)
endif
ifndef CFITSIO
$(error Need to specify CFITSIO if HEALPIX is specified)
endif
ifeq ($(WMAP9COMPILER),gfortran)
WMAP_COMPILE_FLAGS = -D COSMO_WMAP9 -D WMAP9_GFORT
else
WMAP_COMPILE_FLAGS = -D COSMO_WMAP9 -D WMAP9_IFORT
endif
WMAP_INCLUDE_FLAGS = -I $(CFITSIO)/include
WMAP_LIB_FLAGS1 = -L $(CFITSIO)/LIB -L $(WMAP9DIR)
WMAP_LIB_FLAGS2 = $(WMAP9LIBFLAGS) -lcfitsio
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

LFLAGS1 = $(HEALPIX_LIB_FLAGS1) $(LAPACKPP_LIB_FLAGS1) $(CLASS_LIB_FLAGS1) $(MINUIT_LIB_FLAGS1) $(MULTINEST_LIB_FLAGS1) $(PLANCK_LIB_FLAGS1) $(WMAP_LIB_FLAGS1) $(ANN_LIB_FLAGS1) -fopenmp

LFLAGS2 = -lstdc++ $(HEALPIX_LIB_FLAGS2) $(LAPACKPP_LIB_FLAGS2) $(CLASS_LIB_FLAGS2) $(MINUIT_LIB_FLAGS2) $(MULTINEST_LIB_FLAGS2) $(PLANCK_LIB_FLAGS2) $(WMAP_LIB_FLAGS2) $(ANN_LIB_FLAGS2)

MYFLAGS = $(MY_COMPILE_FLAGS) $(HEALPIX_COMPILE_FLAGS) $(LAPACKPP_COMPILE_FLAGS) $(CLASS_COMPILE_FLAGS) $(MINUIT_COMPILE_FLAGS) $(MULTINEST_COMPILE_FLAGS) $(PLANCK_COMPILE_FLAGS) $(WMAP_COMPILE_FLAGS) $(ANN_COMPILE_FLAGS) $(MPI_COMPILE_FLAGS) $(OMP_COMPILE_FLAGS)

INCLUDE_FLAGS = -I include $(HEALPIX_INCLUDE_FLAGS) $(LAPACKPP_INCLUDE_FLAGS) $(CLASS_INCLUDE_FLAGS) $(MINUIT_INCLUDE_FLAGS) $(MULTINEST_INCLUDE_FLAGS) $(PLANCK_INCLUDE_FLAGS) $(WMAP_INCLUDE_FLAGS) $(ANN_INCLUDE_FLAGS)

CFLAGS = $(MYFLAGS) -c $(COMPILER_FLAGS) $(INCLUDE_FLAGS) $(OMP_FLAG)

#Header file tree
MACROS_HPP = include/macros.hpp
EXCEPTION_HANDLER_HPP = include/exception_handler.hpp
COSMO_MPI_HPP = include/cosmo_mpi.hpp include/macros.hpp
PROGRESS_METER_HPP = include/progress_meter.hpp $(MACROS_HPP)
TIMER_HPP = include/timer.hpp $(MACROS_HPP)
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
LEGENDRE_HPP = include/legendre.hpp
SPHERICAL_HARMONICS_HPP = include/spherical_harmonics.hpp $(COMPLEX_TYPES_HPP) $(MATH_CONSTANTS_HPP)
TABLE_FUNCTION_HPP = include/table_function.hpp $(MACROS_HPP) $(FUNCTION_HPP)
CUBIC_SPLINE_HPP = include/cubic_spline.hpp $(MACROS_HPP) $(FUNCTION_HPP)
GAUSS_SMOOTH_HPP = include/gauss_smooth.hpp $(MACROS_HPP) $(FUNCTION_HPP)
THREE_VECTOR_HPP = include/three_vector.hpp $(ANGULAR_COORDINATES_HPP) $(MATH_CONSTANTS_HPP)
THREE_ROTATION_HPP = include/three_rotation.hpp $(MACROS_HPP) $(THREE_VECTOR_HPP) $(NUMERICS_HPP)
MCMC_HPP = include/mcmc.hpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(LIKELIHOOD_FUNCTION_HPP) $(RANDOM_HPP)
MN_SCANNER_HPP = include/mn_scanner.hpp $(LIKELIHOOD_FUNCTION_HPP) $(TABLE_FUNCTION_HPP)
WIGNER_3J_HPP = include/wigner_3j.hpp $(MACROS_HPP)
CONJUGATE_GRADIENT_HPP = include/conjugate_gradient.hpp $(MACROS_HPP)
MARKOV_CHAIN_HPP = include/markov_chain.hpp $(FUNCTION_HPP) $(TABLE_FUNCTION_HPP) $(RANDOM_HPP)
K_NEAREST_NEIGHBORS_HPP = include/k_nearest_neighbors.hpp
FAST_APPROXIMATOR_HPP = include/fast_approximator.hpp $(MACROS_HPP) $(K_NEAREST_NEIGHBORS_HPP) $(TIMER_HPP) $(PROGRESS_METER_HPP)
FAST_APPROXIMATOR_ERROR_HPP = include/fast_approximator_error.hpp $(FUNCTION_HPP) $(FAST_APPROXIMATOR_HPP) $(MARKOV_CHAIN_HPP)
LEARN_AS_YOU_GO_HPP = include/learn_as_you_go.hpp $(MACROS_HPP) $(FUNCTION_HPP) $(RANDOM_HPP) $(FAST_APPROXIMATOR_HPP) $(FAST_APPROXIMATOR_ERROR_HPP)
GAUSSIAN_PROCESS_HPP = include/gaussian_process.hpp

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
MATTER_LIKELIHOOD_HPP = include/matter_likelihood.hpp $(FUNCTION_HPP) $(COSMOLOGICAL_PARAMS_HPP)

TEST_UNIT_CONVERSIONS_HPP = include/test_unit_conversions.hpp $(TEST_FRAMEWORK_HPP)
TEST_INT_OPERATIONS_HPP = include/test_int_operations.hpp $(TEST_FRAMEWORK_HPP)
TEST_INTEGRAL_HPP = include/test_integral.hpp $(TEST_FRAMEWORK_HPP)
TEST_CONJUGATE_GRADIENT_HPP = include/test_conjugate_gradient.hpp $(TEST_FRAMEWORK_HPP)
TEST_POLYNOMIAL_HPP = include/test_polynomial.hpp $(TEST_FRAMEWORK_HPP)
TEST_LEGENDRE_HPP = include/test_legendre.hpp $(TEST_FRAMEWORK_HPP)
TEST_SPHERICAL_HARMONICS_HPP = include/test_spherical_harmonics.hpp $(TEST_FRAMEWORK_HPP)
TEST_MCMC_HPP = include/test_mcmc.hpp $(TEST_FRAMEWORK_HPP)
TEST_MULTINEST_HPP = include/test_multinest.hpp $(TEST_FRAMEWORK_HPP)
TEST_MCMC_PLANCK_HPP = include/test_mcmc_planck.hpp $(TEST_FRAMEWORK_HPP)
TEST_MULTINEST_PLANCK_HPP = include/test_multinest_planck.hpp $(TEST_FRAMEWORK_HPP)
TEST_CMB_HPP = include/test_cmb.hpp $(TEST_FRAMEWORK_HPP)
TEST_CMB_GIBBS_HPP = include/test_cmb_gibbs.hpp $(TEST_FRAMEWORK_HPP)
TEST_FIT_HPP = include/test_fit.hpp $(TEST_FRAMEWORK_HPP)
TEST_PLANCK_LIKE_HPP = include/test_planck_like.hpp $(TEST_FRAMEWORK_HPP)
TEST_WMAP9_LIKE_HPP = include/test_wmap9_like.hpp $(TEST_FRAMEWORK_HPP)
TEST_LIKE_HIGH_HPP = include/test_like_high.hpp $(TEST_FRAMEWORK_HPP)
TEST_LIKE_LOW_HPP = include/test_like_low.hpp $(TEST_FRAMEWORK_HPP)
TEST_WIGNER_3J_HPP = include/test_wigner_3j.hpp $(TEST_FRAMEWORK_HPP)
TEST_TABLE_FUNCTION_HPP = include/test_table_function.hpp $(TEST_FRAMEWORK_HPP)
TEST_CUBIC_SPLINE_HPP = include/test_cubic_spline.hpp $(TEST_FRAMEWORK_HPP)
TEST_THREE_ROTATION_HPP = include/test_three_rotation.hpp $(TEST_FRAMEWORK_HPP)
TEST_MASK_APODIZER_HPP = include/test_mask_apodizer.hpp $(TEST_FRAMEWORK_HPP)
TEST_MATTER_LIKELIHOOD_HPP = include/test_matter_likelihood.hpp $(TEST_FRAMEWORK_HPP)
TEST_K_NEAREST_NEIGHBORS_HPP = include/test_k_nearest_neighbors.hpp $(TEST_FRAMEWORK_HPP)
TEST_FAST_APPROXIMATOR_HPP = include/test_fast_approximator.hpp $(TEST_FRAMEWORK_HPP)
TEST_FAST_APPROXIMATOR_ERROR_HPP = include/test_fast_approximator_error.hpp $(TEST_FRAMEWORK_HPP)
TEST_GAUSSIAN_PROCESS_HPP = include/test_gaussian_process.hpp $(TEST_FRAMEWORKH_HPP)

GENERAL_TARGET =

all: lib/libcosmopp.a bin/test $(GENERAL_TARGET) $(HEALPIX_AND_LAPACKPP_TARGET) $(HEALPIX_TARGET) $(CLASS_TARGET) $(PLANCK_AND_CLASS_TARGET) $(PLANCK_AND_MULTINEST_AND_CLASS_TARGET) $(LAPACKPP_TARGET)

GENERAL_OBJ = obj/macros.o obj/test_framework.o obj/whole_matrix.o obj/scale_factor.o obj/markov_chain.o

OBJ_LIBRARY = $(GENERAL_OBJ) $(HEALPIX_OBJ) $(LAPACKPP_OBJ) $(CLASS_AND_LAPACK_OBJ) $(HEALPIX_AND_LAPACKPP_OBJ) $(CLASS_OBJ) $(MULTINEST_OBJ) $(PLANCK_OBJ) $(PLANCK_AND_CLASS_OBJ) $(WMAP_OBJ) $(WMAP_AND_CLASS_OBJ) $(ANN_OBJ) $(ANN_AND_LAPACK_OBJ) $(MINUIT_AND_LAPACKPP_OBJ)
lib/libcosmopp.a: $(OBJ_LIBRARY)
	ar rcs $@ $(OBJ_LIBRARY)

GENERAL_TEST_OBJ = obj/test.o obj/test_unit_conversions.o obj/test_int_operations.o obj/test_integral.o obj/test_conjugate_gradient.o obj/test_polynomial.o obj/test_legendre.o obj/test_spherical_harmonics.o obj/test_wigner_3j.o obj/test_table_function.o obj/test_cubic_spline.o obj/test_three_rotation.o

OBJ_TEST = $(OBJ_LIBRARY) $(GENERAL_TEST_OBJ) $(MULTINEST_TEST_OBJ) $(PLANCK_TEST_OBJ) $(PLANCK_AND_MULTINEST_AND_CLASS_TEST_OBJ) $(CLASS_TEST_OBJ) $(CLASS_AND_LAPACK_TEST_OBJ) $(PLANCK_AND_CLASS_TEST_OBJ) $(LAPACKPP_TEST_OBJ) $(HEALPIX_TEST_OBJ) $(HEALPIX_AND_LAPACK_AND_CLASS_TEST_OBJ) $(MINUIT_TEST_OBJ) $(WMAP_AND_CLASS_TEST_OBJ) $(PLANCK_AND_CLASS_AND_LAPACK_TEST_OBJ) $(ANN_TEST_OBJ) $(ANN_AND_LAPACK_TEST_OBJ) $(MINUIT_AND_LAPACKPP_TEST_OBJ)
bin/test: $(OBJ_TEST)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_TEST) $(LFLAGS2)

ifdef LAPACKPP_TARGET
OBJ_EXAMPLE_METROPOLIS_HASTINGS = obj/example_metropolis_hastings.o $(OBJ_LIBRARY)
bin/example_metropolis_hastings: $(OBJ_EXAMPLE_METROPOLIS_HASTINGS)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_EXAMPLE_METROPOLIS_HASTINGS) $(LFLAGS2)
endif

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

ifdef CLASS_TARGET
OBJ_EXAMPLE_CL = obj/example_cl.o $(OBJ_LIBRARY)
bin/example_cl: $(OBJ_EXAMPLE_CL)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_EXAMPLE_CL) $(LFLAGS2)
endif

ifdef PLANCK_AND_CLASS_TARGET
OBJ_EXAMPLE_PLANCK = obj/example_planck.o obj/macros.o obj/planck_like.o obj/cmb.o
bin/example_planck: $(OBJ_EXAMPLE_PLANCK)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_EXAMPLE_PLANCK) $(LFLAGS2)
endif

ifdef PLANCK_AND_MULTINEST_AND_CLASS_TARGET
OBJ_EXAMPLE_MN_PLANCK = obj/example_mn_planck.o obj/macros.o obj/planck_like.o obj/cmb.o obj/mn_scanner.o
bin/example_mn_planck: $(OBJ_EXAMPLE_MN_PLANCK)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_EXAMPLE_MN_PLANCK) $(LFLAGS2)
endif

obj/macros.o: source/macros.cpp $(MACROS_HPP)
	$(CC) $(CFLAGS) source/macros.cpp -o $@

obj/c_matrix.o: source/c_matrix.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(UTILS_HPP) $(C_MATRIX_HPP)
	$(CC) $(CFLAGS) source/c_matrix.cpp -o $@

obj/c_matrix_generator.o: source/c_matrix_generator.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(ANGULAR_COORDINATES_HPP) $(THREE_ROTATION_HPP) $(PROGRESS_METER_HPP) $(C_MATRIX_GENERATOR_HPP) $(UTILS_HPP) $(RANDOM_HPP) $(LEGENDRE_HPP)
	$(CC) $(CFLAGS) source/c_matrix_generator.cpp -o $@

obj/cmb.o: source/cmb.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(CMB_HPP) $(TIMER_HPP)
	$(CC) $(CFLAGS) source/cmb.cpp -o $@

obj/likelihood.o: source/likelihood.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(PROGRESS_METER_HPP) $(THREE_ROTATION_HPP) $(LIKELIHOOD_HPP) $(UTILS_HPP)
	$(CC) $(CFLAGS) source/likelihood.cpp -o $@

obj/master.o: source/master.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(PROGRESS_METER_HPP) $(WIGNER_3J_HPP) $(MASTER_HPP)
	$(CC) $(CFLAGS) source/master.cpp -o $@

ifdef CLASS_TARGET
obj/example_cl.o: examples/example_cl.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(CMB_HPP) $(COSMOLOGICAL_PARAMS_HPP)
	$(CC) $(CFLAGS) examples/example_cl.cpp -o $@
endif

ifdef PLANCK_AND_MULTINEST_AND_CLASS_TARGET
obj/example_mn_planck.o: examples/example_mn_planck.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(LIKELIHOOD_FUNCTION_HPP) $(PLANCK_LIKE_HPP) $(MN_SCANNER_HPP)
	$(CC) $(CFLAGS) examples/example_mn_planck.cpp -o $@
endif

ifdef MULTINEST
obj/mn_scanner.o: source/mn_scanner.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(MN_SCANNER_HPP)
	$(CC) $(CFLAGS) source/mn_scanner.cpp -o $@
endif

obj/mode_directions.o: source/mode_directions.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(ANGULAR_COORDINATES_HPP) $(THREE_ROTATION_HPP) $(PROGRESS_METER_HPP) $(MATH_CONSTANTS_HPP) $(NUMERICS_HPP) $(MODE_DIRECTIONS_HPP)
	$(CC) $(CFLAGS) source/mode_directions.cpp -o $@

ifdef PLANCKDIR
obj/planck_like.o: source/planck_like.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(PLANCK_LIKE_HPP) $(TIMER_HPP)
	$(CC) $(CFLAGS) source/planck_like.cpp -o $@
endif

ifdef WMAP9DIR
obj/wmap9_like.o: source/wmap9_like.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(WMAP9_LIKE_HPP)
	$(CC) $(CFLAGS) source/wmap9_like.cpp -o $@
endif

obj/scale_factor.o: source/scale_factor.cpp $(UNIT_CONVERSIONS_HPP) $(TABLE_FUNCTION_HPP) $(SCALE_FACTOR_HPP)
	$(CC) $(CFLAGS) source/scale_factor.cpp -o $@

obj/simulate.o: source/simulate.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(WHOLE_MATRIX_HPP) $(RANDOM_HPP) $(SIMULATE_HPP)
	$(CC) $(CFLAGS) source/simulate.cpp -o $@

obj/test.o: source/test.cpp $(COSMO_MPI_HPP) $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(TEST_FRAMEWORK_HPP) $(TEST_UNIT_CONVERSIONS_HPP) $(TEST_INT_OPERATIONS_HPP) $(TEST_INTEGRAL_HPP) $(TEST_CONJUGATE_GRADIENT_HPP) $(TEST_POLYNOMIAL_HPP) $(TEST_LEGENDRE_HPP) $(TEST_SPHERICAL_HARMONICS_HPP) $(TEST_MCMC_HPP) $(TEST_MULTINEST_HPP) $(TEST_MCMC_PLANCK_HPP) $(TEST_CMB_HPP) $(TEST_CMB_GIBBS_HPP) $(TEST_FIT_HPP) $(TEST_PLANCK_LIKE_HPP) $(TEST_WMAP9_LIKE_HPP) $(TEST_LIKE_HIGH_HPP) $(TEST_LIKE_LOW_HPP) $(TEST_WIGNER_3J_HPP) $(TEST_TABLE_FUNCTION_HPP) $(TEST_CUBIC_SPLINE_HPP) $(TEST_THREE_ROTATION_HPP) $(TEST_MASK_APODIZER_HPP) $(TEST_MATTER_LIKELIHOOD_HPP) $(TEST_K_NEAREST_NEIGHBORS_HPP) $(TEST_FAST_APPROXIMATOR_HPP) $(TEST_FAST_APPROXIMATOR_ERROR_HPP) $(TEST_GAUSSIAN_PROCESS_HPP)
	$(CC) $(CFLAGS) source/test.cpp -o $@

obj/test_framework.o: source/test_framework.cpp $(COSMO_MPI_HPP) $(MACROS_HPP) $(NUMERICS_HPP) $(TEST_FRAMEWORK_HPP)
	$(CC) $(CFLAGS) source/test_framework.cpp -o $@

ifdef PLANCK_AND_CLASS_TARGET
obj/example_planck.o: examples/example_planck.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(TIMER_HPP) $(PLANCK_LIKE_HPP)
	$(CC) $(CFLAGS) examples/example_planck.cpp -o $@
endif

obj/utils.o: source/utils.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(THREE_VECTOR_HPP) $(UTILS_HPP)
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

ifdef CLASS_AND_LAPACK_OBJ
obj/matter_likelihood.o: source/matter_likelihood.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(SCALE_FACTOR_HPP) $(UNIT_CONVERSIONS_HPP) $(MATTER_LIKELIHOOD_HPP)
	$(CC) $(CFLAGS) source/matter_likelihood.cpp -o $@
endif

obj/mcmc.o: source/mcmc.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MCMC_HPP)
	$(CC) $(CFLAGS) source/mcmc.cpp -o $@

obj/markov_chain.o: source/markov_chain.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(CUBIC_SPLINE_HPP) $(GAUSS_SMOOTH_HPP) $(PROGRESS_METER_HPP) $(MARKOV_CHAIN_HPP)
	$(CC) $(CFLAGS) source/markov_chain.cpp -o $@

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

obj/test_legendre.o: source/test_legendre.cpp $(MACROS_HPP) $(TEST_LEGENDRE_HPP) $(LEGENDRE_HPP)
	$(CC) $(CFLAGS) source/test_legendre.cpp -o $@

obj/test_spherical_harmonics.o: source/test_spherical_harmonics.cpp $(MACROS_HPP) $(TEST_SPHERICAL_HARMONICS_HPP) $(SPHERICAL_HARMONICS_HPP)
	$(CC) $(CFLAGS) source/test_spherical_harmonics.cpp -o $@

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

obj/test_planck_like.o: source/test_planck_like.cpp $(PLANCK_LIKE_HPP) $(TEST_PLANCK_LIKE_HPP)
	$(CC) $(CFLAGS) source/test_planck_like.cpp -o $@

obj/test_wmap9_like.o: source/test_wmap9_like.cpp $(WMAP9_LIKE_HPP) $(NUMERICS_HPP) $(TEST_WMAP9_LIKE_HPP)
	$(CC) $(CFLAGS) source/test_wmap9_like.cpp -o $@

obj/test_like_high.o: source/test_like_high.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(UTILS_HPP) $(MASK_APODIZER_HPP) $(RANDOM_HPP) $(HISTOGRAM_HPP) $(LIKELIHOOD_HPP) $(CMB_HPP) $(MASTER_HPP) $(SIMULATE_HPP) $(PROGRESS_METER_HPP) $(CHI_SQUARED_HPP) $(TEST_LIKE_HIGH_HPP)
	$(CC) $(CFLAGS) source/test_like_high.cpp -o $@

obj/test_like_low.o: source/test_like_low.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(UTILS_HPP) $(RANDOM_HPP) $(HISTOGRAM_HPP) $(C_MATRIX_GENERATOR_HPP) $(LIKELIHOOD_HPP) $(CMB_HPP) $(SIMULATE_HPP) $(PROGRESS_METER_HPP) $(CHI_SQUARED_HPP) $(TEST_LIKE_LOW_HPP)
	$(CC) $(CFLAGS) source/test_like_low.cpp -o $@

obj/test_wigner_3j.o: source/test_wigner_3j.cpp $(WIGNER_3J_HPP) $(TEST_WIGNER_3J_HPP)
	$(CC) $(CFLAGS) source/test_wigner_3j.cpp -o $@

obj/test_table_function.o: source/test_table_function.cpp $(TABLE_FUNCTION_HPP) $(TEST_TABLE_FUNCTION_HPP)
	$(CC) $(CFLAGS) source/test_table_function.cpp -o $@

obj/test_cubic_spline.o: source/test_cubic_spline.cpp $(CUBIC_SPLINE_HPP) $(TEST_CUBIC_SPLINE_HPP)
	$(CC) $(CFLAGS) source/test_cubic_spline.cpp -o $@

obj/test_three_rotation.o: source/test_three_rotation.cpp $(THREE_ROTATION_HPP) $(MATH_CONSTANTS_HPP) $(TEST_THREE_ROTATION_HPP)
	$(CC) $(CFLAGS) source/test_three_rotation.cpp -o $@

obj/test_mask_apodizer.o: source/test_mask_apodizer.cpp $(MACROS_HPP) $(MASK_APODIZER_HPP) $(MATH_CONSTANTS_HPP) $(TEST_MASK_APODIZER_HPP)
	$(CC) $(CFLAGS) source/test_mask_apodizer.cpp -o $@

obj/test_matter_likelihood.o: source/test_matter_likelihood.cpp $(MACROS_HPP) $(CMB_HPP) $(MATTER_LIKELIHOOD_HPP) $(TEST_MATTER_LIKELIHOOD_HPP)
	$(CC) $(CFLAGS) source/test_matter_likelihood.cpp -o $@

obj/example_metropolis_hastings.o: examples/example_metropolis_hastings.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MCMC_HPP) $(MARKOV_CHAIN_HPP)
	$(CC) $(CFLAGS) examples/example_metropolis_hastings.cpp -o $@

obj/k_nearest_neighbors.o: source/k_nearest_neighbors.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(TIMER_HPP) $(K_NEAREST_NEIGHBORS_HPP)
	$(CC) $(CFLAGS) source/k_nearest_neighbors.cpp -o $@

obj/test_k_nearest_neighbors.o: source/test_k_nearest_neighbors.cpp $(MACROS_HPP) $(K_NEAREST_NEIGHBORS_HPP) $(NUMERICS_HPP) $(TEST_K_NEAREST_NEIGHBORS_HPP)
	$(CC) $(CFLAGS) source/test_k_nearest_neighbors.cpp -o $@

obj/fast_approximator.o: source/fast_approximator.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(GAUSSIAN_PROCESS_HPP) $(FAST_APPROXIMATOR_HPP)
	$(CC) $(CFLAGS) source/fast_approximator.cpp -o $@

obj/fast_approximator_error.o: source/fast_approximator.cpp $(FAST_APPROXIMATOR_ERROR_HPP)
	$(CC) $(CFLAGS) source/fast_approximator_error.cpp -o $@

obj/learn_as_you_go.o: source/learn_as_you_go.cpp $(COSMO_MPI_HPP) $(LEARN_AS_YOU_GO_HPP)
	$(CC) $(CFLAGS) source/learn_as_you_go.cpp -o $@

obj/test_fast_approximator.o: source/test_fast_approximator.cpp $(RANDOM_HPP) $(FAST_APPROXIMATOR_HPP) $(TEST_FAST_APPROXIMATOR_HPP)
	$(CC) $(CFLAGS) source/test_fast_approximator.cpp -o $@

obj/test_fast_approximator_error.o: source/test_fast_approximator_error.cpp $(RANDOM_HPP) $(FAST_APPROXIMATOR_ERROR_HPP) $(TEST_FAST_APPROXIMATOR_ERROR_HPP)
	$(CC) $(CFLAGS) source/test_fast_approximator_error.cpp -o $@

obj/gaussian_process.o: source/gaussian_process.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(TIMER_HPP) $(GAUSSIAN_PROCESS_HPP)
	$(CC) $(CFLAGS) source/gaussian_process.cpp -o $@

obj/test_gaussian_process.o: source/test_gaussian_process.cpp $(MACROS_HPP) $(GAUSSIAN_PROCESS_HPP) $(RANDOM_HPP) $(MARKOV_CHAIN_HPP) $(TEST_GAUSSIAN_PROCESS_HPP)
	$(CC) $(CFLAGS) source/test_gaussian_process.cpp -o $@

clean:
	rm obj/* bin/* lib/*

test: all
	./bin/test fast
