include make.inc

CC=gcc

PLANCK_AND_MULTINEST_TARGET = bin/mn_scan_planck

ifdef MINUIT
MINUIT_INCLUDE_FLAGS = -I $(MINUIT)/include
MINUIT_LIB_FLAGS1 = -L $(MINUIT)/lib
MINUIT_LIB_FLAGS2 = -lMinuit2
else
MINUIT_INCLUDE_FLAGS =
MINUIT_LIB_FLAGS1 =
MINUIT_LIB_FLAGS2 =
endif

ifdef MULTINEST
MULTINEST_INCLUDE_FLAGS = -I $(MULTINEST)/includes
MULTINEST_LIB_FLAGS1 = -L $(MULTINEST)/lib
MULTINEST_LIB_FLAGS2 = -lmultinest_mpi
MULTINEST_OBJ = obj/mn_scanner.o
else
MULTINEST_INCLUDE_FLAGS =
MULTINEST_LIB_FLAGS1 =
MULTINEST_LIB_FLAGS2 =
MULTINEST_OBJ =
PLANCK_AND_MULTINEST_TARGET =
endif

ifdef PLANCKDIR
PLANCK_COMPILE_FLAGS = -D LAPACK_CLIK -D PLANCK_DATA_DIR=$(PLANCKDATADIR)
PLANCK_INCLUDE_FLAGS = -I $(PLANCKDIR)/include
PLANCK_LIB_FLAGS1 = -L $(PLANCKDIR)/lib
PLANCK_LIB_FLAGS2 = -lclik
PLANCK_OBJ = obj/planck_like.o
PLANCK_TARGET = bin/test_planck
else
PLANCK_COMPILE_FLAGS =
PLANCK_INCLUDE_FLAGS =
PLANCK_LIB_FLAGS1 =
PLANCK_LIB_FLAGS2 =
PLANCK_OBJ =
PlANCK_TARGET =
PLANCK_AND_MULTINEST_TARGET =
endif

LFLAGS1 = -L $(CFITSIO)/lib -L $(HEALPIX)/lib -L $(HEALPIXPP)/lib -L $(LAPACKPPLIBDIR) -fopenmp -L $(CLASSDIR) $(MINUIT_LIB_FLAGS1) $(MULTINEST_LIB_FLAGS1) $(PLANCK_LIB_FLAGS1)

LFLAGS2 = -lstdc++ -lchealpix -lhealpix_cxx -lcxxsupport -lfftpack -lc_utils -lcfitsio -llapackpp -lclass $(MINUIT_LIB_FLAGS2) $(MULTINEST_LIB_FLAGS2) $(PLANCK_LIB_FLAGS2)

MYFLAGS = $(MY_COMPILE_FLAGS) -D HEALPIX_DATA_DIR=$(HEALPIX)/data $(PLANCK_COMPILE_FLAGS)

CFLAGS=$(MYFLAGS) -c -g -O2 -fpic -I include -I $(CFITSIO)/include -I $(HEALPIX)/include -I $(HEALPIXPP)/include -I $(BOOSTDIR) -I $(LAPACKPPINCDIR) -I $(CLASSDIR)/include $(MINUIT_INCLUDE_FLAGS) $(MULTINEST_INCLUDE_FLAGS) $(PLANCK_INCLUDE_FLAGS) -fopenmp

#Header file tree
MACROS_HPP = include/macros.hpp
EXCEPTION_HANDLER_HPP = include/exception_handler.hpp
PROGRESS_METER_HPP = include/progress_meter.hpp $(MACROS_HPP)

ANGULAR_COORDINATES_HPP = include/angular_coordinates.hpp
COMPLEX_TYPES_HPP = include/complex_types.hpp
MATH_CONSTANTS_HPP = include/math_constants.hpp
PHYS_CONSTANTS_HPP = include/phys_constants.hpp $(MATH_CONSTANTS_HPP)
UNIT_CONVERSIONS_HPP = include/unit_conversions.hpp $(PHYS_CONSTANTS_HPP)
FUNCTION_HPP = include/function.hpp
HISTOGRAM_HPP = include/histogram.hpp
INT_OPERATION_HPP = include/int_operations.hpp
LIKELIHOOD_FUNCTION_HPP = include/likelihood_function.hpp
NUMERICS_HPP = include/numerics.hpp
INTEGRAL_HPP = include/integral.hpp $(MACROS_HPP) $(FUNCTION_HPP)
PARAMETRIC_FUNCTION_HPP = include/parametric_function.hpp $(MACROS_HPP) $(FUNCTION_HPP)
FIT_HPP = include/fit.hpp $(MACROS_HPP) $(PARAMETRIC_FUNCTION_HPP)
POLYNOMIAL_HPP = include/polynomial.hpp $(PARAMETRIC_FUNCTION_HPP)
TABLE_FUNCTION_HPP = include/table_function.hpp $(MACROS_HPP) $(FUNCTION_HPP)
CUBIC_SPLINE_HPP = include/cubic_spline.hpp $(MACROS_HPP) $(FUNCTION_HPP)
THREE_VECTOR_HPP = include/three_vector.hpp $(ANGULAR_COORDINATES_HPP) $(MATH_CONSTANTS_HPP)
THREE_ROTATION_HPP = include/three_rotation.hpp $(MACROS_HPP) $(THREE_VECTOR_HPP) $(NUMERICS_HPP)
MCMC_HPP = include/mcmc.hpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(LIKELIHOOD_FUNCTION_HPP)
MN_SCANNER_HPP = include/mn_scanner.hpp $(LIKELIHOOD_FUNCTION_HPP) $(TABLE_FUNCTION_HPP)
WIGNER_3J_HPP = include/wigner_3j.hpp $(MACROS_HPP)
RANDOM_HPP = include/random.hpp
SPIN_SPHERICAL_HARMONICS_HPP = include/spin_spherical_harmonics.hpp $(MACROS_HPP) $(ANGULAR_COORDINATES_HPP) $(COMPLEX_TYPES_HPP) $(MATH_CONSTANTS_HPP)
CONJUGATE_GRADIENT_HPP = include/conjugate_gradient.hpp $(MACROS_HPP)

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
SCALE_FACTOR_HPP = include/scale_factor.hpp $(MACROS_HPP) $(TABLE_FUNCTION_HPP) $(COSMOLOGICAL_PARAMS_HPP)
CMG_GIBBS_HPP = include/cmb_gibbs.hpp
MASK_APODIZER_HPP = include/mask_apodizer.hpp


all: lib/libcosmopp.a bin/sort_chain bin/test bin/generate_white_noise bin/apodize_mask $(PLANCK_TARGET) $(PLANCK_AND_MULTINEST_TARGET)


OBJ_LIBRARY = obj/whole_matrix.o obj/utils.o obj/c_matrix.o obj/c_matrix_generator.o obj/simulate.o obj/likelihood.o obj/master.o obj/mode_directions.o obj/scale_factor.o obj/cmb.o obj/cmb_gibbs.o obj/mask_apodizer.o $(MULTINEST_OBJ) $(PLANCK_OBJ) 
lib/libcosmopp.a: $(OBJ_LIBRARY)
	ar rcs $@ $(OBJ_LIBRARY)

OBJ_TEST = obj/test.o obj/scale_factor.o obj/cmb.o
bin/test: $(OBJ_TEST)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_TEST) $(LFLAGS2)

OBJ_GENERATE_WHITE_NOISE = obj/generate_white_noise.o obj/simulate.o obj/whole_matrix.o
bin/generate_white_noise: $(OBJ_GENERATE_WHITE_NOISE)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_GENERATE_WHITE_NOISE) $(LFLAGS2)

OBJ_APODIZE_MASK = obj/apodize_mask.o obj/mask_apodizer.o
bin/apodize_mask: $(OBJ_APODIZE_MASK)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_APODIZE_MASK) $(LFLAGS2)

ifdef PLANCKDIR
OBJ_TEST_PLANCK = obj/test_planck.o obj/planck_like.o obj/cmb.o
bin/test_planck: $(OBJ_TEST_PLANCK)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_TEST_PLANCK) $(LFLAGS2)
endif

ifdef PLANCK_AND_MULTINEST_TARGET
OBJ_MN_SCAN_PLANCK = obj/mn_scan_planck.o obj/planck_like.o obj/cmb.o obj/mn_scanner.o
bin/mn_scan_planck: $(OBJ_MN_SCAN_PLANCK)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_MN_SCAN_PLANCK) $(LFLAGS2)
endif

OBJ_SORT_CHAIN = obj/sort_chain.o
bin/sort_chain: $(OBJ_SORT_CHAIN)
	$(CC) $(LFLAGS1) -o $@ $(OBJ_SORT_CHAIN) $(LFLAGS2)



obj/c_matrix.o: source/c_matrix.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(UTILS_HPP) $(C_MATRIX_HPP)
	$(CC) $(CFLAGS) source/c_matrix.cpp -o $@

obj/c_matrix_generator.o: source/c_matrix_generator.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(ANGULAR_COORDINATES_HPP) $(THREE_ROTATION_HPP) $(PROGRESS_METER_HPP) $(C_MATRIX_GENERATOR_HPP) $(UTILS_HPP)
	$(CC) $(CFLAGS) source/c_matrix_generator.cpp -o $@

obj/cmb.o: source/cmb.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(CMB_HPP)
	$(CC) $(CFLAGS) source/cmb.cpp -o $@

obj/likelihood.o: source/likelihood.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(PROGRESS_METER_HPP) $(THREE_ROTATION_HPP) $(LIKELIHOOD_HPP) $(UTILS_HPP)
	$(CC) $(CFLAGS) source/likelihood.cpp -o $@

obj/master.o: source/master.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(PROGRESS_METER_HPP) $(WIGNER_3J_HPP) $(MASTER_HPP)
	$(CC) $(CFLAGS) source/master.cpp -o $@

ifdef PLANCK_AND_MULTINEST_TARGET
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
	$(CC) $(CFLAGS) source/planck_like.cpp -o $@
endif

obj/scale_factor.o: source/scale_factor.cpp $(UNIT_CONVERSIONS_HPP) $(TABLE_FUNCTION_HPP) $(SCALE_FACTOR_HPP)
	$(CC) $(CFLAGS) source/scale_factor.cpp -o $@

obj/simulate.o: source/simulate.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(WHOLE_MATRIX_HPP) $(SIMULATE_HPP)
	$(CC) $(CFLAGS) source/simulate.cpp -o $@

obj/sort_chain.o: source/sort_chain.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP)
	$(CC) $(CFLAGS) source/sort_chain.cpp -o $@

obj/test.o: source/test.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(PHYS_CONSTANTS_HPP) $(SCALE_FACTOR_HPP) $(COSMOLOGICAL_PARAMS_HPP) $(CMB_HPP)
	$(CC) $(CFLAGS) source/test.cpp -o $@

ifdef PLANCKDIR
obj/test_planck.o: source/test_planck.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(PLANCK_LIKE_HPP)
	$(CC) $(CFLAGS) source/test_planck.cpp -o $@
endif

obj/utils.o: source/utils.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(UTILS_HPP)
	$(CC) $(CFLAGS) source/utils.cpp -o $@

obj/whole_matrix.o: source/whole_matrix.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(NUMERICS_HPP) $(THREE_ROTATION_HPP) $(WHOLE_MATRIX_HPP)
	$(CC) $(CFLAGS) source/whole_matrix.cpp -o $@

obj/cmb_gibbs.o: source/cmb_gibbs.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MATH_CONSTANTS_HPP) $(PROGRESS_METER_HPP) $(UTILS_HPP) $(CONJUGATE_GRADIENT_HPP) $(NUMERICS_HPP) $(CMB_GIBBS_HPP)
	$(CC) $(CFLAGS) source/cmb_gibbs.cpp -o $@

obj/generate_white_noise.o: source/generate_white_noise.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(SIMULATE_HPP)
	$(CC) $(CFLAGS) source/generate_white_noise.cpp -o $@

obj/mask_apodizer.o: source/mask_apodizer.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(PROGRESS_METER_HPP) $(THREE_VECTOR_HPP) $(MATH_CONSTANTS_HPP) $(MASK_APODIZER_HPP)
	$(CC) $(CFLAGS) source/mask_apodizer.cpp -o $@

obj/apodize_mask.o: source/apodize_mask.cpp $(MACROS_HPP) $(EXCEPTION_HANDLER_HPP) $(MASK_APODIZER_HPP)
	$(CC) $(CFLAGS) source/apodize_mask.cpp -o $@

clean:
	rm obj/* bin/* lib/*
