#ifndef COSMO_PP_MACROS_HPP
#define COSMO_PP_MACROS_HPP

// Possible pre-defined macros from the compiler (e.g. gcc -D MACRO)
// CHECKS_ON - turn on the checks (assertions)
// SLOW_CHECKS_ON turn on all the checks, including the slow ones
// CREATE_LOG create the log file log.txt and output the stuff from output_log(...)
// VERBOSE turn on verbose mode, i.e. print on the screen stuff from output_screen(...)

#ifdef SLOW_CHECKS_ON
#define CHECKS_ON
#endif

#ifdef CREATE_LOG
#include <fstream>

class LOG984375987 : public std::ofstream
{
private:
	inline LOG984375987();
	inline ~LOG984375987() {}
public:
	inline static LOG984375987& create();
};

LOG984375987& LOG984375987::create()
{
	static LOG984375987 object;
	return object;
}

LOG984375987::LOG984375987() : std::ofstream("log.txt")
{
}

#define output_log(A) LOG984375987::create() << A; LOG984375987::create().close(); LOG984375987::create().open("log.txt", std::ios::app);
#else
#define output_log(A)
#endif


#ifdef CHECKS_ON

#include <iostream>
#include <cstdlib>
#include <exception_handler.hpp>

#define check(a, b) {if(!(a)){ \
std::cout << "CHECK FAILED in " << __FILE__ << " on line " << __LINE__ << "!!! " << b << std::endl << "Stopping execution!!!" << std::endl; \
output_log("CHECK FAILED in " << __FILE__ << " on line " << __LINE__ << "!!! " << b << std::endl << "Stopping execution!!!" << std::endl); \
StandardException exc; \
std::string exceptionStr = "CHECK FAILED"; \
exc.set(exceptionStr); \
throw exc;}}

#else
#define check(a, b)
#endif

#ifdef SLOW_CHECKS_ON
#define slow_check(a, b) check(a, b)
#else
#define slow_check(a, b)
#endif

#ifdef VERBOSE
#include <iostream>
#define output_screen(A) std::cout << A; \
output_log(A);
#else
#define output_screen(A) output_log(A);
#endif

#endif
