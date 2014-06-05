#ifndef COSMO_PP_MACROS_HPP
#define COSMO_PP_MACROS_HPP

// Possible pre-defined macros from the compiler (e.g. gcc -D MACRO)
// CHECKS_ON - turn on the checks (assertions)
// SLOW_CHECKS_ON turn on all the checks, including the slow ones
// CREATE_LOG create the log file log.txt and output the stuff from output_log(...)
// VERBOSE turn on verbose mode, i.e. print on the screen stuff from output_screen(...)

#ifdef SLOW_CHECKS_ON
#ifndef CHECKS_ON
#define CHECKS_ON
#endif
#endif

#ifdef VERBOSE2
#ifndef VERBOSE1
#define VERBOSE1
#endif
#endif

#ifdef VERBOSE1
#ifndef VERBOSE
#define VERBOSE
#endif
#endif

#ifdef CREATE_LOG2
#ifndef CREATE_LOG1
#define CREATE_LOG1
#endif
#endif

#ifdef CREATE_LOG1
#ifndef CREATE_LOG
#define CREATE_LOG
#endif
#endif

bool IS_PARALLEL();
int CURRENT_PROCESS();
int NUM_PROCESSES();

int TOTAL_NUM_THREADS();
int CURRENT_THREAD_NUM();

#ifdef CREATE_LOG
#include <fstream>
#include <string>
#include <sstream>

#ifdef COSMO_OMP
#include <omp.h>
#endif

class LOG984375987 : public std::ofstream
{
private:
	inline LOG984375987();
	inline ~LOG984375987()
    {
#ifdef COSMO_OMP
        omp_destroy_lock(&lock_);
#endif
    }
public:
	inline static LOG984375987& create();
    inline static std::string fileName() { std::stringstream str; str << "log"; if(IS_PARALLEL()) str << "_" << CURRENT_PROCESS() << "_of_" << NUM_PROCESSES(); str << ".txt"; return str.str(); }

#ifdef COSMO_OMP
    inline void setLock() { omp_set_lock(&lock_); }
    inline void unsetLock() { omp_unset_lock(&lock_); }
private:
    omp_lock_t lock_;
#else
    inline void setLock() {}
    inline void unsetLock() {}
#endif
};

LOG984375987& LOG984375987::create()
{
	static LOG984375987 object;
	return object;
}

LOG984375987::LOG984375987() : std::ofstream(fileName().c_str())
{
#ifdef COSMO_OMP
    omp_init_lock(&lock_);
#endif
    close();
}

#define output_log(A) { const int tnt = TOTAL_NUM_THREADS(); if(tnt > 1) {  LOG984375987::create().setLock(); } LOG984375987::create().open(LOG984375987::fileName().c_str(), std::ios::app); if(tnt > 1) LOG984375987::create() << "Thread " << CURRENT_THREAD_NUM() << " / " << tnt << ":\t"; LOG984375987::create() << A; LOG984375987::create().close(); if(tnt > 1) { LOG984375987::create().unsetLock(); }}
#else
#define output_log(A)
#endif

#ifdef CREATE_LOG1
#define output_log1(A) output_log(A)
#else
#define output_log1(A)
#endif

#ifdef CREATE_LOG2
#define output_log2(A) output_log(A)
#else
#define output_log2(A)
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
#ifdef COSMO_OMP
#include <omp.h>
#endif
class OS23459783987
{
private:
    inline OS23459783987()
    {
#ifdef COSMO_OMP
        omp_init_lock(&lock_);
#endif
    }

    inline ~OS23459783987()
    {
#ifdef COSMO_OMP
        omp_destroy_lock(&lock_);
#endif
    }
    
public:
    inline static OS23459783987& create()
    {
        static OS23459783987 obj;
        return obj;
    }

#ifdef COSMO_OMP
    inline void setLock() { omp_set_lock(&lock_); }
    inline void unsetLock() { omp_unset_lock(&lock_); }
private:
    omp_lock_t lock_;
#else
    inline void setLock() {}
    inline void unsetLock() {}
#endif
};

#define output_screen_clean(A) { const int tnt = TOTAL_NUM_THREADS(); if(tnt > 1) { OS23459783987::create().setLock(); } std::cout << A; if(tnt > 1) { OS23459783987::create().unsetLock(); } }

#define output_screen(A) { const int tnt = TOTAL_NUM_THREADS(); if(tnt > 1) { OS23459783987::create().setLock(); } if(IS_PARALLEL()) std::cout << "MPI process " << CURRENT_PROCESS() << " / " << NUM_PROCESSES() << ":\t"; if(tnt > 1) std::cout << "OMP thread " << CURRENT_THREAD_NUM() << " / " << tnt << ":\t"; std::cout << A; if(tnt > 1) { OS23459783987::create().unsetLock(); }}
#else
#define output_screen_clean(A)
#define output_screen(A)
#endif

#ifdef VERBOSE1
#define output_screen_clean1(A) output_screen_clean(A)
#define output_screen1(A)
#else
#define output_screen_clean1(A) output_log1(A)
#define output_screen1(A)
#endif

#ifdef VERBOSE2
#define output_screen_clean2(A) output_screen_clean(A)
#define output_screen2(A) output_screen(A)
#else
#define output_screen_clean2(A)
#define output_screen2(A)
#endif

#endif
