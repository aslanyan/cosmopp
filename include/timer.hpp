#ifndef COSMO_PP_TIMER_HPP
#define COSMO_PP_TIMER_HPP

#include <string>
#include <chrono>

#include <macros.hpp>

class Timer
{
public:
    inline Timer(const std::string& name = "") : name_(name) 
#ifdef CHECKS_ON
        , started_(false)
#endif
    {}

    inline void start()
    {
        check(!started_, "the timer has already been started");
#ifdef CHECKS_ON
        started_ = true;
#endif
        start_ = std::chrono::steady_clock::now();
        output_screen1("Starting the timer " << name_ << "..." << std::endl);
        //output_log("Starting the timer " << name_ << "..." << std::endl);
    }

    inline unsigned long end()
    {
        check(started_, "the timer has not been started");
#ifdef CHECKS_ON
        started_ = false;
#endif
        std::chrono::steady_clock::time_point end_ = std::chrono::steady_clock::now();
        const unsigned long c = std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count();
        output_screen1("The timer " << name_ << " has ended." << std::endl << "Duration: " << c << " microseconds." << std::endl);
        //output_log("The timer " << name_ << " has ended." << std::endl << "Duration: " << c << " microseconds." << std::endl);
        return c;
    }

private:
    std::string name_;
    std::chrono::steady_clock::time_point start_;
#ifdef CHECKS_ON
    bool started_;
#endif
};

#endif

