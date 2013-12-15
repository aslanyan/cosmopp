#ifndef COSMO_PP_PROGRESS_METER_HPP
#define COSMO_PP_PROGRESS_METER_HPP

#include <macros.hpp>

/// Progress Meter class.

/// This class helps monitor the progress of certain operations and outputs the percentage completed on the screen.
class ProgressMeter
{
public:
    /// Constructor.

    /// Constructs the progress meter.
    /// \param total The total number of operations to be performed.
	ProgressMeter(unsigned long total) : total_(total), completed_(0), previous_(0) { check(total_ > 0, ""); output_screen(previous_ << '%' << std::flush); }

    /// Advace the meter.

    /// This function advances the progress meter by a given amount of step (default is 1).
    /// \param delta The step size to advance the meter.
	void advance(unsigned long delta = 1)
	{
		completed_ += delta;
		check(completed_ <= total_, "called too many times");
		unsigned progress = 100 * completed_ / total_;
		if(progress == previous_)
			return;
		previous_ = progress;
		output_screen('\r' << previous_ << '%' << std::flush);
		
		if(previous_ == 100)
			output_screen(std::endl);
	}
	
private:
	unsigned long total_;
	unsigned long completed_;
	unsigned previous_;
};

#endif
