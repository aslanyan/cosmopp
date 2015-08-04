#include <macros.hpp>
#include <exception_handler.hpp>
#include <markov_chain.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 2)
        {
            std::string exceptionString = "Need to specify the chain file name.";
            exc.set(exceptionString);
            throw exc;
        }

        MarkovChain chain(argv[1]);

        output_screen("95\% confidence limits:" << std::endl);

        for(int i = 0; i < 5; ++i)
        {
            std::stringstream fileName;
            fileName << "post_v_" << i << ".txt";
            Posterior1D *p = chain.posterior(i + 4, Posterior1D::GAUSSIAN_SMOOTHING);

            p->writeIntoFile(fileName.str().c_str());

            const double median = p->median();
            double lower, upper;
            p->get2SigmaTwoSided(lower, upper);
            const double sigma = (upper - lower) / 2.0;
            output_screen("v" << i << " = " << median << " + " << upper - median << " - " << median - lower << std::endl);
        }

        return 0;
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
