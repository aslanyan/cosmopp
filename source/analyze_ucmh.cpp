#include <memory>

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

        output_screen("68\% confidence limits:" << std::endl);

        for(int i = 0; i < 9; ++i)
        {
            if(i == 4 || i == 5)
                continue;
            std::string paramName;
            switch(i)
            {
            case 0:
                paramName = "epsilon";
                break;
            case 1:
                paramName = "eta";
                break;
            case 2:
                paramName = "xi";
                break;
            case 3:
                paramName = "omega";
                break;
            case 6:
                paramName = "ns";
                break;
            case 7:
                paramName = "alpha_s";
                break;
            case 8:
                paramName = "alpha_s_prime";
                break;
            default:
                check(false, "");
                break;
            }

            std::unique_ptr<Posterior1D> p(chain.posterior(i + 4, Posterior1D::GAUSSIAN_SMOOTHING));

            const double median = p->median();
            double lower, upper;
            p->get1SigmaTwoSided(lower, upper);
            const double sigma = (upper - lower) / 2.0;
            output_screen(paramName << " = " << median << " + " << upper - median << " - " << median - lower << std::endl);
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
