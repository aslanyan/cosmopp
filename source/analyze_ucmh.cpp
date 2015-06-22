#include <fstream>
#include <sstream>
#include <memory>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 3)
        {
            exc.set("Need to specify the chain file and the output contour file.");
            throw exc;
        }
        MarkovChain chain(argv[1]);
        Posterior2D* cont = chain.posterior(5, 6);
        cont->writeIntoFile(argv[2]);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
