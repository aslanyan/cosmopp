#include <fstream>
#include <map>
#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 4)
        {
            std::string exceptionStr = "Input chain, output 1 sigma file, and output 2 sigma file (excludes 1 sigma) must be specified.";
            exc.set(exceptionStr);
            throw exc;
        }

        std::map<double, std::string> chain;
        unsigned long count;

        std::ifstream in(argv[1]);
        if(!in)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Input chain file " << argv[1] << " cannot be read.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        while(!in.eof())
        {
            std::string s;
            std::getline(in, s);
            if(s == "")
                break;

            double dummy, like;
            std::stringstream str(s);
            str >> dummy >> like;
            chain[like] = s;
            ++count;
        }

        in.close();

        const unsigned long count1 = (unsigned long) (0.683 * double(count)), count2 = (unsigned long) (0.955 * double(count));

        check(count2 > count1, "");
        check(count2 < count, "");

        std::map<double, std::string>::const_iterator it = chain.begin();
        unsigned long i = 0;
        
        std::ofstream out(argv[2]);
        if(!out)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot write into output file " << argv[2] << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        while(i < count1)
        {
            out << (*it).second << std::endl;
            ++it;
            ++i;
        }

        out.close();
        out.open(argv[3]);
        if(!out)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot write into output file " << argv[3] << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        while(i < count2)
        {
            out << (*it).second << std::endl;
            ++it;
            ++i;
        }

        out.close();
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

