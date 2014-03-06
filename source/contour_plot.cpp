#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <macros.hpp>
#include <exception_handler.hpp>

struct ChainElement
{
    ChainElement() : p(0), like(0)
    {
    }

    ChainElement(const ChainElement& other) : p(other.p), like(other.like), params(other.params)
    {
    }

    double p;
    double like;
    std::vector<double> params;

    bool operator < (const ChainElement& other) const
    {
        return like < other.like;
    }
};

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 5)
        {
            std::string exceptionStr = "Input chain file, the indices of the two parameters for the contour plot (starting from 1), the output file, and the posterior resolution (i.e. number of points, optional, default = 100) must be specified.";
            exc.set(exceptionStr);
            throw exc;
        }

        std::stringstream paramsStr;
        paramsStr << argv[2] << ' ' << argv[3];
        int p1, p2;
        paramsStr >> p1 >> p2;

        if(p1 <= 0)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid parameter index " << argv[2] << ". Needs to be a positive integer.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(p2 <= 0)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid parameter index " << argv[3] << ". Needs to be a positive integer.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        int res = 100;
        if(argc > 5)
        {
            std::stringstream resStr;
            resStr << argv[5];
            resStr >> res;

            if(res <= 0)
            {
                std::stringstream exceptionStr;
                exceptionStr << "Invalid resolution " << argv[5] << ". Needs to be a positive integer.";
                exc.set(exceptionStr.str());
                throw exc;
            }
        }

        std::vector<ChainElement> chain;

        std::ifstream in(argv[1]);
        if(!in)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Input chain file " << argv[1] << " cannot be read.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        const int n = (p1 > p2 ? p1 : p2);

        output_screen("Reading the chain..." << std::endl);
        double maxP = 0;
        while(!in.eof())
        {
            std::string s;
            std::getline(in, s);
            if(s == "")
                break;

            ChainElement e;
            chain.push_back(e);
            
            std::stringstream str(s);
            str >> e.p >> e.like;
            e.params.resize(n);
            for(int i = 0; i < n; ++i)
                str >> e.params[i];

            if(e.p > maxP)
                maxP = e.p;

            chain.push_back(e);
        }

        in.close();
        output_screen("OK" << std::endl);

        const double minP = maxP * 1e-6;


        output_screen("Creating the posterior distributions..." << std::endl);
        double min1 = 1e100, max1 = -1e100, min2 = 1e100, max2 = -1e100;
        for(unsigned long i = 0; i < chain.size(); ++i)
        {
            const ChainElement& e = chain[i];
            if(e.p < minP)
                continue;

            if(e.params[p1] < min1)
                min1 = e.params[p1];

            if(e.params[p1] > max1)
                max1 = e.params[p1];

            if(e.params[p2] < min2)
                min2 = e.params[p2];

            if(e.params[p2] > max2)
                max2 = e.params[p2];
        }
        const double d1 = (max1 - min1) / res;
        const double d2 = (max2 - min2) / res;
        std::vector<double> x(res), y(res);
        for(int i = 0; i < res; ++i)
        {
            x[i] = min1 + d1 * i + d1 / 2;
            y[i] = min2 + d2 * i + d2 / 2;
        }
        std::vector<std::vector<double> > z(res);
        for(int i = 0; i < res; ++i)
            z[i].resize(res, 0.0);

        double totalP = 0;
        for(unsigned long i = 0; i < chain.size(); ++i)
        {
            const ChainElement& e = chain[i];
            if(e.p < minP)
                continue;

            const double param1 = e.params[p1];
            check(param1 >= min1, "");
            int j1 = (int)std::floor((param1 - min1) / d1);
            check(j1 >= 0, "");
            if(j1 >= res)
                j1 = res - 1;

            const double param2 = e.params[p2];
            check(param2 >= min2, "");
            int j2 = (int)std::floor((param2 - min2) / d2);
            check(j2 >= 0, "");
            if(j2 >= res)
                j2 = res - 1;

            z[j1][j2] += e.p;
            totalP += e.p;
        }
        output_screen("OK" << std::endl);

        output_screen("Writing the output..." << std::endl);
        std::ofstream out(argv[4]);
        if(!out)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot write into the output file " << argv[4] << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        double prob = 0;
        unsigned long i = 0;
        while(prob < 1 - 1e-5 && i < chain.size())
        {
            if(chain[i].p > minP)
            {
                out << chain[i].params[p1 - 1] << ' ' << chain[i].params[p2 - 1] << ' ' << prob << std::endl;
            }
            prob += chain[i].p;
            ++i;
        }
        out.close();
        output_screen("OK" << std::endl);
        
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

