#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <cubic_spline.hpp>
#include <table_function.hpp>

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
};

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 3)
        {
            std::string exceptionStr = "Input chain file, the number of parameters and the posterior resolution (i.e. number of points, optional, default = 100) must be specified.";
            exc.set(exceptionStr);
            throw exc;
        }

        std::stringstream countStr;
        countStr << argv[2];
        int n;
        countStr >> n;

        int res = 100;
        if(argc > 3)
        {
            std::stringstream resStr;
            resStr << argv[3];
            resStr >> res;

            if(res <= 0)
            {
                std::stringstream exceptionStr;
                exceptionStr << "Invalid resolution " << argv[3] << ". Needs to be a positive integer.";
                exc.set(exceptionStr.str());
                throw exc;
            }
        }

        if(n <= 0)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid number of parameters " << argv[2] << ". Needs to be a positive integer.";
            exc.set(exceptionStr.str());
            throw exc;
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
        for(int i = 0; i < n; ++i)
        {
            double min = 1e100, max = -1e100;
            for(unsigned long j = 0; j < chain.size(); ++j)
            {
                const ChainElement& e = chain[j];
                if(e.p < minP)
                    continue;

                if(e.params[i] < min)
                    min = e.params[i];
                if(e.params[i] > max)
                    max = e.params[i];
            }

            const double d = (max - min) / res;
            std::vector<double> x(res), y(res, 0);
            for(int j = 0; j < res; ++j)
                x[j] = min + d * j + d / 2;

            for(unsigned long j = 0; j < chain.size(); ++j)
            {
                const ChainElement& e = chain[j];
                if(e.p < minP)
                    continue;

                const double p = e.params[i];
                check(p >= min, "");
                int k = (int)std::floor((p - min) / d);
                check(k >= 0, "");
                if(k >= res)
                    k = res - 1;

                y[k] += e.p;
            }

            Math::CubicSpline cs(x, y);
            std::stringstream fileName;
            fileName << "posterior_param_" << i + 1 << ".txt";
            std::ofstream out(fileName.str().c_str());
            if(!out)
            {
                std::stringstream exceptionStr;
                exceptionStr << "Cannot write into output file " << fileName.str() << ".";
                exc.set(exceptionStr.str());
                throw exc;
            }


            const int N = 10000;
            const double delta = (x[x.size() - 1] - x[0]) / N;
            double maxVal = 0, maxX;
            Math::TableFunction<double, double> cumul, cumulInv;
            double sum = 0;
            for(int j = 0; j <= N; ++j)
            {
                double v = x[0] + j * delta;
                
                if(j == N)
                    v = x[x.size() - 1];

                check(v <= x[x.size() - 1], "");

                double y = cs.evaluate(v);
                if(y < 0)
                    y = 0;

                if(y > maxVal)
                {
                    maxVal = y;
                    maxX = v;
                }

                sum += y * delta;
                cumul[sum] = v;
                cumulInv[v] = sum;
            }

            for(int j = 0; j <= N; ++j)
            {
                double v = x[0] + j * delta;
                
                if(j == N)
                    v = x[x.size() - 1];

                check(v <= x[x.size() - 1], "");
                double y = cs.evaluate(v);
                if(y < 0)
                    y = 0;
                out << v << '\t' << y / sum << std::endl;
            }

            out.close();

            const double cumulBest = cumulInv.evaluate(maxX);
            const double cumul1Min = cumulBest - 0.5 * 0.683 * sum, cumul1Max = cumulBest + 0.5 * 0.683 * sum;
            const double cumul2Min = cumulBest - 0.5 * 0.955 * sum, cumul2Max = cumulBest + 0.5 * 0.955 * sum;
            const double u1 = 0.683 * sum, u2 = 0.955 * sum;
            const double l1 = sum - u1, l2 = sum - u2;

            double min1, max1, min2, max2, upper1, upper2, lower1, lower2;
            Math::TableFunction<double, double>::const_iterator it = cumul.lower_bound(cumul1Min);
            check(it != cumul.end(), "");
            min1 = (*it).second;

            it = cumul.lower_bound(cumul1Max);
            if(it == cumul.end())
                --it;
            check(it != cumul.end(), "");
            max1 = (*it).second;

            it = cumul.lower_bound(cumul2Min);
            check(it != cumul.end(), "");
            min2 = (*it).second;

            it = cumul.lower_bound(cumul2Max);
            if(it == cumul.end())
                --it;
            check(it != cumul.end(), "");
            max2 = (*it).second;

            it = cumul.lower_bound(u1);
            check(it != cumul.end(), "");
            upper1 = (*it).second;

            it = cumul.lower_bound(u2);
            check(it != cumul.end(), "");
            upper2 = (*it).second;

            it = cumul.lower_bound(l1);
            check(it != cumul.end(), "");
            lower1 = (*it).second;

            it = cumul.lower_bound(l2);
            check(it != cumul.end(), "");
            lower2 = (*it).second;

            output_screen("Parameter " << i + 1 << ":\t" << maxX << "\t(" << min1 << ", " << max1 << ")\t(" << min2 << ", " << max2 << ")\n   Upper bounds:\t" << upper1 << "\t" << upper2 << "\n   Lower bounds:\t" << lower1 << "\t" << lower2 << std::endl);
        }
        output_screen("OK" << std::endl);
        
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}

