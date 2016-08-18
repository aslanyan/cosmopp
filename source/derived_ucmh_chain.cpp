#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <progress_meter.hpp>
#include <taylor_pk.hpp>
#include <modecode.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 3)
        {
            std::string exceptionString = "Need to specify the input and output chain files.";
            exc.set(exceptionString);
            throw exc;
        }

        std::ifstream in(argv[1]);
        if(!in)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot read input file " << argv[1] << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        std::ofstream out(argv[2]);
        if(!out)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot write into output file " << argv[2] << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        const double kPiv = 0.05;
        const double kPivR = 0.002;
        const double kMin = 5e-6;
        const double kMax = 1.2;
        const int kPerDecade = 500;
        TaylorPk taylor(kPiv, kMin, kMax, kPerDecade);
        taylor.addKValue(10, 0, 1e10, 0, 1e10);
        taylor.addKValue(1e3, 0, 1e10, 0, 1e10);
        taylor.addKValue(1e6, 0, 1e10, 0, 1e10);
        taylor.addKValue(1e9, 0, 1e10, 0, 1e10);

        ModeCode::initialize(12, kPiv, 55, false, false, false, true, kMin, kMax, 500);

        std::vector<double> v(17);

        int line = 0;

        double totalWeight = 0;

        // initial scan
        while(!in.eof())
        {
            std::string s;
            std::getline(in, s);
            if(s == "")
                break;

            ++line;
            std::stringstream str(s);
            double w;
            str >> w;

            totalWeight += w;
        }

        in.close();
        in.open(argv[1]);
        check(in, "problem reading the same file for the second time");

        ProgressMeter meter(line);
        line = 0;

        const double deltaLogK = 0.05;
        const double deltaLogK1 = 0.5;

        double badWeight = 0;

        while(!in.eof())
        {
            std::string s;
            std::getline(in, s);
            if(s == "")
                break;

            ++line;
            std::stringstream str(s);
            for(int i = 0; i < 12; ++i)
                str >> v[i];

            std::vector<double> vParams(v.cbegin() + 6, v.cbegin() + 11);
            vParams[4] += std::log(vParams[0]) / std::log(10.0); // correction

            //const bool res = taylor.calculate(vParams);
            const bool res = ModeCode::calculate(vParams);
            if(!res)
            {
                /*
                std::stringstream exceptionStr;
                exceptionStr << "Problem calculating pk for line " << line << ". v parameters are:" << std::endl << vParams[0];
                for(int i = 1; i < 5; ++i)
                    exceptionStr << " " << vParams[i];
                exceptionStr << std::endl;
                exc.set(exceptionStr.str());
                throw exc;
                */
                badWeight += v[0];
            }
            else
            {
                //const Math::TableFunction<double, double>& ps = taylor.getScalarPs();
                const Math::TableFunction<double, double>& ps = ModeCode::getScalarPs();
                const Math::TableFunction<double, double>& psTensor = ModeCode::getTensorPs();

                /*
                std::ofstream out_ps("ps_class.txt");
                for(auto it = ps.cbegin(); it != ps.cend(); ++it)
                    out_ps << it->first << '\t' << it->second << std::endl;
                out_ps.close();
                */

                const double pk = std::log(ps.evaluate(kPiv));
                double pkNext = std::log(ps.evaluate(std::exp(std::log(kPiv) + deltaLogK)));
                double pkPrev = std::log(ps.evaluate(std::exp(std::log(kPiv) - deltaLogK)));
                const double ns = 1 + (pkNext - pkPrev) / (2 * deltaLogK);

                pkNext = std::log(ps.evaluate(std::exp(std::log(kPiv) + deltaLogK1)));
                pkPrev = std::log(ps.evaluate(std::exp(std::log(kPiv) - deltaLogK1)));
                const double pkNext2 = std::log(ps.evaluate(std::exp(std::log(kPiv) + 2 * deltaLogK1)));
                const double pkPrev2 = std::log(ps.evaluate(std::exp(std::log(kPiv) - 2 * deltaLogK1)));

                const double run = (pkNext - 2 * pk + pkPrev) / (deltaLogK1 * deltaLogK1);
                const double runRun = (pkNext2 - 2 * pkNext + 2 * pkPrev - pkPrev2) / (2 * deltaLogK1 * deltaLogK1 * deltaLogK1);

                const double r = psTensor.evaluate(kPivR) / ps.evaluate(kPivR);

                const double nPiv = ModeCode::getNPivot();

                v[12] = ns;
                v[13] = run;
                v[14] = runRun;
                v[15] = r;
                v[16] = nPiv;
                auto it = v.cbegin();
                out << *(it++);
                for(; it != v.cend(); ++it)
                    out << '\t' << *it;
                out << std::endl;
            }
            meter.advance();
        }
        out.close();

        output_screen("Had problem running for " << badWeight / totalWeight * 100 << "\% of the points." << std::endl);

        return 0;
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
