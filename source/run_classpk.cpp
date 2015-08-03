#include <vector>
#include <fstream>

#include <macros.hpp>
#include <taylor_pk.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 6)
        {
            std::string exceptionString = "Need to specify the 5 potential params.";
            exc.set(exceptionString);
            throw exc;
        }

        std::stringstream args;
        args << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5];
        std::vector<double> v(5);
        args >> v[0] >> v[1] >> v[2] >> v[3] >> v[4];

        const double kMin = 5e-6, kMax = 1.2;
        const double kPerDecade = 10;
        const double kPivot = 0.002;

        TaylorPk tpk(kPivot, kMin, kMax, kPerDecade);

        tpk.addKValue(10, 0, 1e10, 0, 1e10);
        tpk.addKValue(1e2, 0, 1e10, 0, 1e10);
        tpk.addKValue(1e3, 0, 1e10, 0, 1e10);
        tpk.addKValue(1e4, 0, 1e10, 0, 1e10);
        tpk.addKValue(1e5, 0, 1e10, 0, 1e10);
        tpk.addKValue(1e6, 0, 1e10, 0, 1e10);
        tpk.addKValue(1e7, 0, 1e10, 0, 1e10);
        tpk.addKValue(1e8, 0, 1e10, 0, 1e10);
        tpk.addKValue(1e9, 0, 1e10, 0, 1e10);

        if(tpk.calculate(v))
        {
            output_screen("SUCCESS" << std::endl);
            std::ofstream out("ps_class.txt");
            const Math::TableFunction<double, double>& scalar = tpk.getScalarPs();
            const Math::TableFunction<double, double>& tensor = tpk.getTensorPs();

            for(Math::TableFunction<double, double>::const_iterator it = scalar.begin(); it != scalar.end(); ++it)
            {
                const double k = it->first;
                const double s = it->second;
                const double t = tensor.evaluate(k);
                out << k << ' ' << s << ' ' << t << std::endl;
            }
            out.close();
        }
        else {
            output_screen("FAIL" << std::endl);
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
