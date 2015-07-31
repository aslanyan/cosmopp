#include <fstream>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <modecode.hpp>

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
        const double NPivot = 70;
        ModeCode::initialize(12, 0.002, NPivot, false, false, false, true, 5e-6, 1.2, 500);

        /*
        ModeCode::addKValue(10, 0, 1e10, 0, 1e10);
        ModeCode::addKValue(1e2, 0, 1e10, 0, 1e10);
        ModeCode::addKValue(1e3, 0, 1e10, 0, 1e10);
        ModeCode::addKValue(1e4, 0, 1e10, 0, 1e10);
        ModeCode::addKValue(1e5, 0, 1e10, 0, 1e10);
        ModeCode::addKValue(1e6, 0, 1e10, 0, 1e10);
        ModeCode::addKValue(1e7, 0, 1e10, 0, 1e10);
        ModeCode::addKValue(1e8, 0, 1e10, 0, 1e10);
        ModeCode::addKValue(1e9, 0, 1e10, 0, 1e10);
        */

        //std::vector<double>  v1{0.063324347453767582117, 0.087259287071610569519, -0.07146682502014208016, -0.012870094586320449115,    -8.6053603888697871582};

        //output_screen("Now: " << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << v[3] << ' ' << v[4] << std::endl);
        //ModeCode::calculate(v1);

        const bool res = ModeCode::calculate(v);


        if(!res)
        {
            output_screen("ModeCode failed!" << std::endl);
            return 1;
        }
        
        const Math::TableFunction<double, double>& scalarPs = ModeCode::getScalarPs();
        const Math::TableFunction<double, double>& tensorPs = ModeCode::getTensorPs();

        std::ofstream out("ps.txt");
        for(Math::TableFunction<double, double>::const_iterator it = scalarPs.begin(); it != scalarPs.end(); ++it)
        {
            const double k = it->first;
            const double s = it->second;
            check(tensorPs.find(k) != tensorPs.end(), "");
            const double t = tensorPs.evaluate(k);
            out << k << ' ' << s << ' ' << t << std::endl;
        }

        out.close();
        return 0;

    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
