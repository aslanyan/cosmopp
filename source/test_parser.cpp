#include <macros.hpp>
#include <exception_handler.hpp>
#include <parser.hpp>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 2)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Need to pass at least one parameters file as an argument!";
            exc.set(exceptionStr.str());
            throw(exc);
        }
        
        Parser p(argv[1]);
        for(int i = 2; i < argc; ++i)
            p.readFile(argv[i]);

        p.dump();
        double d = p.getDouble("d");
        output_screen("d = " << d << std::endl);

        int x= p.getInt("x", 15);
        output_screen("x = " << x << std::endl);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
