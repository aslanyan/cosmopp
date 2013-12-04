#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <simulate.hpp>

#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 2)
        {
            std::string exceptionStr = "Nside, noise (optional, default = 1.0), output file name(optional, white_noise_map.fits by default), and the seed (optional) must be specified.";
            exc.set(exceptionStr);
            throw exc;
        }
        
        long nSide;
        std::stringstream nSideStr;
        nSideStr << argv[1];
        nSideStr >> nSide;
        
        double noise = 1.0;
        if(argc > 2)
        {
            std::stringstream noiseStr;
            noiseStr << argv[2];
            noiseStr >> noise;
        }

        std::string outFileName("white_noise_map.fits");
        if(argc > 3)
            outFileName = argv[3];

        time_t seed = 0;
        if(argc > 4)
        {
            std::stringstream seedStr;
            seedStr << argv[4];
            seedStr >> seed;
        }
        
        Healpix_Map<double> map;
        map.SetNside(nSide, NEST);
        Simulate::simulateWhiteNoise(map, noise, seed);

        fitshandle outh;
        outh.create(outFileName);
        write_Healpix_map_to_fits(outh, map, PLANCK_FLOAT64);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
