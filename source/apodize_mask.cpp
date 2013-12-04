#include <string>
#include <sstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <mask_apodizer.hpp>

#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 5)
        {
            std::string exceptionStr = "Input mask file, apodization type (cosine or gaussian), apodization angle, and output file must be specified.";
            exc.set(exceptionStr);
            throw exc;
        }
        
        fitshandle outh;

        MaskApodizer::ApodizationType type = MaskApodizer::APODIZATION_TYPE_MAX;
        
        if(std::string(argv[2]) == std::string("cosine"))
            type = MaskApodizer::COSINE_APODIZATION;

        if(std::string(argv[2]) == std::string("gaussian"))
            type = MaskApodizer::GAUSSIAN_APODIZATION;

        if(type == MaskApodizer::APODIZATION_TYPE_MAX)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid apodization type " << type << ", must be cosine or gaussian.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        std::stringstream inStr;
        inStr << argv[3];
        double angle;
        inStr >> angle;

        if(angle <= 0)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid apodization angle " << angle << ", it must be positive.";
            exc.set(exceptionStr.str());
            throw exc;
        }
        
        Healpix_Map<double> mask;
        output_screen("Reading the input mask..." << std::endl);
        read_Healpix_map_from_fits(std::string(argv[1]), mask);
        output_screen("OK" << std::endl);

        MaskApodizer ap(mask);
        Healpix_Map<double> result;
        ap.apodize(type, angle, result);

        output_screen("Writing the output mask..." << std::endl);
        outh.create(std::string(argv[4]));
        write_Healpix_map_to_fits(outh, result, PLANCK_FLOAT64);
        output_screen("OK" << std::endl);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
