#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <utils.hpp>

#include "chealpix.h"
#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include "fitsio.h"

#define MY_STRINGIZE1(P) #P
#define MY_STRINGIZE(P) MY_STRINGIZE1(P)
#define HEALPIX_DATA_DIR_STR MY_STRINGIZE(HEALPIX_DATA_DIR)

void
Utils::readMask(const char* maskFileName, long& nSide, std::vector<int>& goodPixels)
{
    Healpix_Map<double> mask;
    read_Healpix_map_from_fits(std::string(maskFileName), mask);
    
    if(mask.Scheme() != NEST)
    {
        StandardException exc;
        std::string exceptionStr = "The mask must have nested ordering.";
        exc.set(exceptionStr);
        throw exc;
    }

    nSide = mask.Nside();
    
    const int nPix = (int)mask.Npix();
    
    goodPixels.clear();
    
    for(int i = 0; i < nPix; ++i)
    {
        if(mask[i] > 0.5)
        {
            goodPixels.push_back(i);
        }
    }
}

double
Utils::beamFunction(int l, double fwhm)
{
    if(fwhm == 0)
        return 1.0;
    
    check(fwhm > 0, "invalid fwhm");
    
    const double sigma = std::sqrt(8 * std::log(2.0)) / (fwhm * Math::pi / 180);
    return std::exp(-l * (l + 1) / (2 * sigma * sigma));
}

void
Utils::readPixelWindowFunction(std::vector<double>& f, long nSide, int lMax, double fwhm, bool polarization)
{
    StandardException exc;

    check(nSide > 0, "invalid nSide = " << nSide);

    int numDigits = 1;
    long nSideCopy = nSide;
    while((nSideCopy /= 10) != 0)
        ++numDigits;

    check(numDigits <= 4, "nSide = " << nSide << " is too big");
    
    f.resize(lMax + 1);
    
    std::string healpixDataDir = HEALPIX_DATA_DIR_STR;

    std::stringstream pixelWindowFileName;
    pixelWindowFileName << healpixDataDir << "/pixel_window_n";
    for(int i = 0; i < 4 - numDigits; ++i)
        pixelWindowFileName << 0;
    pixelWindowFileName << nSide << ".fits";

	fitsfile* fptr;
	char card[FLEN_CARD];
	int status = 0, nkeys, i;
	
	fits_open_file(&fptr, pixelWindowFileName.str().c_str(), READONLY, &status);
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);

	int hdunum;
	
	fits_get_num_hdus(fptr, &hdunum, &status);

    if(hdunum < 2)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Invalid format of the pixel windows function file " << pixelWindowFileName.str() << ". Number of HDUs is " << hdunum << ", needs to be at least 2.";
        exc.set(exceptionStr.str());
        throw exc;
    }

	int hdutype;
	fits_movabs_hdu(fptr, 2, &hdutype, &status);
	
	int ncols;
	long nrows;
	fits_get_num_cols(fptr, &ncols, &status);
	fits_get_num_rows(fptr, &nrows, &status);

    if(ncols < 2)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Invalid format of the pixel windows function file " << pixelWindowFileName.str() << ". Number of columns is " << ncols << ", needs to be at least 2.";
        exc.set(exceptionStr.str());
        throw exc;
    }

    if(nrows < lMax + 1)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The pixel windows function file " << pixelWindowFileName.str() << " contains values only up to l = " << nrows - 1 << ". Cannot read up to lMax = " << lMax << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    i = (polarization ? 2 : 1);

    std::stringstream s;
    s << i;
    char colname[100];
    char colnumstr[10];
    std::strcpy(colnumstr, s.str().c_str());
    int colnum, type, dispwidth, anynul;
    long repeat, width;
    fits_get_colname(fptr, CASEINSEN, colnumstr, colname, &colnum, &status);
    fits_get_coltype(fptr, colnum, &type, &repeat, &width, &status);
    fits_get_col_display_width(fptr, colnum, &dispwidth, &status);
    check(colnum == i, "");
    check(repeat == 1, "");
    
    char** column = new char*[nrows];
    for(int j = 0; j < nrows; ++j)
        column[j] = new char[dispwidth + 1];
    
    fits_read_col(fptr, TSTRING, colnum, 1, 1, (LONGLONG) nrows, 0, column, &anynul, &status);
    
    for(int l = 0; l <= lMax; ++l)
    {
        std::stringstream str;
        str << column[l];
        str >> f[l];
        f[l] *= Utils::beamFunction(l, fwhm);
    }

    for(int j = 0; j < nrows; ++j)
        delete column[j];
    delete column;

	fits_close_file(fptr, &status);
	
	if(status)
		fits_report_error(stderr, status);
}

void
Utils::readClFromFile(const char* fileName, std::vector<double>& cl, bool hasL, bool isDl)
{
    StandardException exc;
    std::ifstream in(fileName);
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot open the input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    cl.clear();
    int l = 0;
    while(!in.eof())
    {
        std::string s;
        std::getline(in, s);
        if(s == "")
            break;

        std::stringstream str(s);

        int dummyL;
        double val;
        if(hasL)
        {
            str >> dummyL;
            if(dummyL != l)
            {
                std::stringstream exceptionStr;
                exceptionStr << "Invalid format of the input file " << fileName << ". Expected to read l = " << l << " but found l = " << dummyL << ".";
                exc.set(exceptionStr.str());
                throw exc;
            }
        }

        str >> val;
        if(isDl && l)
            val *= (2 * Math::pi / (l * (l + 1)));

        cl.push_back(val);
        ++l;
    }
    in.close();
}
