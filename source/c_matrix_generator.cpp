#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <limits>
#include <iomanip>
#include <cmath>
#include <ctime>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <numerics.hpp>
#include <angular_coordinates.hpp>
#include <three_rotation.hpp>
#include <progress_meter.hpp>
#include <utils.hpp>
#include <random.hpp>
#include <legendre.hpp>
#include <c_matrix_generator.hpp>

#include "healpix_base.h"
#include "alm.h"
#include "alm_healpix_tools.h"
#include "xcomplex.h"
#include "healpix_map.h"
#include "rotmatrix.h"
#include "alm_powspec_tools.h"
#include "chealpix.h"

LegendrePolynomialContainer::LegendrePolynomialContainer(int lMax, long nSide, const std::vector<int>* goodPixels)
{
    check(lMax >= 0, "");
    data_.resize(lMax + 1);
    const int nPix = (goodPixels ? goodPixels->size() : (int)nside2npix(nSide));
    
    std::vector<Math::ThreeVectorDouble> pixels;
    
    for(int i = 0; i < nPix; ++i)
    {
        double theta, phi;
        const int index = (goodPixels ? (*goodPixels)[i] : i);
        pix2ang_nest(nSide, index, &theta, &phi);
        const Math::ThreeVectorDouble pix(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
        pixels.push_back(pix);
    }
    
    Math::Legendre legendre;
    ProgressMeter meter((lMax + 1) * nPix * (nPix + 1) / 2);
    
    for(int l = 0; l <= lMax; ++l)
    {
        data_[l].resize(nPix);
        for(int j = 0; j < nPix; ++j)
        {
            data_[l][j].resize(j + 1);
            for(int i = 0; i <= j; ++i)
            {
                double dot = pixels[i] * pixels[j];
                if(dot > 1)
                {
                    check(Math::areEqual(dot, 1.0, 0.0001), "");
                    dot = 1;
                }
                if(dot < -1)
                {
                    check(Math::areEqual(dot, -1.0, 0.0001), "");
                    dot = -1;
                }
                
                data_[l][j][i] = legendre.calculate(l, dot);
                
                meter.advance();
            }
        }
    }
}

double
LegendrePolynomialContainer::value(int l, int j, int i) const
{
    check(l < data_.size() && l >= 0, "invalid l = " << l);
    check(j < data_[l].size() && j >= 0, "invalid j");
    check(i <= j && i >= 0, "invalid i");
    check(data_[l][j].size() == j + 1, "");
    
    return data_[l][j][i];
}

LegendrePolynomialContainer::LegendrePolynomialContainer(const char* fileName)
{
    std::ifstream in;
    in.open(fileName, std::ios::in | std::ios::binary);
    
    StandardException exc;
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot open input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    int lMax;
    in.read((char*)(&lMax), sizeof(int));
    if(lMax < 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Invalid lMax = " << lMax << " read from file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    int nPix;
    in.read((char*)(&nPix), sizeof(int));
    if(nPix < 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Invalid nPix = " << nPix << " read from file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    data_.resize(lMax + 1);
    
    for(int l = 0; l <= lMax; ++l)
    {
        data_[l].resize(nPix);
        for(int j = 0; j < nPix; ++j)
        {
            data_[l][j].resize(j + 1);
            in.read((char*)(&(data_[l][j][0])), (j + 1) * sizeof(double));
        }
    }
    in.close();
}

void
LegendrePolynomialContainer::writeIntoFile(const char* fileName) const
{
    std::ofstream out(fileName, std::ios::binary | std::ios::out);
    if(!out)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    const int lMax = data_.size() - 1;
    check(lMax >= 0, "");
    out.write((char*)(&lMax), sizeof(int));
    const int nPix = data_[0].size();
    out.write((char*)(&nPix), sizeof(int));
    for(int l = 0; l <= lMax; ++l)
    {
        for(int j = 0; j < nPix; ++j)
        {
            out.write((char*)(&(data_[l][j][0])), (j + 1) * sizeof(double));
        }
    }
    out.close();
}

CMatrix*
CMatrixGenerator::clToCMatrix(const std::vector<double>& cl, long nSide, double fwhm, const std::vector<int>* goodPixels, const LegendrePolynomialContainer* lp)
{
    check(!cl.empty(), "");
    
    const int lMax = cl.size() - 1;
    const int nPix = (goodPixels ? goodPixels->size() : (int)nside2npix(nSide));
    
    std::vector<double> beam;
    Utils::readPixelWindowFunction(beam, nSide, lMax, fwhm);
    
    std::vector<Math::ThreeVectorDouble> pixels;
    
    if(!lp)
        for(int i = 0; i < nPix; ++i)
        {
            double theta, phi;
            const int index = (goodPixels ? (*goodPixels)[i] : i);
            pix2ang_nest(nSide, index, &theta, &phi);
            const Math::ThreeVectorDouble pix(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
            pixels.push_back(pix);
        }
    
    CMatrix* cMat = new CMatrix(nPix);
    std::vector<double> clCopy(cl.size());
    
    for(int l = 2; l <= lMax; ++l)
    {
        clCopy[l] = cl[l] * (2 * l + 1) / (4 * Math::pi);
    }
    
    Math::Legendre legendre;
    ProgressMeter meter(nPix * (nPix + 1) / 2);
    
    for(int j = 0; j < nPix; ++j)
    {
        for(int i = 0; i <= j; ++i)
        {
            double dot = 0;
            if(!lp)
            {
                dot = pixels[i] * pixels[j];
                if(dot > 1)
                {
                    check(Math::areEqual(dot, 1.0, 0.0001), "");
                    dot = 1;
                }
                if(dot < -1)
                {
                    check(Math::areEqual(dot, -1.0, 0.0001), "");
                    dot = -1;
                }
            }
            
            double element = 0;
            for(int l = 2; l <= lMax; ++l)
            {
                const double leg = (lp ? lp->value(l, j, i) : legendre.calculate(l, dot));
                element += clCopy[l] * leg * beam[l] * beam[l];
            }
            
            cMat->element(i, j) = element;
            cMat->element(j, i) = element; // just in case if the implementation of CMatrix changes
            
            meter.advance();
        }
    }
    return cMat;
}

CMatrix*
CMatrixGenerator::clToCMatrix(const char* clFileName, long nSide, int lMax, double fwhm, const std::vector<int>* goodPixels, const LegendrePolynomialContainer* lp)
{
    std::vector<double> cl;
    Utils::readClFromFile(clFileName, cl);
    return clToCMatrix(cl, nSide, fwhm, goodPixels, lp);
}

void
CMatrixGenerator::clToWholeMatrix(const char* clFileName, WholeMatrix& tt, WholeMatrix& te, WholeMatrix& ee)
{
    const int lMin = tt.getLMin(), lMax = tt.getLMax();
    check(lMin == te.getLMin(), "");
    check(lMax == te.getLMax(), "");
    check(lMin == ee.getLMin(), "");
    check(lMax == ee.getLMax(), "");
    
#ifdef CHECKS_ON
    for(int l1 = lMin; l1 <= lMax; ++l1)
    {
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            for(int l = lMin; l <= lMax; ++l)
            {
                for(int m = -l; m <= l; ++m)
                {
                    check(tt.element(l1, m1, l, m) == 0, "");
                    check(te.element(l1, m1, l, m) == 0, "");
                    check(ee.element(l1, m1, l, m) == 0, "");
                }
            }
        }
    }
#endif
    
    std::ifstream inCl(clFileName);
    if(!inCl)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read the file " << clFileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    for(int l = 2; l <= lMax; ++l)
    {
        int dummyL;
        double ttL, eeL, teL;
        inCl >> dummyL >> ttL >> eeL >> teL;
        if(dummyL != l)
        {
            StandardException exc;
            std::stringstream exceptionStr;
            exceptionStr << "Invalid format of file " << clFileName << ". Expected to read l = " << l << " but reading " << dummyL << " instead.";
            exc.set(exceptionStr.str());
            throw exc;
        }
        
        if(l < lMin)
            continue;
        
        for(int m = -l; m <= l; ++m)
        {
            tt.element(l, m, l, m) = ttL * 2 * Math::pi / (l * (l + 1));
            te.element(l, m, l, m) = teL * 2 * Math::pi / (l * (l + 1));
            ee.element(l, m, l, m) = eeL * 2 * Math::pi / (l * (l + 1));
        }
    }
    inCl.close();
}

CMatrix*
CMatrixGenerator::wholeMatrixToCMatrix(const WholeMatrix& wholeMatrix, long nSide, double fwhm, double phi, double theta, double psi)
{
    StandardException exc;
    
    const int nPix = (int) nside2npix(nSide);
    
    CMatrix* mat = new CMatrix(nPix);
    
    const Math::ThreeRotationMatrix rot(phi, theta, psi);
    const rotmatrix rotationMatrix(rot[0][0], rot[1][0], rot[2][0], rot[0][1], rot[1][1], rot[2][1], rot[0][2], rot[1][2], rot[2][2]);
    
    const int lMin = wholeMatrix.getLMin(), lMax = wholeMatrix.getLMax();
    
    check(lMin < lMax, "");
    
    //first index is l', second index is l' + m'
    std::vector<std::vector<Alm<xcomplex<double> > > > re(lMax + 1), im(lMax + 1);
    std::vector<std::vector<Healpix_Map<double> > > reMap(lMax + 1), imMap(lMax + 1);
    
    for(int l1 = 0; l1 <= lMax; ++l1)
    {
        re[l1].resize(2 * l1 + 1);
        im[l1].resize(2 * l1 + 1);
        reMap[l1].resize(2 * l1 + 1);
        imMap[l1].resize(2 * l1 + 1);
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            re[l1][l1 + m1].Set(lMax, lMax);
            im[l1][l1 + m1].Set(lMax, lMax);
            reMap[l1][l1 + m1].SetNside(nSide, RING);
            imMap[l1][l1 + m1].SetNside(nSide, RING);
        }
    }
    
    std::vector<double> beam;
    Utils::readPixelWindowFunction(beam, nSide, lMax, fwhm);
    
    //output_screen("Calculating real and imaginary parts of first conversion..." << std::endl);
    for(int l1 = 0; l1 <= lMax; ++l1)
    {
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            for(int l = 0; l <= lMax; ++l)
            {
                for(int m = 0; m <= l; ++m)
                {
                    if(l1 < lMin || l < lMin)
                    {
                        re[l1][l1 + m1](l, m) = xcomplex<double>(0, 0);
                        im[l1][l1 + m1](l, m) = xcomplex<double>(0, 0);
                        continue;
                    }
                    const int minus1m = (m % 2 ? -1 : 1);
                    re[l1][l1 + m1](l, m) = xcomplex<double>((wholeMatrix.element(l1, m1, l, m) + minus1m * wholeMatrix.element(l1, m1, l, -m)) / 2, 0);
                    im[l1][l1 + m1](l, m) = xcomplex<double>(0, -(wholeMatrix.element(l1, m1, l, m) - minus1m * wholeMatrix.element(l1, m1, l, -m)) / 2);
                    
                    // multiply by the beam
                    re[l1][l1 + m1](l, m) *= (beam[l] * beam[l1]);
                    im[l1][l1 + m1](l, m) *= (beam[l] * beam[l1]);
                }
            }
        }
    }
    //output_screen("OK" << std::endl);
    
    //output_screen("Converting to pixel space for l and m..." << std::endl);
    //ProgressMeter meter1((lMax + 1) * (lMax + 1));
    for(int l1 = 0; l1 <= lMax; ++l1)
    {
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            rotate_alm(re[l1][l1 + m1], rotationMatrix);
            rotate_alm(im[l1][l1 + m1], rotationMatrix);
            
            alm2map(re[l1][l1 + m1], reMap[l1][l1 + m1]);
            alm2map(im[l1][l1 + m1], imMap[l1][l1 + m1]);
            reMap[l1][l1 + m1].swap_scheme();
            imMap[l1][l1 + m1].swap_scheme();
            
            //meter1.advance();
        }
    }
    //output_screen("OK" << std::endl);
    
    //output_screen("Calculating real and imaginary parts of second conversion" << std::endl);
    std::vector<Alm<xcomplex<double> > > rePix(nPix);
#ifdef CHECKS_ON
    std::vector<Alm<xcomplex<double> > > imPix(nPix);
#endif
    for(int i = 0; i < nPix; ++i)
    {
        rePix[i].Set(lMax, lMax);
#ifdef CHECKS_ON
        imPix[i].Set(lMax, lMax);
#endif
        for(int l1 = 0; l1 <= lMax; ++l1)
        {
            for(int m1 = 0; m1 <= l1; ++m1)
            {
                const xcomplex<double> al1m1 = xcomplex<double>(reMap[l1][l1 + m1][i], -imMap[l1][l1 + m1][i]);
                const xcomplex<double> al1minusm1 = xcomplex<double>(reMap[l1][l1 - m1][i], -imMap[l1][l1 - m1][i]);
                
                const double minus1m1 = (m1 % 2 ? -1.0 : 1.0);

                rePix[i](l1, m1) = 0.5 * (al1m1 + minus1m1 * al1minusm1.conj());
#ifdef CHECKS_ON
                imPix[i](l1, m1) = xcomplex<double>(0, -0.5) * (al1m1 - minus1m1 * al1minusm1.conj());
#endif
            }
        }
    }
    //output_screen("OK" << std::endl);
    
    //output_screen("Converting to pixel space for l' and m'..." << std::endl);
    std::vector<Healpix_Map<double> > rePixMap(nPix);
#ifdef CHECKS_ON
    std::vector<Healpix_Map<double> > imPixMap(nPix);
#endif
    //ProgressMeter meter2(nPix);
    for(int i = 0; i < nPix; ++i)
    {
        rePixMap[i].SetNside(nSide, RING);
        rotate_alm(rePix[i], rotationMatrix);
        alm2map(rePix[i], rePixMap[i]);
        rePixMap[i].swap_scheme();
#ifdef CHECKS_ON
        imPixMap[i].SetNside(nSide, RING);
        rotate_alm(imPix[i], rotationMatrix);
        alm2map(imPix[i], imPixMap[i]);
        imPixMap[i].swap_scheme();
#endif
        //meter2.advance();
    }
    //output_screen("OK" << std::endl);
    
    for(int j = 0; j < nPix; ++j)
    {
        for(int i = 0; i <= j; ++i)
        {
            mat->element(i, j) = rePixMap[i][j];
            mat->element(j, i) = rePixMap[i][j];
            check(Math::areEqual(imPixMap[i][j], 0.0, 1e-15), "");
        }
    }
    
    return mat;
}

CMatrix*
CMatrixGenerator::polarizationEEWholeMatrixToCMatrix(const WholeMatrix& ee, long nSide, double fwhm, double phi, double theta, double psi)
{
    StandardException exc;
    
    const int nPix = (int) nside2npix(nSide);
    
    CMatrix* mat = new CMatrix(2 * nPix);
    mat->comment() = "polarization";
    
    const Math::ThreeRotationMatrix rot(phi, theta, psi);
    const rotmatrix rotationMatrix(rot[0][0], rot[1][0], rot[2][0], rot[0][1], rot[1][1], rot[2][1], rot[0][2], rot[1][2], rot[2][2]);
    const int lMin = ee.getLMin(), lMax = ee.getLMax();
    
    check(lMin < lMax, "");
    
    //first index is l', second index is l' + m'
    std::vector<std::vector<Alm<xcomplex<double> > > > t(lMax + 1), reE(lMax + 1), imE(lMax + 1), b(lMax + 1);
    std::vector<std::vector<Healpix_Map<double> > > tMap(lMax + 1), reQMap(lMax + 1), imQMap(lMax + 1), reUMap(lMax + 1), imUMap(lMax + 1);
    
    for(int l1 = 0; l1 <= lMax; ++l1)
    {
        t[l1].resize(2 * l1 + 1);
        reE[l1].resize(2 * l1 + 1);
        imE[l1].resize(2 * l1 + 1);
        b[l1].resize(2 * l1 + 1);
        
        tMap[l1].resize(2 * l1 + 1);
        reQMap[l1].resize(2 * l1 + 1);
        imQMap[l1].resize(2 * l1 + 1);
        reUMap[l1].resize(2 * l1 + 1);
        imUMap[l1].resize(2 * l1 + 1);
        
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            t[l1][l1 + m1].Set(lMax, lMax);
            reE[l1][l1 + m1].Set(lMax, lMax);
            imE[l1][l1 + m1].Set(lMax, lMax);
            b[l1][l1 + m1].Set(lMax, lMax);
            
            tMap[l1][l1 + m1].SetNside(nSide, RING);
            reQMap[l1][l1 + m1].SetNside(nSide, RING);
            imQMap[l1][l1 + m1].SetNside(nSide, RING);
            reUMap[l1][l1 + m1].SetNside(nSide, RING);
            imUMap[l1][l1 + m1].SetNside(nSide, RING);
        }
    }
    
    std::vector<double> beam;
    Utils::readPixelWindowFunction(beam, nSide, lMax, fwhm);
    
    //output_screen("Calculating real and imaginary parts of first conversion..." << std::endl);
    for(int l1 = 0; l1 <= lMax; ++l1)
    {
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            for(int l = 0; l <= lMax; ++l)
            {
                for(int m = 0; m <= l; ++m)
                {
                    if(l1 < lMin || l < lMin)
                    {
                        t[l1][l1 + m1](l, m) = xcomplex<double>(0, 0);
                        reE[l1][l1 + m1](l, m) = xcomplex<double>(0, 0);
                        imE[l1][l1 + m1](l, m) = xcomplex<double>(0, 0);
                        b[l1][l1 + m1](l, m) = xcomplex<double>(0, 0);
                        continue;
                    }
                    const int minus1m = (m % 2 ? -1 : 1);
                    t[l1][l1 + m1](l, m) = xcomplex<double>(0, 0);
                    b[l1][l1 + m1](l, m) = xcomplex<double>(0, 0);
                    reE[l1][l1 + m1](l, m) = xcomplex<double>((ee.element(l1, m1, l, m) + minus1m * ee.element(l1, m1, l, -m)) / 2, 0);
                    imE[l1][l1 + m1](l, m) = xcomplex<double>(0, -(ee.element(l1, m1, l, m) - minus1m * ee.element(l1, m1, l, -m)) / 2);
                    
                    // multiply by the beam
                    reE[l1][l1 + m1](l, m) *= (beam[l] * beam[l1]);
                    imE[l1][l1 + m1](l, m) *= (beam[l] * beam[l1]);
                }
            }
        }
    }
    //output_screen("OK" << std::endl);
    
    //output_screen("Converting to pixel space for l and m..." << std::endl);
    //ProgressMeter meter1((lMax + 1) * (lMax + 1));
    for(int l1 = 0; l1 <= lMax; ++l1)
    {
        for(int m1 = -l1; m1 <= l1; ++m1)
        {
            rotate_alm(t[l1][l1 + m1], reE[l1][l1 + m1], b[l1][l1 + m1], rotationMatrix);
            rotate_alm(t[l1][l1 + m1], imE[l1][l1 + m1], b[l1][l1 + m1], rotationMatrix);
            
            alm2map_pol(t[l1][l1 + m1], reE[l1][l1 + m1], b[l1][l1 + m1], tMap[l1][l1 + m1], reQMap[l1][l1 + m1], reUMap[l1][l1 + m1]);
            alm2map_pol(t[l1][l1 + m1], imE[l1][l1 + m1], b[l1][l1 + m1], tMap[l1][l1 + m1], imQMap[l1][l1 + m1], imUMap[l1][l1 + m1]);
            reQMap[l1][l1 + m1].swap_scheme();
            imQMap[l1][l1 + m1].swap_scheme();
            reUMap[l1][l1 + m1].swap_scheme();
            imUMap[l1][l1 + m1].swap_scheme();
            
            //meter1.advance();
        }
    }
    //output_screen("OK" << std::endl);
    
    //output_screen("Calculating real and imaginary parts of second conversion" << std::endl);
    std::vector<Alm<xcomplex<double> > > tPix(nPix), reQPix(nPix), reUPix(nPix), bPix(nPix);
#ifdef CHECKS_ON
    std::vector<Alm<xcomplex<double> > > imQPix(nPix), imUPix(nPix);
#endif
    for(int i = 0; i < nPix; ++i)
    {
        tPix[i].Set(lMax, lMax);
        reQPix[i].Set(lMax, lMax);
        reUPix[i].Set(lMax, lMax);
        bPix[i].Set(lMax, lMax);
#ifdef CHECKS_ON
        imQPix[i].Set(lMax, lMax);
        imUPix[i].Set(lMax, lMax);
#endif
        for(int l1 = 0; l1 <= lMax; ++l1)
        {
            for(int m1 = 0; m1 <= l1; ++m1)
            {
                const xcomplex<double> al1m1Q = xcomplex<double>(reQMap[l1][l1 + m1][i], -imQMap[l1][l1 + m1][i]);
                const xcomplex<double> al1minusm1Q = xcomplex<double>(reQMap[l1][l1 - m1][i], -imQMap[l1][l1 - m1][i]);
                const xcomplex<double> al1m1U = xcomplex<double>(reUMap[l1][l1 + m1][i], -imUMap[l1][l1 + m1][i]);
                const xcomplex<double> al1minusm1U = xcomplex<double>(reUMap[l1][l1 - m1][i], -imUMap[l1][l1 - m1][i]);
                
                const double minus1m1 = (m1 % 2 ? -1.0 : 1.0);
                
                tPix[i](l1, m1) = xcomplex<double>(0, 0);
                reQPix[i](l1, m1) = 0.5 * (al1m1Q + minus1m1 * al1minusm1Q.conj());
                reUPix[i](l1, m1) = 0.5 * (al1m1U + minus1m1 * al1minusm1U.conj());
                bPix[i](l1, m1) = xcomplex<double>(0, 0);
#ifdef CHECKS_ON
                imQPix[i](l1, m1) = xcomplex<double>(0, -0.5) * (al1m1Q - minus1m1 * al1minusm1Q.conj());
                imUPix[i](l1, m1) = xcomplex<double>(0, -0.5) * (al1m1U - minus1m1 * al1minusm1U.conj());
#endif
            }
        }
    }
    //output_screen("OK" << std::endl);
    
    //output_screen("Converting to pixel space for l' and m'..." << std::endl);
    std::vector<Healpix_Map<double> > tPixMap(nPix), reQQPixMap(nPix), reQUPixMap(nPix), reUQPixMap(nPix), reUUPixMap(nPix);
#ifdef CHECKS_ON
    std::vector<Healpix_Map<double> > imQQPixMap(nPix), imQUPixMap(nPix), imUQPixMap(nPix), imUUPixMap(nPix);
#endif
    //ProgressMeter meter2(nPix);
    for(int i = 0; i < nPix; ++i)
    {
        tPixMap[i].SetNside(nSide, RING);
        reQQPixMap[i].SetNside(nSide, RING);
        reQUPixMap[i].SetNside(nSide, RING);
        reUQPixMap[i].SetNside(nSide, RING);
        reUUPixMap[i].SetNside(nSide, RING);
        
        rotate_alm(tPix[i], reQPix[i], bPix[i], rotationMatrix);
        rotate_alm(tPix[i], reUPix[i], bPix[i], rotationMatrix);
        
        alm2map_pol(tPix[i], reQPix[i], bPix[i], tPixMap[i], reQQPixMap[i], reQUPixMap[i]);
        alm2map_pol(tPix[i], reUPix[i], bPix[i], tPixMap[i], reUQPixMap[i], reUUPixMap[i]);
        
        reQQPixMap[i].swap_scheme();
        reQUPixMap[i].swap_scheme();
        reUQPixMap[i].swap_scheme();
        reUUPixMap[i].swap_scheme();
        
#ifdef CHECKS_ON
        imQQPixMap[i].SetNside(nSide, RING);
        imQUPixMap[i].SetNside(nSide, RING);
        imUQPixMap[i].SetNside(nSide, RING);
        imUUPixMap[i].SetNside(nSide, RING);
        
        rotate_alm(tPix[i], imQPix[i], bPix[i], rotationMatrix);
        rotate_alm(tPix[i], imUPix[i], bPix[i], rotationMatrix);
        
        alm2map_pol(tPix[i], imQPix[i], bPix[i], tPixMap[i], imQQPixMap[i], imQUPixMap[i]);
        alm2map_pol(tPix[i], imUPix[i], bPix[i], tPixMap[i], imUQPixMap[i], imUUPixMap[i]);
        
        imQQPixMap[i].swap_scheme();
        imQUPixMap[i].swap_scheme();
        imUQPixMap[i].swap_scheme();
        imUUPixMap[i].swap_scheme();
#endif
        //meter2.advance();
    }
    //output_screen("OK" << std::endl);
    
    for(int j = 0; j < nPix; ++j)
    {
        for(int i = 0; i < nPix; ++i)
        {
            mat->element(i, j) = reQQPixMap[i][j];
            mat->element(i, nPix + j) = reQUPixMap[i][j];
            mat->element(nPix + i, j) = reUQPixMap[i][j];
            mat->element(nPix + i, nPix + j) = reUUPixMap[i][j];
            
            check(Math::areEqual(reQQPixMap[i][j], reQQPixMap[j][i], 1e-7), reQQPixMap[i][j] << ' ' << reQQPixMap[j][i]);
            check(Math::areEqual(reQUPixMap[i][j], reUQPixMap[j][i], 1e-7), "");
            check(Math::areEqual(reUUPixMap[i][j], reUUPixMap[j][i], 1e-7), "");
            
            check(Math::areEqual(imQQPixMap[i][j], 0.0, 1e-15), "");
            check(Math::areEqual(imQUPixMap[i][j], 0.0, 1e-15), "");
            check(Math::areEqual(imUQPixMap[i][j], 0.0, 1e-15), "");
            check(Math::areEqual(imUUPixMap[i][j], 0.0, 1e-15), "");
        }
    }
    
    return mat;
}

CMatrix*
CMatrixGenerator::getFiducialMatrix(const char* clFileName, long nSide, int lMax, double fwhm, const std::vector<int>* goodPixels, const LegendrePolynomialContainer* lp)
{
    std::vector<double> cl;
    Utils::readClFromFile(clFileName, cl);
    return getFiducialMatrix(cl, nSide, lMax, fwhm, goodPixels, lp);
}

CMatrix*
CMatrixGenerator::getFiducialMatrix(const std::vector<double>& cl, long nSide, int lMax, double fwhm, const std::vector<int>* goodPixels, const LegendrePolynomialContainer* lp)
{
    const int lMaxMax = 4 * nSide;
    check(cl.size() >= lMaxMax + 1, "");
    
    const int nPix = (goodPixels ? goodPixels->size() : (int)nside2npix(nSide));
    
    std::vector<double> beam;
    Utils::readPixelWindowFunction(beam, nSide, lMaxMax, fwhm);
    
    std::vector<Math::ThreeVectorDouble> pixels;
    
    for(int i = 0; i < nPix; ++i)
    {        
        double theta, phi;
        const int index = (goodPixels ? (*goodPixels)[i] : i);
        pix2ang_nest(nSide, index, &theta, &phi);
        const Math::ThreeVectorDouble pix(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
        pixels.push_back(pix);
    }
    
    CMatrix* fiducialMat = new CMatrix(nPix);
    fiducialMat->comment() = "fiducial matrix";
    
    Math::Legendre legendre;
    ProgressMeter meter(nPix * (nPix + 1) / 2);
    
    for(int j = 0; j < nPix; ++j)
    {
        for(int i = 0; i <= j; ++i)
        {
            double dot;
            if(lp)
                dot = lp->value(1, j, i);
            else
            {
                dot = pixels[i] * pixels[j];
                if(dot > 1)
                {
                    check(Math::areEqual(dot, 1.0, 0.0001), "");
                    dot = 1;
                }
                if(dot < -1)
                {
                    check(Math::areEqual(dot, -1.0, 0.0001), "");
                    dot = -1;
                }
            }
            
            double element = 0;
            for(int l = lMax + 1; l <= lMaxMax; ++l)
            {
                const double leg = (lp ? lp->value(l, j, i) : legendre.calculate(l, dot));
                element += cl[l] * ((2 * l + 1) / (4 * Math::pi)) * leg * beam[l] * beam[l];
            }
            
            //marginalize monopole and dipole
            element += 100 * cl[2] * (1 + dot) * beam[2] * beam[2];
            
            fiducialMat->element(i, j) = element;
            fiducialMat->element(j, i) = element; // just in case if the implementation of CMatrix changes
            
            meter.advance();
        }
    }
    return fiducialMat;
}

CMatrix*
CMatrixGenerator::generateNoiseMatrix(long nSide, double noise)
{
    const int nPix = (int)nside2npix(nSide);
    
    CMatrix* mat = new CMatrix(nPix);
    mat->comment() = "noise matrix";
    
    for(int i = 0; i < nPix; ++i)
    {
        mat->element(i, i) = noise * noise;
    }
    return mat;
}

CMatrix*
CMatrixGenerator::calculateNoiseMatrix(const char* maskFileName, const char* noiseDataFileName, double sigma0, double fwhm, long nSideOriginal, double fwhmOriginal)
{
    StandardException exc;
    
    long nSide;
    std::vector<int> goodPixels;
    Utils::readMask(maskFileName, nSide, goodPixels);
    
    const int nPix = (int)nside2npix(nSide);
    
    const int goodPixelsSize = goodPixels.size();
    
    const long nPixOriginal = nside2npix(nSideOriginal);
    
    Healpix_Map<double> noiseMap;
    noiseMap.SetNside(nSideOriginal, NEST);
    
    CMatrix* noiseMatrix = new CMatrix(goodPixelsSize);
    noiseMatrix->comment() = "noise matrix";
    
    std::ifstream inNoise(noiseDataFileName);
    if(!inNoise)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read the noise input file " << noiseDataFileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    for(int i = 0; i < nPixOriginal; ++i)
    {
        long index;
        double dummy, nObs;
        inNoise >> index >> dummy >> nObs;
        if(index != i)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid format of input file. Expected to read pixel " << i <<  " but found " << index << " instead.";
            exc.set(exceptionStr.str());
            throw exc;
        }
        
        noiseMap[i] = sigma0 / std::sqrt(nObs);
    }
    
    const int N = 1000;
    
    std::vector<Healpix_Map<double> > noiseSimLowRes(N);
    
    arr<double> weight(2 * nSideOriginal);
    
    for(int i = 0; i < 2 * nSideOriginal; ++i)
        weight[i] = 1;
    
    Math::GaussianGenerator generator(std::time(0), 0, 1);
    
    const int lMax = 2 * nSide;
    output_screen("Reading the pixel window functions..." << std::endl);
    std::vector<double> beam, originalBeam;
    Utils::readPixelWindowFunction(beam, nSide, lMax, fwhm);
    Utils::readPixelWindowFunction(originalBeam, nSideOriginal, lMax, fwhmOriginal);
    
    output_screen("OK" << std::endl);
    
    output_screen("Simulating " << N << " instances of noise..." << std::endl);
    ProgressMeter meter(N);
    for(int i = 0; i < N; ++i)
    {
        Healpix_Map<double> noiseSim;
        noiseSim.SetNside(nSideOriginal, NEST);
        for(int j = 0; j < nPixOriginal; ++j)
        {
            noiseSim[j] = generator.generate() * noiseMap[j];
        }
        noiseSim.swap_scheme();
        Alm<xcomplex<double> > alm(lMax, lMax);
        map2alm_iter(noiseSim, alm, 10, weight);
        
        for(int l = 0; l <= lMax; ++l)
        {
            for(int m = 0; m <= l; ++m)
            {
                alm(l, m) *= (beam[l] / originalBeam[l]);
            }
        }
        
        noiseSimLowRes[i].SetNside(nSide, RING);
        alm2map(alm, noiseSimLowRes[i]);
        noiseSimLowRes[i].swap_scheme();
        meter.advance();
    }
    output_screen("OK" << std::endl);
    
    double averageNoise = 0;
    double maxAvgPerPixel = 0;
    output_screen("Averaging out the simulations..." << std::endl);
    
    for(int i = 0; i < goodPixelsSize; ++i)
    {
        for(int j = i; j < goodPixelsSize; ++j)
        {
            double element = 0;
            for(int k = 0; k < N; ++k)
                element += noiseSimLowRes[k][goodPixels[i]] * noiseSimLowRes[k][goodPixels[j]];
            
            element /= N;
            
            noiseMatrix->element(i, j) = element;
            noiseMatrix->element(j, i) = element;
        }
    }
    for(int j = 0; j < goodPixelsSize; ++j)
    {
        double avg = 0;
        for(int i = 0; i < N; ++i)
            avg += noiseSimLowRes[i][goodPixels[j]];
        
        avg /= N;
        if(std::abs(avg) > maxAvgPerPixel)
            maxAvgPerPixel = std::abs(avg);
        
        /*double sigma = 0;
        for(int i = 0; i < N; ++i)
        {
            const double diff = noiseSimLowRes[i][goodPixels[j]] - avg;
            sigma += diff * diff;
        }
        
        sigma /= (N - 1);
        
        averageNoise += std::sqrt(sigma);
        noiseMatrix->element(j, j) = sigma;*/
    }
    averageNoise /= goodPixelsSize;
    output_screen("Average noise per pixel is " << averageNoise << std::endl);
    
    output_screen("Maximum average of simulations is " << maxAvgPerPixel << std::endl);
    
    return noiseMatrix;
}
