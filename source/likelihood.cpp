#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <numerics.hpp>
#include <progress_meter.hpp>
#include <three_rotation.hpp>
#include <likelihood.hpp>
#include <utils.hpp>

#include <gmd.h>
#include <lavd.h>
#include <laslv.h>
#include <lavli.h>
#include <blas2pp.h>
#include <blas3pp.h>

/*
#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
*/

#include "healpix_base.h"
#include "alm.h"
#include "alm_healpix_tools.h"
#include "xcomplex.h"
#include "healpix_map.h"
#include "rotmatrix.h"
#include "alm_powspec_tools.h"
#include "chealpix.h"
#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include <omp.h>

Likelihood::Likelihood(const CMatrix& cMatrix, const CMatrix& fiducialMatrix, const CMatrix& noiseMatrix, const char* maskFileName, const char* foregroundFileName)
{
    construct(cMatrix, fiducialMatrix, noiseMatrix, maskFileName, foregroundFileName);
}

Likelihood::Likelihood(const CMatrix& cMatrix, const CMatrix& fiducialMatrix, const CMatrix& noiseMatrix, const std::vector<int>& goodPixels, const LaVectorDouble& foreground)
{
    construct(cMatrix, fiducialMatrix, noiseMatrix, goodPixels, foreground);
}

void
Likelihood::construct(const CMatrix& cMatrix, const CMatrix& fiducialMatrix, const CMatrix& noiseMatrix, const char* maskFileName, const char* foregroundFileName)
{
    long nSideMask, nSideFore;
    std::vector<int> goodPixels;
    LaVectorDouble f;
    
    Utils::readMask(maskFileName, nSideMask, goodPixels);
    if(foregroundFileName != NULL)
    {
        readForeground(foregroundFileName, goodPixels, nSideFore, f);
    
        if(nSideMask != nSideFore)
        {
            StandardException exc;
            std::stringstream exceptionStr;
            exceptionStr << "Mask file " << maskFileName << " has nSide = " << nSideMask << " while the foreground file " << foregroundFileName << " has nSide = " << nSideFore << ". They need to be the same.";
            exc.set(exceptionStr.str());
            throw exc;
        }
    }
    
    construct(cMatrix, fiducialMatrix, noiseMatrix, goodPixels, f);
}

void
Likelihood::construct(const CMatrix& cMatrix, const CMatrix& fiducialMatrix, const CMatrix& noiseMatrix, const std::vector<int>& goodPixels, const LaVectorDouble& foreground)
{
    StandardException exc;
    
    goodPixels_ = goodPixels;
    f_ = foreground;
    
    if(cMatrix.getNPix() != goodPixels_.size())
    {
        std::stringstream exceptionStr;
        exceptionStr << "There are " << goodPixels_.size() << " unmasked pixels, however the covariance matrix in c.dat corresponds to " << cMatrix.getNPix() << ". Please generate the covariance matrix with the same mask.";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    if(fiducialMatrix.getNPix() != goodPixels_.size())
    {
        std::stringstream exceptionStr;
        exceptionStr << "There are " << goodPixels_.size() << " unmasked pixels, however the fiducial covariance matrix in c_fiducial.dat corresponds to " << fiducialMatrix.getNPix() << ". Please generate the fiducial covariance matrix with the same mask.";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    if(noiseMatrix.getNPix() != goodPixels_.size())
    {
        std::stringstream exceptionStr;
        exceptionStr << "There are " << goodPixels_.size() << " unmasked pixels, however the noise covariance matrix in c_noise.dat corresponds to " << noiseMatrix.getNPix() << ". Please generate the noise covariance matrix with the same mask.";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    cInv_.resize(goodPixels_.size(), goodPixels_.size());
    
#pragma omp parallel for default(shared)
    for(int i = 0; i < goodPixels_.size(); ++i)
    {
        for(int j = 0; j < goodPixels_.size(); ++j)
        {
            cInv_(i, j) = cMatrix.element(i, j) + fiducialMatrix.element(i, j) + noiseMatrix.element(i, j);
        }
    }
    
    //output_screen("LU factorization of the covariance matrix..." << std::endl);
    LaVectorLongInt pivotC(goodPixels_.size());
    LUFactorizeIP(cInv_, pivotC);
    //output_screen("OK" << std::endl);
    
    //output_screen("Calculating the determinant of c..." << std::endl);
    logDet_ = 0;
    int signC = 1;
    for(int i = 0; i < goodPixels_.size(); ++i)
    {
        if(Math::areEqual(cInv_(i, i), 0.0, 1e-15))
        {
            std::string exceptionStr = "The determinant of the covariance matrix is 0. The covariance matrix must be positive definite.";
            exc.set(exceptionStr);
            throw exc;
        }
        signC *= (cInv_(i, i) < 0 ? -1 : 1);
        logDet_ += std::log(std::abs(cInv_(i, i)));
    }
    /*if(signC != 1)
     {
     std::string exceptionStr = "The determinant of the covariance matrix is not positive. The covariance matrix must be positive definite.";
     exc.set(exceptionStr);
     throw exc;
     }*/
    
    const double detOffset = -29677.0566; // to be done better
    logDet_ -= detOffset;
    //output_screen("OK" << std::endl);
    
    //output_screen("Taking the inverse of the covariance matrix from LU factorization..." << std::endl);
    LaLUInverseIP(cInv_, pivotC);
    //output_screen("OK" << std::endl);
}

double
Likelihood::vmv(int n, const LaVectorDouble& a, const LaGenMatDouble& matrix, const LaVectorDouble& b) const
{
    std::vector<double> res(n, 0);
    
#pragma omp parallel for default(shared)
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
			res[i] += a(i) * matrix(i, j) * b(j);
	}
    
    double result = 0;
    for(int i = 0; i < n; ++i)
        result += res[i];
	return result;
}

double
Likelihood::calculate(const char* mapName, const char* noiseMapName, double& chi2, double& logDet) const
{
    LaVectorDouble t;
    long nSide;
    readMapAndNoise(mapName, noiseMapName, goodPixels_, nSide, t);
    return calculate(t, chi2, logDet);
}

double
Likelihood::calculate(const LaVectorDouble& t, double& chi2, double& logDet) const
{
    check(t.size() == goodPixels_.size(), "");
    
    chi2 = vmv(goodPixels_.size(), t, cInv_, t);
    logDet = logDet_;
    
    if(f_.size() > 0)
    {
        const double tCinvf = vmv(goodPixels_.size(), t, cInv_, f_);
        const double fCinvf = vmv(goodPixels_.size(), f_, cInv_, f_);
        logDet += std::log(fCinvf / goodPixels_.size());
        chi2 -= tCinvf * tCinvf / fCinvf;
    }
    
    return chi2 + logDet;
}

void
Likelihood::readMapAndNoise(const char* mapName, const char* noiseMapName, const std::vector<int>& goodPixels, long& nSide, LaVectorDouble& t)
{
    StandardException exc;
    
    Healpix_Map<double> map, noise;
    read_Healpix_map_from_fits(std::string(mapName), map);
    
    if(map.Scheme() != NEST)
    {
        std::string exceptionStr = "The map must have nested ordering.";
        exc.set(exceptionStr);
        throw exc;
    }

    read_Healpix_map_from_fits(std::string(noiseMapName), noise);
    
    if(noise.Scheme() != NEST)
    {
        std::string exceptionStr = "The noise map must have nested ordering.";
        exc.set(exceptionStr);
        throw exc;
    }
    
    if(map.Nside() != noise.Nside())
    {
        std::stringstream exceptionStr;
        exceptionStr << "Map and noise must have the same NSide, for the map it is " << map.Nside() << " and for the noise it is " << noise.Nside() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    nSide = map.Nside();
    
    t.resize(goodPixels.size());
    for(int i = 0; i < goodPixels.size(); ++i)
    {
        int index = goodPixels[i];
        t(i) = (double)(map[index] + noise[index]);
    }
}

void
Likelihood::readForeground(const char* foregroundFileName, const std::vector<int>& goodPixels, long& nSide, LaVectorDouble& f)
{
    Healpix_Map<double> fore;
    read_Healpix_map_from_fits(std::string(foregroundFileName), fore);
    
    if(fore.Scheme() != NEST)
    {
        StandardException exc;
        std::string exceptionStr = "The foreground map must have nested ordering.";
        exc.set(exceptionStr);
        throw exc;
    }

    nSide = fore.Nside();
    
    f.resize(goodPixels.size());
    for(int i = 0; i < goodPixels.size(); ++i)
    {
        int index = goodPixels[i];
        f(i) = (double)(fore[index]);
    }
}

void
Likelihood::readInput(const char* inputListName, const std::vector<int>& goodPixels, std::vector<LaVectorDouble>& t, std::vector<std::string>& mapNames)
{
    StandardException exc;
    std::ifstream inList(inputListName);
    if(!inList)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read the map list file " << inputListName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    int numOfMaps;
    inList >> numOfMaps;
    if(numOfMaps < 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The number of maps cannot be negative. It is " << numOfMaps << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    t.resize(numOfMaps);
    mapNames.resize(numOfMaps);
    
    long nSide;
    for(int i = 0; i < numOfMaps; ++i)
    {
        std::string mapName, noiseName;
        inList >> mapName >> noiseName;
        mapNames[i] = mapName;
        readMapAndNoise(mapName.c_str(), noiseName.c_str(), goodPixels, nSide, t[i]);
    }
}

void
Likelihood::calculateAll(const std::vector<LaVectorDouble>& t, const std::vector<std::string>& mapNames, std::vector<LikelihoodResult>& results) const
{
    StandardException exc;
    omp_lock_t lock;
    //output_screen("Calculating likelihood for all the maps..." << std::endl);
    const int numOfMaps = t.size();
    check(mapNames.size() == numOfMaps, "");
    std::vector<double> chi2Results(numOfMaps), logDetResults(numOfMaps);
    
    //ProgressMeter meter(numOfMaps);
    
    omp_init_lock(&lock);
    
//#pragma omp parallel for default(shared)
    for(int i = 0; i < numOfMaps; ++i)
    {
        double chi2, logDet;
        calculate(t[i], chi2, logDet);
        
        omp_set_lock(&lock);
        logDetResults[i] = logDet;
        chi2Results[i] = chi2;
        //meter.advance();
        omp_unset_lock(&lock);
    }
    //output_screen("OK" << std::endl);
    
    LikelihoodResult res;
    
    for(int i = 0; i < numOfMaps; ++i)
    {
        res.mapName = mapNames[i];
        res.logDet = logDetResults[i];
        res.chi2 = chi2Results[i];
        res.like = res.logDet + res.chi2;
        results.push_back(res);
    }
}

void
Likelihood::calculateAll(const char* inputListName, std::vector<LikelihoodResult>& results) const
{
    std::vector<std::string> mapNames;
    std::vector<LaVectorDouble> t;
    readInput(inputListName, goodPixels_, t, mapNames);
    check(mapNames.size() == t.size(), "");
    calculateAll(t, mapNames, results);
}

//LikelihoodPolarization
LikelihoodPolarization::LikelihoodPolarization(const CMatrix& cMatrix, long nSide, const std::vector<int>& goodPixels, int lMax, const WholeMatrix& etttInverse, double phi, double theta, double psi) : goodPixels_(goodPixels), nSide_(nSide), lMax_(lMax), etttInverse_(etttInverse), phi_(phi), theta_(theta), psi_(psi)
{
    StandardException exc;
    
    const int size = cMatrix.getNPix(), goodSize = goodPixels_.size();
    check(size % 2 == 0, "polarization c matrix includes q and u parts, so size must be even");
    check(size / 2 >= goodSize, "");
    
    std::vector<std::vector<double> > nInv(size);
    
    for(int i = 0; i < size; ++i)
        nInv[i].resize(size);
    
    std::ifstream in("n_inv.txt");
    if(!in)
    {
        std::string exceptionStr = "Cannot read the input file n_inv.txt.";
        exc.set(exceptionStr);
        throw exc;
    }
    
    for(int i = 0; i < size; ++i)
    {
        for(int j = 0; j < size; ++j)
        {
            in >> nInv[i][j];
        }
    }
    in.close();
    
    nInv_.resize(2 * goodSize, 2 * goodSize);
    LaGenMatDouble cMat(2 * goodSize, 2 * goodSize), nInvCMat(2 * goodSize, 2 * goodSize);
    cInv_.resize(2 * goodSize, 2 * goodSize);
    
    for(int i = 0; i < goodSize; ++i)
    {
        for(int j = 0; j < goodSize; ++j)
        {
            nInv_(i, j) = nInv[goodPixels_[i]][goodPixels_[j]];
            nInv_(i, goodSize + j) = nInv[goodPixels_[i]][size / 2 + goodPixels_[j]];
            nInv_(goodSize + i, j) = nInv[size / 2 + goodPixels_[i]][goodPixels_[j]];
            nInv_(goodSize + i, goodSize + j) = nInv[size / 2 + goodPixels_[i]][size / 2 + goodPixels_[j]];
            
            cMat(i, j) = cMatrix.element(goodPixels_[i], goodPixels_[j]);
            cMat(i, goodSize + j) = cMatrix.element(goodPixels_[i], size / 2 + goodPixels_[j]);
            cMat(goodSize + i, j) = cMatrix.element(size / 2 + goodPixels_[i], goodPixels_[j]);
            cMat(goodSize + i, goodSize + j) = cMatrix.element(size / 2 + goodPixels_[i], size / 2 + goodPixels_[j]);
        }
    }
    
    cInv_ = nInv_;
    Blas_Mat_Mat_Mult(nInv_, cMat, nInvCMat);
    Blas_Mat_Mat_Mult(nInvCMat, nInv_, cInv_, 1.0, 1.0);
    
    LaVectorLongInt pivotC(2 * goodSize);
    LUFactorizeIP(cInv_, pivotC);
    
    //output_screen("Calculating the determinant of c..." << std::endl);
    logDet_ = 0;
    int signC = 1;
    for(int i = 0; i < 2 * goodSize; ++i)
    {
        if(Math::areEqual(cInv_(i, i), 0.0, 1e-15))
        {
            std::string exceptionStr = "The determinant of the covariance matrix is 0. The covariance matrix must be positive definite.";
            exc.set(exceptionStr);
            throw exc;
        }
        signC *= (cInv_(i, i) < 0 ? -1 : 1);
        logDet_ += std::log(std::abs(cInv_(i, i)));
    }
    
    const double detOffset = 16078.083180; // to be done better
    logDet_ -= detOffset;
    LaLUInverseIP(cInv_, pivotC);
    
    Utils::readPixelWindowFunction(beam_, nSide, lMax, 0);
}

void
LikelihoodPolarization::readMaps(const char* qMapName, const char* uMapName, const char* almTTFileName, int lMax, const std::vector<int>& goodPixels, std::vector<double>& v, AlmType& alm)
{
    StandardException exc;
    
    long qNSide, uNSide;
    char coordSys[20], ordering[20];
    
    float* qMap = read_healpix_map(qMapName, &qNSide, coordSys, ordering);
    
    if(ordering[0] != 'N')
    {
        std::string exceptionStr = "The map must have nested ordering.";
        exc.set(exceptionStr);
        throw exc;
    }
    
    float* uMap = read_healpix_map(uMapName, &uNSide, coordSys, ordering);
    
    if(ordering[0] != 'N')
    {
        std::string exceptionStr = "The map must have nested ordering.";
        exc.set(exceptionStr);
        throw exc;
    }
    
    if(qNSide != uNSide)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Q map has NSide = " << qNSide << " while the U map has NSide = " << uNSide << ". They need to be the same.";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    v.clear();
    v.resize(2 * goodPixels.size());
    for(int i = 0; i < goodPixels.size(); ++i)
    {
        v[i] = qMap[goodPixels[i]];
        v[goodPixels.size() + i] = uMap[goodPixels[i]];
    }
    
    delete qMap;
    delete uMap;
    
    std::ifstream inTT(almTTFileName);
    if(!inTT)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read the input file " << almTTFileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    alm.Set(lMax, lMax);
    alm.SetToZero();
    for(int l = 0; l <= lMax; ++l)
    {
        for(int m = 0; m <= l; ++m)
        {
            ComplexDouble x;
            inTT >> x;
            alm(l, m) = x;
        }
    }
    
    inTT.close();
}

void
LikelihoodPolarization::readInput(const char* inputListName, const std::vector<int>& goodPixels, std::vector<LaVectorDouble>& t, std::vector<std::string>& mapNames, const std::vector<int>& goodPixelsPol, std::vector<std::vector<double> > & v, std::vector<AlmType>& alm, int lMaxPol)
{
    StandardException exc;
    std::ifstream inList(inputListName);
    if(!inList)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot read the map list file " << inputListName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    int numOfMaps;
    inList >> numOfMaps;
    if(numOfMaps < 0)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The number of maps cannot be negative. It is " << numOfMaps << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    t.resize(numOfMaps);
    mapNames.resize(numOfMaps);
    v.clear();
    v.resize(numOfMaps);
    alm.clear();
    alm.resize(numOfMaps);
    
    long nSide;
    for(int i = 0; i < numOfMaps; ++i)
    {
        std::string mapName, noiseName, almName, qName, uName;
        inList >> mapName >> noiseName >> almName >> qName >> uName;
        mapNames[i] = mapName;
        Likelihood::readMapAndNoise(mapName.c_str(), noiseName.c_str(), goodPixels, nSide, t[i]);
        readMaps(qName.c_str(), uName.c_str(), almName.c_str(), lMaxPol, goodPixelsPol, v[i], alm[i]);
    }
}

void
LikelihoodPolarization::calculateAll(const std::vector<std::vector<double> >& v, const std::vector<AlmType>& alm, const std::vector<std::string>& mapNames, std::vector<LikelihoodResult>& results) const
{
    int n = v.size();
    check(n == alm.size(), "");
    check(n == mapNames.size(), "");
    
    LikelihoodResult res;
    for(int i = 0; i < n; ++i)
    {
        res.mapName = mapNames[i];
        res.like = calculate(v[i], alm[i], res.chi2, res.logDet);
        check(res.like == res.chi2 + res.logDet, "");
        results.push_back(res);
    }
}

double
LikelihoodPolarization::calculate(const std::vector<double>& v, const AlmType& alm, double& chi2, double& logDet) const
{
    const Math::ThreeRotationMatrix rot(phi_, theta_, psi_);
    const rotmatrix rotationMatrix(rot[0][0], rot[0][1], rot[0][2], rot[1][0], rot[1][1], rot[1][2], rot[2][0], rot[2][1], rot[2][2]);
    
    AlmType almCopy(alm);
    rotate_alm(almCopy, rotationMatrix);
    
    std::vector<double> vCopy(v);
    
    check(vCopy.size() == 2 * goodPixels_.size(), "");
    
    AlmType almE(lMax_, lMax_);
    almE.SetToZero();
    
    for(int l = etttInverse_.getLMin(); l <= lMax_; ++l)
    {
#pragma omp parallel for default(shared)
        for(int m = 0; m <= l; ++m)
        {
            for(int l1 = etttInverse_.getLMin(); l1 <= lMax_; ++l1)
            {
                for(int m1 = -l1; m1 <= l1; ++m1)
                {
                    const int m1Abs = (m1 >= 0 ? m1 : -m1);
                    xcomplex<double> a = almCopy(l1, m1Abs);
                    if(m1 < 0)
                    {
                        a = a.conj();
                        if(m1Abs % 2)
                            a *= -1.0;
                    }
                    almE(l, m) += etttInverse_.element(l, m, l1, m1) * a;
                }
            }
        }
    }
    
    const rotmatrix rotationMatrixInv(rot[0][0], rot[1][0], rot[2][0], rot[0][1], rot[1][1], rot[2][1], rot[0][2], rot[1][2], rot[2][2]);
    AlmType alm1(lMax_, lMax_), alm2(lMax_, lMax_);
    alm1.SetToZero();
    alm2.SetToZero();
    
    rotate_alm(alm1, almE, alm2, rotationMatrixInv);
    
    Healpix_Map<double> mapT, mapQ, mapU;
    mapT.SetNside(nSide_, RING);
    mapQ.SetNside(nSide_, RING);
    mapU.SetNside(nSide_, RING);
    alm2map_pol(alm1, almE, alm2, mapT, mapQ, mapU);
    
    mapQ.swap_scheme();
    mapU.swap_scheme();
    
    std::vector<double> subtractMap;
    for(int i = 0; i < goodPixels_.size(); ++i)
        subtractMap.push_back(mapQ[goodPixels_[i]]);
    for(int i = 0; i < goodPixels_.size(); ++i)
        subtractMap.push_back(mapU[goodPixels_[i]]);
    
    for(int i = 0; i < 2 * goodPixels_.size(); ++i)
    {
        for(int j = 0; j < 2 * goodPixels_.size(); ++j)
            vCopy[i] -= nInv_(i, j) * subtractMap[j];
    }
    
    chi2 = 0;
    
    for(int i = 0; i < 2 * goodPixels_.size(); ++i)
    {
        for(int j = 0; j < 2 * goodPixels_.size(); ++j)
            chi2 += vCopy[i] * cInv_(i, j) * vCopy[j];
    }
    
    logDet = logDet_;
    return chi2 + logDet;
}

double
LikelihoodPolarization::calculate(const char* qMapName, const char* uMapName, const char* almTTFileName, double& chi2, double& logDet) const
{
    StandardException exc;
    
    std::vector<double> v;
    AlmType alm;
    readMaps(qMapName, uMapName, almTTFileName, lMax_, goodPixels_, v, alm);
    
    return calculate(v, alm, chi2, logDet);
}

double
LikelihoodPolarization::getNInv(int i, int j) const
{
    check(i >= 0 && i < 2 * goodPixels_.size(), "");
    check(j >= 0 && j < 2 * goodPixels_.size(), "");
    
    return nInv_(i, j);
}

int myIndex(int l, int m, int lMin)
{
    return l * (l + 1) + m - (lMin * (lMin + 1) - lMin);
}

void
LikelihoodPolarization::combineWholeMatrices(const WholeMatrix& tt, const WholeMatrix& te, const WholeMatrix& ee, WholeMatrix& combined, WholeMatrix& etttInverse)
{
    StandardException exc;
    
    const int lMin = combined.getLMin(), lMax = combined.getLMax();
    check(etttInverse.getLMin() == lMin, "");
    check(etttInverse.getLMax() == lMax, "");
    
    //inverting the tt matrix
    int size = myIndex(lMax, lMax, lMin) + 1;
    LaGenMatDouble ttMat(size, size);
    for(int l1 = lMin; l1 <= lMax; ++l1)
        for(int m1 = -l1; m1 <= l1; ++m1)
            for(int l = lMin; l <= lMax; ++l)
                for(int m = -l; m <= l; ++m)
                {
                    int i = myIndex(l, m, lMin);
                    int j = myIndex(l1, m1, lMin);
                    check(i < size, "");
                    check(j < size, "");
                    ttMat(i, j) = tt.element(l1, m1, l, m);
                    
                    check(Math::areEqual(tt.element(l1, m1, l, m), tt.element(l, m, l1, m1), 1e-10), l1 << ' ' << m1 << ' ' << l << ' ' << m << ' ' << tt.element(l1, m1, l, m) << ' ' << tt.element(l, m, l1, m1));
                }
    
    LaVectorLongInt pivot(size);
    LUFactorizeIP(ttMat, pivot);
    LaLUInverseIP(ttMat, pivot);
    
    WholeMatrix ttInverse(lMin, lMax);
    for(int l1 = lMin; l1 <= lMax; ++l1)
        for(int m1 = -l1; m1 <= l1; ++m1)
            for(int l = lMin; l <= lMax; ++l)
                for(int m = -l; m <= l; ++m)
                {
                    int i = myIndex(l, m, lMin);
                    int j = myIndex(l1, m1, lMin);
                    check(i < size, "");
                    check(j < size, "");
                    ttInverse.element(l1, m1, l, m) = ttMat(i, j);
                }
    
    //Multiplying matrices
    for(int l1 = lMin; l1 <= lMax; ++l1)
        for(int m1 = -l1; m1 <= l1; ++m1)
            for(int l = lMin; l <= lMax; ++l)
                for(int m = -l; m <= l; ++m)
                {
                    etttInverse.element(l1, m1, l, m) = 0;
                    for(int l2 = lMin; l2 <= lMax; ++l2)
                        for(int m2 = -l2; m2 <= l2; ++m2)
                            etttInverse.element(l1, m1, l, m) += te.element(l2, m2, l1, m1) * ttInverse.element(l2, m2, l, m);
                    
                    check(Math::areEqual(te.element(l1, m1, l, m), te.element(l1, -m1, l, -m), 1e-10), l1 << ' ' << m1 << ' ' << l << ' ' << m << ' ' << te.element(l1, m1, l, m) << ' ' << te.element(l1, -m1, l, -m));
                    check(Math::areEqual(ttInverse.element(l1, m1, l, m), ttInverse.element(l1, -m1, l, -m), 1e-10), l1 << ' ' << m1 << ' ' << l << ' ' << m << ' ' << ttInverse.element(l1, m1, l, m) << ' ' << ttInverse.element(l1, -m1, l, -m));
                    check(Math::areEqual(ttInverse.element(l1, m1, l, m), ttInverse.element(l, m, l1, m1), 1e-10), l1 << ' ' << m1 << ' ' << l << ' ' << m << ' ' << ttInverse.element(l1, m1, l, m) << ' ' << ttInverse.element(l, m, l1, m1));
                }
    
    for(int l1 = lMin; l1 <= lMax; ++l1)
        for(int m1 = -l1; m1 <= l1; ++m1)
            for(int l = lMin; l <= lMax; ++l)
                for(int m = -l; m <= l; ++m)
                {
                    combined.element(l1, m1, l, m) = ee.element(l1, m1, l, m);
                    for(int l2 = lMin; l2 <= lMax; ++l2)
                        for(int m2 = -l2; m2 <= l2; ++m2)
                            combined.element(l1, m1, l, m) -= etttInverse.element(l1, m1, l2, m2) * te.element(l2, m2, l, m);
                    
                    check(Math::areEqual(etttInverse.element(l1, m1, l, m), etttInverse.element(l1, -m1, l, -m), 1e-8), l1 << ' ' << m1 << ' ' << l << ' ' << m << ' ' << etttInverse.element(l1, m1, l, m) << ' ' << etttInverse.element(l1, -m1, l, -m));
                }
    
#ifdef CHECKS_ON
    WholeMatrix etttInversett(lMin, lMax);
    for(int l1 = lMin; l1 <= lMax; ++l1)
        for(int m1 = -l1; m1 <= l1; ++m1)
            for(int l = lMin; l <= lMax; ++l)
                for(int m = -l; m <= l; ++m)
                {
                    etttInversett.element(l1, m1, l, m) = 0;
                    for(int l2 = lMin; l2 <= lMax; ++l2)
                        for(int m2 = -l2; m2 <= l2; ++m2)
                            etttInversett.element(l1, m1, l, m) += etttInverse.element(l1, m1, l2, m2) * tt.element(l2, m2, l, m);
                    
                    check(Math::areEqual(combined.element(l1, m1, l, m), combined.element(l, m, l1, m1), 1e-7), l1 << ' ' << m1 << ' ' << l << ' ' << m << ' ' << combined.element(l1, m1, l, m) << ' ' << combined.element(l, m, l1, m1));
                    
                    check(Math::areEqual(etttInversett.element(l1, m1, l, m), te.element(l, m, l1, m1), 1e-7), l1 << ' ' << m1 << ' ' << l << ' ' << m << ' ' << etttInversett.element(l1, m1, l, m) << ' ' << te.element(l, m, l1, m1));
                }
#endif
}

void readClFromFile(const char* fileName, std::vector<double>& cl, int lMax)
{
    StandardException exc;
    std::ifstream in;
    in.open(fileName);
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot open input file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    
    cl.resize(lMax + 1, 0);
    
    while(!in.eof())
    {
        std::string s;
        std::getline(in, s);
        if(s == "")
            break;
        
        if(s[0] == '#')
            continue;
        
        int l;
        double c;
        std::stringstream str(s);
        str >> l >> c;
        if(l < 0)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Invalid l = " << l << ". It must be non-negative.";
            exc.set(exceptionStr.str());
            throw exc;
        }
        
        if(l <= lMax)
            cl[l] = c;
    }
    in.close();

}

void
LikelihoodHigh::construct(const char* couplingKernelFileName)
{
    StandardException exc;
    std::ifstream in(couplingKernelFileName);
    if(!in)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot open input file " << couplingKernelFileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    int lMaxCK;
    in >> lMaxCK;
    if(lMaxCK < lMax_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "The coupling kernel has lMax = " << lMaxCK << ", needed up to lMax = " << lMax_ << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    coupling_.resize(lMaxCK + 1);
    for(int l1 = 0; l1 <= lMaxCK; ++l1)
    {
        coupling_[l1].resize(lMaxCK + 1);
        for(int l2 = 0; l2 <= lMaxCK; ++l2)
            in >> coupling_[l1][l2];
    }
    in.close();
}

LikelihoodHigh::LikelihoodHigh(const std::vector<double>& cl, const std::vector<double>& nl, const char* couplingKernelFileName, int lMin, int lMax) : lMin_(lMin), lMax_(lMax), offset_((lMax_ - lMin_ + 1) * std::log(2 * Math::pi)), cl_(cl), nl_(nl)
{
    check(lMin_ >= 0, "invalid lMin");
    check(lMax_ >= lMin_, "invalid lMax");

    check(cl_.size() - 1 >= lMax_, "");
    check(nl_.size() - 1 >= lMax_, "");

    construct(couplingKernelFileName);
}

LikelihoodHigh::LikelihoodHigh(const char* dataClFileName, const char* noiseClFileName, const char* couplingKernelFileName, int lMin, int lMax) : lMin_(lMin), lMax_(lMax), offset_((lMax_ - lMin_ + 1) * std::log(2 * Math::pi))
{
    check(lMin_ >= 0, "invalid lMin");
    check(lMax_ >= lMin_, "invalid lMax");

    readClFromFile(dataClFileName, cl_, lMax_);
    readClFromFile(noiseClFileName, nl_, lMax_);

    construct(couplingKernelFileName);
}

double
LikelihoodHigh::calculate(const char* clFileName) const
{
    std::vector<double> cl;
    readClFromFile(clFileName, cl, lMax_);
    return calculate(cl);
}

double
LikelihoodHigh::calculate(const std::vector<double>& cl) const
{
    check(cl.size() >= lMax_ + 1, "");

    double res = 0.0;
    std::vector<double> resVec(lMax_ + 1, 0);

#pragma omp parallel for default(shared)
    for(int l1 = lMin_; l1 <= lMax_; ++l1)
    {
        const double clTot1 = cl[l1] + nl_[l1];
        check(clTot1 != 0, "");
        const double deltaCl1 = cl[l1] - cl_[l1];
        //const double deltaLn1 = (cl_[l1] + nl_[l1] > 0 ? std::log(cl_[l1] + nl_[l1]) - std::log(clTot1) : 0.0);
        for(int l2 = lMin_; l2 <= lMax_; ++l2)
        {
            const double clTot2 = cl[l2] + nl_[l2];
            check(clTot2 != 0, "");
            const double deltaCl2 = cl[l2] - cl_[l2];
            //const double deltaLn2 = (cl_[l2] + nl_[l2] > 0 ? std::log(cl_[l2] + nl_[l2]) - std::log(clTot2) : 0.0);

            const double f = double(2 * l1 + 1) * coupling_[l1][l2] / (2.0 * clTot1 * clTot2);
            //const double lnF = f * clTot1 * clTot2;
            //const double alpha = 1.0; //(deltaLn1 != 0 && deltaLn2 != 0 ? 1.0 / 3.0 : 1.0);
            //chi2 += alpha * deltaCl1 * deltaCl2 * f + (1 - alpha) * deltaLn1 * deltaLn2 * lnF;
            resVec[l1] += deltaCl1 * deltaCl2 * f;
            //output_screen(l1 << ' ' << l2 << ' ' << ' ' << coupling_[l1][l2] << ' ' << clTot1 << ' ' << clTot2 << ' ' << ' ' << deltaCl1 << ' ' << deltaCl2 << ' ' << chi2 << std::endl);
        }
    }

    for(int l1 = lMin_; l1 <= lMax_; ++l1)
        res += resVec[l1];

    return res;
}

