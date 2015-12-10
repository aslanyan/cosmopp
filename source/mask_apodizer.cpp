#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <progress_meter.hpp>
#include <three_vector.hpp>
#include <math_constants.hpp>
#include <mask_apodizer.hpp>

#include <healpix_map.h>
#include <chealpix.h>
#include <healpix_base.h>

#include <fitshandle.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#ifdef COSMO_OMP
#include <omp.h>
#endif

double
MaskApodizer::correctTheta(double theta) const
{
    if(theta < 0)
        return -theta;
    if(theta > Math::pi)
        return 2 * Math::pi - theta;
    return theta;
}

double
MaskApodizer::cosApodization(double sigma, double x) const
{
   if(x > sigma)
       return 1;
   return 1 - std::cos(Math::pi * x / (2 * sigma));
}

double
MaskApodizer::gaussApodization(double sigma, double x) const
{
   if(x > sigma)
       return 1;
    check(sigma > 0, "invalid sigma = " << sigma);
    return 1 - std::exp(- 9 * x * x / (2 * sigma * sigma));
}

/*
double
MaskApodizer::distance(const Healpix_Map<double> &map, long i, long j) const
{
    check(map.Scheme() == NEST, "");
    double theta1, phi1, theta2, phi2;
    pix2ang_nest(map.Nside(), i, &theta1, &phi1);
    pix2ang_nest(map.Nside(), j, &theta2, &phi2);
    const double x1 = std::sin(theta1) * std::cos(phi1), y1 = std::sin(theta1) * std::sin(phi1), z1 = std::cos(theta1);
    const double x2 = std::sin(theta2) * std::cos(phi2), y2 = std::sin(theta2) * std::sin(phi2), z2 = std::cos(theta2);

    return std::acos(x1 * x2 + y1 * y2 + z1 * z2);
}

void
MaskApodizer::findNeighbors(const Healpix_Base2 &base2, const Healpix_Map<double> &map, long i, std::unordered_set<long> *region, double maxAngle) const
{
    region->clear();

    std::queue<long> q;
    region->insert(i);
    q.push(i);

    while(!q.empty())
    {
        const long j = q.front();
        q.pop();
        fix_arr<int64, 8> neighbors;
        base2.neighbors(j, neighbors);
        for(int k = 0; k < neighbors.size(); ++k)
        {
            if(region->count(neighbors[k]) == 1)
                continue;
            if(distance(map, i, neighbors[k]) > maxAngle)
                continue;
            q.push(neighbors[k]);
            region->insert(neighbors[k]);
        }
    }
}
*/

void
MaskApodizer::apodize(ApodizationType type, double angle, Healpix_Map<double>& result) const
{
    check(type >= COSINE_APODIZATION && type < APODIZATION_TYPE_MAX, "invalid apodization type");
    check(angle > 0, "invalid angle = " << angle);

    const int edgeNSide = std::min(128, mask_.Nside());
    //const int edgeNSide = mask_.Nside();

    Healpix_Map<double> edgeMap;
    edgeMap.SetNside(edgeNSide, NEST);

    for(int i = 0; i < edgeMap.Npix(); ++i)
        edgeMap[i] = 0;

    result.SetNside(mask_.Nside(), mask_.Scheme());
    result.Import(mask_);

    /*
    bool swapScheme = false;
    if(result.Scheme() == RING)
    {
        result.swap_scheme();
        swapScheme = true;
    }
    */

    const long nPix = result.Npix();
    const double pixSize = std::sqrt(4 * Math::pi / nPix);

    Healpix_Base2 base2(result.Nside(), result.Scheme(), SET_NSIDE);

    output_screen("Input mask has " <<nPix << " pixels." << std::endl);
    
    output_screen("Finding the edge pixels..." << std::endl);
    std::vector<Math::ThreeVectorDouble> edge;
    ProgressMeter meter(nPix);
    for(long i = 0; i < nPix; ++i)
    {
        if(result[i] == 0)
        {
            meter.advance();
            continue;
        }
        bool isOnEdge = false;
        fix_arr<int64, 8> neighbors;
        base2.neighbors(i, neighbors);
        const long s = neighbors.size();
        for(int j = 0; j < s; ++j)
        {
            if(result[neighbors[j]] == 0)
                isOnEdge = true;
        }
        double theta, phi;
        result.Scheme() == NEST ? pix2ang_nest(result.Nside(), i, &theta, &phi) : pix2ang_ring(result.Nside(), i, &theta, &phi);
        long index;
        /*
        result.Scheme() == NEST ? ang2pix_nest(result.Nside(), correctTheta(theta + pixSize), phi, &index) : ang2pix_ring(result.Nside(), correctTheta(theta + pixSize), phi, &index);
        if(result[index] == 0)
            isOnEdge = true;
        result.Scheme() == NEST ? ang2pix_nest(result.Nside(), correctTheta(theta - pixSize), phi, &index) : ang2pix_ring(result.Nside(), correctTheta(theta - pixSize), phi, &index);
        if(result[index] == 0)
            isOnEdge = true;
        result.Scheme() == NEST ? ang2pix_nest(result.Nside(), theta, phi + pixSize, &index) : ang2pix_ring(result.Nside(), theta, phi + pixSize, &index);
        if(result[index] == 0)
            isOnEdge = true;
        result.Scheme() == NEST ? ang2pix_nest(result.Nside(), theta, phi - pixSize, &index) : ang2pix_ring(result.Nside(), theta, phi - pixSize, &index);
        if(result[index] == 0)
            isOnEdge = true;
        */
        if(isOnEdge)
        {
            edge.push_back(Math::ThreeVectorDouble(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)));
            ang2pix_nest(edgeNSide, theta, phi, &index);
            edgeMap[index] = 1;
            //edgeMap[i] = 1;
        }
        meter.advance();
    }
    output_screen("OK" << std::endl);

    output_screen("Found " << edge.size() << " edge pixels." << std::endl);

    fitshandle outh;
    outh.create("edge.fits");
    write_Healpix_map_to_fits(outh, edgeMap, PLANCK_FLOAT64);

    output_screen("Apodizing mask..." << std::endl);
    ProgressMeter meter1(nPix);
#ifdef COSMO_OMP
    omp_lock_t lock;
    omp_init_lock(&lock);
#endif
#pragma omp parallel for default(shared)
    for(long i = 0; i < nPix; ++i)
    {
        if(result[i] == 0)
        {
#ifdef COSMO_OMP
            omp_set_lock(&lock);
#endif
            meter1.advance();
#ifdef COSMO_OMP
            omp_unset_lock(&lock);
#endif
            continue;
        }

        double theta, phi;
        result.Scheme() == NEST ? pix2ang_nest(result.Nside(), i, &theta, &phi) : pix2ang_ring(result.Nside(), i, &theta, &phi);
        Math::ThreeVectorDouble current(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
        double minDistance = 2 * Math::pi;
        for(int j = 0; j < edge.size(); ++j)
        {
            const double distance = std::acos(current * edge[j]);
            if(distance < minDistance)
                minDistance = distance;
        }
        result[i] = (type == COSINE_APODIZATION ? cosApodization(angle, minDistance) : gaussApodization(angle, minDistance));

#ifdef COSMO_OMP
        omp_set_lock(&lock);
#endif
        meter1.advance();
#ifdef COSMO_OMP
        omp_unset_lock(&lock);
#endif
    }
}
