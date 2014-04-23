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

#include <omp.h>

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

void
MaskApodizer::apodize(ApodizationType type, double angle, Healpix_Map<double>& result) const
{
    check(type >= COSINE_APODIZATION && type < APODIZATION_TYPE_MAX, "invalid apodization type");
    check(angle > 0, "invalid angle = " << angle);

    Healpix_Map<double> edgeMap;
    edgeMap.SetNside(256, NEST);

    result.SetNside(mask_.Nside(), mask_.Scheme());
    result.Import(mask_);

    const long nPix = result.Npix();
    const double pixSize = std::sqrt(4 * Math::pi / nPix);

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
        double theta, phi;
        result.Scheme() == NEST ? pix2ang_nest(result.Nside(), i, &theta, &phi) : pix2ang_ring(result.Nside(), i, &theta, &phi);
        long index;
        result.Scheme() == NEST ? ang2pix_nest(result.Nside(), correctTheta(theta + pixSize), phi, &index) : ang2pix_nest(result.Nside(), correctTheta(theta + pixSize), phi, &index);
        if(result[index] == 0)
            isOnEdge = true;
        result.Scheme() == NEST ? ang2pix_nest(result.Nside(), correctTheta(theta - pixSize), phi, &index) : ang2pix_nest(result.Nside(), correctTheta(theta - pixSize), phi, &index);
        if(result[index] == 0)
            isOnEdge = true;
        result.Scheme() == NEST ? ang2pix_nest(result.Nside(), theta, phi + pixSize, &index) : ang2pix_nest(result.Nside(), theta, phi + pixSize, &index);
        if(result[index] == 0)
            isOnEdge = true;
        result.Scheme() == NEST ? ang2pix_nest(result.Nside(), theta, phi - pixSize, &index) : ang2pix_nest(result.Nside(), theta, phi - pixSize, &index);
        if(result[index] == 0)
            isOnEdge = true;
        if(isOnEdge)
        {
            edge.push_back(Math::ThreeVectorDouble(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)));
            ang2pix_nest(256, theta, phi, &index);
            edgeMap[index] = 1;
        }
        meter.advance();
    }
    output_screen("OK" << std::endl);

    output_screen("Found " << edge.size() << " edge pixels." << std::endl);

    output_screen("Apodizing mask..." << std::endl);
    ProgressMeter meter1(nPix);
    omp_lock_t lock;
    omp_init_lock(&lock);
#pragma omp parallel for default(shared)
    for(long i = 0; i < nPix; ++i)
    {
        if(result[i] == 0)
        {
            omp_set_lock(&lock);
            meter1.advance();
            omp_unset_lock(&lock);
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
        omp_set_lock(&lock);
        meter1.advance();
        omp_unset_lock(&lock);
    }
}
