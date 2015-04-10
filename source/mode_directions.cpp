#include <macros.hpp>
#include <exception_handler.hpp>
#include <angular_coordinates.hpp>
#include <three_rotation.hpp>
#include <progress_meter.hpp>
#include <math_constants.hpp>
#include <numerics.hpp>
#include <mode_directions.hpp>

#include <healpix_base.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <xcomplex.h>
#include <healpix_map.h>
#include <rotmatrix.h>
#include <alm_powspec_tools.h>
#include <chealpix.h>

#ifdef COSMO_OMP
#include <omp.h>
#endif

ModeDirections::ModeDirections(const Alm<xcomplex<double> >& alm) : alm_(alm.Lmax(), alm.Mmax())
{
    check(alm.Lmax() == alm.Mmax(), "");
    
    for(int l = 0; l <= alm.Lmax(); ++l)
    {
        for(int m = 0; m <= l; ++m)
            alm_(l, m) = alm(l, m);
    }
}

double
ModeDirections::calculateAngularMomentumDispersion(int l, double theta, double phi) const
{
    check(l >= 0 && l <= alm_.Lmax(), "invalid l");
    
    Alm<xcomplex<double> > almCopy(alm_.Lmax(), alm_.Mmax());
    check(alm_.Lmax() == alm_.Mmax(), "");
    for(int l1 = 0; l1 <= alm_.Lmax(); ++l1)
    {
        for(int m1 = 0; m1 <= l1; ++m1)
            almCopy(l1, m1) = alm_(l1, m1);
    }
    
    const Math::ThreeRotationMatrix rot(phi + Math::pi / 2, theta, 0);
    //const rotmatrix rotationMatrix(rot[0][0], rot[1][0], rot[2][0], rot[0][1], rot[1][1], rot[2][1], rot[0][2], rot[1][2], rot[2][2]);
    const rotmatrix rotationMatrix(rot[0][0], rot[0][1], rot[0][2], rot[1][0], rot[1][1], rot[1][2], rot[2][0], rot[2][1], rot[2][2]);
    
    rotate_alm(almCopy, rotationMatrix);
    
    double res = 0;
    for(int m = 1; m <= l; ++m)
    {
        xcomplex<double> x = almCopy(l, m);
        res += m * m * (x.real() * x.real() + x.imag() * x.imag());
    }
    
    return res;
}

void
ModeDirections::maximizeAngularMomentumDispersion(int l, long nSide, double& theta, double& phi, Healpix_Map<double>* map) const
{
    long nPix = nside2npix(nSide);
    double max = 0;
    output_screen("Maximizing over directions..." << std::endl);
    
    if(map)
        map->SetNside(nSide, NEST);
    
    ProgressMeter meter(nPix);
    
#ifdef COSMO_OMP
    omp_lock_t lock;
    omp_init_lock(&lock);
#endif

#pragma omp parallel for default(shared)
    for(long i = 0; i < nPix; ++i)
    {
        double t, p;
        pix2ang_nest(nSide, i, &t, &p);
        double disp = calculateAngularMomentumDispersion(l, t, p);
        
#ifdef COSMO_OMP
        omp_set_lock(&lock);
#endif
        if(map)
            (*map)[i] = disp;
        
        if(disp > max)
        {
            max = disp;
            theta = t;
            phi = p;
        }
        meter.advance();
#ifdef COSMO_OMP
        omp_unset_lock(&lock);
#endif
    }
    
    if(theta < Math::pi / 2)
    {
        theta = Math::pi - theta;
        phi += Math::pi;
        if(phi > 2 * Math::pi)
            phi -= 2 * Math::pi;
        
        check(Math::areEqual(max, calculateAngularMomentumDispersion(l, theta, phi), 1e-7), "");
    }
    output_screen("OK" << std::endl);
}
