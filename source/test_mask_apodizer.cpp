#include <cmath>

#include <macros.hpp>
#include <mask_apodizer.hpp>
#include <math_constants.hpp>
#include <test_mask_apodizer.hpp>

#include <chealpix.h>

std::string
TestMaskApodizer::name() const
{
    return std::string("CMB TESTER");
}

unsigned int
TestMaskApodizer::numberOfSubtests() const
{
    return 2;
}

void
TestMaskApodizer::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 2, "invalid index " << i);

    using namespace Math;

    const unsigned long nSide = 2048;

    Healpix_Map<double> mask;
    mask.SetNside(nSide, RING);
    for(unsigned long i = 0; i < mask.Npix(); ++i)
        mask[i] = 1;

    for(unsigned long i = 0; i < mask.Npix(); ++i)
    {
        double theta, phi;
        pix2ang_ring(nSide, i, &theta, &phi);

        if(theta < pi / 2 + pi / 3 && theta > pi / 2 - pi / 3)
            mask[i] = 0;
    }

    MaskApodizer ap(mask);
    const double apAngle = pi / 10;
    Healpix_Map<double> apodizedMask;
    ap.apodize((i == 0 ? MaskApodizer::COSINE_APODIZATION : MaskApodizer::GAUSSIAN_APODIZATION), apAngle, apodizedMask);

    double theta1 = pi / 2 + pi / 3 + apAngle / 2;
    long index;
    ang2pix_ring(nSide, theta1, 2.0, &index);
    res = apodizedMask[index];
    expected = 1.0 - (i == 0 ? std::cos(pi / 4) : std::exp(-9.0 / 8.0));

    if(i == 0)
        subTestName = "cos";
    else
        subTestName = "gauss";
}
