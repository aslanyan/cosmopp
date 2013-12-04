#ifndef COSMO_CPP_MASK_APODIZER_HPP
#define COSMO_CPP_MASK_APODIZER_HPP

#include <healpix_map.h>

/// A mask apodizing class.
class MaskApodizer
{
public:
    /// The apodization type (cosine or gaussian).
    enum ApodizationType { COSINE_APODIZATION = 0, GAUSSIAN_APODIZATION, APODIZATION_TYPE_MAX };

    /// Constructor.
    /// \param originalMask The original mask to be apodized.
    MaskApodizer(const Healpix_Map<double>& originalMask) : mask_(originalMask) {}

    /// Perform the apodization.
    /// \param type The apodization type.
    /// \param angle The apodization angle.
    /// \param result The resulting apodized mask will be written here.
    void apodize(ApodizationType type, double angle, Healpix_Map<double>& result) const;

private:
    double correctTheta(double theta) const;
    double cosApodization(double sigma, double x) const;
    double gaussApodization(double sigma, double x) const;

private:
    const Healpix_Map<double>& mask_;
};

#endif

