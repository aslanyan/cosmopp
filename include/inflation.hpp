#ifndef COSMO_PP_INFLATION_HPP
#define COSMO_PP_INFLATION_HPP

#include <function.hpp>
#include <table_function.hpp>

/// Single Field Inflationary Models.

/// A class for numerically calculating the scalar and tensor power spectra for single field inflationary models.
class SingleFieldInflation
{
public:
    /// Constructor.
    /// \param potential The inflationary potential. Both the field and the potential are unitless, assuming the reduced Planck mass is 1.
    /// \param potentialDeriv The first derivative of the inflationary potential. Both the field and the potential are unitless, assuming the reduced Planck mass is 1.
    /// \param phi0 The starting value of the field.
    /// \param phiDot0 The starting value of the derivative of the field.
    /// \param phiEnd The value of the field at which inflation ends.
    /// \param kPivot The pivot scale in inverse Megaparsecs.
    /// \param NPivot The number of e-folds for the pivot scale. If this is 0, 55 e-folds will be assumed by default.
    /// \param a0 The starting value of the scale factor of the universe. This does not affect the results unless the universe is not flat.
    /// \param curvature The curvature of space. Should be 0 or +/-1.
    SingleFieldInflation(const Math::RealFunction& potential, const Math::RealFunction& potentialDeriv, double phi0, double phiDot0, double phiEnd, double kPivot = 0.05, double NPivot = 0, double a0 = 1, int curvature = 0);

    /// The field value as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& phi() const { return phi_; }

    /// The field derivative as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& phiDot() const { return phiDot_; }

    /// The Hubble parameter as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& H() const { return H_; }

    /// The scale factor as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& a() const { return a_; }

    /// The natural log of the scale factor as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& logA() const { return logA_; }

    /// The Mukhanov z variable as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& z() const { return z_; }

    /// The first derivative of z with respect to the conformal time as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& zPrime() const { return zPrime_; }

    /// The second derivative of z with respect to the conformal time as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& zDoublePrime() const { return zDoublePrime_; }

    /// The first derivative of a with respect to the conformal time as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& aPrime() const { return aPrime_; }

    /// The second derivative of a with respect to the conformal time as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& aDoublePrime() const { return aDoublePrime_; }

    /// The conformal time as a function of time. Inflation starts at t = 0. t is unitless, assuming the reduced Planck mass is 1.
    const Math::TableFunction<double, double>& tau() const { return tau_; }

    /// Calculate the scalar power spectrum for a given scale.
    /// \param k The scale in inverse Megaparsecs.
    /// \return The scalar power spectrum value.
    double scalarP(double k) const;

    /// Calculate the scalar power spectrum in a given range.
    /// \param kMin The lower end of the range in inverse Megaparsecs.
    /// \param kMax The upper end of the range in inverse Megaparsecs.
    /// \param nPoints The number of points in the range for evaluating the power spectrum. Points will be uniformly separated in log space.
    /// \return A pointer to the TableFunction containing the scalar power spectrum. Must be deleted after use.
    Math::TableFunction<double, double>* scalarPs(double kMin, double kMax, unsigned int nPoints = 100) const;

    /// Calculate the tensor power spectrum for a given scale.
    /// \param k The scale in inverse Megaparsecs.
    /// \return The tensor power spectrum value.
    double tensorP(double k) const;

    /// Calculate the tensor power spectrum in a given range.
    /// \param kMin The lower end of the range in inverse Megaparsecs.
    /// \param kMax The upper end of the range in inverse Megaparsecs.
    /// \param nPoints The number of points in the range for evaluating the power spectrum. Points will be uniformly separated in log space.
    /// \return A pointer to the TableFunction containing the tensor power spectrum. Must be deleted after use.
    Math::TableFunction<double, double>* tensorPs(double kMin, double kMax, unsigned int nPoints = 100) const;

    /// The scale factor at the end of inflation.
    double aEnd() const { return aEnd_; }

    /// The Hubble parameter at the end of inflation.
    double HEnd() const { return HEnd_; }

    /// The current value of the scale factor. This is determined from the starting value, the pivot scale, and the number of e-folds for the pivot scale.
    double aNow() const { return aNow_; }

    /// The field value at the pivot scale.
    double phiPivot() const { return phiPivot_; }

private:
    double kUnitlessFromMpcInv(double k) const;

private:
    Math::TableFunction<double, double> phi_, phiDot_, H_, a_, logA_, z_, zPrime_, zDoublePrime_, tau_, zRatioTau_, aHInverse_, aPrime_, aDoublePrime_, aRatioTau_, logAInv_;
    int K_;
    double tauStart_, tStart_;
    double aHStart_, aHEnd_, HEnd_, aEnd_, aNow_;
    double phiPivot_;
};

#endif

