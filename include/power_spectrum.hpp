#ifndef COSMO_PP_POWER_SPECTRUM_HPP
#define COSMO_PP_POWER_SPECTRUM_HPP

#include <vector>
#include <cmath>

#include <function.hpp>
#include <table_function.hpp>
#include <cubic_spline.hpp>

/// A standard amplitude-tilt-running scalar power spectrum class.
class StandardPowerSpectrum : public Math::RealFunction
{
public:
    /// Constructor.
    /// \param as The scalar amplitude.
    /// \param ns The scalar tilt.
    /// \param pivot The pivot point in Mpc^(-1).
    /// \param run The scalar running, 0 by default.
    StandardPowerSpectrum(double as, double ns, double pivot, double run = 0, double runRun = 0) : as_(as), ns_(ns), pivot_(pivot), run_(run), runRun2_(runRun * runRun) {}

    /// Get the scalar tilt.
    /// \return The scalar tilt.
    double getNs() const { return ns_; }

    /// Set the scalar tilt.
    /// \param ns The scalar tilt.
    void setNs(double ns) { ns_ = ns; }

    /// Get the scalar amplitude.
    /// \return The scalar amplitude.
    double getAs() const { return as_; }

    /// Set the scalar amplitude.
    /// \param as The scalar amplitude.
    void setAs(double as) { as_ = as; }

    /// Get the pivot point.
    /// \return The pivot point in Mpc^(-1).
    double getPivot() const { return pivot_; }

    /// Get the scalar running.
    /// \return The scalar running.
    double getRun() const { return run_; }

    /// Set the scalar running.
    /// \param run The scalar running.
    void setRun(double run) { run_ = run; }

    /// Get the running of the running.
    /// \return The running of the running.
    double getRunRun() const { return std::sqrt(runRun2_); }
    
    /// Set the running of the running.
    /// \param runRun The running of the running.
    void setRunRun(double runRun) { runRun2_ = runRun * runRun; }

    /// Calculate the scalar power spectrum at a given point.
    /// \param k The k value in Mpc^(-1).
    /// \return The scalar power spectrum value.
    virtual double evaluate(double k) const
    {
        const double lkPiv = std::log(k / pivot_);
        const double p = ns_ - 1.0 + 0.5 * run_ * lkPiv + 1.0 / 6.0 * runRun2_ * lkPiv * lkPiv;
        return as_ * std::pow(k / pivot_, p);
    }

private:
    double ns_;
    double as_;
    double pivot_;
    double run_;
    double runRun2_;
};

class CutoffPowerSpectrum : public StandardPowerSpectrum
{
public:
    CutoffPowerSpectrum(double kCut, double as, double ns, double pivot, double run = 0) : StandardPowerSpectrum(as, ns, pivot, run), kCut_(kCut) { check(kCut >= 0, "invalid kCut"); }

    double getKCut() const { return kCut_; }

    virtual double evaluate(double k) const
    {
        if(k < kCut_)
            return 1e-100;

        return StandardPowerSpectrum::evaluate(k);
    }

private:
    double kCut_;
};

/// A standard amplitude-tilt tensor power spectrum class.
class StandardPowerSpectrumTensor : public Math::RealFunction
{
public:
    /// Constructor.
    /// \param at The tensor amplitude (equal to r times the scalar amplitude).
    /// \param nt The tensor tilt.
    /// \param pivot The pivot point in Mpc^(-1).
    /// \param run The tensor running, 0 by default.
    StandardPowerSpectrumTensor(double at, double nt, double pivot, double run = 0) : at_(at), nt_(nt), pivot_(pivot), run_(run) {}

    /// Constructor.
    /// \param scalarPs A constant reference to the scalar power spectrum. This is used to set the amplitude from the tensor-to-scalar ratio.
    /// \param r The tensor-to-scalar ratio.
    /// \param nt The tensor tilt.
    /// \param pivot The pivot point in Mpc^(-1).
    /// \param run The tensor running, 0 by default.
    StandardPowerSpectrumTensor(const Math::RealFunction& scalarPs, double r, double nt, double pivot, double run = 0)
    {
        set(scalarPs, r, nt, pivot, run);
    }

    /// Get the tensor tilt.
    /// \return The tensor tilt.
    double getNt() const { return nt_; }

    /// Get the tensor amplitude.
    /// \return The tensor amplitude.
    double getAt() const { return at_; }

    /// Get the pivot point.
    /// \return The pivot point in Mpc^(-1).
    double getPivot() const { return pivot_; }

    /// Set new parameters.
    /// \param at The tensor amplitude (equal to r times the scalar amplitude).
    /// \param nt The tensor tilt.
    /// \param pivot The pivot point in Mpc^(-1).
    /// \param run The tensor running, 0 by default.
    void set(double at, double nt, double pivot, double run = 0)
    {
        at_ = at;
        nt_ = nt;
        pivot_ = pivot;
        run_ = run;
    }

    /// Set new parameters.
    /// \param scalarPs A constant reference to the scalar power spectrum. This is used to set the amplitude from the tensor-to-scalar ratio.
    /// \param r The tensor-to-scalar ratio.
    /// \param nt The tensor tilt.
    /// \param pivot The pivot point in Mpc^(-1).
    /// \param run The tensor running, 0 by default.
    void set(const Math::RealFunction& scalarPs, double r, double nt, double pivot, double run = 0)
    {
        at_ = 1.0;
        nt_ = nt;
        pivot_ = pivot;
        run_ = run;

        // fix at
        const double s = scalarPs.evaluate(pivot);
        const double t = evaluate(pivot);
        at_ = r * s / t;
    }

    /// Calculate the tensor power spectrum at a given point.
    /// \param k The k value in Mpc^(-1).
    /// \return The tensor power spectrum value.
    virtual double evaluate(double k) const
    {
        const double p = nt_ + 0.5 * run_ * std::log(k / pivot_);
        return at_ * std::pow(k / pivot_, p);
    }

private:
    double nt_;
    double at_;
    double pivot_;
    double run_;
};

/// A scalar power spectrum function with binning with a linear spline.
class LinearSplinePowerSpectrum : public Math::RealFunction
{
public:
    /// Constructor.
    /// \param kVals A vector containig the points in k space in Mpc^(-1). At least two points are needed (the endpoints), and the calculation should always be done between these endpoints.
    /// \param amplitudes A vector containing the amplitudes corresponding to the k values. The size of the vector needs to be the same as kVals.
    LinearSplinePowerSpectrum(const std::vector<double>& kVals, const std::vector<double>& amplitudes)
    {
        check(kVals.size() >= 2, "at least two points needed");
        check(kVals.size() == amplitudes.size(), "");

        for(int i = 0; i < kVals.size(); ++i)
            tf_[std::log(kVals[i])] = std::log(amplitudes[i]);

        ns_ = 1; // just whatever
        as_ = 2e-9; // just whatever
        pivot_ = 0.05; // just whatever
    }

    /// Calculate the scalar power spectrum at a given point.
    /// \param k The k value in Mpc^(-1).
    /// \return The scalar power spectrum value.
    virtual double evaluate(double k) const
    {
        return std::exp(tf_.evaluate(std::log(k)));
    }

    /// This function should not be used for this class, just written for compatibility.
    double getNs() const { return ns_; }

    /// This function should not be used for this class, just written for compatibility.
    double getAs() const { return as_; }

    /// This function should not be used for this class, just written for compatibility.
    double getPivot() const { return pivot_; }

private:
    Math::TableFunction<double, double> tf_;

    double ns_;
    double as_;
    double pivot_;
};

/// A scalar power spectrum function with binning with a cubic spline.
class CubicSplinePowerSpectrum : public Math::RealFunction
{
public:
    /// Constructor.
    /// \param kVals A vector containig the points in k space in Mpc^(-1). At least two points are needed (the endpoints), and the calculation should always be done between these endpoints.
    /// \param amplitudes A vector containing the amplitudes corresponding to the k values. The size of the vector needs to be the same as kVals.
    CubicSplinePowerSpectrum(const std::vector<double>& kVals, const std::vector<double>& amplitudes)
    {
        check(kVals.size() >= 2, "at least two points needed");
        check(kVals.size() == amplitudes.size(), "");

        std::vector<double> logK, logA;
        for(int i = 0; i < kVals.size(); ++i)
        {
            logK.push_back(std::log(kVals[i]));
            logA.push_back(std::log(amplitudes[i]));
        }

        cs_ = new Math::CubicSpline(logK, logA);

        ns_ = 1; // just whatever
        as_ = 2e-9; // just whatever
        pivot_ = 0.05; // just whatever
    }

    /// Destructor.
    ~CubicSplinePowerSpectrum()
    {
        delete cs_;
    }

    /// Calculate the scalar power spectrum at a given point.
    /// \param k The k value in Mpc^(-1).
    /// \return The scalar power spectrum value.
    virtual double evaluate(double k) const
    {
        return std::exp(cs_->evaluate(std::log(k)));
    }

    /// This function should not be used for this class, just written for compatibility.
    double getNs() const { return ns_; }

    /// This function should not be used for this class, just written for compatibility.
    double getAs() const { return as_; }

    /// This function should not be used for this class, just written for compatibility.
    double getPivot() const { return pivot_; }

private:
    Math::CubicSpline* cs_;

    double ns_;
    double as_;
    double pivot_;
};

#endif

