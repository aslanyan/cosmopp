#ifndef COSMO_PP_MODECODE_HPP
#define COSMO_PP_MODECODE_HPP

#include <vector>

#include <macros.hpp>
#include <table_function.hpp>
#include <cosmological_params.hpp>

class ModeCode
{
private:
    ModeCode() { check(false, "should not be called"); }
    ~ModeCode() { check(false, "should not be called"); }

public:
    static void initialize(int potentialChoice = 1, double kPivot = 0.05, double NPivot = 55, bool instantReheating = true, bool slowRollEnd = true, double kMin = 1e-6, double kMax = 1.0, int nPoints = 100);

    static void setNPivot(double NPivot);

    static bool calculate(const std::vector<double>& vParams);

    inline static const Math::TableFunction<double, double>& getScalarPs() { return scalarPs_; }
    inline static const Math::TableFunction<double, double>& getTensorPs() { return tensorPs_; }

private:
    static double* vParams_;

    static Math::TableFunction<double, double> scalarPs_;
    static Math::TableFunction<double, double> tensorPs_;
};

class ModeCodeCosmologicalParams : public LambdaCDMParams
{
public:
    ModeCodeCosmologicalParams(double omBH2, double omCH2, double h, double tau, double kPivot, double NPivot, int potentialChoice, bool instReheat = true, bool slowRollEnd = true, double kMin = 1e-6, double kMax = 1.0, int nPoints = 500) : LambdaCDMParams(omBH2, omCH2, h, tau, 1.0, 1.0, kPivot)
    {
        ModeCode::initialize(potentialChoice, kPivot, NPivot, instReheat, slowRollEnd, kMin, kMax, nPoints);
    }

    ~ModeCodeCosmologicalParams()
    {
    }

    void setBaseParams(double omBH2, double omCH2, double h, double tau)
    {
        omBH2_ = omBH2;
        omCH2_ = omCH2;
        h_ = h;
        tau_ = tau;
    }

    void setNPivot(double NPivot) { ModeCode::setNPivot(NPivot); }
    bool setVParams(const std::vector<double>& vParams) { return ModeCode::calculate(vParams); }

    virtual const Math::RealFunction& powerSpectrum() const { return ModeCode::getScalarPs(); }
    virtual const Math::RealFunction& powerSpectrumTensor() const { return ModeCode::getTensorPs(); }
};

#endif

