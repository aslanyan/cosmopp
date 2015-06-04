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
    static void initialize(int potentialChoice = 1, double kPivot = 0.05, double NPivot = 55, bool instantReheating = true, bool physicalPriors = true, bool slowRollEnd = true, bool eternalInflOK = true, double kMin = 1e-6, double kMax = 1.0, int nPoints = 100);

    static void setNPivot(double NPivot);
    static double getNPivot();

    static void addKValue(double k, double sMin = 0, double sMax = 1, double tMin = 0, double tMax = 1);

    static bool calculate(const std::vector<double>& vParams);

    static int getNumVParams() { return nVPar_; }

    inline static const Math::TableFunction<double, double>& getScalarPs() { return scalarPs_; }
    inline static const Math::TableFunction<double, double>& getTensorPs() { return tensorPs_; }

private:
    static int nVPar_;
    static double* vParams_;

    static Math::TableFunction<double, double> scalarPs_;
    static Math::TableFunction<double, double> tensorPs_;

    static std::map<double, double> scalarLower_, scalarUpper_, tensorLower_, tensorUpper_;
};

class ModeCodeCosmologicalParams : public LambdaCDMParams
{
public:
    ModeCodeCosmologicalParams(double omBH2, double omCH2, double h, double tau, double kPivot, double NPivot, int potentialChoice, bool instReheat = true, bool physicalPriors = true, bool slowRollEnd = true, bool eternalInflOK = true, double kMin = 8e-7, double kMax = 1.2, int nPoints = 500) : LambdaCDMParams(omBH2, omCH2, h, tau, 1.0, 1.0, kPivot)
    {
        ModeCode::initialize(potentialChoice, kPivot, NPivot, instReheat, physicalPriors, slowRollEnd, eternalInflOK, kMin, kMax, nPoints);
        vParams_.resize(ModeCode::getNumVParams());
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

    virtual void getAllParameters(std::vector<double>& v) const
    {
        check(vParams_.size() == ModeCode::getNumVParams(), "");
        v.resize(5 + vParams_.size());
        v[0] = omBH2_;
        v[1] = omCH2_;
        v[2] = h_;
        v[3] = tau_;
        v[4] = ModeCode::getNPivot();
        for(int i = 0; i < vParams_.size(); ++i)
            v[5 + i] = vParams_[i];
    }

    virtual bool setAllParameters(const std::vector<double>& v)
    {
        check(v.size() == 5 + ModeCode::getNumVParams(), "");

        //output_screen_clean("Param values:");
        //for(int i = 0; i < v.size(); ++i)
            //output_screen_clean("\t" << v[i]);
        //output_screen_clean(std::endl);

        setBaseParams(v[0], v[1], v[2], v[3]);
        setNPivot(v[4]);

        check(vParams_.size() == ModeCode::getNumVParams(), "");
        for(int i = 0; i < vParams_.size(); ++i)
            vParams_[i] = v[5 + i];

        return setVParams(vParams_);
    }

private:
    std::vector<double> vParams_;
};

#endif

