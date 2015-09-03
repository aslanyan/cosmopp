#include <memory>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <polychord.hpp>
#include <planck_like.hpp>
#include <power_spectrum.hpp>
#include <cosmological_params.hpp>

namespace
{

class SplineWithDegenerateNeutrinosParams : public LambdaCDMParams
{
public:
    // the ns, as, and pivot values for the base class should not be used and are set to arbitrary values
    SplineWithDegenerateNeutrinosParams(bool isLinear, double omBH2, double omCH2, double h, double tau, const std::vector<double>& kVals, const std::vector<double>& amplitudes,  double nEff, int nMassive, double sumMNu, bool varyNEff, bool varySumMNu) : LambdaCDMParams(omBH2, omCH2, h, tau, 1.0, 2e-9, 0.05), isLinear_(isLinear), nEff_(nEff), nMassive_(nMassive), sumMNu_(sumMNu), kVals_(kVals), amplitudes_(amplitudes), lcdmParams_(6), varyNEff_(varyNEff), varySumMNu_(varySumMNu)
    {
        check(nEff > 0, "invalid nEff = " << nEff);
        check(sumMNu >= 0, "invalid sumMNu = " << sumMNu);
        check(nMassive >= 0, "number of massive neutrinos is negative: " << nMassive);
        check(nEff > nMassive, "nEff needs to be more than the number of massive neutrinos");

        check(kVals_.size() >= 2, "");
        check(kVals_.size() == amplitudes.size(), "");

        check(!varyNEff_ || nMassive_ > 0, "");

        resetPS();
    }

    ~SplineWithDegenerateNeutrinosParams() {}

    virtual double getNEff() const { return nEff_ - nMassive_; }
    virtual int getNumNCDM() const { return nMassive_; }
    virtual double getNCDMParticleMass(int i) const
    {
        check(i >= 0 && i < nMassive_, "invalid index = " << i);
        return sumMNu_ / nMassive_;
    }

    virtual double getNCDMParticleTemp(int i) const
    {
        check(i >= 0 && i < nMassive_, "invalid index = " << i);
        //return 0.715985;
        return 0.713765855506013;
    }

    virtual const Math::RealFunction& powerSpectrum() const { return *ps_; }

    virtual void getAllParameters(std::vector<double>& v) const
    {
        check(kVals_.size() >= 2, "");
        check(kVals_.size() == amplitudes_.size(), "");

        int nPar = 4 + 2 * kVals_.size() - 2;

        if(varyNEff_)
            ++nPar;
        if(varySumMNu_)
            ++nPar;

        v.resize(nPar);

        std::vector<double>::iterator it = v.begin();
        *(it++) = getOmBH2();
        *(it++) = getOmCH2();
        *(it++) = getH();
        *(it++) = getTau();

        if(varyNEff_)
            *(it++) = nEff_;

        if(varySumMNu_)
            *(it++) = sumMNu_;

        for(int i = 1; i < kVals_.size() - 1; ++i)
        {
            check(it < v.end(), "");
            *(it++) = std::log(kVals_[i]);
        }

        for(int i = 0; i < amplitudes_.size(); ++i)
        {
            check(it < v.end(), "");
            *(it++) = std::log(amplitudes_[i] * 1e10);
        }

        check(it == v.end(), "");
    }

    virtual bool setAllParameters(const std::vector<double>& v, double *badLike = NULL)
    {
        int nPar = 4 + 2 * kVals_.size() - 2;

        if(varyNEff_)
            ++nPar;
        if(varySumMNu_)
            ++nPar;

        check(v.size() == nPar, "");
        check(kVals_.size() >= 2, "");
        check(kVals_.size() == amplitudes_.size(), "");

        check(lcdmParams_.size() == 6, "");

        std::vector<double>::const_iterator it = v.begin();

        for(int i = 0; i < 4; ++i)
            lcdmParams_[i] = *(it++);

        lcdmParams_[4] = 1.0; // arbitrary ns, doesn't matter
        lcdmParams_[5] = 2e-9; // arbitrary as, doesn't matter

        if(!LambdaCDMParams::setAllParameters(lcdmParams_, badLike))
            return false;

        if(varyNEff_)
            nEff_ = *(it++);

        if(varySumMNu_)
            sumMNu_ = *(it++);

        for(int i = 1; i < kVals_.size() - 1; ++i)
        {
            check(it < v.end(), "");
            kVals_[i] = std::exp(*(it++));
        }

        for(int i = 0; i < amplitudes_.size(); ++i)
        {
            check(it < v.end(), "");
            amplitudes_[i] = std::exp(*(it++)) / 1e10;
        }

        resetPS();

        return true;
    }

private:
    void resetPS()
    {
        if(isLinear_)
            ps_.reset(new LinearSplinePowerSpectrum(kVals_, amplitudes_));
        else
            ps_.reset(new CubicSplinePowerSpectrum(kVals_, amplitudes_));
    }

private:
    const bool isLinear_;

    double nEff_;
    int nMassive_;
    double sumMNu_;

    const bool varyNEff_;
    const bool varySumMNu_;

    std::unique_ptr<Math::RealFunction> ps_;
    std::vector<double> kVals_, amplitudes_;
    std::vector<double> lcdmParams_;
};

} // namespace

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 3)
        {
            std::string exceptionString = "Need to specify the spline type (linear or cubic) and the number of (internal) knots.";
            exc.set(exceptionString);
            throw exc;
        }

        bool isLinear = true;
        if(argv[1][0] == 'c')
            isLinear = false;

        const int nKnots = std::atoi(argv[2]);

        if(nKnots < 0 || nKnots > 100)
        {
            std::stringstream exceptionStr;
            exceptionStr << "The number of knots " << nKnots << " is invalid. Needs to be positive and not larger than 100.";
            exc.set(exceptionStr.str());
            throw exc;
        }

        bool varyNEff = false;
        bool varySumMNu = false;

        for(int i = 3; i < argc; ++i)
        {
            if(std::string(argv[i]) == "neff")
                varyNEff = true;
            if(std::string(argv[i]) == "sum_mnu")
                varySumMNu = true;
        }

        int nPar = 4 + 2 * (nKnots + 2) - 2;

        if(varyNEff)
        {
            ++nPar;
            output_screen("Varying N_eff!" << std::endl);
        }
        else
        {
            output_screen("N_eff not being varied. To vary it give \"neff\" as an extra argument." << std::endl);
        }
        
        if(varySumMNu)
        {
            ++nPar;
            output_screen("Varying sum_mnu!" << std::endl);
        }
        else
        {
            output_screen("sum_mnu not being varied. To vary it give \"sum_mnu\" as an extra argument." << std::endl);
        }

        // for A_planck
        ++nPar;

        const double kMin = 0.8e-6;
        const double kMax = 1.2;
        const double aMin = -2;
        const double aMax = 4;

        std::vector<double> kVals(nKnots + 2);
        std::vector<double> amplitudes(nKnots + 2, 2e-9);

        kVals[0] = kMin;
        kVals.back() = kMax;

        const double deltaLogK = (std::log(kMax) - std::log(kMin)) / (nKnots + 1);

        for(int i = 1; i < kVals.size() - 1; ++i)
            kVals[i] = std::exp(std::log(kMin) + i * deltaLogK);

        SplineWithDegenerateNeutrinosParams params(isLinear, 0.022, 0.12, 0.7, 0.1, kVals, amplitudes, 3.046, 0, 0, varyNEff, varySumMNu);

#ifdef COSMO_PLANCK_15
        PlanckLikelihood planckLike(true, true, true, false, true, false, false, false, 100);
#else
        ERROR NOT IMPLEMENTED;
#endif
        planckLike.setModelCosmoParams(&params);

        std::stringstream root;
        root << "knotted_";
        if(isLinear)
            root << "linear_";
        else
            root << "cubic_";
        root << nKnots;
        if(varyNEff)
            root << "_neff";
        if(varySumMNu)
            root << "_summnu";
        PolyChord pc(nPar, planckLike, 300, root.str(), 3 * (4 + (varyNEff ? 1 : 0) + (varySumMNu ? 1 : 0)));

        int paramIndex = 0;

        pc.setParam(paramIndex++, "ombh2", 0.02, 0.025, 1);
        pc.setParam(paramIndex++, "omch2", 0.1, 0.2, 1);
        pc.setParam(paramIndex++, "h", 0.55, 0.85, 1);
        pc.setParam(paramIndex++, "tau", 0.02, 0.20, 1);
        if(varyNEff)
            pc.setParam(paramIndex++, "n_eff", 2.0, 5.0, 1);
        if(varySumMNu)
            pc.setParam(paramIndex++, "sum_mnu", 0.001, 3.0, 1);

        for(int i = 1; i < kVals.size() - 1; ++i)
        {
            std::stringstream paramName;
            paramName << "k_" << i;
            pc.setParamSortedUniform(paramIndex++, paramName.str(), std::log(kMin), std::log(kMax), 0, 2);
        }

        for(int i = 0; i < amplitudes.size(); ++i)
        {
            std::stringstream paramName;
            paramName << "a_" << i;
            pc.setParam(paramIndex++, paramName.str(), aMin, aMax, 2);
        }
#ifdef COSMO_PLANCK_15
        pc.setParamGauss(paramIndex++, "A_planck", 1, 0.0025, 3);
#else
        ERROR NOT IMPLEMENTED;
#endif
        check(paramIndex == nPar, "");

        const std::vector<double> fracs{0.7, 0.25, 0.05};
        pc.setParameterHierarchy(fracs);

        pc.run(true);
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
