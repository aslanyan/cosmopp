#include <cosmo_mpi.hpp>

#include <fstream>
#include <sstream>
#include <memory>
#include <iomanip>
#include <memory>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <planck_like.hpp>
#include <mn_scanner.hpp>
#include <polychord.hpp>
#include <mcmc.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>
#include <modecode.hpp>
#include <taylor_pk.hpp>
#include <ucmh_likelihood.hpp>
#include <progress_meter.hpp>

namespace
{

class TaylorParamsUCMH : public LambdaCDMParams
{
public:
    TaylorParamsUCMH(double omBH2, double omCH2, double h, double tau, double kPivot, double NPivot, int potentialChoice, bool slowRollEnd, bool eternalInflOK, double kMin = 8e-7, double kMax = 1.2, int nPoints = 500, bool useClass = false) : LambdaCDMParams(omBH2, omCH2, h, tau, 1.0, 1.0, kPivot), useClass_(useClass), potentialChoice_(potentialChoice)
    {
        check(!useClass_ || potentialChoice_ == 12, "");

        if(useClass)
        {
            taylor_.reset(new TaylorPk(kPivot, kMin, kMax, 10));
            vParams_.resize(5);
        }
        else
        {
            ModeCode::initialize(potentialChoice, kPivot, NPivot, false, false, slowRollEnd, eternalInflOK, kMin, kMax, nPoints);
            vParams_.resize(ModeCode::getNumVParams());
        }

        output_screen("Using " << vParams_.size() << " modecode parameters." << std::endl);
    }

    ~TaylorParamsUCMH()
    {
    }

    void addKValue(double k, double sMin = 0, double sMax = 1, double tMin = 0, double tMax = 1)
    {
        if(useClass_)
            taylor_->addKValue(k, sMin, sMax, tMin, tMax);
        else
            ModeCode::addKValue(k, sMin, sMax, tMin, tMax);
    }

    void setBaseParams(double omBH2, double omCH2, double h, double tau)
    {
        omBH2_ = omBH2;
        omCH2_ = omCH2;
        h_ = h;
        tau_ = tau;
    }

    bool setVParams(const std::vector<double>& vParams, double *badLike)
    {
        if(useClass_)
            return taylor_->calculate(vParams, badLike);

        return ModeCode::calculate(vParams, badLike);
    }

    virtual const Math::RealFunction& powerSpectrum() const
    {
        if(useClass_) return taylor_->getScalarPs();

        return ModeCode::getScalarPs();
    }

    virtual const Math::RealFunction& powerSpectrumTensor() const
    {
        if(useClass_) return taylor_->getTensorPs();

        return ModeCode::getTensorPs();
    }

    virtual void getAllParameters(std::vector<double>& v) const
    {
        check(useClass_ || vParams_.size() == ModeCode::getNumVParams(), "");
        v.resize(4 + vParams_.size());
        v[0] = omBH2_;
        v[1] = omCH2_;
        v[2] = h_;
        v[3] = tau_;
        for(int i = 0; i < vParams_.size(); ++i)
            v[4 + i] = vParams_[i];
    }

    virtual bool setAllParameters(const std::vector<double>& v, double *badLike)
    {
        check(useClass_ || v.size() == 4 + ModeCode::getNumVParams(), "");

        output_screen1("Param values:");
        for(int i = 0; i < v.size(); ++i)
        {
            output_screen_clean1(std::setprecision(20) << "\t" << v[i]);
        }
        output_screen_clean1(std::endl);

        setBaseParams(v[0], v[1], v[2], v[3]);

        check(useClass_ || vParams_.size() == ModeCode::getNumVParams(), "");
        for(int i = 0; i < vParams_.size(); ++i)
            vParams_[i] = v[4 + i];

        check(vParams_[0] != 0, "");

        //hack
        if(potentialChoice_ == 12)
        {
            // last param is log_10(V0 / eps), need to convert to log_10(V0)
            vParams_[4] += std::log(vParams_[0]) / std::log(10.0);
        }

        if(potentialChoice_ == 14)
        {
            for(int i = 1; i < vParams_.size(); ++i)
                vParams_[i] = std::pow(10.0, vParams_[i]);
        }

        const bool res = setVParams(vParams_, badLike);
        if(!useClass_)
        {
            output_screen1("N_piv = " << ModeCode::getNPivot() << std::endl);
        }
        output_screen1("Result = " << res << std::endl);
        return res;
    }

private:
    const bool useClass_;
    const int potentialChoice_;
    std::unique_ptr<TaylorPk> taylor_;
    std::vector<double> vParams_;
};

class CombinedLikelihood : public Math::LikelihoodFunction
{
public:
    CombinedLikelihood(PlanckLikelihood& planck, CosmologicalParams *params, bool ucmhLim, bool noGamma, bool use200, bool use500, bool useWeak, bool lateDec) : planck_(planck), params_(params)
    {
        if(ucmhLim)
        {
            std::stringstream gammaFileName, pulsarFileName;
            gammaFileName << "data/ucmh_gamma_";
            pulsarFileName << "data/ucmh_pulsar_";

            if(useWeak)
            {
                gammaFileName << "weakened";
                pulsarFileName << "weakened";
            }
            else if(use500)
            {
                gammaFileName << "500";
                pulsarFileName << "500";
            }
            else if(use200)
            {
                gammaFileName << "200";
                pulsarFileName << "200";
            }
            else
            {
                gammaFileName << "1000";
                pulsarFileName << "1000";
            }

            gammaFileName << ".txt";
            pulsarFileName << ".txt";
            if(!noGamma)
                gamma_.reset(new UCMHLikelihood(gammaFileName.str().c_str(), lateDec));

            pulsar_.reset(new UCMHLikelihood(pulsarFileName.str().c_str(), lateDec));
        }
    }

    virtual double calculate(double* params, int nParams)
    {
        double l = planck_.calculate(params, nParams);
        if(l <= 1e8)
        {
            if(gamma_)
            {
                const double gammaLike = gamma_->calculate(params_->powerSpectrum());
                if(gammaLike != 0)
                {
                    output_screen("NONZERO GAMMA LIKE: " << gammaLike << std::endl);
                    l += gammaLike;
                }

            }
            if(pulsar_)
            {
                const double pulsarLike = pulsar_->calculate(params_->powerSpectrum());
                if(pulsarLike != 0)
                {
                    output_screen("NONZERO PULSAR LIKE: " << pulsarLike << std::endl);
                    l += pulsarLike;
                }
            }
        }
        return l;
    }
private:
    PlanckLikelihood& planck_;
    CosmologicalParams *params_;
    std::unique_ptr<UCMHLikelihood> gamma_;
    std::unique_ptr<UCMHLikelihood> pulsar_;
};

} // namespace

int main(int argc, char *argv[])
{
    try {
        bool ucmhLim = false;
        bool useClass = false;
        bool usePoly = false;
        bool noGamma = false;
        bool use200 = false;
        bool use500 = false;
        bool useWeak = false;
        bool lateDecoupling = false;

        bool pbhLimits = false;

        bool knottedPotential = false;

        bool slowRollEnd = false;

        bool highP = false;
        bool foreground = false;

        for(int i = 1; i < argc; ++i)
        {
            if(std::string(argv[i]) == std::string("ucmh"))
                ucmhLim = true;

            if(std::string(argv[i]) == std::string("class"))
                useClass = true;

            if(std::string(argv[i]) == std::string("pc"))
                usePoly = true;

            if(std::string(argv[i]) == std::string("no_gamma"))
                noGamma = true;

            if(std::string(argv[i]) == std::string("ucmh_200"))
                use200 = true;

            if(std::string(argv[i]) == std::string("ucmh_500"))
                use500 = true;

            if(std::string(argv[i]) == std::string("ucmh_weak"))
                useWeak = true;

            if(std::string(argv[i]) == std::string("ucmh_late_dec"))
                lateDecoupling = true;

            if(std::string(argv[i]) == std::string("pbh"))
                pbhLimits = true;

            if(std::string(argv[i]) == std::string("knotted_potential"))
                knottedPotential = true;

            if(std::string(argv[i]) == std::string("slow_roll_end"))
                slowRollEnd = true;

            if(std::string(argv[i]) == std::string("highp"))
                highP = true;

            if(std::string(argv[i]) == std::string("foreground"))
                foreground = true;
        }

        const bool eternalInflOK = !slowRollEnd;

        const bool isMaster = CosmoMPI::create().isMaster();

        if(knottedPotential && useClass)
        {
            if(isMaster)
            {
                output_screen("WARNING: cannot do knotted potential using class. Switching back to ModeCode!" << std::endl);
            }
            useClass = false;
        }

        if(isMaster)
        {
            if(knottedPotential && !usePoly)
            {
                output_screen("WARNING: the knotted potential has a large number of parameters. You're using MultiNest which may take forever to converge. RECOMMENDED to use PolyChord in this case by passing \"pc\" as an argument." << std::endl);
            }
            if(foreground && !usePoly)
            {
                output_screen("WARNING: You're varying the foreground parameters. You're using MultiNest which may take forever to converge. RECOMMENDED to use PolyChord in this case by passing \"pc\" as an argument." << std::endl);
            }
            if(useClass)
            {
                output_screen("Using CLASS for calculating pk." << std::endl);
            }
            else
            {
                output_screen("Using Modecode for calculating pk. To use CLASS instead specify \"class\" as an argument." << std::endl);
            }

            if(usePoly)
            {
                output_screen("Using Polychord sampler." << std::endl);
            }
            else
            {
                output_screen("Using MultiNest sampler. To use Polychord instead specify \"pc\" as an argument." << std::endl);
            }

            if(slowRollEnd)
            {
                output_screen("Requiring slow roll to end!" << std::endl);
            }
            else
            {
                output_screen("NOT requiring slow roll to end!" << std::endl);
            }

            if(ucmhLim)
            {
                output_screen("Using the UCMH limits." << std::endl);
                if(noGamma)
                {
                    output_screen("The gamma-ray ucmh limits will NOT be included." << std::endl);
                }
                else
                {
                    output_screen("The gamma-ray ucmh limits are included. To not include those specify \"no_gamma\" as an argument." << std::endl);
                }

                if(useWeak)
                {
                    output_screen("The weak ucmh limits will be used." << std::endl);
                }
                else if(use500)
                {
                    output_screen("z_c = 500 ucmh limits will be used. To use the weak ones specify \"ucmh_weak\" as an argument instead of \"ucmh_500\"." << std::endl);
                }
                else if(use200)
                {
                    output_screen("z_c = 200 ucmh limits will be used. To use the weak ones specify \"ucmh_weak\" as an argument instead of \"ucmh_200\"." << std::endl);
                }
                else
                {
                    output_screen("z_c = 1000 ucmh limits will be used. To use the z_c = 200 instead specify \"ucmh_200\" as an argument. If you want the weak ucmh limits instead specify \"ucmh_weak\" as an argument." << std::endl);
                }

                if(lateDecoupling)
                {
                    output_screen("Using LATE kinetic decoupling for ucmh." << std::endl);
                }
                else
                {
                    output_screen("Using EARLY kinetic decoupling for ucmh. To use late decoupling instead specify \"ucmh_late_dec\" as an argument." << std::endl);
                }
            }
            else
            {
                output_screen("Not using the UCMH limits. To use those specify \"ucmh\" as an argument." << std::endl);
            }

            if(highP)
            {
                output_screen("Including the high-l polarization likelihood." << std::endl);
            }
            else
            {
                output_screen("The high-l polarization is not included. To include it specify \"highp\" as an argument." << std::endl);
            }

            if(foreground)
            {
                output_screen("Varying the foreground parameters." << std::endl);
            }
            else
            {
                output_screen("Not varying the foreground parameters. To vary them specify \"foreground\" as an argument." << std::endl);
            }
        }

        std::string root;

        if(usePoly)
            root = "slow_test_files/pc_ucmh";
        else
            root = "slow_test_files/mn_ucmh";

        const int nKnots = 20;

        PlanckLikelihood planck(true, true, true, highP, !foreground, false, false, true, 500);
        int nPar = 4;
        if(knottedPotential)
            nPar += (nKnots + 1);
        else
            nPar += 5;

        // foreground
        nPar += 1;

        if(foreground)
        {
            nPar += 15;
            if(highP)
                nPar += 18;
        }

        const double kPivot = 0.05;

        int potentialChoice = 12;
        if(knottedPotential)
            potentialChoice = 14;

        TaylorParamsUCMH modelParams(0.02, 0.1, 0.7, 0.1, kPivot, 55, potentialChoice, slowRollEnd, eternalInflOK, 5e-6, 1.2, 500, useClass);

        if(pbhLimits)
        {
            output_screen("Including PBH limits from the file data/PBH_limits.dat" << std::endl);
            std::ifstream inPBH("data/PBH_limits.dat");
            if(!inPBH)
            {
                StandardException exc;
                std::string exceptionStr = "Cannot read the file data/PBH_limits.dat";
                exc.set(exceptionStr);
                throw exc;
            }
            while(true)
            {
                std::string s;
                std::getline(inPBH, s);
                if(s == "")
                    break;
                if(s[0] == '#')
                    continue;

                std::stringstream str(s);
                double k, lim;
                str >> k >> lim;
                lim = std::pow(10.0, lim);

                if(useClass && k > 1e9)
                    continue;

                //output_screen("PBH limit:\t" << k << '\t' << lim << std::endl);
                modelParams.addKValue(k, 0, lim, 0, 1e10);
            }
        }
        else if(ucmhLim)
        {
            modelParams.addKValue(1e3, 0, 1e10, 0, 1e10);
            modelParams.addKValue(3e3, 0, 1e10, 0, 1e10);
            modelParams.addKValue(1e4, 0, 1e10, 0, 1e10);
            modelParams.addKValue(3e4, 0, 1e10, 0, 1e10);
            modelParams.addKValue(1e5, 0, 1e10, 0, 1e10);
            modelParams.addKValue(3e5, 0, 1e10, 0, 1e10);
            modelParams.addKValue(1e6, 0, 1e10, 0, 1e10);
            modelParams.addKValue(3e6, 0, 1e10, 0, 1e10);
            modelParams.addKValue(1e7, 0, 1e10, 0, 1e10);
            modelParams.addKValue(3e7, 0, 1e10, 0, 1e10);
            modelParams.addKValue(1e8, 0, 1e10, 0, 1e10);
            modelParams.addKValue(3e8, 0, 1e10, 0, 1e10);
            modelParams.addKValue(1e9, 0, 1e10, 0, 1e10);
        }

        planck.setModelCosmoParams(&modelParams);

        CombinedLikelihood like(planck, &modelParams, ucmhLim, noGamma, use200, use500, useWeak, lateDecoupling);

        std::unique_ptr<MnScanner> mn;
        std::unique_ptr<PolyChord> pc;
        if(usePoly)
            pc.reset(new PolyChord(nPar, like, 50, root, 10));
        else
            mn.reset(new MnScanner(nPar, like, (pbhLimits ? 2000 : 500), root));


        int nChains;
        unsigned long burnin;
        unsigned int thin;

        const double vMin = -12, vMax = -4;
        const double sigmaVMin = 0, sigmaVMax = 1;

        if(usePoly)
        {
            int paramIndex = 0;
            pc->setParam(paramIndex++, "ombh2", 0.02, 0.025, 1);
            pc->setParam(paramIndex++, "omch2", 0.1, 0.2, 1);
            pc->setParam(paramIndex++, "h", 0.55, 0.85, 1);
            pc->setParam(paramIndex++, "tau", 0.02, 0.20, 1);
            if(knottedPotential)
            {
                pc->setParam(paramIndex++, "sigma_v", sigmaVMin, sigmaVMax, 2);
                for(int i = 0; i < nKnots; ++i)
                {
                    std::stringstream paramName;
                    paramName << "V_" << i;
                    pc->setParam(paramIndex++, paramName.str().c_str(), vMin, vMax, 2);
                }
            }
            else
            {
                pc->setParam(paramIndex++, "v_1", 0, 0.1, 2);
                pc->setParam(paramIndex++, "v_2", -0.1, 0.1, 2);
                pc->setParam(paramIndex++, "v_3", -0.1, 0.1, 2);
                pc->setParam(paramIndex++, "v_4", -0.1, 0.1, 2);
                pc->setParam(paramIndex++, "v_5", -10, -4, 2);
            }

            pc->setParamGauss(paramIndex++, "A_planck", 1.0, 0.0025, 3);

            if(foreground)
            {
                pc->setParam(paramIndex++, "A_cib_217", 0, 200, 3);
                pc->setParamFixed(paramIndex++, "cib_index", -1.3);
                pc->setParam(paramIndex++, "xi_sz_cib", 0, 1, 3);
                pc->setParam(paramIndex++, "A_sz", 0, 10, 3);
                pc->setParam(paramIndex++, "ps_A_100_100", 0, 400, 3);
                pc->setParam(paramIndex++, "ps_A_143_143", 0, 400, 3);
                pc->setParam(paramIndex++, "ps_A_143_217", 0, 400, 3);
                pc->setParam(paramIndex++, "ps_A_217_217", 0, 400, 3);
                pc->setParam(paramIndex++, "k_sz", 0, 10, 3);
                pc->setParamGauss(paramIndex++, "gal545_A_100", 7, 2, 3);
                pc->setParamGauss(paramIndex++, "gal545_A_143", 9, 2, 3);
                pc->setParamGauss(paramIndex++, "gal545_A_143_217", 21, 8.5, 3);
                pc->setParamGauss(paramIndex++, "gal545_A_217", 80, 20, 3);
                if(highP)
                {
                    pc->setParamGauss(paramIndex++, "galf_EE_A_100", 0.06, 0.012, 3);
                    pc->setParamGauss(paramIndex++, "galf_EE_A_100_143", 0.05, 0.015, 3);
                    pc->setParamGauss(paramIndex++, "galf_EE_A_100_217", 0.11, 0.033, 3);
                    pc->setParamGauss(paramIndex++, "galf_EE_A_143", 0.1, 0.02, 3);
                    pc->setParamGauss(paramIndex++, "galf_EE_A_143_217", 0.24, 0.048, 3);
                    pc->setParamGauss(paramIndex++, "galf_EE_A_217", 0.72, 0.14, 3);
                    pc->setParamFixed(paramIndex++, "galf_EE_index", -2.4);
                    pc->setParamGauss(paramIndex++, "galf_TE_A_100", 0.14, 0.042, 3);
                    pc->setParamGauss(paramIndex++, "galf_TE_A_100_143", 0.12, 0.036, 3);
                    pc->setParamGauss(paramIndex++, "galf_TE_A_100_217", 0.3, 0.09, 3);
                    pc->setParamGauss(paramIndex++, "galf_TE_A_143", 0.24, 0.072, 3);
                    pc->setParamGauss(paramIndex++, "galf_TE_A_143_217", 0.6, 0.18, 3);
                    pc->setParamGauss(paramIndex++, "galf_TE_A_217", 1.8, 0.54, 3);
                    pc->setParamFixed(paramIndex++, "galf_TE_index", -2.4);
                }
                pc->setParamGauss(paramIndex++, "calib_100T", 0.999, 0.001, 3);
                pc->setParamGauss(paramIndex++, "calib_217T", 0.995, 0.002, 3);
                if(highP)
                {
                    pc->setParamFixed(paramIndex++, "calib_100P", 1.0);
                    pc->setParamFixed(paramIndex++, "calib_143P", 1.0);
                    pc->setParamFixed(paramIndex++, "calib_217P", 1.0);
                    pc->setParamFixed(paramIndex++, "A_pol", 1.0);
                }
            }

            check(paramIndex == nPar, "");

            std::vector<double> fracs{0.05, 0.5, 0.45};
            pc->setParameterHierarchy(fracs);

            pc->run(true);

            nChains = 1;
            burnin = 0;
            thin = 1;
        }
        else
        {
            int paramIndex = 0;
            mn->setParam(paramIndex++, "ombh2", 0.02, 0.025);
            mn->setParam(paramIndex++, "omch2", 0.1, 0.2);
            mn->setParam(paramIndex++, "h", 0.55, 0.85);
            mn->setParam(paramIndex++, "tau", 0.02, 0.20);
            if(knottedPotential)
            {
                mn->setParam(paramIndex++, "sigma_v", sigmaVMin, sigmaVMax);
                for(int i = 0; i < nKnots; ++i)
                {
                    std::stringstream paramName;
                    paramName << "V_" << i;
                    mn->setParam(paramIndex++, paramName.str().c_str(), vMin, vMax);
                }
            }
            else
            {
                mn->setParam(paramIndex++, "v_1", 0, 0.1);
                mn->setParam(paramIndex++, "v_2", -0.1, 0.1);
                mn->setParam(paramIndex++, "v_3", -0.1, 0.1);
                mn->setParam(paramIndex++, "v_4", -0.1, 0.1);
                mn->setParam(paramIndex++, "v_5", -10, -4);
            }

            mn->setParamGauss(paramIndex++, "A_planck", 1.0, 0.0025);
            if(foreground)
            {
                mn->setParam(paramIndex++, "A_cib_217", 0, 200);
                mn->setParamFixed(paramIndex++, "cib_index", -1.3);
                mn->setParam(paramIndex++, "xi_sz_cib", 0, 1);
                mn->setParam(paramIndex++, "A_sz", 0, 10);
                mn->setParam(paramIndex++, "ps_A_100_100", 0, 400);
                mn->setParam(paramIndex++, "ps_A_143_143", 0, 400);
                mn->setParam(paramIndex++, "ps_A_143_217", 0, 400);
                mn->setParam(paramIndex++, "ps_A_217_217", 0, 400);
                mn->setParam(paramIndex++, "k_sz", 0, 10);
                mn->setParamGauss(paramIndex++, "gal545_A_100", 7, 2);
                mn->setParamGauss(paramIndex++, "gal545_A_143", 9, 2);
                mn->setParamGauss(paramIndex++, "gal545_A_143_217", 21, 8.5);
                mn->setParamGauss(paramIndex++, "gal545_A_217", 80, 20);
                if(highP)
                {
                    mn->setParamGauss(paramIndex++, "galf_EE_A_100", 0.06, 0.012);
                    mn->setParamGauss(paramIndex++, "galf_EE_A_100_143", 0.05, 0.015);
                    mn->setParamGauss(paramIndex++, "galf_EE_A_100_217", 0.11, 0.033);
                    mn->setParamGauss(paramIndex++, "galf_EE_A_143", 0.1, 0.02);
                    mn->setParamGauss(paramIndex++, "galf_EE_A_143_217", 0.24, 0.048);
                    mn->setParamGauss(paramIndex++, "galf_EE_A_217", 0.72, 0.14);
                    mn->setParamFixed(paramIndex++, "galf_EE_index", -2.4);
                    mn->setParamGauss(paramIndex++, "galf_TE_A_100", 0.14, 0.042);
                    mn->setParamGauss(paramIndex++, "galf_TE_A_100_143", 0.12, 0.036);
                    mn->setParamGauss(paramIndex++, "galf_TE_A_100_217", 0.3, 0.09);
                    mn->setParamGauss(paramIndex++, "galf_TE_A_143", 0.24, 0.072);
                    mn->setParamGauss(paramIndex++, "galf_TE_A_143_217", 0.6, 0.18);
                    mn->setParamGauss(paramIndex++, "galf_TE_A_217", 1.8, 0.54);
                    mn->setParamFixed(paramIndex++, "galf_TE_index", -2.4);
                }
                mn->setParamGauss(paramIndex++, "calib_100T", 0.999, 0.001);
                mn->setParamGauss(paramIndex++, "calib_217T", 0.995, 0.002);
                if(highP)
                {
                    mn->setParamFixed(paramIndex++, "calib_100P", 1.0);
                    mn->setParamFixed(paramIndex++, "calib_143P", 1.0);
                    mn->setParamFixed(paramIndex++, "calib_217P", 1.0);
                    mn->setParamFixed(paramIndex++, "A_pol", 1.0);
                }
            }

            check(paramIndex == nPar, "");

            mn->run(true);

            nChains = 1;
            burnin = 0;
            thin = 1;
        }
        
        if(!CosmoMPI::create().isMaster())
            return 0;

        MarkovChain chain(nChains, root.c_str(), burnin, thin);

        std::vector<MarkovChain::Element*> container;
        chain.getRange(container, 1.0, 0.0);

        std::ofstream outParamLimits("slow_test_files/ucmh_param_limits.txt");
        for(int i = 0; i < nPar; ++i)
        {
            std::string paramName = (usePoly ? pc->getParamName(i) : mn->getParamName(i));

            std::stringstream fileName;
            fileName << "slow_test_files/";
            if(usePoly)
                fileName << "pc_";
            else
                fileName << "mn_";
            fileName << "ucmh_" << paramName << ".txt";
            std::unique_ptr<Posterior1D> p(chain.posterior(i));
            p->writeIntoFile(fileName.str().c_str());

            const double median = p->median();
            double lower, upper;
            p->get1SigmaTwoSided(lower, upper);
            const double sigma = (upper - lower) / 2.0;

            outParamLimits << paramName << " = " << median << "+-" << sigma << std::endl;
        }
        outParamLimits.close();
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
