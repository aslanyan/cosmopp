#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <cosmological_params.hpp>
#include <math_constants.hpp>
#include <taylor_pk.hpp>

#include <class.h>

TaylorPk::TaylorPk(double kPivot, double kMin, double kMax, int kPerDecade)
{
    check(kPivot > 0, "invalid k_pivot = " << kPivot);
    check(kMin > 0, "invalid k_min = " << kMin);
    check(kMax >= kMin, "invalid k_max = " << kMax);
    check(kPerDecade > 0, "invalid kPerDecade = " << kPerDecade << ", needs to be positive");

    const int nDecades = (int)std::ceil((std::log(kMax) - std::log(kMin)) / std::log(10.0));
    check(nDecades > 0, "");
    const int nPoints = kPerDecade * nDecades;
    const double deltaLogK = (std::log(kMax) - std::log(kMin)) / nPoints;

    for(int i = 0; i <= nPoints; ++i)
    {
        const double k = (i == nPoints ? kMax : std::exp(std::log(kMin) + i * deltaLogK));
        scalarLower_[k] = 0;
        scalarUpper_[k] = 1e-6;
        tensorLower_[k] = 0;
        tensorUpper_[k] = 1e-6;
    }
    
    pr_ = new precision;
    br_ = new background;
    th_ = new thermo;
    pt_ = new perturbs;
    tr_ = new transfers;
    pm_ = new primordial;
    sp_ = new spectra;
    nl_ = new nonlinear;
    le_ = new lensing;
    op_ = new output;

    StandardException exc;
    if(input_default_precision(pr_) == _FAILURE_)
    {
        std::string exceptionStr = "CLASS: input_default_precision failed!";
        exc.set(exceptionStr);
        throw exc;
    }

    pr_->k_per_decade_primordial = kPerDecade;

    if(input_default_params(br_, th_, pt_, tr_, pm_, sp_, nl_, le_, op_) == _FAILURE_)
    {
        std::string exceptionStr = "CLASS: input_default_params failed!";
        exc.set(exceptionStr);
        throw exc;
    }

    pr_->evolver = rk;

    const double h = 0.6704;
    const double omBH2 = 0.022032;
    const double omCH2 = 0.12038;
    const double tau = 0.0925;
    const double ns = 0.9619;
    const double as = 2.2154e-9;
    const double pivot = 0.05;
    LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

    br_->H0 = params.getH() * 1e5 / _c_;
    br_->h = params.getH();
    br_->T_cmb = params.getTemperature();
    br_->Omega0_b = params.getOmB();
    br_->Omega0_cdm = params.getOmC();
    br_->Omega0_k = params.getOmK();
    br_->Omega0_fld = 0;
    br_->Omega0_g = params.getOmG();
    br_->Omega0_ur = params.getOmNeutrino();
    br_->Omega0_ncdm_tot = 0.0;

    br_->Omega0_lambda = 1 + br_->Omega0_k - br_->Omega0_b - br_->Omega0_cdm - br_->Omega0_g - br_->Omega0_ur - br_->Omega0_ncdm_tot;

    if(background_init(pr_, br_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: background_init failed!" << std::endl << br_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    th_->reio_parametrization = reio_camb;
    th_->reio_z_or_tau = reio_tau;
    th_->tau_reio = params.getTau();
    th_->YHe = _BBN_;
    th_->recombination = recfast;

    if(thermodynamics_init(pr_, br_, th_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: thermodynamics_init failed!" << std::endl << th_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    const int lMax = 3000;

    pt_->has_scalars = true;
    pt_->has_vectors = false;
    pt_->has_tensors = true;
    pt_->has_cls = true;
    pt_->l_scalar_max = lMax;
    pt_->l_tensor_max = lMax;

    pt_->has_perturbations = true;
    pt_->has_cl_cmb_temperature = true;
    pt_->has_cl_cmb_polarization = true;
    pt_->has_cl_cmb_lensing_potential = true;

    if(perturb_init(pr_, br_, th_, pt_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: perturb_init failed!" << std::endl << pt_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    check(kMax <= 1e9, "");

    pt_->k_min = kMin;
    pt_->k_max = kMax;

    if(kMax > 1e7)
        pr_->primordial_inflation_bg_stepsize *= 3;
    if(kMax > 1e8)
        pr_->primordial_inflation_bg_stepsize *= 10;

    pm_->primordial_spec_type = inflation_V;
    pm_->potential = polynomial;
    pm_->k_pivot = kPivot;
}

TaylorPk::~TaylorPk()
{
    StandardException exc;

    if(perturb_free(pt_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: perturb_free failed!" << std::endl << pt_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    if(thermodynamics_free(th_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: thermodynamics_free failed!" << std::endl << th_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    if(background_free(br_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: background_free failed!" << std::endl << br_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }
    delete pr_;
    delete br_;
    delete th_;
    delete pt_;
    delete tr_;
    delete pm_;
    delete sp_;
    delete nl_;
    delete le_;
    delete op_;
}

void
TaylorPk::addKValue(double k, double sMin, double sMax, double tMin, double tMax)
{
    check(k > 0, "");
    check(sMax >= sMin, "");
    check(sMin >= 0, "");
    check(sMax > 0, "");
    check(tMax >= tMin, "");
    check(tMin >= 0, "");
    check(tMax > 0, "");
    scalarLower_[k] = sMin;
    scalarUpper_[k] = sMax;
    tensorLower_[k] = tMin;
    tensorUpper_[k] = tMax;

    check(k >= pt_->k_min, "");

    if(k > pt_->k_max)
    {
        check(k <= 1e9, "");

        if(k > 1e7 && pt_->k_max <= 1e7)
            pr_->primordial_inflation_bg_stepsize *= 3;
        if(k > 1e8 && pt_->k_max <= 1e8)
            pr_->primordial_inflation_bg_stepsize *= 10;

        pt_->k_max = k;
    }
}


bool
TaylorPk::calculate(const std::vector<double>& v, double *badLike)
{
    check(v.size() == 5, "");

    double bad1 = 0;
    for(int i = 0; i < v.size(); ++i)
    {
        bad1 += v[i] * v[i];
    }
    
    if(bad1 > 0)
        while(bad1 < 1) bad1 *= 10;

    check(bad1 >= 0, "");

    pm_->V0 = std::pow(10.0, v[4]) / (64 * Math::pi * Math::pi);
    pm_->V1 = -pm_->V0 * std::sqrt(2.0 * v[0] * 8 * Math::pi);
    pm_->V2 = pm_->V0 * v[1] * 8 * Math::pi;
    pm_->V3 = pm_->V0 * pm_->V0 / (pm_->V1) * v[2] * (64 * Math::pi * Math::pi);
    pm_->V4 = pm_->V0 * pm_->V0 * pm_->V0 / (pm_->V1 * pm_->V1) * v[3] * (512 * Math::pi * Math::pi * Math::pi);


    StandardException exc;
    if(primordial_init(pr_, pt_, pm_) == _FAILURE_)
    {
        output_screen1("Taylor pk failure: " << pm_->error_message << std::endl);
        /*
        if(primordial_free(pm_) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: primordial_free failed!" << std::endl << pm_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }
        */

        if(badLike)
            *badLike = 1e10 * bad1;

        return false;
    }

    scalarPs_.clear();
    tensorPs_.clear();
    for(int i = 0; i < pm_->lnk_size; ++i)
    {
        const double k = std::exp(pm_->lnk[i]);
        const double s = std::exp(pm_->lnpk[0][i]);
        const double t = std::exp(pm_->lnpk[1][i]);

        scalarPs_[k] = s;
        tensorPs_[k] = t;
    }

    if(primordial_free(pm_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: primordial_free failed!" << std::endl << pm_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    double bad = 0;

    for(Math::TableFunction<double, double>::const_iterator it = scalarLower_.begin(); it != scalarLower_.end(); ++it)
    {
        const double k = it->first;
        check(scalarUpper_.find(k) != scalarUpper_.end(), "");
        check(tensorLower_.find(k) != tensorLower_.end(), "");
        check(tensorUpper_.find(k) != tensorUpper_.end(), "");
        const double sMin = scalarLower_[k], sMax = scalarUpper_[k], tMin = tensorLower_[k], tMax = tensorUpper_[k];

        const double s = scalarPs_.evaluate(k);
        const double t = tensorPs_.evaluate(k);

        if(s < sMin)
            bad += (sMin - s) / sMin;
        if(s > sMax)
            bad += (s - sMax) / sMax;
        if(t < tMin)
            bad += (tMin - t) / tMin; 
        if(t > tMax)
            bad += (t - tMax) / tMax;
    }

    check(bad >= 0, "");

    if(badLike)
        *badLike = bad * 1e10;

    return (bad == 0);
}

