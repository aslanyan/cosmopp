#include <vector>
#include <fstream>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <cosmological_params.hpp>
#include <math_constants.hpp>

#include <class.h>

int main(int argc, char *argv[])
{
    try {
        StandardException exc;
        if(argc < 6)
        {
            std::string exceptionString = "Need to specify the 5 potential params.";
            exc.set(exceptionString);
            throw exc;
        }

        std::stringstream args;
        args << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5];
        std::vector<double> v(5);
        args >> v[0] >> v[1] >> v[2] >> v[3] >> v[4];

        const double kMin = 5e-6, kMax = 1.2;
        const double kPerDecade = 500;
        const double kPivot = 0.002;

        const int lMax = 3000;

        precision pr;
        background br;
        thermo th;
        perturbs pt;
        transfers tr;
        primordial pm;
        spectra sp;
        nonlinear nl;
        lensing le;
        output op;

        if(input_default_precision(&pr) == _FAILURE_)
        {
            std::string exceptionStr = "CLASS: input_default_precision failed!";
            exc.set(exceptionStr);
            throw exc;
        }

        pr.k_per_decade_primordial = kPerDecade;

        if(input_default_params(&br, &th, &pt, &tr, &pm, &sp, &nl, &le, &op) == _FAILURE_)
        {
            std::string exceptionStr = "CLASS: input_default_params failed!";
            exc.set(exceptionStr);
            throw exc;
        }

        pr.evolver = rk;

        const double h = 0.6704;
        const double omBH2 = 0.022032;
        const double omCH2 = 0.12038;
        const double tau = 0.0925;
        const double ns = 0.9619;
        const double as = 2.2154e-9;
        const double pivot = 0.05;
        LambdaCDMParams params(omBH2, omCH2, h, tau, ns, as, pivot);

        br.H0 = params.getH() * 1e5 / _c_;
        br.h = params.getH();
        br.T_cmb = params.getTemperature();
        br.Omega0_b = params.getOmB();
        br.Omega0_cdm = params.getOmC();
        br.Omega0_k = params.getOmK();
        br.Omega0_fld = 0;
        br.Omega0_g = params.getOmG();
        br.Omega0_ur = params.getOmNeutrino();
        br.Omega0_ncdm_tot = 0.0;

        br.Omega0_lambda = 1 + br.Omega0_k - br.Omega0_b - br.Omega0_cdm - br.Omega0_g - br.Omega0_ur - br.Omega0_ncdm_tot;

        if(background_init(&pr, &br) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: background_init failed!" << std::endl << br.error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        th.reio_parametrization = reio_camb;
        th.reio_z_or_tau = reio_tau;
        th.tau_reio = params.getTau();
        th.YHe = _BBN_;
        th.recombination = recfast;

        if(thermodynamics_init(&pr, &br, &th) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: thermodynamics_init failed!" << std::endl << th.error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        pt.has_scalars = true;
        pt.has_vectors = false;
        pt.has_tensors = true;
        pt.has_cls = true;
        pt.l_scalar_max = lMax;
        pt.l_tensor_max = lMax;

        pt.has_perturbations = true;
        pt.has_cl_cmb_temperature = true;
        pt.has_cl_cmb_polarization = true;
        pt.has_cl_cmb_lensing_potential = true;

        if(perturb_init(&pr, &br, &th, &pt) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: perturb_init failed!" << std::endl << pt.error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        pm.primordial_spec_type = inflation_V;
        pm.potential = polynomial;
        pm.V0 = std::pow(10.0, v[4]) / (64 * Math::pi * Math::pi);
        pm.V1 = -pm.V0 * std::sqrt(2.0 * v[0] * 8 * Math::pi);
        pm.V2 = pm.V0 * pm.V0 * v[1] * 8 * Math::pi;
        pm.V3 = pm.V0 * pm.V0 / (pm.V1) * v[2] * (64 * Math::pi * Math::pi);
        pm.V4 = pm.V0 * pm.V0 * pm.V0 / (pm.V1 * pm.V1) * v[3] * (512 * Math::pi * Math::pi * Math::pi);

        pm.k_pivot = kPivot;

        if(primordial_init(&pr, &pt, &pm) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: primordial_init failed!" << std::endl << pm.error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        std::ofstream out("ps_class.txt");
        for(int i = 0; i < pm.lnk_size; ++i)
            out << std::exp(pm.lnk[i]) << ' ' << std::exp(pm.lnpk[0][i]) << ' ' << std::exp(pm.lnpk[1][i]) << std::endl;

        out.close();

        if(primordial_free(&pm) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: primordial_free failed!" << std::endl << pm.error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(perturb_free(&pt) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: perturb_free failed!" << std::endl << pt.error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(thermodynamics_free(&th) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: thermodynamics_free failed!" << std::endl << th.error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(background_free(&br) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: background_free failed!" << std::endl << br.error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        return 0;
    } catch (std::exception& e)
    {
        output_screen("EXCEPTION CAUGHT!!! " << std::endl << e.what() << std::endl);
        output_screen("Terminating!" << std::endl);
        return 1;
    }
    return 0;
}
