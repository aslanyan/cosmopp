#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <cmb.hpp>

#include <class.h>

void CMB::allocate()
{
    pr_ = new precision;
    br_ = new background;
    th_ = new thermo;
    pt_ = new perturbs;
    //bs_ = new bessels;
    tr_ = new transfers;
    pm_ = new primordial;
    sp_ = new spectra;
    nl_ = new nonlinear;
    le_ = new lensing;
    op_ = new output;

    for(int i = 0; i < 1000; ++i)
    {
        dum1_[i] = new double[1000];
        dum2_[i] = new double[1000];
    }
}

void
CMB::deAllocate()
{
    delete pr_;
    delete br_;
    delete th_;
    delete pt_;
    //delete bs_;
    delete tr_;
    delete pm_;
    delete sp_;
    delete nl_;
    delete le_;
    delete op_;

    for(int i = 0; i < 1000; ++i)
    {
        delete dum1_[i];
        delete dum2_[i];
    }
}

void
CMB::preInitialize(int lMax, bool wantAllL, bool primordialInitialize, bool includeTensors, int lMaxTensors, double kPerDecade, double kMin, double kMax)
{
    StandardException exc;
    check(lMax >= 2, "invalid lMax = " << lMax);
    lMax_ = lMax;

    if(preInit_)
        preClean();

    if(input_default_precision(pr_) == _FAILURE_)
    {
        std::string exceptionStr = "CLASS: input_default_precision failed!";
        exc.set(exceptionStr);
        throw exc;
    }

    primordialInitialize_ = primordialInitialize;

    if(primordialInitialize)
    {
        check(kPerDecade >= 1, "");
        check(kMin > 0, "");
        check(kMax > kMin, "");

        pr_->k_per_decade_primordial = kPerDecade;
        kPerDecade_ = kPerDecade;
        kMin_ = kMin;
        kMax_ = kMax;
    }

    if(wantAllL)
    {
        pr_->l_logstep = std::log(double(lMax + 1)) / std::log(double(lMax));
        //pr_->l_linstep = 1;
    }

    if(input_default_params(br_, th_, pt_, tr_, pm_, sp_, nl_, le_, op_) == _FAILURE_)
    {
        std::string exceptionStr = "CLASS: input_default_params failed!";
        exc.set(exceptionStr);
        throw exc;
    }

    /*
    bs_->x_max = lMax * pr_->k_scalar_max_tau0_over_l_max;
    bs_->x_step = pr_->bessel_x_step;
    bs_->x_max = ((int)(bs_->x_max * 1.01 / bs_->x_step) + 1) * bs_->x_step;

    bs_->l_max = lMax;

    if(bessel_init(pr_, bs_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: bessel_init failed!" << std::endl << bs_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }
    */

    includeTensors_ = includeTensors;
    if(includeTensors_)
    {
        check(lMaxTensors >= 2, "invalid lMaxTensors = " << lMaxTensors);
        lMaxTensors_ = lMaxTensors;
    }

    preInit_ = true;
}

void
CMB::initialize(const CosmologicalParams& params, bool wantT, bool wantPol, bool wantLensing, bool wantMatterPs, double zMaxPk)
{
    StandardException exc;

    check(preInit_, "need to preInitialize first");

    params_ = &params;

    if(init_)
        clean();

    pr_->k_per_decade_for_pk = 100;
    pr_->k_per_decade_for_bao = 100;

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

    const double ncdm = params.getNumNCDM();
    check(ncdm >= 0, "invalid ncdm = " << ncdm);
    if(ncdm != 0)
    {
        br_->N_ncdm = ncdm;

        br_->T_ncdm = (double*) malloc(ncdm * sizeof(double));
        br_->ksi_ncdm = (double*) malloc(ncdm * sizeof(double));
        br_->deg_ncdm = (double*) malloc(ncdm * sizeof(double));
        br_->M_ncdm = (double*) malloc(ncdm * sizeof(double));
        br_->Omega0_ncdm = (double*) malloc(ncdm * sizeof(double));
        br_->m_ncdm_in_eV = (double*) malloc(ncdm * sizeof(double));

        br_->got_files = (int*) malloc(ncdm * sizeof(int));

        br_->ncdm_psd_files = NULL;
        br_->ncdm_psd_parameters = NULL;

        for(int i = 0; i < ncdm; ++i)
        {
            br_->T_ncdm[i] = params.getNCDMParticleTemp(i);
            br_->ksi_ncdm[i] = br_->ksi_ncdm_default;
            br_->deg_ncdm[i] = br_->deg_ncdm_default;
            br_->M_ncdm[i] = 0.0;
            br_->Omega0_ncdm[i] = 0.0;
            br_->m_ncdm_in_eV[i] = params.getNCDMParticleMass(i);

            br_->got_files[i] = 0;
        }

        if(background_ncdm_init(pr_, br_) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: background_ncdm_init failed!" << std::endl << br_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        for(int i = 0; i < ncdm; ++i)
        {
            double rhoNCDM;
            br_->M_ncdm[i] = br_->m_ncdm_in_eV[i] * Phys::eCharge / Phys::kB / (br_->T_ncdm[i] * br_->T_cmb);
            if(background_ncdm_momenta(br_->q_ncdm_bg[i], br_->w_ncdm_bg[i], br_->q_size_ncdm_bg[i], br_->M_ncdm[i], br_->factor_ncdm[i], 0.0, NULL, &rhoNCDM, NULL, NULL, NULL) == _FAILURE_)
            {
                std::stringstream exceptionStr;
                exceptionStr << "CLASS: background_ncdm_momenta failed!" << std::endl << br_->error_message;
                exc.set(exceptionStr.str());
                throw exc;
            }
            br_->Omega0_ncdm[i] = rhoNCDM / (br_->H0 * br_->H0);
            br_->Omega0_ncdm_tot += br_->Omega0_ncdm[i];
        }
    }

    br_->Omega0_lambda = 1 + br_->Omega0_k - br_->Omega0_b - br_->Omega0_cdm - br_->Omega0_g - br_->Omega0_ur - br_->Omega0_ncdm_tot;

    int backgroundInitFailures = 0;
    while(background_init(pr_, br_) == _FAILURE_)
    {
        ++backgroundInitFailures;
        if(backgroundInitFailures >= 100)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: background_init failed!" << std::endl << br_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        output_screen("background_init failed, trying again!" << std::endl);
        std::ofstream outLog("cmb_error_log.txt", std::ios::app);
        outLog << "background_init failed for:" << std::endl;
        outLog << std::setprecision(30) << "h = " << br_->h << std::endl;
        outLog << "Omega0_b = " << br_->Omega0_b << std::endl;
        outLog << "Omega0_cdm = " << br_->Omega0_cdm << std::endl;
        outLog << "Omega0_g = " << br_->Omega0_g << std::endl;
        outLog << "Omega0_ur = " << br_->Omega0_ur << std::endl;
        outLog << "tau = " << th_->tau_reio << std::endl;
        outLog.close();

        br_->h += 1e-5;
        br_->H0 = br_->h * 1e5 / _c_;
    }

    //pr_->k_min_tau0 = kMin_ * br_->conformal_age;
    //pr_->k_max_tau0_over_l_max = kMax_ * br_->conformal_age / lMax_;

    th_->reio_parametrization = reio_camb;
    th_->reio_z_or_tau = reio_tau;
    th_->tau_reio = params.getTau();
    
    const double yHe = params.getYHe();
    if(yHe == 0.0)
        th_->YHe = _BBN_;
    else
    {
        check(yHe > 0 && yHe < 1, "invalid yHe = " << yHe);
        th_->YHe = yHe;
    }

    th_->recombination = recfast;
    
    int thermoInitFailures = 0;
    while(thermodynamics_init(pr_, br_, th_) == _FAILURE_)
    {
        ++thermoInitFailures;
        if(thermoInitFailures >= 100)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: thermodynamics_init failed!" << std::endl << th_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        output_screen("thermodynamics_init failed, trying again!" << std::endl);
        std::ofstream outLog("cmb_error_log.txt", std::ios::app);
        outLog << "thermodynamics_init failed for:" << std::endl;
        outLog << std::setprecision(30) << "h = " << br_->h << std::endl;
        outLog << "Omega0_b = " << br_->Omega0_b << std::endl;
        outLog << "Omega0_cdm = " << br_->Omega0_cdm << std::endl;
        outLog << "Omega0_g = " << br_->Omega0_g << std::endl;
        outLog << "Omega0_ur = " << br_->Omega0_ur << std::endl;
        outLog << "tau = " << th_->tau_reio << std::endl;
        outLog.close();

        background_free(br_);
        br_->Omega0_b -= 1e-4;
        background_init(pr_, br_);
    }

    if(wantMatterPs)
    {
        check(zMaxPk >= 0, "invalid zMaxPk = " << zMaxPk);
        pt_->has_pk_matter = true;
        pt_->has_density_transfers = true;
        sp_->z_max_pk = zMaxPk;
        //nl_->method = nl_halofit;
        //pr_->halofit_dz = 0.1;
        //pr_->halofit_min_k_nonlinear = 0.0035;
        //pr_->halofit_sigma_precision = 0.05;
    }

    pt_->has_scalars = true;
    pt_->has_vectors = false;
    pt_->has_tensors = includeTensors_;
    //pt_->has_cmb = true;
    pt_->has_cls = true;
    pt_->l_scalar_max = lMax_;
    if(includeTensors_)
        pt_->l_tensor_max = lMaxTensors_;

    pt_->has_perturbations = true;
    pt_->has_cl_cmb_temperature = wantT;
    pt_->has_cl_cmb_polarization = wantPol;
    pt_->has_cl_cmb_lensing_potential = wantLensing;
    //pt_->has_source_t = wantT;
    //pt_->has_source_p = wantPol;
    //pt_->has_source_g = wantLensing;
    
    int perturbInitFailures = 0;
    while(perturb_init(pr_, br_, th_, pt_) == _FAILURE_)
    {
        ++perturbInitFailures;
        if(perturbInitFailures >= 100)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: perturb_init failed!" << std::endl << pt_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        output_screen("perturb_init failed, trying again!" << std::endl);
        std::ofstream outLog("cmb_error_log.txt", std::ios::app);
        outLog << "perturb_init failed for:" << std::endl;
        outLog << std::setprecision(30) << "h = " << br_->h << std::endl;
        outLog << "Omega0_b = " << br_->Omega0_b << std::endl;
        outLog << "Omega0_cdm = " << br_->Omega0_cdm << std::endl;
        outLog << "Omega0_g = " << br_->Omega0_g << std::endl;
        outLog << "Omega0_ur = " << br_->Omega0_ur << std::endl;
        outLog << "tau = " << th_->tau_reio << std::endl;
        outLog << "Message: " << pt_->error_message;
        outLog.close();

        thermodynamics_free(th_);
        background_free(br_);
        br_->h += 1e-5;
        br_->H0 = br_->h * 1e5 / _c_;

        if(ncdm != 0)
        {
            br_->N_ncdm = ncdm;

            br_->T_ncdm = (double*) malloc(ncdm * sizeof(double));
            br_->ksi_ncdm = (double*) malloc(ncdm * sizeof(double));
            br_->deg_ncdm = (double*) malloc(ncdm * sizeof(double));
            br_->M_ncdm = (double*) malloc(ncdm * sizeof(double));
            br_->Omega0_ncdm = (double*) malloc(ncdm * sizeof(double));
            br_->m_ncdm_in_eV = (double*) malloc(ncdm * sizeof(double));

            br_->got_files = (int*) malloc(ncdm * sizeof(int));

            br_->ncdm_psd_files = NULL;
            br_->ncdm_psd_parameters = NULL;

            for(int i = 0; i < ncdm; ++i)
            {
                br_->T_ncdm[i] = params.getNCDMParticleTemp(i);
                br_->ksi_ncdm[i] = br_->ksi_ncdm_default;
                br_->deg_ncdm[i] = br_->deg_ncdm_default;
                br_->M_ncdm[i] = 0.0;
                br_->Omega0_ncdm[i] = 0.0;
                br_->m_ncdm_in_eV[i] = params.getNCDMParticleMass(i);

                br_->got_files[i] = 0;
            }

            if(background_ncdm_init(pr_, br_) == _FAILURE_)
            {
                std::stringstream exceptionStr;
                exceptionStr << "CLASS: background_ncdm_init failed!" << std::endl << br_->error_message;
                exc.set(exceptionStr.str());
                throw exc;
            }

            for(int i = 0; i < ncdm; ++i)
            {
                double rhoNCDM;
                br_->M_ncdm[i] = br_->m_ncdm_in_eV[i] * Phys::eCharge / Phys::kB / (br_->T_ncdm[i] * br_->T_cmb);
                if(background_ncdm_momenta(br_->q_ncdm_bg[i], br_->w_ncdm_bg[i], br_->q_size_ncdm_bg[i], br_->M_ncdm[i], br_->factor_ncdm[i], 0.0, NULL, &rhoNCDM, NULL, NULL, NULL) == _FAILURE_)
                {
                    std::stringstream exceptionStr;
                    exceptionStr << "CLASS: background_ncdm_momenta failed!" << std::endl << br_->error_message;
                    exc.set(exceptionStr.str());
                    throw exc;
                }
                br_->Omega0_ncdm[i] = rhoNCDM / (br_->H0 * br_->H0);
                br_->Omega0_ncdm_tot += br_->Omega0_ncdm[i];
            }
        }


        background_init(pr_, br_);
        thermodynamics_init(pr_, br_, th_);

    }

    pm_->n_s = params.getNs();
    pm_->A_s = params.getAs();
    pm_->k_pivot = params.getPivot();
    if(includeTensors_)
    {
        pm_->r = params.getR();
        pm_->n_t = params.getNt();
    }

    if(primordialInitialize_)
    {
        std::stringstream pkFileNameStr;
        pkFileNameStr << "cosmo_pk";
#ifdef COSMO_MPI
        int mpif;
        MPI_Initialized(&mpif);
        if(mpif)
        {
            int n;
            MPI_Comm_size(MPI_COMM_WORLD, &n);
            if(n > 1)
            {
                int p;
                MPI_Comm_rank(MPI_COMM_WORLD, &p);

                pkFileNameStr << '_' << p;
            }
        }
#endif
        pkFileNameStr << ".txt";
        std::string pkFileName = pkFileNameStr.str();

        std::ofstream outPk(pkFileName.c_str());
        if(!outPk)
        {
            std::stringstream exceptionStr;
            exceptionStr << "Cannot output the primordial power spectrum into the file " << pkFileName << ".";
            exc.set(exceptionStr.str());
            throw exc;
        }

        const double decades = std::log(kMax_ / kMin_) / std::log(10.0);
        check(decades > 0, "");
        const int nPoints = int(decades * kPerDecade_);
        check(nPoints > 0, "");
        const double kDelta = (std::log(kMax_) - std::log(kMin_)) / nPoints;

        for(int i = -2; i <= nPoints + 2; ++i)
        {
            //const double k = (i == nPoints ? kMax_ : (i == 0 ? kMin_ : std::exp(std::log(kMin_) + i * kDelta)));
            const double k = std::exp(std::log(kMin_) + i * kDelta);
            outPk << k << ' ' << params.powerSpectrum().evaluate(k);
            
            if(includeTensors_)
                outPk << ' ' << params.powerSpectrumTensor().evaluate(k);
            outPk << std::endl;
        }
        outPk.close();

        std::stringstream classCommand;
        classCommand << "cat " << pkFileName;

        pm_->primordial_spec_type = external_Pk;
        pm_->command = (char *) malloc(classCommand.str().size() + 1);
        std::strcpy(pm_->command, classCommand.str().c_str());
    }

    if(primordial_init(pr_, pt_, pm_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: primordial_init failed!" << std::endl << pm_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    /*
    if(primordialInitialize_)
    {
        // Altering the primordial power spectrum
        const int iMdIndex = pt_->index_md_scalars;
        for(int i = 0; i < pm_->lnk_size; ++i)
        {
            const double lnK = pm_->lnk[i];
            const double k = std::exp(lnK);
            const double pk = params.powerSpectrum().evaluate(k);
            const double lnPk = std::log(pk);
            pm_->lnpk[iMdIndex][i * pm_->ic_ic_size[iMdIndex] + pt_->index_ic_ad] = lnPk;
        }

        if(array_spline_table_lines(pm_->lnk, pm_->lnk_size, pm_->lnpk[iMdIndex], pm_->ic_ic_size[iMdIndex], pm_->ddlnpk[iMdIndex], _SPLINE_EST_DERIV_, pm_->error_message) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: array_spline_table_lines failed!" << std::endl << pm_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(includeTensors_)
        {
            const int iMdIndex = pt_->index_md_tensors;
            for(int i = 0; i < pm_->lnk_size; ++i)
            {
                const double lnK = pm_->lnk[i];
                const double k = std::exp(lnK);
                const double pk = params.powerSpectrumTensor().evaluate(k);
                const double lnPk = std::log(pk);
                pm_->lnpk[iMdIndex][i * pm_->ic_ic_size[iMdIndex] + pt_->index_ic_ad] = lnPk;
            }

            if(array_spline_table_lines(pm_->lnk, pm_->lnk_size, pm_->lnpk[iMdIndex], pm_->ic_ic_size[iMdIndex], pm_->ddlnpk[iMdIndex], _SPLINE_EST_DERIV_, pm_->error_message) == _FAILURE_)
            {
                std::stringstream exceptionStr;
                exceptionStr << "CLASS: array_spline_table_lines failed!" << std::endl << pm_->error_message;
                exc.set(exceptionStr.str());
                throw exc;
            }
        }
    }
    */

    if(nonlinear_init(pr_, br_, th_, pt_, pm_, nl_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: nonlinear_init failed!" << std::endl << nl_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }


    if(transfer_init(pr_, br_, th_, pt_, nl_, tr_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: transfer_init failed!" << std::endl << tr_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }


    sp_->has_tt = wantT;
    sp_->has_ee = wantPol;
    sp_->has_te = (wantT && wantPol);
    sp_->has_bb = wantPol;
    sp_->has_pp = wantLensing;
    sp_->has_tp = (wantT && wantLensing);
    sp_->has_ep = (wantPol && wantLensing);

    if(spectra_init(pr_, br_, pt_, pm_, nl_, tr_, sp_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: spectra_init failed!" << std::endl << sp_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    lensing_ = wantLensing;

    if(wantLensing)
    {
        le_->has_lensed_cls = true;
        le_->has_tt = wantT;
        le_->has_ee = wantPol;
        le_->has_te = (wantT && wantPol);
        le_->has_bb = wantPol;

        if(lensing_init(pr_, pt_, sp_, nl_, le_)  == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: lensing_int failed!" << std::endl << le_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }
    }

    init_ = true;
}

void
CMB::getCl(std::vector<double>* clTT, std::vector<double>* clEE, std::vector<double>* clTE, std::vector<double>* clPP, std::vector<double>* clTP, std::vector<double>* clEP, std::vector<double>* clBB)
{
    StandardException exc;

    check(init_, "need to initialize first");

    check(sp_->has_tt || clTT == NULL, "T has not been requested earlier");
    check(sp_->has_ee || clEE == NULL, "Pol has not been requested earlier");
    check(sp_->has_te || clTE == NULL, "T and Pol have not been requested earlier");
    check(sp_->has_pp || clPP == NULL, "Lensing has not been requested earlier");
    check(sp_->has_tp || clTP == NULL, "T and Lensing have not been requested earlier");
    check(sp_->has_ep || clEP == NULL, "Pol and Lensing have not been requested earlier");
    check(sp_->has_bb || clBB == NULL, "Pol has not been requested earlier");

    const int lMax = lMax_;

    if(clTT)
    {
        clTT->clear();
        clTT->resize(lMax + 1, 0);
    }
    if(clEE)
    {
        clEE->clear();
        clEE->resize(lMax + 1, 0);
    }
    if(clTE)
    {
        clTE->clear();
        clTE->resize(lMax + 1, 0);
    }
    if(clPP)
    {
        clPP->clear();
        clPP->resize(lMax + 1, 0);
    }
    if(clTP)
    {
        clTP->clear();
        clTP->resize(lMax + 1, 0);
    }
    if(clEP)
    {
        clEP->clear();
        clEP->resize(lMax + 1, 0);
    }
    if(clBB)
    {
        clBB->clear();
        clBB->resize(lMax + 1, 0);
    }

    for(int l = 2; l <= lMax; ++l)
    {
        double clTot[1000];
        if(spectra_cl_at_l(sp_, l, clTot, dum1_, dum2_) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: spectra_cl_at_l failed!" << std::endl << sp_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(clTT)
            (*clTT)[l] = br_->T_cmb * br_->T_cmb * 1e12 * clTot[sp_->index_ct_tt];
        if(clEE)
            (*clEE)[l] = br_->T_cmb * br_->T_cmb * 1e12 * clTot[sp_->index_ct_ee];
        if(clTE)
            (*clTE)[l] = br_->T_cmb * br_->T_cmb * 1e12 * clTot[sp_->index_ct_te];
        if(clPP)
            (*clPP)[l] = clTot[sp_->index_ct_pp];
        if(clTP)
            (*clTP)[l] = br_->T_cmb * 1e6 * clTot[sp_->index_ct_tp];
        if(clEP)
            (*clEP)[l] = br_->T_cmb * 1e6 * clTot[sp_->index_ct_ep];
        if(clBB)
            (*clBB)[l] = br_->T_cmb * br_->T_cmb * 1e12 * clTot[sp_->index_ct_bb];
    }
}

void
CMB::getLensedCl(std::vector<double>* clTT, std::vector<double>* clEE, std::vector<double>* clTE, std::vector<double>* clBB)
{
    StandardException exc;

    check(init_, "need to initialize first");
    check(lensing_, "lensing not requested");

    check(le_->has_tt || clTT == NULL, "T has not been requested earlier");
    check(le_->has_ee || clEE == NULL, "Pol has not been requested earlier");
    check(le_->has_te || clTE == NULL, "T and Pol have not been requested earlier");
    check(le_->has_bb || clBB == NULL, "Pol has not been requested earlier");

    const int lMax = le_->l_lensed_max;

    if(clTT)
    {
        clTT->clear();
        clTT->resize(lMax + 1, 0);
    }
    if(clEE)
    {
        clEE->clear();
        clEE->resize(lMax + 1, 0);
    }
    if(clTE)
    {
        clTE->clear();
        clTE->resize(lMax + 1, 0);
    }
    if(clBB)
    {
        clBB->clear();
        clBB->resize(lMax + 1, 0);
    }

    for(int l = 2; l <= lMax; ++l)
    {
        double clTot[1000];
        if(lensing_cl_at_l(le_, l, clTot) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: lensing_cl_at_l failed!" << std::endl << le_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(clTT)
            (*clTT)[l] = br_->T_cmb * br_->T_cmb * 1e12 * clTot[le_->index_lt_tt];
        if(clEE)
            (*clEE)[l] = br_->T_cmb * br_->T_cmb * 1e12 * clTot[le_->index_lt_ee];
        if(clTE)
            (*clTE)[l] = br_->T_cmb * br_->T_cmb * 1e12 * clTot[le_->index_lt_te];
        if(clBB)
            (*clBB)[l] = br_->T_cmb * br_->T_cmb * 1e12 * clTot[le_->index_lt_bb];
    }
}

void
CMB::preClean()
{
    StandardException exc;

    if(init_)
    {
        check(preInit_, "");
        clean();
    }

    if(!preInit_)
        return;

    /*
    if(bessel_free(bs_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: bessel_free failed!" << std::endl << bs_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }
    */
}

void
CMB::clean()
{
    StandardException exc;

    if(!init_)
        return;

    if(lensing_)
    {
        if(lensing_free(le_) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: lensing_free failed!" << std::endl << le_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }

        if(nonlinear_free(nl_) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: nonlinear_free failed!" << std::endl << nl_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }
    }

    if(spectra_free(sp_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: spectra_free failed!" << std::endl << sp_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    if(transfer_free(tr_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: transfer_free failed!" << std::endl << tr_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    if(primordial_free(pm_) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: primordial_free failed!" << std::endl << pm_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

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
}

void
CMB::getTransfer(int l, Math::TableFunction<double, double>* t, Math::TableFunction<double, double>* e, Math::TableFunction<double, double>* p)
{
    StandardException exc;

    check(init_, "need to initialize first");
    check(l >= 2 && l <= lMax_, "invalid value of l = " <<l);

    check(sp_->has_tt || t == NULL, "T has not been requested earlier");
    check(sp_->has_ee || e == NULL, "Pol has not been requested earlier");
    check(sp_->has_pp || p == NULL, "T and Pol have not been requested earlier");

    const int iMdSc = pt_->index_md_scalars;
    const int iIcAd = pt_->index_ic_ad;

    // find index l
    const int lSize = tr_->l_size[iMdSc];
    int indexL = -1;
    for(int i = 0; i < lSize; ++i)
    {
        if(tr_->l[i] == l)
        {
            indexL = i;
            break;
        }
    }

    if(indexL == -1)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Transfer function not calculated for l = " << l << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    if(t)
        t->clear();

    if(e)
        e->clear();

    if(p)
        p->clear();

    for (int i = 0; i < tr_->q_size; ++i)
    {
        const double k = tr_->q[i];
        if(t)
        {
            double val0, val1, val2, val;
            if(transfer_functions_at_q(tr_, iMdSc, iIcAd, tr_->index_tt_t0, indexL, k, &val0) == _FAILURE_)
            {
                std::stringstream exceptionStr;
                exceptionStr << "CLASS: transfer_functions_at_k failed!" << std::endl << tr_->error_message;
                exc.set(exceptionStr.str());
                throw exc;
            }
            if(transfer_functions_at_q(tr_, iMdSc, iIcAd, tr_->index_tt_t1, indexL, k, &val1) == _FAILURE_)
            {
                std::stringstream exceptionStr;
                exceptionStr << "CLASS: transfer_functions_at_k failed!" << std::endl << tr_->error_message;
                exc.set(exceptionStr.str());
                throw exc;
            }
            if(transfer_functions_at_q(tr_, iMdSc, iIcAd, tr_->index_tt_t2, indexL, k, &val2) == _FAILURE_)
            {
                std::stringstream exceptionStr;
                exceptionStr << "CLASS: transfer_functions_at_k failed!" << std::endl << tr_->error_message;
                exc.set(exceptionStr.str());
                throw exc;
            }
            val = val0 + val1 + val2;
            (*t)[k] = val;
        }

        if(e)
        {
            double val;
            if(transfer_functions_at_q(tr_, iMdSc, iIcAd, tr_->index_tt_e, indexL, k, &val) == _FAILURE_)
            {
                std::stringstream exceptionStr;
                exceptionStr << "CLASS: transfer_functions_at_k failed!" << std::endl << tr_->error_message;
                exc.set(exceptionStr.str());
                throw exc;
            }
            (*e)[k] = val;
        }

        if(p)
        {
            double val;
            if(transfer_functions_at_q(tr_, iMdSc, iIcAd, tr_->index_tt_lcmb, indexL, k, &val) == _FAILURE_)
            {
                std::stringstream exceptionStr;
                exceptionStr << "CLASS: transfer_functions_at_k failed!" << std::endl << tr_->error_message;
                exc.set(exceptionStr.str());
                throw exc;
            }
            (*p)[k] = val;
        }
    }
}

void
CMB::getMatterPs(double z, Math::TableFunction<double, double>* ps)
{
    StandardException exc;
    check(init_, "need to initialize first");
    check(pt_->has_pk_matter, "matter ps not requested");
    check(z >= 0 && z <= sp_->z_max_pk, "invalid z = " << z);
    const int kSize = sp_->ln_k_size;
    check(kSize > 0, "");

    // To be done better
    check(kSize < 10000, "");
    double outTot[100000];
    double outIc[100000];

    if(spectra_pk_at_z(br_, sp_, linear, z, outTot, outIc) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: spectra_pk_at_z failed!" << std::endl << sp_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }

    ps->clear();

    for(int i = 0; i < kSize; ++i)
    {
        const double k = std::exp(sp_->ln_k[i]);
        (*ps)[k] = outTot[i];
    }
}

void
CMB::getMatterTransfer(double z, Math::TableFunction<double, double>* tk)
{
    StandardException exc;
    check(init_, "need to initialize first");
    check(pt_->has_density_transfers, "matter ps not requested");
    check(z >= 0 && z <= sp_->z_max_pk, "invalid z = " << z);
    const int kSize = sp_->ln_k_size;
    check(kSize > 0, "");

    // To be done better
    check(kSize < 10000, "");
    double out[100000];

    tk->clear();

    for(int i = 0; i < kSize; ++i)
    {
        const double k = std::exp(sp_->ln_k[i]);
        if(spectra_tk_at_k_and_z(br_, sp_, k, z, out) == _FAILURE_)
        {
            std::stringstream exceptionStr;
            exceptionStr << "CLASS: spectra_tk_at_k_and_z failed!" << std::endl << sp_->error_message;
            exc.set(exceptionStr.str());
            throw exc;
        }
        (*tk)[k] = out[sp_->index_tr_delta_tot];
    }
}

double
CMB::sigma8()
{
    StandardException exc;
    check(init_, "need to initialize first");
    check(pt_->has_density_transfers, "matter ps not requested");

    // from CLASS
    /*
    double sigma;
    if(spectra_sigma(br_, pm_, sp_, 8.0 / params_->getH(), 0.0, &sigma) == _FAILURE_)
    {
        std::stringstream exceptionStr;
        exceptionStr << "CLASS: spectra_sigma failed!" << std::endl << sp_->error_message;
        exc.set(exceptionStr.str());
        throw exc;
    }
    return sigma;
    */

    Math::TableFunction<double, double> pk;
    getMatterPs(0.0, &pk);
    check(!pk.empty(), "");

    double res = 0;
    int N = 10000;
    double logKMin, logKMax;
    Math::TableFunction<double, double>:: const_iterator it = pk.begin();
    logKMin = std::log((*it).first);
    it = pk.end();
    --it;
    logKMax = std::log((*it).first);
    check(logKMax > logKMin, "");
    const double logKDelta = (logKMax - logKMin) / N;
    double yPrev = 0;
    for(int i = 1; i <= N; ++i)
    {
        const double x = logKMin + i * logKDelta;
        double k = std::exp(x);
        if(i == N)
            k = (*it).first;

        const double v = 8 * k / (params_->getH());
        check(v > 0, "");
        const double w = 3 * (std::sin(v) - v * std::cos(v)) / (v * v * v);

        const double p = pk.evaluate(k);
        const double y = k * k * k * w * w * p / (2 * Math::pi * Math::pi);
        res += (y + yPrev) * logKDelta / 2;
        yPrev = y;
    }

    return std::sqrt(res);
}


/*
extern "C"
{
    void getcl_(double* omb, double* omc, double* h, double* tau, double* ns, double* pivot, int* lmax, double* cl);
}

void
CMB::calculateCl(const CosmologicalParams& params, int lMax, std::vector<double>& cl)
{
    check(lMax >= 2, "invalid lMax " << lMax);
    cl.clear();
    cl.resize(lMax + 1);

    double myomb = params.getOmB();
    double myomc = params.getOmC();
    double myh = params.getH();
    double mytau = params.getTau();
    double myns = params.getNs();
    double mypivot = params.getPivot();

    getcl_(&myomb, &myomc, &myh, &mytau, &myns, &mypivot, &lMax, &(cl[0]));

    for(int l = 2; l <= lMax; ++l)
    {
        cl[l] *= params.getAs() * params.getTemperature() * params.getTemperature() * 1e12; // muK^2
        cl[l] *= (2.0 * Math::pi) / (double(l) * (l + 1));
    }
}
*/

