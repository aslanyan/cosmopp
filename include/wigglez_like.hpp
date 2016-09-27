#pragma once
#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <string>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <macros.hpp>
#include <matrix_impl.hpp>
#include <cubic_spline.hpp>
#include <likelihood_function.hpp>
#include <cmb.hpp>

#include <class.h>

/// WiggleZ Likelihood (written by Nicolas Canac)
class WiggleZLikelihood : public Math::LikelihoodFunction
{
public:
    /// Constructor.
    /// \param path The path to the data directory. It must then contain a WiggleZ directory with the corresponding data files.
    /// \param cmb A reference to a CMB object to be used.
    /// \param redshiftBin The redshift bin. Must be either 'a' or 'b' or 'c' or 'd'.
    /// \param initializeCMBAtEachStep Determines if the cmb object should be initialized every time setCosmoParams is called. Note that the WiggleZ likelihood may be used in combination with other likelihoods, such as Planck, for some kind of parameter space scanning (MCMC, MultiNest). In this case the slowest part is the initialization of CMB. It will make sense to have one CMB object serve multiple likelihoods, so the initialize of CMB will be called outside. If this is the case then this parameter needs to be set to false. The assumption then will be that the CMB object will be initialized with the same cosmological parameters as setCosmoParams.
    WiggleZLikelihood(std::string path, CMB& cmb, char redshiftBin, bool initializeCMBAtEachStep = true) : initializeCMBAtEachStep_(initializeCMBAtEachStep), cmb_(&cmb)
    {
        //bool Q_marge = false;
        //double Q_mid = 4.0;
        //double Q_sigma = 1.5;
        //double Ag = 1.4;
        //bool Use_jennings = false;
        //bool Use_simpledamp = false;
        giggleZ_fidpk_size_ = 500;
        used_region_ = std::vector<bool> {true, true, true, true, true, true, true};

        switch(redshiftBin)
        {
        case 'a':
            d_angular_fid_ = 736.293;
            d_radial_fid_ = 3876.81;
            //double sigmav = 360;
            giggleZ_fidpoly_ = std::vector<double> {4.61900,-13.7787,58.9410,-175.240,284.321,-187.284};
            redshift_ = 0.22;
            break;
        case 'b':
            d_angular_fid_ = 1134.87;
            d_radial_fid_ = 3511.96;
            //double sigmav = 308;
            giggleZ_fidpoly_ = std::vector<double> {4.63079, -12.6293, 42.9265, -91.8068, 97.808, -37.633};
            redshift_ = 0.41;
            break;
        case 'c':
            d_angular_fid_ = 1396.05;
            d_radial_fid_ = 3160.38;
            //double sigmav = 325;
            giggleZ_fidpoly_ = std::vector<double> {4.69659, -12.7287, 42.5681, -89.5578, 96.664, -41.2564};
            redshift_ = 0.60;
            break;
        case 'd':
            d_angular_fid_ = 1558.68;
            d_radial_fid_ = 2852.95;
            //double sigmav = 212;
            giggleZ_fidpoly_ = std::vector<double> {4.6849, -13.4747, 53.7172, -145.832, 216.638, -132.782};
            redshift_ = 0.78;
            break;
        default:
            check(false, "invalid redshift bin " << redshiftBin << "! Must be either a, b, c, or d.");
            break;
        }

        // Number of points and kbands in the input files
        num_mpk_points_full_ = 50;
        num_mpk_kbands_full_ = 100;

        // Decide which bandpowers to use, min to max
        min_mpk_points_use_ = 3;
        max_mpk_points_use_ = 20; // cut off at k/h~0.2
        min_mpk_kbands_use_ = 1;
        max_mpk_kbands_use_ = 100;

        k_size_ = max_mpk_kbands_use_ - min_mpk_kbands_use_ + 1;
        int mu_size = 1;
        k_.resize(k_size_);
        kh_.resize(k_size_);

        // Read in data file containing kbands
        root_ = path + "/WiggleZ/";
        std::ifstream datafile(root_ + "wigglez_nov11_kbands.txt");
        for(int i = 0; i < num_mpk_kbands_full_; ++i)
            if((i+2 > min_mpk_kbands_use_) && (i < max_mpk_kbands_use_))
                datafile >> kh_[i-min_mpk_kbands_use_+1];
        datafile.close();

        double khmax = kh_[k_size_-1];

        // Read in data file containing fiducial power spectrum
        // to determine k_fid_size and ifid_discard. These are
        // basically parameters that select the values of k in
        // the fiducial model that match with the values of k
        // from the data.
        int ifid_discard;
        datafile.open(root_ + "gigglezfiducialmodel_matterpower_" + redshiftBin + ".dat");
        double kdum, dum1;
        datafile >> kdum >> dum1;
        int line_number = 1;
        while(kdum < kh_[0])
        {
            check(line_number <= giggleZ_fidpk_size_, "read too many lines");
            datafile >> kdum >> dum1;
            ++line_number;
        }
        ifid_discard = line_number - 2;
        while(kdum < khmax)
        {
            check(line_number <= giggleZ_fidpk_size_, "read too many lines");
            datafile >> kdum >> dum1;
            ++line_number;
        }
        datafile.close();
        k_fid_size_ = line_number - ifid_discard + 1;
        khmax = kdum;

        khmax = khmax * 2; // If using halofit. Why?

        num_regions_ = used_region_.size();
        num_regions_used_ = 0;
        for(int i = 0; i < num_regions_; ++i)
            if(used_region_[i])
                ++num_regions_used_;
        check(num_regions_used_ > 0, "Number of regions used is not greater than zero.");

        // Read in window functions
        n_size_ = max_mpk_points_use_ - min_mpk_points_use_ + 1;
        window_.resize(num_regions_, Math::Matrix<double>(n_size_, k_size_, 0));
        datafile.open(root_ + "wigglez_nov11" + redshiftBin + "_windows.txt");
        std::string line;
        for(int i_region = 0; i_region < num_regions_; ++i_region)
        {
            if(num_regions_ > 1)
                getline(datafile, line);
            for(int i = 0; i < num_mpk_points_full_; ++i)
            {
                // Read a line and store in  vector
                double dummy;
                std::vector<double> vec;
                getline(datafile, line);
                std::istringstream iss(line);
                while(iss >> dummy)
                    vec.push_back(dummy);
                if(i == 0 && i_region == 0)
                    check(vec[0] >= -1.01 && vec[0] <= -0.99, "Check failed: Reading window functions");
                if((i+2) > min_mpk_points_use_ && i < max_mpk_points_use_)
                {
                    for(int j = 0; j < k_size_; ++j)
                        window_[i_region](i-min_mpk_points_use_+1, j) = vec[j+min_mpk_kbands_use_-1];
                }
            }
        }
        datafile.close();

        // Read in measurements
        P_obs_.resize(num_regions_, n_size_, 0); // P(k) / h^-3Mpc^3
        P_err_.resize(num_regions_, n_size_, 0);
        datafile.open(root_ + "wigglez_nov11" + redshiftBin + "_measurements.txt");
        // Skip first two lines
        for(int i_region = 0; i_region < num_regions_; ++i_region)
        {
            std::getline(datafile, line);
            std::getline(datafile, line);
            for(int i = 0; i < num_mpk_points_full_; ++i)
            {
                std::getline(datafile, line);
                if((i+2 > min_mpk_points_use_) && (i < max_mpk_points_use_))
                {
                    std::istringstream iss(line);
                    double khdum, klodum, khidum, pdum, errdum, dum;
                    iss >> khdum >> klodum >> khidum >> pdum >> errdum >> dum;
                    P_obs_(i_region, i-min_mpk_points_use_+1) = pdum;
                    P_err_(i_region, i-min_mpk_points_use_+1) = errdum;
                }
            }
        }
        datafile.close();

        // Read in covariance matrices
        invcov_.resize(num_regions_, Math::Matrix<double>(n_size_, n_size_, 0));
        Math::Matrix<double> cov(n_size_, n_size_, 0);
        Math::Matrix<double> invcov_tmp(n_size_, n_size_, 0);
        datafile.open(root_ + "wigglez_nov11" + redshiftBin + "_cov.txt");
        for(int i_region = 0; i_region < num_regions_; ++i_region)
        {
            std::getline(datafile, line);
            for(int i = 0; i < num_mpk_points_full_; ++i)
            {
                double dummy;
                std::vector<double> vec;
                std::getline(datafile, line);
                std::istringstream iss(line);
                while(iss >> dummy)
                    vec.push_back(dummy);
                if(i == 0 && i_region == 0 && redshiftBin == 'a')
                    check(vec[0] >= -1.01 && vec[0] <= -0.99, "Check failed: Reading covariance matrices");
                if((i+2 > min_mpk_points_use_) && (i < max_mpk_points_use_))
                {
                    for(int j = 0; j < num_mpk_points_full_; ++j)
                    {
                        if((j+2 > min_mpk_points_use_) && (j < max_mpk_points_use_))
                            cov(i-min_mpk_points_use_+1, j-min_mpk_points_use_+1) = vec[j];
                    }
                }
            }
            cov.getInverse(&invcov_tmp);
            check(invcov_tmp.rows() == n_size_ && invcov_tmp.cols() == n_size_,
                    "Check failed: Wrong number of rows and columns in covariance matrix.")
            for(int i = 0; i < n_size_; ++i)
                for(int j = 0; j < n_size_; ++j)
                    invcov_[i_region](i, j) = invcov_tmp(i, j);
        }
        datafile.close();

        // Read in fiducial model again, but different from first time
        P_fid_.resize(k_fid_size_, 1, 0);
        k_fid_.resize(k_fid_size_, 1, 0);
        datafile.open(root_ + "gigglezfiducialmodel_matterpower_" + redshiftBin + ".dat");
        //Skip first ifid_discard lines
        for(int i = 0; i < ifid_discard; ++i)
            std::getline(datafile, line);
        for(int i = 0; i < k_fid_size_; ++i)
        {
            std::getline(datafile, line);
            std::istringstream iss(line);
            double kdummy, pdummy;
            iss >> kdummy >> pdummy; 
            k_fid_(i, 0) = kdummy;
            P_fid_(i, 0) = pdummy;
        }
        datafile.close();
    }

    /// Destructor.
    ~WiggleZLikelihood() {}

    /// Calculate the likelihood. Note that setModelCosmoParams must be called before this function.
    /// \param params A pointer to the array of parameters.
    /// \param nPar The number of parameters.
    /// \return -2ln(likelihood).
    double calculate(double* params, int nPar)
    {
        // Check to see that modelParams_ is set. This is set by calling setModelCosmoParams().
        check(modelParams_, "model params must be set before calling this function");
        // Check to see that vModel_ is not empty
        check(!vModel_.empty(), "");
        const int nModel = vModel_.size();
        
        // Check that number of parameters passed in is same as number of parameters in nModel
        check(nPar == nModel, "");
    
        // Set all the parameters in vModel_ to the values in params
        for(int i = 0; i < nModel; ++i)
            vModel_[i] = params[i];
    
        // Sets the parameters in modelParams_ to the values in vModel_
        modelParams_->setAllParameters(vModel_);
        // Set the cosmological parameters to modelParams_.
        setCosmoParams(*modelParams_);
    
        return likelihood();
    }

    /// Full likelihood
    /// \return -2ln(likelihood).
    double likelihood()
    {
        double h = params_->getH();

        //output_screen("h: " << h << std::endl);

        double d_angular = cmb_->getAngularDistance(redshift_);

        //output_screen("d_angular: " << d_angular << std::endl);

        std::vector<double> z_array {redshift_};
        std::vector<double> rHz = cmb_->z_of_r(z_array);
        double r = rHz[0];
        double Hz = rHz[1];
        //output_screen("r, Hz: " << r << " " << Hz << std::endl);
        double d_radial = 1.0/Hz;
        //output_screen("d_radial: " << d_radial << std::endl);

        double scaling = pow(pow(d_angular_fid_/d_angular, 2.0)*(d_radial_fid_/d_radial), 1.0/3.0);
        // Rescale k_ and P_obs with correct factors of h
        for(int i = 0; i < k_size_; ++i)
            k_[i] = kh_[i]*h*scaling;

        // Calculate the matter power spectrum
        Math::Matrix<double> P(k_fid_size_, 1, 0);
        //Math::TableFunction<double, double> P_function;
        //cmb_->getMatterPsNL2(redshift_, &P_function);
        for(int i = 0; i < k_fid_size_; ++i)
        {
            //P(i, 0) = P_function.evaluate(k_fid_(i, 0)*h);
            P(i, 0) = cmb_->getPkNLatk(k_fid_(i, 0)*h, redshift_);
            //output_screen(k_fid_(i, 0) << " " << P(i, 0) << std::endl);
            double power = 0;
            for(int j = 0; j < 6; ++j)
                power += giggleZ_fidpoly_[j]*pow(k_fid_(i, 0), j);
            // rescale P by fiducial model convert to (Mpc/h)^3
            P(i, 0) = P(i, 0) * pow(10.0, power) * pow(h/scaling, 3.0) / P_fid_(i, 0);
        }

        // get P_lin by interpolation
        Math::TableFunction<double, double> P_lin_function;
        for(int i = 0; i < k_fid_size_; ++i)
            P_lin_function[k_fid_(i, 0)] = P(i, 0);

        Math::Matrix<double> P_lin(k_size_, 1, 0);
        for(int i = 0; i < k_size_; ++i)
        {
            P_lin(i, 0) = P_lin_function.evaluate(kh_[i]);
        }

        Math::Matrix<double> W_P_th(n_size_, 1, 0);

        // Power spectrum data
        Math::Matrix<double> P_data_large(n_size_*num_regions_used_, 1, 0);
        // Windowed matter power spectrum
        Math::Matrix<double> W_P_th_large(n_size_*num_regions_used_, 1, 0);
        // Covariance matrices
        Math::Matrix<double> cov_dat_large(n_size_*num_regions_used_, 1, 0);
        Math::Matrix<double> cov_th_large(n_size_*num_regions_used_, 1, 0);
        
        Math::Matrix<double> P_th = P_lin;

        for(int i_region = 0; i_region < num_regions_; ++i_region)
        {
            if(used_region_[i_region])
            {
                int imin = i_region * n_size_;
                int imax = (i_region + 1) * n_size_ - 1;
                
                Math::Matrix<double>::multiplyMatrices(window_[i_region], P_th, &W_P_th);
                for(int i = 0; i < n_size_; ++i)
                {
                    P_data_large(imin+i, 0) = P_obs_(i_region, i);
                    W_P_th_large(imin+i, 0) = W_P_th(i, 0);

                    Math::Matrix<double> invcov_row(1, n_size_);
                    Math::Matrix<double> P_obs_row(n_size_, 1);
                    for(int j = 0; j < n_size_; ++j)
                    {
                        invcov_row(0, j) = invcov_[i_region](i, j);
                        P_obs_row(j, 0) = P_obs_(i_region, j);
                    }
                    Math::Matrix<double> tempMat;
                    Math::Matrix<double>::multiplyMatrices(invcov_row, P_obs_row, &tempMat);
                    cov_dat_large(imin+i, 0) = tempMat(0, 0);
                    Math::Matrix<double>::multiplyMatrices(invcov_row, W_P_th, &tempMat);
                    cov_th_large(imin+i, 0) = tempMat(0, 0);
                }
            }
        }

        // Calculate normV
        double normV = 0;
        Math::Matrix<double> tempMat;
        Math::Matrix<double>::multiplyMatrices(W_P_th_large.getTranspose(), cov_th_large, &tempMat);
        normV = normV + tempMat(0, 0);

        // Calculate bias factor
        //double b_out = 0;
        //Math::Matrix<double>::multiplyMatrices(W_P_th.getTranspose(), cov_dat, &tempMat);
        //b_out = b_out + tempMat(0, 0);
        //Math::Matrix<double>::multiplyMatrices(W_P_th.getTranspose(), cov_th, &tempMat);
        //b_out = b_out / tempMat(0, 0);

        double chisq = 0;
        Math::Matrix<double>::multiplyMatrices(P_data_large.getTranspose(), cov_dat_large, &tempMat);
        chisq = tempMat(0, 0);
        Math::Matrix<double>::multiplyMatrices(W_P_th_large.getTranspose(), cov_dat_large, &tempMat); 
        chisq = chisq - pow(tempMat(0, 0), 2.0) / normV;

        // Return -2ln(L)
        return chisq;
    }

    /// Set the cosmological parameters. This sets the parameters for which the WiggleZ likelihoods will be subsequently calculated.
    /// The initializeCMBAtEachStep parameter of the constructor is important for this function. If that parameter was set to true
    /// then the cmb object will be initialized with these parameters. However, if it was set to false then it is expected that the cmb
    /// object will be initialized separately outside this function with the same parameters.
    /// \param params The cosmological parameters object.
    void setCosmoParams(const CosmologicalParams& params)
    {
        params_ = &params;
        // Does wantT need to be true?
        if(initializeCMBAtEachStep_)
            cmb_->initialize(params, true, false, false, true, 1.0);
    }

    /// Set the model cosmological parameters. This function must be called before the calculate function.
    /// The model cosmological parameters will be assumed to be the same model as the parameters passed into calculate.
    void setModelCosmoParams(CosmologicalParams *params)
    {
        modelParams_ = params;
        modelParams_->getAllParameters(vModel_);
    }

private:
    CMB* cmb_;
    const bool initializeCMBAtEachStep_;

    const CosmologicalParams* params_;

    CosmologicalParams* modelParams_;
    std::vector<double> vModel_;

    int num_mpk_points_full_;
    int num_mpk_kbands_full_;

    int min_mpk_points_use_;
    int max_mpk_points_use_;
    int min_mpk_kbands_use_;
    int max_mpk_kbands_use_;

    int k_size_;
    int n_size_;
    int k_fid_size_;
    int num_regions_used_;
    int num_regions_;
    double zerowindowfxnsubdatnorm_;

    // Parameters
    int giggleZ_fidpk_size_;
    std::vector<bool> used_region_;
    double d_angular_fid_;
    double d_radial_fid_;
    std::vector<double> giggleZ_fidpoly_;
    double redshift_;

    // Data vectors
    std::vector<double> kh_;
    std::vector<double> k_;
    std::vector< Math::Matrix<double> > window_;
    Math::Matrix<double> zerowindowfxn_;
    Math::Matrix<double> zerowindowfxnsubdat_;
    //std::vector< std::vector<double> > window_;
    Math::Matrix<double> P_obs_;
    //std::vector<double> P_obs_;
    Math::Matrix<double> P_err_;
    //std::vector<double> P_err_;
    std::vector< Math::Matrix<double> > invcov_;
    //Math::Matrix<double> invcov_;
    //std::vector< std::vector<double> > invcov_;
    Math::Matrix<double> k_fid_;
    Math::Matrix<double> P_fid_;

    // Location of data/model files
    std::string root_;
};
