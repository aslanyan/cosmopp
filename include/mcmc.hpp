#ifndef COSMO_PP_MCMC_HPP
#define COSMO_PP_MCMC_HPP

#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <ctime>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <likelihood_function.hpp>
#include <random.hpp>
#include <matrix_impl.hpp>

namespace Math
{

/// An abstract class for a prior function, used in the MetropolisHastings class.
class PriorFunctionBase
{
public:
    /// Calculate the prior.
    /// \param params The parameters vector (passed as a pointer to the first element).
    /// \param nPar The number of parameters.
    /// \return The prior distribution value.
    virtual double calculate(double* params, int nPar) = 0;
};

/// An abstract class for a proposal distribution, used in the MetropolisHastings class.
class ProposalFunctionBase
{
public:
    /// Generate the next set of parameters.
    /// \param params A vector of all of the parameters of the PREVIOUS sample (passed as a pointer to the first element).
    /// \param nPar The total number of parameters.
    /// \param blockParams A vector to contain the PROPOSED parameters of a given block (passed as a pointer to the first element).
    /// \param i The index of the block to be generated.
    virtual void generate(double* params, int nPar, double* blockParams, int i) = 0;

    /// Calculate the proposal distribution value.
    /// \param params A vector of all of the parameters of the PREVIOUS sample (passed as a pointer to the first element).
    /// \param nPar The total number of parameters.
    /// \param blockParams A vector containing the block of PROPOSED parameters (passed as a pointer to the first element).
    /// \param i The index of the block.
    /// \return The value of the proposal distribution.
    virtual double calculate(double* params, int nPar, double* blockParams, int i) = 0;
    
    /// Tells if the proposal distribution is symmetric for the given block.
    /// \param i The index of the block.
    /// \return true if the proposal distribution is symmetric for block i.
    virtual bool isSymmetric(int i) = 0;
};

/// A Metropolis-Hastings scanner.
class MetropolisHastings
{
private:
    enum PRIOR_MODE { UNIFORM_PRIOR = 0, GAUSSIAN_PRIOR, PRIOR_MODE_MAX };

public:
    enum CONVERGENCE_DIAGNOSTIC { GELMAN_RUBIN = 0, ACCURACY, CONVERGENCE_DIAGNOSTIC_MAX };

    /// Constructor.
    /// \param nPar The number of parameters.
    /// \param like The likelihood function.
    /// \param fileRoot The root for filenames produced by MetropolisHastings.
    /// \param seed A random seed. If set to 0 (the default value), it will be determined from the current time.
    MetropolisHastings(int nPar, LikelihoodFunction& like, std::string fileRoot, time_t seed = 0, bool isLikelihoodApproximate = false);

    /// Destructor.
    ~MetropolisHastings();

    /// Define a given parameter to have a uniform prior. One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param min The minimum value of the parameter (the lower bound for the prior).
    /// \param max The maximum value of the parameter (the upper bound for the prior).
    /// \param starting The starting value of the parameter. If not set, it will be set to the midpoint of the range by default.
    /// \param startingWidth The starting value will be chosen within this width of starting parameter. Useful for multiple chains to start at different places. If not set, by default it will be set to 1/100-th of the width of the range.
    /// \param samplingWidth The sampling width of the parameter (the width of the Gaussian proposal distribution). If not set, by default it will be set to 1/100-th of the width of the range.
    /// \param accuracy The accuracy with which the parameter needs to be determined (used to choose the stopping time). If not set, by default it will be set to 1/10-th of the sampling width.
    void setParam(int i, const std::string& name, double min, double max, double starting = std::numeric_limits<double>::max(), double startingWidth = 0.0, double samplingWidth = 0.0, double accuracy = 0.0);

    /// Define a given parameter to have a gaussian prior. One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param mean The mean of the prior
    /// \param sigma The sigma of the prior
    /// \param starting The starting value of the parameter. If not set, it will be set to the midpoint of the range by default.
    /// \param startingWidth The starting value will be chosen within this width of starting parameter. Useful for multiple chains to start at different places. If not set, by default it will be set to 1/100-th of the width of the range.
    /// \param samplingWidth The sampling width of the parameter (the width of the Gaussian proposal distribution). If not set, by default it will be set to 1/100-th of the width of the range.
    /// \param accuracy The accuracy with which the parameter needs to be determined (used to choose the stopping time). If not set, by default it will be set to 1/10-th of the sampling width.
    void setParamGauss(int i, const std::string& name, double mean, double sigma, double starting = std::numeric_limits<double>::max(), double startingWidth = 0.0, double samplingWidth = 0.0, double accuracy = 0.0);

    /// Get the name of a parameter.
    /// \param i The index of the parameter.
    /// \return The name of the parameter.
    const std::string& getParamName(int i) const { check(i >= 0 && i < n_, "invalid index " << i); return paramNames_[i]; }

    /// Set the blocks in which the parameters are varied. If this function is not called, each paramter will be assigned to a separate block, by default.
    /// \param blocks A vector defining the indices of the parameters in each block. Each element of the vector is the index following the end of the corresponding block. There are as many elements as there are blocks. For example, if all of the parameters are to belong to one block, the vector should contain one element with value equal to the number of the parameters.
    void specifyParameterBlocks(const std::vector<int>& blocks);

    /// Set an external prior function for all of the parameters. The values set by setParam or setParamGauss will then be ignored. 
    /// One of these functions still needs to be called for each parameter to set their names, starting values, sampling widths, and accuracies.
    /// \param prior A pointer to the external prior function.
    void useExternalPrior(PriorFunctionBase* prior) { externalPrior_ = prior; }

    /// Set an external proposal distribution for all of the parameters. The sampling width value set by setParam or setParamGauss will then be ignored.
    /// One of these functions still needs to be called for each parameter to set their names, priors, starting values, and accuracies.
    /// \param proposal A pointer to the external proposal distribution.
    void useExternalProposal(ProposalFunctionBase* proposal) { externalProposal_ = proposal; }

    /// Run the scan. Should be called after all of the other necessary functions have been called to set all of the necessary settings. The resulting chain is written in the file (fileRoot).txt. The first column is the number of repetitions of the element, the second column is -2ln(likelihood), the following columns are the values of all of the parameters.
    /// \param maxChainLength The maximum length of the chain (1000000 by default). The scan will stop when the chain reaches that length, even if the required accuracy for the parameters has not been achieved. If the accuracies are achieved earlier the scan will stop earlier.
    /// \param writeResumeInformationEvery Defines if resume information should be written in a file and how often. This will allow an interrupted run to resume. 0 will mean no resume information will be written. The default setting of 1 is recommended in most cases. However, if the likelihood calculation is very fast, so that the likelihood computing time is faster or comparable to writing out a small binary file, this parameter should be set to higher value. The reason is that it will slow down the scan significantly, and the chance of the resume file being corrupt and useless will be high (this will happen if the code is stopped during writing out the resume file).
    /// \param burnin The burnin length. These elements will still be written out into the chain but will be ignored for determining convergence.
    /// \param cd Convergence diagnostic to be used.
    /// \param convergenceCriterion A number used to determine convergence. For Gelman-Rubin diagnostic this is the number below which (R - 1) absolute values need to be for all the parameters.
    /// \param adaptiveProposal This turns on the usage of the Adaptive Metropolis algorithm (optional, true by default). The proposal distribution will be continuously updated during the run based on the covariance of the existing elements. This typically speeds up the run by about 1 order of magnitude! HIBHLY RECOMMENDED to keep this argument true.
    /// \return The number of chains generated.
    int run(unsigned long maxChainLength = 1000000, int writeResumeInformationEvery = 1, unsigned long burnin = 0, CONVERGENCE_DIAGNOSTIC cd = ACCURACY, double convergenceCriterion = 0.01, bool adaptiveProposal = true);

private:
    void useAdaptiveProposal();

    inline double uniformPrior(double min, double max, double x) const;
    inline double gaussPrior(double mean, double sigma, double x) const;
    inline double calculatePrior();
    inline void calculateStoppingData();
    inline bool stop();
    inline bool checkStoppingCrit();
    inline double generateNewPoint(int i) const { return current_[i] + generator_->generate() * samplingWidth_[i]; }
    inline void openOut(bool append);
    inline void closeOut() { out_.close(); }
    inline void writeChainElement();
    inline void update();
    inline void writeResumeInfo() const;
    inline bool readResumeInfo();

    inline void writeCommInfo(std::ofstream& out) const;
    inline void readCommInfo(std::ifstream& in);
    bool synchronizeCommInfo();

    inline void logProgress() const;

    inline bool isMaster() const { return currentChainI_ == 0; }
    void communicate();
    void sendHaveStopped();

    inline void calculateMeanVar(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end, double& mean, double& var);

    struct BadResumeInfo
    {
        BadResumeInfo(int a) : n(a) {}
        int n;
    };
    
private:
    int n_;
    LikelihoodFunction* like_;
    LikelihoodFunction* spareLike_;
    std::string fileRoot_, resumeFileName_;
    std::vector<double> param1_, param2_, starting_, samplingWidth_, accuracy_;
    std::vector<PRIOR_MODE> priorMods_;
    std::vector<std::string> paramNames_;
    std::vector<double> paramSum_, paramSquaredSum_, corSum_, reachedSigma_;
    PriorFunctionBase* externalPrior_;
    ProposalFunctionBase* externalProposal_;
    std::vector<int> blocks_;

    unsigned long covarianceElementsNum_;
    Math::SymmetricMatrix<double> covariance_;
    std::vector<double> generatedVec_, rotatedVec_;
    Math::SymmetricMatrix<double> cholesky_;

    bool covarianceReady_;
    std::vector<double> paramMean_;
    std::vector<double> paramMeanNew_;
    const double covEpsilon_;
    const double covFactor_;
    bool adapt_;

    std::vector<double> rGelmanRubin_;

    CONVERGENCE_DIAGNOSTIC cd_;
    double cc_;

    time_t seed_;
    Math::UniformRealGenerator* uniformGen_;
    Math::GaussianGenerator* generator_;

    std::vector<double> prev_, current_;
    unsigned long maxChainLength_;
    unsigned long iteration_;
    double currentLike_;
    double currentPrior_;

    struct CommunicationInfo
    {
        CommunicationInfo(int n = 0) : sums(n), sqSums(n), stdMean(n) {}
        CommunicationInfo(const CommunicationInfo& other) : sums(other.sums), sqSums(other.sqSums), stdMean(other.stdMean), iter(other.iter) {}

        std::vector<double> sums, sqSums, stdMean;
        double iter;

        inline void writeIntoFile(std::ofstream& out) const
        {
            const int n = sums.size();
            check(sqSums.size() == n, "");
            check(stdMean.size() == n, "");

            out.write((char*)(&iter), sizeof(iter));
            out.write((char*)(&(sums[0])), n * sizeof(double));
            out.write((char*)(&(sqSums[0])), n * sizeof(double));
            out.write((char*)(&(stdMean[0])), n * sizeof(double));
        };

        inline void readFromFile(std::ifstream& in)
        {
            const int n = sums.size();
            check(sqSums.size() == n, "");
            check(stdMean.size() == n, "");

            in.read((char*)(&iter), sizeof(iter));
            in.read((char*)(&(sums[0])), n * sizeof(double));
            in.read((char*)(&(sqSums[0])), n * sizeof(double));
            in.read((char*)(&(stdMean[0])), n * sizeof(double));
        }
    };
    // first index is chain index, second index is the communication number, third index is parameter index;
    std::vector<std::list<CommunicationInfo> > commInfo_;


    std::vector<std::vector<double> > communicationBuff_;
    std::vector<double> myStdMean_;

    std::vector<std::vector<double> > sendComBuff_;

    std::vector<std::vector<double> > covUpdateBuff_;

    const int resumeCode_;

    std::ofstream out_;
    int nChains_, currentChainI_;
    double burnin_;

    bool stop_;
    int stopRequestMessage_;
    int stopRequestTag_;
    bool stopRequestSent_;

    void* receiveStopRequest_;
    void* sendStopRequest_;

    bool stopMessageRequested_;
    int stopMessageBuff_;

    int haveStoppedMessage_;
    std::vector<int> haveStoppedBuff_;
    int haveStoppedMessageTag_;

    void* haveStoppedMesReq_;
    std::vector<void*> haveStoppedReceiveReq_;

    int updateReqTag_;
    std::vector<void*> updateRequests_;
    bool firstUpdateRequested_;
    std::vector<void*> updateReceiveReq_;

    std::vector<void*> covarianceUpdateRequests_;
    int covUpdateReqTag_;

    bool firstCovUpdateRequested_;
    void* receiveCovUpdateRequest_;
    std::vector<double> eigenUpdateBuff_;

    struct CovarianceMatrixUpdateInfo
    {
        CovarianceMatrixUpdateInfo(int dim)
        {
            check(dim > 0, "");

            n = 0;
            paramSum.resize(dim, 0);
            matrixSum.resize(dim);
            for(int i = 0; i < dim; ++i)
                matrixSum[i].resize(dim, 0);
        }

        unsigned long n;
        std::vector<double> paramSum;
        std::vector<std::vector<double> > matrixSum;

        inline void writeIntoFile(std::ofstream& out) const
        {
            const int dim = paramSum.size();
            check(matrixSum.size() == dim, "");
            out.write((char*)(&n), sizeof(n));
            out.write((char*)(&(paramSum[0])), dim * sizeof(double));
            for(int i = 0; i < dim; ++i)
            {
                check(matrixSum[i].size() == dim, "");
                out.write((char*)(&(matrixSum[i][0])), dim * sizeof(double));
            }
        }

        inline void readFromFile(std::ifstream& in)
        {
            const int dim = paramSum.size();
            check(matrixSum.size() == dim, "");
            in.read((char*)(&n), sizeof(n));
            in.read((char*)(&(paramSum[0])), dim * sizeof(double));
            for(int i = 0; i < dim; ++i)
            {
                check(matrixSum[i].size() == dim, "");
                in.read((char*)(&(matrixSum[i][0])), dim * sizeof(double));
            }
        }

        inline void flush()
        {
            n = 0;
            const int dim = paramSum.size();
            check(matrixSum.size() == dim, "");
            for(int i = 0; i < dim; ++i)
            {
                paramSum[i] = 0;
                check(matrixSum[i].size() == dim, "");
                for(int j = 0; j < dim; ++j)
                    matrixSum[i][j] = 0;
            }
        }
    };

    inline void updateCovarianceMatrix(const CovarianceMatrixUpdateInfo& info);

    CovarianceMatrixUpdateInfo myCovUpdateInfo_;
    CovarianceMatrixUpdateInfo tempCovUpdateInfo_;

    bool likelihoodApproximate_;
};

double
MetropolisHastings::uniformPrior(double min, double max, double x) const
{
    check(max > min, "");
    if(x >= min && x <= max)
        return 1.0 / (max - min);

    return 0.0;
}

double
MetropolisHastings::gaussPrior(double mean, double sigma, double x) const
{
    check(sigma > 0, "");
    const double norm = 1.0 / (std::sqrt(2 * Math::pi) * sigma);
    return norm * std::exp(-(x - mean) * (x - mean) / (2 * sigma * sigma));
}

double
MetropolisHastings::calculatePrior()
{
    if(externalPrior_)
        return externalPrior_->calculate(&(current_[0]), n_);

    double result = 1.0;
    for(int i = 0; i < n_; ++i)
    {
        switch(priorMods_[i])
        {
        case UNIFORM_PRIOR:
            result *= uniformPrior(param1_[i], param2_[i], current_[i]);
            break;

        case GAUSSIAN_PRIOR:
            result *= gaussPrior(param1_[i], param2_[i], current_[i]);
            break;

        default:
            check(false, "invalid prior mode");
            break;
        }
    }

    return result;
}

void
MetropolisHastings::logProgress() const
{
    check(isMaster(), "");

    output_log("MCMC iteration " << iteration_ << ":" << std::endl);
    switch(cd_)
    {
    case GELMAN_RUBIN:
        for(int i = 0; i < n_; ++i)
        {
            output_log("MCMC parameter " << i << ": Gelman-Rubin estimated potential scale reduction = " << rGelmanRubin_[i] << std::endl);
        }
        break;
    case ACCURACY:
        for(int i = 0; i < n_; ++i)
        {
            output_log("MCMC parameter " << i << ": reached accuracy = " << reachedSigma_[i] << " , expected accuracy = " << accuracy_[i] << std::endl);
        }
        break;
    default:
        check(false ,"");
        break;
    }
}

bool
MetropolisHastings::stop()
{
    check(iteration_ >= 0, "");

    if(!isMaster())
        return stop_;

    if(iteration_ < burnin_ + 100)
        return false;

    if(iteration_ >= maxChainLength_)
        return true;

    return checkStoppingCrit();
}

bool
MetropolisHastings::checkStoppingCrit()
{
    check(isMaster(), "");

    unsigned long chainSize, firstHalfEnd, secondHalfStart;

    bool doStop = true;
    double s, x, mean1, mean2, var1, var2;
    double total, totalMean;
    double B, W;
    
    unsigned long minChainNum;

    std::vector<double> means(nChains_);

    const bool synchronized = synchronizeCommInfo();
    std::list<CommunicationInfo>::const_iterator firstIt;

    const int nParamsToCheck = n_;
        
    switch(cd_)
    {
    case GELMAN_RUBIN:
        check(nChains_ > 1, "");
        if(!synchronized)
        {
            for(int i = 0; i < n_; ++i)
                rGelmanRubin_[i] = 100;

            doStop = false;
            break;
        }

        check(!(commInfo_[0].empty()), "");
        firstIt = commInfo_[0].begin();
        total = (*firstIt).iter;
#ifdef CHECKS_ON
        for(int i = 0; i < nChains_; ++i)
        {
            firstIt = commInfo_[i].begin();
            check((*firstIt).iter == total, "");
        }
#endif

        if(total < 100)
            return false;

        for(int i = 0; i < nParamsToCheck; ++i)
        {
            totalMean = 0;
            for(int j = 0; j < nChains_; ++j)
            {
                const CommunicationInfo& ci = *(commInfo_[j].begin());
                means[j] = ci.sums[i] / total;
                totalMean += means[j];
            }
            totalMean /= nChains_;
            B = 0;
            W = 0;
            for(int j = 0; j < nChains_; ++j)
            {
                const CommunicationInfo& ci = *(commInfo_[j].begin());
                const double diff = means[j] - totalMean;
                B += diff * diff;
                const double s2 = (ci.sqSums[i] - 2 * means[j] * ci.sums[i] + total * means[j] * means[j]) / (total - 1);
                check(s2 >= 0, "i = " << i << ", j = " << j << ", total = " << total << ", squared sum = " << ci.sqSums[i] << ", sum = " << ci.sums[i] << ", mean = " << means[j]);
                W += s2;
            }
            B *= total;
            B /= (nChains_ - 1);
            W /= nChains_;

            check(B >= 0, "");
            check(W >= 0, "");

            const double var = (total - 1) * W / total + B / total;
            check(var >= 0, "");

            if(W == 0)
                rGelmanRubin_[i] = (var == 0 ? 1.0 : 100.0);
            else
                rGelmanRubin_[i] = std::sqrt(var / W);

            if(std::abs(rGelmanRubin_[i] - 1) > cc_)
                doStop = false;
        }
        break;
    case ACCURACY:
        for(int i = 0; i < nParamsToCheck; ++i)
        {
            s = 0;
            for(int j = 0; j < nChains_; ++j)
            {
                const CommunicationInfo& ci = *(commInfo_[j].begin());
                x = ci.stdMean[i];
                if(x == -1)
                    return false;

                s += x * x;
            }

            reachedSigma_[i] = std::sqrt(s) / nChains_;
            if(reachedSigma_[i] > accuracy_[i])
                doStop = false;
        }
        break;
    default:
        check(false, "");
        break;
    }
    return doStop;
}

void
MetropolisHastings::calculateMeanVar(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end, double& mean, double& var)
{
    mean = 0;
    double meanSq = 0;

    unsigned long count = 0;
    for(std::vector<double>::const_iterator it = begin; it != end; ++it)
    {
        ++count;
        const double& val = (*it);
        mean += val;
        meanSq += val * val;
    }

    check(count > 0, "");

    mean /= count;
    meanSq /= count;

    var = meanSq - mean;
}

void
MetropolisHastings::calculateStoppingData()
{
    check(iteration_ > burnin_, "");

    for(int i = 0; i < n_; ++i)
    {
        const double mean = paramSum_[i] / (iteration_ - burnin_);
        const double meanSq = paramSquaredSum_[i] / (iteration_ - burnin_);

        const double stdev = std::sqrt(meanSq - mean * mean);
        double stdMean = stdev / std::sqrt(double(iteration_ - burnin_));

        const double cor = (corSum_[i] / (iteration_ - burnin_) - mean * mean) / (stdev * stdev);
        if(cor < 1 && cor > -1)
            stdMean *= std::sqrt((1 + cor) / (1 - cor));
        myStdMean_[i] = stdMean;
    }
}

void
MetropolisHastings::openOut(bool append)
{
    std::stringstream fileName;
    fileName << fileRoot_;
    if(nChains_ > 1)
        fileName << '_' << currentChainI_;
    fileName << ".txt";

    if(append)
        out_.open(fileName.str().c_str(), std::ios::app);
    else
        out_.open(fileName.str().c_str());
    
    if(!out_)
    {
        StandardException exc;
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into output file " << fileName.str() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
}

void
MetropolisHastings::writeChainElement()
{
    check(out_, "");
    out_ << 1 << "   " << currentLike_;
    for(int i = 0; i < n_; ++i)
        out_ << "   " << current_[i];
    out_ << std::endl;
}

void
MetropolisHastings::update()
{
    if(iteration_ > burnin_)
    {
        for(int i = 0; i < n_; ++i)
        {
            paramSum_[i] += current_[i];
            paramSquaredSum_[i] += current_[i] * current_[i];
            corSum_[i] += current_[i] * prev_[i];
        }
    }

    if(adapt_)
    {
        ++myCovUpdateInfo_.n;
        check(myCovUpdateInfo_.paramSum.size() == n_, "");
        check(myCovUpdateInfo_.matrixSum.size() == n_, "");
        for(int i = 0; i < n_; ++i)
        {
            myCovUpdateInfo_.paramSum[i] += current_[i];
            check(myCovUpdateInfo_.matrixSum[i].size() == n_, "");
            for(int j = 0; j < n_; ++j)
                myCovUpdateInfo_.matrixSum[i][j] += current_[i] * current_[j];
        }
    }

    prev_ = current_;
}

void
MetropolisHastings::updateCovarianceMatrix(const CovarianceMatrixUpdateInfo& info)
{
    check(isMaster(), "");
    check(adapt_, "");

    check(info.n >= 2, "at least 2 elements needed for updating");
    check(info.paramSum.size() == n_, "");
    for(int i = 0; i < n_; ++i)
        paramMeanNew_[i] = double(covarianceElementsNum_) / double(covarianceElementsNum_ + info.n) * paramMean_[i] + info.paramSum[i] / double(covarianceElementsNum_ + info.n);

    check(info.matrixSum.size() == n_, "");
    for(int i = 0; i < n_; ++i)
    {
        check(info.matrixSum[i].size() == n_, "");
        for(int j = i; j < n_; ++j)
        {
            covariance_(i, j) = double(covarianceElementsNum_ - 1) / double(covarianceElementsNum_ + info.n - 1) * covariance_(i, j) + covFactor_ / double(covarianceElementsNum_ + info.n - 1) * (double(covarianceElementsNum_) * paramMean_[i] * paramMean_[j] - double(covarianceElementsNum_ + info.n) * paramMeanNew_[i] * paramMeanNew_[j] + info.matrixSum[i][j] + info.n * (i == j ? covEpsilon_ : 0.0));
        }
    }

    for(int i = 0; i < n_; ++i)
        paramMean_[i] = paramMeanNew_[i];

    covarianceElementsNum_ += info.n;

    for(int i = 0; i < n_; ++i)
    {
        for(int j = 0; j < n_; ++j)
            cholesky_(i, j) = covariance_(i, j);
    }

    cholesky_.choleskyFactorize();

    // to be removed
    covariance_.writeIntoTextFile("mcmc_covariance_matrix.txt");
    cholesky_.writeIntoTextFile("mcmc_cholesky_matrix.txt");

    if(covarianceElementsNum_ > 100)
        covarianceReady_ = true;
}

void
MetropolisHastings::writeCommInfo(std::ofstream& out) const
{
    const int n = commInfo_.size();
    check(n == nChains_, "");

    out.write((char*)(&n), sizeof(n));

    std::vector<int> sizes(n);
    
    for(int i = 0; i < n; ++i)
        sizes[i] = commInfo_[i].size();

    out.write((char*)(&(sizes[0])), n * sizeof(int));
    for(int i = 0; i < n; ++i)
        for(std::list<CommunicationInfo>::const_iterator it = commInfo_[i].begin(); it != commInfo_[i].end(); ++it)
            (*it).writeIntoFile(out);
}

void
MetropolisHastings::readCommInfo(std::ifstream& in)
{
    commInfo_.clear();
    int n;
    in.read((char*)(&n), sizeof(n));
    if(n != nChains_)
    {
        BadResumeInfo bri(n);
        throw bri;
    }
    std::vector<int> sizes(n);
    in.read((char*)(&(sizes[0])), n * sizeof(int));
    commInfo_.resize(n);
    for(int i = 0; i < n; ++i)
    {
        if(sizes[i] < 0 || sizes[i] > 1000000)
        {
            BadResumeInfo bri(sizes[i]);
            throw bri;
        }
    }

    for(int i = 0; i < n; ++i)
        for(int j = 0; j < sizes[i]; ++j)
        {
            CommunicationInfo info(n_);
            commInfo_[i].push_back(info);
            
            std::list<CommunicationInfo>::iterator it = commInfo_[i].end();
            --it;

            (*it).readFromFile(in);
        }
}

void
MetropolisHastings::writeResumeInfo() const
{
    std::ofstream out(resumeFileName_.c_str(), std::ios::binary | std::ios::out);
    if(!out)
        return;

    out.write((char*)(&maxChainLength_), sizeof(unsigned long));
    out.write((char*)(&iteration_), sizeof(unsigned long));
    out.write((char*)(&currentLike_), sizeof(double));
    out.write((char*)(&currentPrior_), sizeof(double));
    out.write((char*)(&(current_[0])), n_ * sizeof(double));
    out.write((char*)(&(prev_[0])), n_ * sizeof(double));
    out.write((char*)(&(paramSum_[0])), n_ * sizeof(double));
    out.write((char*)(&(paramSquaredSum_[0])), n_ * sizeof(double));
    out.write((char*)(&(corSum_[0])), n_ * sizeof(double));

    if(adapt_)
    {
        myCovUpdateInfo_.writeIntoFile(out);

        out.write((char*)(&covarianceReady_), sizeof(covarianceReady_));

        check(covariance_.size() == n_, "");
        check(cholesky_.size() == n_, "");

        out.write((char*)(&covarianceElementsNum_), sizeof(covarianceElementsNum_));
        for(int i = 0; i < n_; ++i)
        {
            for(int j = 0; j < n_; ++j)
            {
                double x = covariance_(i, j);
                out.write((char*)(&x), sizeof(double));
                x = cholesky_(i, j);
                out.write((char*)(&x), sizeof(double));
            }
        }

        check(paramMean_.size() == n_, "");
        out.write((char*)(&(paramMean_[0])), n_ * sizeof(double));
    }

    if(isMaster())
        writeCommInfo(out);

    out.write((char*)(&resumeCode_), sizeof(int));

    out.close();
}

bool
MetropolisHastings::readResumeInfo()
{
    std::ifstream in(resumeFileName_.c_str(), std::ios::binary | std::ios::in);
    if(!in)
        return false;

    in.read((char*)(&maxChainLength_), sizeof(unsigned long));
    in.read((char*)(&iteration_), sizeof(unsigned long));
    in.read((char*)(&currentLike_), sizeof(double));
    in.read((char*)(&currentPrior_), sizeof(double));
    in.read((char*)(&(current_[0])), n_ * sizeof(double));
    in.read((char*)(&(prev_[0])), n_ * sizeof(double));
    in.read((char*)(&(paramSum_[0])), n_ * sizeof(double));
    in.read((char*)(&(paramSquaredSum_[0])), n_ * sizeof(double));
    in.read((char*)(&(corSum_[0])), n_ * sizeof(double));

    if(adapt_)
    {
        myCovUpdateInfo_.readFromFile(in);

        in.read((char*)(&covarianceReady_), sizeof(covarianceReady_));

        check(covariance_.size() == n_, "");
        check(cholesky_.size() == n_, "");

        in.read((char*)(&covarianceElementsNum_), sizeof(covarianceElementsNum_));
        for(int i = 0; i < n_; ++i)
        {
            for(int j = 0; j < n_; ++j)
            {
                double x;
                in.read((char*)(&x), sizeof(double));
                covariance_(i, j) = x;
                in.read((char*)(&x), sizeof(double));
                cholesky_(i, j) = x;
            }
        }

        check(paramMean_.size() == n_, "");
        in.read((char*)(&(paramMean_[0])), n_ * sizeof(double));
    }

    if(isMaster())
    {
        try{
            readCommInfo(in);
        } catch (BadResumeInfo& bri)
        {
            if(bri.n > 0)
            {
                output_screen1("Seems like the resume info is for a different number of chains. Currently running " << nChains_ << " resume file indicates " << bri.n << std::endl);
            }
            else
            {
                output_screen1("Problem in the resume file, read a size of" << bri.n << std::endl);
            }
            return false;
        }
    }

    int code = 0;

    in.read((char*)(&code), sizeof(int));

    in.close();

    if(code != resumeCode_)
    {
        output_screen("Resume file is corrupt or not complete!" << std::endl);
        return false;
    }
    return true;
}
} // namespace Math

#endif

