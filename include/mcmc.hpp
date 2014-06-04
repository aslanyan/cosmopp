#ifndef COSMO_PP_MCMC_HPP
#define COSMO_PP_MCMC_HPP

#include <fstream>
#include <vector>
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

    /// Constructor.
    /// \param nPar The number of parameters.
    /// \param like The likelihood function.
    /// \param fileRoot The root for filenames produced by MetropolisHastings.
    /// \param seed A random seed. If set to 0 (the default value), it will be determined from the current time.
    MetropolisHastings(int nPar, LikelihoodFunction& like, std::string fileRoot, time_t seed = 0);

    /// Destructor.
    ~MetropolisHastings();

    /// Define a given parameter to have a uniform prior. One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param min The minimum value of the parameter (the lower bound for the prior).
    /// \param max The maximum value of the parameter (the upper bound for the prior).
    /// \param starting The starting value of the parameter. If not set, it will be set to the midpoint of the range by default.
    /// \param samplingWidth The sampling width of the parameter (the width of the Gaussian proposal distribution). If not set, by default it will be set to 1/100-th of the width of the range.
    /// \param accuracy The accuracy with which the parameter needs to be determined (used to choose the stopping time). If not set, by default it will be set to 1/10-th of the sampling width.
    void setParam(int i, const std::string& name, double min, double max, double starting = std::numeric_limits<double>::max(), double samplingWidth = 0.0, double accuracy = 0.0);

    /// Define a given parameter to have a gaussian prior. One of the parameter setting functions must be called for each parameter before the run.
    /// \param i The index of the parameter, 0 <= i < number of parameters.
    /// \param name The name of the parameter.
    /// \param mean The mean of the prior
    /// \param sigma The sigma of the prior
    /// \param starting The starting value of the parameter. If not set, it will be set to the midpoint of the range by default.
    /// \param samplingWidth The sampling width of the parameter (the width of the Gaussian proposal distribution). If not set, by default it will be set to 1/100-th of the width of the range.
    /// \param accuracy The accuracy with which the parameter needs to be determined (used to choose the stopping time). If not set, by default it will be set to 1/10-th of the sampling width.
    void setParamGauss(int i, const std::string& name, double mean, double sigma, double starting = std::numeric_limits<double>::max(), double samplingWidth = 0.0, double accuracy = 0.0);

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
    /// \return The number of chains generated.
    int run(unsigned long maxChainLength = 1000000, int writeResumeInformationEvery = 1, unsigned long burnin = 0);

private:
    inline double uniformPrior(double min, double max, double x) const;
    inline double gaussPrior(double mean, double sigma, double x) const;
    inline double calculatePrior();
    inline void calculateStoppingData();
    inline bool stop() const;
    inline bool checkStoppingCrit() const;
    inline double generateNewPoint(int i) const { return current_[i] + generator_->generate() * samplingWidth_[i]; }
    inline void openOut(bool append);
    inline void closeOut() { out_.close(); }
    inline void writeChainElement();
    inline void update();
    inline void writeResumeInfo() const;
    inline bool readResumeInfo();

    inline bool isMaster() const { return currentChainI_ == 0; }
    void communicate();
    void sendHaveStopped();
    
private:
    int n_;
    LikelihoodFunction& like_;
    std::string fileRoot_, resumeFileName_;
    std::vector<double> param1_, param2_, starting_, samplingWidth_, accuracy_;
    std::vector<PRIOR_MODE> priorMods_;
    std::vector<std::string> paramNames_;
    std::vector<double> paramSum_, paramSquaredSum_, corSum_;
    PriorFunctionBase* externalPrior_;
    ProposalFunctionBase* externalProposal_;
    std::vector<int> blocks_;

    time_t seed_;
    Math::UniformRealGenerator* uniformGen_;
    Math::GaussianGenerator* generator_;

    std::vector<double> prev_, current_;
    unsigned long maxChainLength_;
    unsigned long iteration_;
    double currentLike_;
    double currentPrior_;

    // first index is chain index, second index is param index;
    std::vector<std::vector<double> > stdMean_;
    std::vector<std::vector<double> > stdMeanBuff_;
    std::vector<double> myStdMean_;

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

bool
MetropolisHastings::stop() const
{
    check(iteration_ >= 0, "");

    if(iteration_ < burnin_ + 100)
        return false;

    if(iteration_ >= maxChainLength_)
        return true;

    if(!isMaster())
        return stop_;

    return checkStoppingCrit();
}

bool
MetropolisHastings::checkStoppingCrit() const
{
    check(isMaster(), "");
    
    for(int i = 0; i < n_; ++i)
    {
        double s = 0;
        for(int j = 0; j < nChains_; ++j)
        {
            const double x = stdMean_[j][i];
            if(x == -1)
                return false;

            s += x * x;
        }

        const double sigma = std::sqrt(s) / nChains_;
        if(sigma > accuracy_[i])
            return false;
    }

    return true;
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
        exceptionStr << "Cannot write into output file " << fileName << ".";
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
    if(iteration_ >= burnin_)
    {
        for(int i = 0; i < n_; ++i)
        {
            paramSum_[i] += current_[i];
            paramSquaredSum_[i] += current_[i] * current_[i];
            corSum_[i] += current_[i] * prev_[i];
        }
    }
    prev_ = current_;
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

