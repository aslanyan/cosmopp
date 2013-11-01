#ifndef COSMO_CPP_MCMC_HPP
#define COSMO_CPP_MCMC_HPP

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstdlib>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <math_constants.hpp>
#include <likelihood_function.hpp>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/normal.hpp>

namespace Math
{

class PriorFunctionBase
{
public:
    virtual double calculate(double* params, int nPar) = 0;
};

class ProposalFunctionBase
{
public:
    virtual void generate(double* params, int nPar, double* blockParams, int i) = 0;
    virtual double calculate(double* params, int nPar, double* blockParams, int i) = 0;
    virtual bool isSymmetric(int i) = 0;
};

class MetropolisHastings
{
private:
    enum PRIOR_MODE { UNIFORM_PRIOR = 0, GAUSSIAN_PRIOR, PRIOR_MODE_MAX };
public:
    MetropolisHastings(int nPar, LikelihoodFunction& like, std::string fileRoot, time_t seed = 0) : n_(nPar), like_(like), fileRoot_(fileRoot), paramNames_(nPar), param1_(nPar, 0), param2_(nPar, 0), starting_(nPar, std::numeric_limits<double>::max()), current_(nPar), prev_(nPar), samplingWidth_(nPar, 0), accuracy_(nPar, 0), paramSum_(nPar, 0), paramSquaredSum_(nPar, 0), corSum_(nPar, 0), priorMods_(nPar, PRIOR_MODE_MAX), externalPrior_(NULL), externalProposal_(NULL), resumeCode_(123456)
    {
        check(nPar > 0, "");
        for(int i = 1; i <= nPar; ++i)
            blocks_.push_back(i);

        if(seed == 0)
            seed_ = std::time(0);
        else
            seed_ = seed;

        std::srand((unsigned) seed_);
        generator_ = new boost::variate_generator<boost::mt19937, boost::normal_distribution<> >(boost::mt19937(seed_), boost::normal_distribution<>(0, 1));

        std::stringstream resFileName;
        resFileName << fileRoot_ << "resume.dat";
        resumeFileName_ = resFileName.str();
    }

    ~MetropolisHastings() { delete generator_; }

    inline void setParam(int i, const std::string& name, double min, double max, double starting = std::numeric_limits<double>::max(), double samplingWidth = 0.0, double accuracy = 0.0);
    inline void setParamGauss(int i, const std::string& name, double mean, double sigma, double starting = std::numeric_limits<double>::max(), double samplingWidth = 0.0, double accuracy = 0.0);

    inline void specifyParameterBlocks(const std::vector<int>& blocks);
    void useExternalPrior(PriorFunctionBase* prior) { externalPrior_ = prior; }

    void useExternalProposal(ProposalFunctionBase* proposal) { externalProposal_ = proposal; }
    inline void run(unsigned long maxChainLength = 1000000, bool writeResumeInformation = true);

private:
    double uniformPrior(double min, double max, double x) const
    {
        check(max > min, "");
        if(x >= min && x <= max)
            return 1.0 / (max - min);

        return 0.0;
    }

    double gaussPrior(double mean, double sigma, double x) const
    {
        check(sigma > 0, "");
        const double norm = 1.0 / (std::sqrt(2 * Math::pi) * sigma);
        return norm * std::exp(-(x - mean) * (x - mean) / (2 * sigma * sigma));
    }

    double calculatePrior()
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

    bool stop() const
    {
        check(iteration_ >= 0, "");

        if(iteration_ < 100)
            return false;

        if(iteration_ >= maxChainLength_)
            return true;

        for(int i = 0; i < n_; ++i)
        {
            const double mean = paramSum_[i] / iteration_;
            const double meanSq = paramSquaredSum_[i] / iteration_;

            const double stdev = std::sqrt(meanSq - mean * mean);
            double stdMean = stdev / std::sqrt(double(iteration_));

            const double cor = (corSum_[i] / iteration_ - mean * mean) / (stdev * stdev);
            if(cor < 1 && cor > -1)
                stdMean *= std::sqrt((1 + cor) / (1 - cor));
            if(stdMean > accuracy_[i])
                return false;
        }

        return true;
    }

    double generateNewPoint(int i) const
    {
        return current_[i] + (*generator_)() * samplingWidth_[i];
    }

    void openOut(bool append)
    {
        std::stringstream fileName;
        fileName << fileRoot_ << ".txt";

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

    void closeOut()
    {
        out_.close();
    }

    void writeChainElement()
    {
        check(out_, "");
        out_ << 1 << "   " << currentLike_ / 2.0;
        for(int i = 0; i < n_; ++i)
            out_ << "   " << current_[i];
        out_ << std::endl;
    }

    void update()
    {
        for(int i = 0; i < n_; ++i)
        {
            paramSum_[i] += current_[i];
            paramSquaredSum_[i] += current_[i] * current_[i];
            corSum_[i] += current_[i] * prev_[i];
        }
        prev_ = current_;
    }

    void writeResumeInfo() const
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

    bool readResumeInfo()
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
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> >* generator_;

    std::vector<double> prev_, current_;
    unsigned long maxChainLength_;
    unsigned long iteration_;
    double currentLike_;
    double currentPrior_;

    const int resumeCode_;

    std::ofstream out_;
};

void
MetropolisHastings::setParam(int i, const std::string& name, double min, double max, double starting, double samplingWidth, double accuracy)
{
    check(i >= 0 && i < n_, "invalid index = " << i);
    check(max > min, "max = " << max << ", min = " << min << ". Need max > min.")

    paramNames_[i] = name;
    param1_[i] = min;
    param2_[i] = max;
    priorMods_[i] = UNIFORM_PRIOR;

    if(starting == std::numeric_limits<double>::max())
        starting_[i] = (max + min) / 2.0;
    else
        starting_[i] = starting;

    check(samplingWidth >= 0, "invalid sampling width " << samplingWidth);
    if(samplingWidth == 0.0)
        samplingWidth_[i] = (max - min) / 100;
    else
        samplingWidth_[i] = samplingWidth;

    check(accuracy >= 0, "invalid accuracy = " << accuracy);
    if(accuracy == 0.0)
        accuracy_[i] = samplingWidth_[i] / 10;
    else
        accuracy_[i] = accuracy;
}

void
MetropolisHastings::setParamGauss(int i, const std::string& name, double mean, double sigma, double starting, double samplingWidth, double accuracy)
{
    check(i >= 0 && i < n_, "invalid index = " << i);
    check(sigma > 0, "invalid sigma = " << sigma);

    paramNames_[i] = name;
    param1_[i] = mean;
    param2_[i] = sigma;
    priorMods_[i] = GAUSSIAN_PRIOR;

    if(starting == std::numeric_limits<double>::max())
        starting_[i] = mean;
    else
        starting_[i] = starting;

    check(samplingWidth >= 0, "invalid sampling width " << samplingWidth);
    if(samplingWidth == 0.0)
        samplingWidth_[i] = sigma / 100;
    else
        samplingWidth_[i] = samplingWidth;

    check(accuracy >= 0, "invalid accuracy = " << accuracy);
    if(accuracy == 0.0)
        accuracy_[i] = samplingWidth_[i] / 10;
    else
        accuracy_[i] = accuracy;
}

void
MetropolisHastings::specifyParameterBlocks(const std::vector<int>& blocks)
{
    check(!blocks.empty(), "");
#ifdef CHECKS_ON
    for(int i = 1; i < blocks.size(); ++i)
    {
        check(blocks[i] > blocks[i - 1], "");
        check(blocks[i] <= n_, "");
    }
#endif

    blocks_ = blocks;
}

void
MetropolisHastings::run(unsigned long maxChainLength, bool writeResumeInformation)
{
    check(maxChainLength > 0, "invalid maxChainLength = " << maxChainLength);
    check(!blocks_.empty(), "");

    // Creating the paramnames file
    StandardException exc;
    std::stringstream paramNamesFileName;
    paramNamesFileName << fileRoot_ << ".paramnames";
    std::ofstream outPar(paramNamesFileName.str().c_str());

    if(!outPar)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into paramnames file " << paramNamesFileName.str() << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    for(int i = 0; i < n_; ++i)
    {
        outPar << paramNames_[i] << '\t' << paramNames_[i] << std::endl;
    }
    outPar.close();

    if(readResumeInfo())
    {
        output_screen("Resuming from previous run, already have " << iteration_ << " iterations." << std::endl);
        openOut(true);
    }
    else
    {
        output_screen("No resume file found (or the resume file is not complete), starting from scratch." << std::endl);

        maxChainLength_ = maxChainLength;

        current_ = starting_;
        currentLike_ = like_.calculate(&(current_[0]), n_);
        currentPrior_ = calculatePrior();
        prev_ = current_;
        iteration_ = 0;

        for(int i = 0; i < n_; ++i)
        {
            paramSum_[i] = 0;
            paramSquaredSum_[i] = 0;
            corSum_[i] = 0;
        }

        openOut(false);
    }

    std::vector<unsigned long> accepted(blocks_.size(), 0);
    while(!stop())
    {
        int blockBegin = 0;
        for(int i = 0; i < blocks_.size(); ++i)
        {
            int blockEnd = blocks_[i];

            std::vector<double> block(blockEnd - blockBegin);

            if(externalProposal_)
                externalProposal_->generate(&(current_[0]), n_, &(block[0]), i);
            else
            {
                for(int j = blockBegin; j < blockEnd; ++j)
                    block[j - blockBegin] = generateNewPoint(j);
            }

            std::vector<double> currentOld = current_;
            for(int j = blockBegin; j < blockEnd; ++j)
                current_[j] = block[j - blockBegin];

            const double oldLike = currentLike_;
            currentLike_ = like_.calculate(&(current_[0]), n_);
            const double newPrior = calculatePrior();

            double p = newPrior / currentPrior_;
            const double deltaLike = currentLike_ - oldLike;
            p *= std::exp(-deltaLike / 2.0);

            if(externalProposal_ && !externalProposal_->isSymmetric(i))
            {
                std::vector<double> oldBlock(blockEnd - blockBegin);
                for(int j = blockBegin; j < blockEnd; ++j)
                    oldBlock[j - blockBegin] = currentOld[j];

                p *= externalProposal_->calculate(&(current_[0]), n_, &(oldBlock[0]), i);
                p /= externalProposal_->calculate(&(currentOld[0]), n_, &(block[0]), i);
            }
            if(p > 1)
                p = 1;
            
            const int r = std::rand();
            const double q = ((double)r) / RAND_MAX;

            if(q <= p)
            {
                currentPrior_ = newPrior;
                ++accepted[i];
            }
            else
            {
                current_ = currentOld;
                currentLike_ = oldLike;
            }

            blockBegin = blockEnd;
        }

        writeChainElement();
        ++iteration_;
        update();

        if(writeResumeInformation)
            writeResumeInfo();

        if(iteration_ % 1000 == 0)
        {
            closeOut();
            openOut(true);

            output_screen(std::endl << std::endl << "Total iterations: " << iteration_ << std::endl);
            for(int i = 0; i < blocks_.size(); ++i)
            {
                output_screen("Acceptance rate for parameter block " << i << " = " << double(accepted[i]) / double(iteration_) << std::endl);
            }
        }
    }

    closeOut();

    if(iteration_ >= maxChainLength_)
    {
        output_screen("Maximum number of iterations (" << maxChainLength_ << ") reached, stopping!" << std::endl);
    }
    else
    {
        output_screen("The chain has converged to the requested accuracy after " << iteration_ << " iterations, stopping!" << std::endl);
    }

    for(int i = 0; i < blocks_.size(); ++i)
    {
        output_screen("Acceptance rate for parameter block " << i << " = " << double(accepted[i]) / double(iteration_) << std::endl);
    }
}

} // namespace Math

#endif

