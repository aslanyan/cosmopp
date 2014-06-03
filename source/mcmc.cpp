#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <macros.hpp>
#include <exception_handler.hpp>
#include <mcmc.hpp>

namespace Math
{

MetropolisHastings::MetropolisHastings(int nPar, LikelihoodFunction& like, std::string fileRoot, time_t seed) : n_(nPar), like_(like), fileRoot_(fileRoot), paramNames_(nPar), param1_(nPar, 0), param2_(nPar, 0), starting_(nPar, std::numeric_limits<double>::max()), current_(nPar), prev_(nPar), samplingWidth_(nPar, 0), accuracy_(nPar, 0), paramSum_(nPar, 0), paramSquaredSum_(nPar, 0), corSum_(nPar, 0), priorMods_(nPar, PRIOR_MODE_MAX), externalPrior_(NULL), externalProposal_(NULL), resumeCode_(123456), nChains_(1), currentChainI_(0)
{
#ifdef COSMO_MPI
    int hasMpiInitialized;
    MPI_Initialized(&hasMpiInitialized);
    if(!hasMpiInitialized)
        MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nChains_);
    MPI_Comm_rank(MPI_COMM_WORLD, &currentChainI_);
#endif

    check(nPar > 0, "");
    for(int i = 1; i <= nPar; ++i)
        blocks_.push_back(i);

    if(seed == 0)
        seed_ = std::time(0);
    else
        seed_ = seed;

    UniformRealGenerator temp(seed_, 0, 1000000);

    for(int i = 0; i < 2 * currentChainI_; ++i)
        temp.generate();

    uniformGen_ = new UniformRealGenerator(int(temp.generate()), 0, 1);
    generator_ = new GaussianGenerator(int(temp.generate()), 0, 1);

    std::stringstream resFileName;
    resFileName << fileRoot_ << "resume";
    if(nChains_ > 1)
        resFileName << '_' << currentChainI_;
    resFileName << ".dat";
    resumeFileName_ = resFileName.str();
}

MetropolisHastings::~MetropolisHastings()
{
    delete uniformGen_;
    delete generator_;
}

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
    {
        check(starting >= min && starting <= max, "invalid starting value " << starting << ", needs to be between " << min << " and " << max);
        starting_[i] = starting;
    }

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

int
MetropolisHastings::run(unsigned long maxChainLength, int writeResumeInformationEvery, unsigned long burnin)
{
    check(maxChainLength > 0, "invalid maxChainLength = " << maxChainLength);
    check(!blocks_.empty(), "");
    
    burnin_ = burnin;

    bool mainProcess = true;

#ifdef COSMO_MPI
    if(currentChainI_ == 0)
    {
        output_screen("Running the MPI version of MetropolisHastings with " << nChains_ << " tasks!!!" << std::endl << std::endl);
    }
    else
        mainProcess = false;
#endif

    StandardException exc;

    // Creating the paramnames file
    if(mainProcess)
    {
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
    }

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
            
            const double q = uniformGen_->generate(); 

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

        if(writeResumeInformationEvery && iteration_ % writeResumeInformationEvery == 0)
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

    return nChains_;
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

} // namespace Math
