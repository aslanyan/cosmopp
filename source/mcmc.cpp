#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <macros.hpp>
#include <exception_handler.hpp>
#include <mcmc.hpp>

namespace Math
{

MetropolisHastings::MetropolisHastings(int nPar, LikelihoodFunction& like, std::string fileRoot, time_t seed) : n_(nPar), like_(like), fileRoot_(fileRoot), paramNames_(nPar), param1_(nPar, 0), param2_(nPar, 0), starting_(nPar, std::numeric_limits<double>::max()), current_(nPar), prev_(nPar), samplingWidth_(nPar, 0), accuracy_(nPar, 0), paramSum_(nPar, 0), paramSquaredSum_(nPar, 0), corSum_(nPar, 0), priorMods_(nPar, PRIOR_MODE_MAX), externalPrior_(NULL), externalProposal_(NULL), resumeCode_(123456), nChains_(1), currentChainI_(0), stop_(false), stopRequestMessage_(111222), stopRequestTag_(0), stopRequestSent_(false), stopMessageRequested_(false), haveStoppedMessage_(476901), haveStoppedMessageTag_(100000), updateReqTag_(200000), firstUpdateRequested_(false)
{
#ifdef COSMO_MPI
    int hasMpiInitialized;
    MPI_Initialized(&hasMpiInitialized);
    if(!hasMpiInitialized)
        MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nChains_);
    MPI_Comm_rank(MPI_COMM_WORLD, &currentChainI_);

    sendStopRequest_ = new MPI_Request;
    receiveStopRequest_ = new MPI_Request;
    haveStoppedMesReq_ = new MPI_Request;

    if(isMaster())
    {
        haveStoppedBuff_.resize(nChains_);
        haveStoppedReceiveReq_.resize(nChains_);
        updateReceiveReq_.resize(nChains_);
        for(int i = 0; i < nChains_; ++i)
        {
            haveStoppedReceiveReq_[i] = new MPI_Request;
            updateReceiveReq_[i] = new MPI_Request;
        }
    }
#endif

    stdMean_.resize(nChains_);
    stdMeanBuff_.resize(nChains_);
    for(int i = 0; i < nChains_; ++i)
    {
        stdMean_[i].resize(n_, -1);
        stdMeanBuff_[i].resize(n_, -1);
    }

    myStdMean_.resize(n_, -1);

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

#ifdef COSMO_MPI
    delete (MPI_Request*) sendStopRequest_;
    delete (MPI_Request*) receiveStopRequest_;

    if(!isMaster())
    {
        MPI_Status st;
        MPI_Wait((MPI_Request*) haveStoppedMesReq_, &st);

        for(int i = 0; i < updateRequests_.size(); ++i)
        {
            MPI_Status updateSt;
            MPI_Wait((MPI_Request*) updateRequests_[i], &updateSt);
            delete (MPI_Request*) updateRequests_[i];
        }
    }
    delete (MPI_Request*) haveStoppedMesReq_;

    if(isMaster())
    {
        for(int i = 0; i < nChains_; ++i)
        {
            delete (MPI_Request*) haveStoppedReceiveReq_[i];
            delete (MPI_Request*) updateReceiveReq_[i];
        }
    }
#endif
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

void
MetropolisHastings::communicate()
{
    if(iteration_ >= burnin_ + 100)
        calculateStoppingData();

    if(isMaster())
    {
        for(int i = 0; i < n_; ++i)
            stdMean_[0][i] = myStdMean_[i];
    }

#ifdef COSMO_MPI
    if(!isMaster() && !stop_)
    {
        output_screen("CHAIN " << currentChainI_ << ": sending updates about progress to master." << std::endl);
        MPI_Request* updateReq = new MPI_Request;
        updateRequests_.push_back(updateReq);
        MPI_Isend(&(myStdMean_[0]), n_, MPI_DOUBLE, 0, updateReqTag_ + currentChainI_, MPI_COMM_WORLD, updateReq);
    }

    if(isMaster())
    {
        if(!firstUpdateRequested_)
        {
            for(int i = 1; i < nChains_; ++i)
            {
                MPI_Irecv(&(stdMeanBuff_[i][0]), n_, MPI_DOUBLE, i, updateReqTag_ + i, MPI_COMM_WORLD, (MPI_Request*) updateReceiveReq_[i]);
            }

            firstUpdateRequested_ = true;
        }

        if(!stopRequestSent_)
        {
            for(int i = 1; i < nChains_; ++i)
            {
                int updateFlag = 0;
                MPI_Status updateSt;
                MPI_Test((MPI_Request*) updateReceiveReq_[i], &updateFlag, &updateSt);
                if(updateFlag)
                {
                    output_screen("CHAIN " << currentChainI_ << ": Received an update from chain " << i << "." << std::endl);
                    for(int j = 0; j < n_; ++j)
                        stdMean_[i][j] = stdMeanBuff_[i][j];

                    MPI_Irecv(&(stdMeanBuff_[i][0]), n_, MPI_DOUBLE, i, updateReqTag_ + i, MPI_COMM_WORLD, (MPI_Request*) updateReceiveReq_[i]);
                }
            }
        }
    }

    if(isMaster())
    {
        if(stop_ && !stopRequestSent_)
        {
            for(int i = 1; i < nChains_; ++i)
            {
                output_screen("CHAIN " << currentChainI_ << ": Sending stop request to chain " << i << "." << std::endl);
                MPI_Isend(&stopRequestMessage_, 1, MPI_INT, i, stopRequestTag_ + i, MPI_COMM_WORLD, (MPI_Request*) sendStopRequest_);

                MPI_Irecv(&(haveStoppedBuff_[i]), 1, MPI_INT, i, haveStoppedMessageTag_ + i, MPI_COMM_WORLD, (MPI_Request*) haveStoppedReceiveReq_[i]);
            }

            stopRequestSent_ = true;
        }

        return;
    }

    if(!stopMessageRequested_)
    {
        stopMessageBuff_ = 0;
        MPI_Irecv(&stopMessageBuff_, 1, MPI_INT, 0, stopRequestTag_ + currentChainI_, MPI_COMM_WORLD, (MPI_Request*) receiveStopRequest_);
        stopMessageRequested_ = true;
    }

    if(!stop_)
    {
        int stopFlag = 0;
        MPI_Status stopSt;
        MPI_Test((MPI_Request*) receiveStopRequest_, &stopFlag, &stopSt);
        if(stopFlag)
        {
            check(stopMessageBuff_ == stopRequestMessage_, "wrong message received");
            output_screen("CHAIN " << currentChainI_ << ": Received stop request." << std::endl);
            stop_ = true;
        }
    }
#endif
}

void
MetropolisHastings::sendHaveStopped()
{
#ifdef COSMO_MPI
    check(!isMaster(), "");

    output_screen("CHAIN " << currentChainI_ << ": Informing master that I have stopped." << std::endl);
    MPI_Isend(&haveStoppedMessage_, 1, MPI_INT, 0, haveStoppedMessageTag_ + currentChainI_, MPI_COMM_WORLD, (MPI_Request*) haveStoppedMesReq_);
#endif
}

int
MetropolisHastings::run(unsigned long maxChainLength, int writeResumeInformationEvery, unsigned long burnin)
{
    check(maxChainLength > 0, "invalid maxChainLength = " << maxChainLength);
    check(!blocks_.empty(), "");
    
    burnin_ = burnin;

#ifdef COSMO_MPI
    if(currentChainI_ == 0)
    {
        output_screen("Running the MPI version of MetropolisHastings with " << nChains_ << " tasks!!!" << std::endl << std::endl);
    }
#endif

    StandardException exc;

    // Creating the paramnames file
    if(isMaster())
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

            const double newPrior = calculatePrior();
            const double oldLike = currentLike_;
            if(newPrior != 0)
                currentLike_ = like_.calculate(&(current_[0]), n_);

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

        if(iteration_ % 100)
            communicate();

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

    if(isMaster())
        stop_ = true;

    communicate();

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

    if(!isMaster())
        sendHaveStopped();
    else
    {
#ifdef COSMO_MPI
        for(int i = 1; i < nChains_; ++i)
        {
            MPI_Status st;
            MPI_Wait((MPI_Request*) haveStoppedReceiveReq_[i], &st);
            int f;
            MPI_Test((MPI_Request*) haveStoppedReceiveReq_[i], &f, &st);
            check(f, "");
            check(haveStoppedBuff_[i] == haveStoppedMessage_, "wrong message received");
            output_screen("CHAIN " << currentChainI_ << ": heard from chain " << i << " that it has stopped." << std::endl);

            //MPI_Wait((MPI_Request*) updateReceiveReq_[i], &st);
        }
#endif
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
