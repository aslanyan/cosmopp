#ifdef COSMO_MPI
#include <mpi.h>
#endif

#include <blas2pp.h>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <mcmc.hpp>

namespace Math
{

MetropolisHastings::MetropolisHastings(int nPar, LikelihoodFunction& like, std::string fileRoot, time_t seed) : n_(nPar), like_(&like), spareLike_(NULL), fileRoot_(fileRoot), paramNames_(nPar), param1_(nPar, 0), param2_(nPar, 0), starting_(nPar, std::numeric_limits<double>::max()), current_(nPar), prev_(nPar), samplingWidth_(nPar, 0), accuracy_(nPar, 0), paramSum_(nPar, 0), paramSquaredSum_(nPar, 0), corSum_(nPar, 0), priorMods_(nPar, PRIOR_MODE_MAX), externalPrior_(NULL), externalProposal_(NULL), resumeCode_(123456), nChains_(1), currentChainI_(0), stop_(false), stopRequestMessage_(111222), stopRequestTag_(0), stopRequestSent_(false), stopMessageRequested_(false), haveStoppedMessage_(476901), haveStoppedMessageTag_(100000), updateReqTag_(200000), covUpdateReqTag_(300000), firstUpdateRequested_(false), reachedSigma_(nPar, -1), rGelmanRubin_(nPar, -1), adapt_(false), covEpsilon_(1e-7), covFactor_(2.4 * 2.4 / nPar), myCovUpdateInfo_(nPar), tempCovUpdateInfo_(nPar), covarianceReady_(false), firstCovUpdateRequested_(false)
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

    receiveCovUpdateRequest_ = new MPI_Request;

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

    commInfo_.resize(nChains_);

    communicationBuff_.resize(nChains_);
    for(int i = 0; i < nChains_; ++i)
    {
        communicationBuff_[i].resize(1 + 3 * n_, -1);
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

    delete (MPI_Request*) receiveCovUpdateRequest_;

    if(!isMaster())
    {
        //MPI_Status st;
        //MPI_Wait((MPI_Request*) haveStoppedMesReq_, &st);

        for(int i = 0; i < updateRequests_.size(); ++i)
        {
            //MPI_Status updateSt;
            //MPI_Wait((MPI_Request*) updateRequests_[i], &updateSt);
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

        for(int i = 0; i < covarianceUpdateRequests_.size(); ++i)
            delete (MPI_Request*) covarianceUpdateRequests_[i];
    }
#endif
}

void
MetropolisHastings::useAdaptiveProposal()
{
    adapt_ = true;
    covarianceElementsNum_ = 0;
    covariance_.resize(n_, n_);

    paramMean_.resize(n_, 0);
    paramMeanNew_.resize(n_);
    generatedVec_.resize(n_);
    rotatedVec_.resize(n_);

    choleskyMat_.resize(n_, n_);


    cholesky_.resize(n_, 2 * n_ + 1);

    for(int i = 0; i < n_; ++i)
        for(int j = 0; j < n_; ++j)
        {
            covariance_(i, j) = 0;
            choleskyMat_(i, j) = 0;
        }

    for(int i = 0; i < nChains_; ++i)
    {
        communicationBuff_[i].resize(n_ * n_ + n_ + 1 + 1 + 3 * n_, -1);
    }
}

void
MetropolisHastings::setParam(int i, const std::string& name, double min, double max, double starting, double startingWidth, double samplingWidth, double accuracy)
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

    check(startingWidth >= 0 && startingWidth <= (max - min), "invalid starting width " << startingWidth);
    if(startingWidth == 0)
        startingWidth = (max - min) / 100;

    double thisStarting = starting_[i];
    do
    {
        thisStarting = starting_[i] + startingWidth * generator_->generate();
    }
    while(thisStarting < min || thisStarting > max);

    starting_[i] = thisStarting;

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
MetropolisHastings::setParamGauss(int i, const std::string& name, double mean, double sigma, double starting, double startingWidth, double samplingWidth, double accuracy)
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

    check(startingWidth >= 0, "invalid starting width " << startingWidth);
    
    starting_[i] += startingWidth * generator_->generate();

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
        CommunicationInfo info;
        commInfo_[0].push_back(info);
        std::list<CommunicationInfo>::iterator it = commInfo_[0].end();
        --it;
        (*it).sums = paramSum_;
        (*it).sqSums = paramSquaredSum_;
        (*it).stdMean = myStdMean_;
        (*it).iter = double(iteration_ - burnin_);

        if(adapt_ && !stop_)
        {
            updateCovarianceMatrix(myCovUpdateInfo_);
            myCovUpdateInfo_.flush();
        }
    }

    int covarianceUpdateSize = 0;
    if(adapt_ && !stop_)
        covarianceUpdateSize = n_ * n_ + n_ + 1;

#ifdef COSMO_MPI
    if(!isMaster() && !stop_)
    {
        output_screen1("Sending updates about progress to master." << std::endl);
        MPI_Request* updateReq = new MPI_Request;
        updateRequests_.push_back(updateReq);
        sendComBuff_.resize(sendComBuff_.size() + 1); 
        std::vector<double>& currentCom = sendComBuff_[sendComBuff_.size() - 1];

        currentCom.resize(covarianceUpdateSize + 1 + 3 * n_);

        if(adapt_ && !stop_)
        {
            check(myCovUpdateInfo_.matrixSum.size() == n_, "");
            for(int i = 0; i < n_; ++i)
            {
                check(myCovUpdateInfo_.matrixSum[i].size() == n_, "");
                for(int j = 0; j < n_; ++j)
                    currentCom[i * n_ + j] = myCovUpdateInfo_.matrixSum[i][j];
            }

            check(myCovUpdateInfo_.paramSum.size() == n_, "");
            for(int i = 0; i < n_; ++i)
                currentCom[n_ * n_ + i] = myCovUpdateInfo_.paramSum[i];

            currentCom[n_ * n_ + n_] = double(myCovUpdateInfo_.n);

            myCovUpdateInfo_.flush();
        }
        for(int i = 0; i < n_; ++i)
        {
            currentCom[covarianceUpdateSize + i] = paramSum_[i];
            currentCom[covarianceUpdateSize + n_ + i] = paramSquaredSum_[i];
            currentCom[covarianceUpdateSize + 2 * n_ + i] = myStdMean_[i];
        }
        currentCom[covarianceUpdateSize + 3 * n_] = double(iteration_ - burnin_);
        MPI_Isend(&(currentCom[0]), covarianceUpdateSize + 1 + 3 * n_, MPI_DOUBLE, 0, updateReqTag_ + currentChainI_, MPI_COMM_WORLD, updateReq);
    }

    if(isMaster())
    {
        if(!firstUpdateRequested_)
        {
            for(int i = 1; i < nChains_; ++i)
            {
                check(communicationBuff_[i].size() == covarianceUpdateSize + 1 + 3 * n_, "");
                MPI_Irecv(&(communicationBuff_[i][0]), covarianceUpdateSize + 1 + 3 * n_, MPI_DOUBLE, i, updateReqTag_ + i, MPI_COMM_WORLD, (MPI_Request*) updateReceiveReq_[i]);
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
                    output_screen1("Received an update from chain " << i << "." << std::endl);
                    CommunicationInfo temp(n_);
                    commInfo_[i].push_back(temp);
                    std::list<CommunicationInfo>::iterator it = commInfo_[i].end();
                    --it;
                    for(int j = 0; j < n_; ++j)
                    {
                        (*it).sums[j] = communicationBuff_[i][covarianceUpdateSize + j];
                        (*it).sqSums[j] = communicationBuff_[i][covarianceUpdateSize + n_ + j];
                        (*it).stdMean[j] = communicationBuff_[i][covarianceUpdateSize + 2 * n_ + j];
                    }
                    (*it).iter = communicationBuff_[i][covarianceUpdateSize + 3 * n_];

                    if(adapt_ && !stop_)
                    {
                        check(tempCovUpdateInfo_.matrixSum.size() == n_, "");
                        for(int k = 0; k < n_; ++k)
                        {
                            check(tempCovUpdateInfo_.matrixSum[k].size() == n_, "");
                            for(int j = 0; j < n_; ++j)
                                tempCovUpdateInfo_.matrixSum[k][j] = communicationBuff_[i][k * n_ + j];
                        }

                        check(tempCovUpdateInfo_.paramSum.size() == n_, "");
                        for(int j = 0; j < n_; ++j)
                            tempCovUpdateInfo_.paramSum[j] = communicationBuff_[i][n_ * n_ + j];

                        tempCovUpdateInfo_.n = (unsigned long)(communicationBuff_[i][n_ * n_ + n_]);

                        updateCovarianceMatrix(tempCovUpdateInfo_);
                    }

                    MPI_Irecv(&(communicationBuff_[i][0]), covarianceUpdateSize + 1 + 3 * n_, MPI_DOUBLE, i, updateReqTag_ + i, MPI_COMM_WORLD, (MPI_Request*) updateReceiveReq_[i]);
                }
            }

            // update everybody's covariance matrix
            if(adapt_ && !stop_ && covarianceReady_)
            {
                for(int i = 1; i < nChains_; ++i)
                {
                    output_screen1("Sending the covariance matrix update to chain " << i << "." << std::endl);

                    MPI_Request* covUpdateReq = new MPI_Request;
                    covarianceUpdateRequests_.push_back(covUpdateReq);
                    covUpdateBuff_.resize(covUpdateBuff_.size() + 1); 
                    std::vector<double>& currentCom = covUpdateBuff_[covUpdateBuff_.size() - 1];

                    currentCom.resize(n_ * n_);
                    for(int j = 0; j < n_; ++j)
                    {
                        for(int k = 0; k < n_; ++k)
                            currentCom[j * n_ + k] = choleskyMat_(j, k);
                    }

                    MPI_Isend(&(currentCom[0]), n_ * n_, MPI_DOUBLE, i, covUpdateReqTag_ + i, MPI_COMM_WORLD, covUpdateReq);
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
                output_screen1("Sending stop request to chain " << i << "." << std::endl);
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

    if(adapt_ && !firstCovUpdateRequested_)
    {
        eigenUpdateBuff_.resize(n_ * n_);
        MPI_Irecv(&(eigenUpdateBuff_[0]), n_ * n_, MPI_DOUBLE, 0, covUpdateReqTag_ + currentChainI_, MPI_COMM_WORLD, (MPI_Request*) receiveCovUpdateRequest_);
        firstCovUpdateRequested_ = true;
    }

    if(adapt_)
    {
        int covUpdateFlag = 0;
        MPI_Status covUpdateSt;
        MPI_Test((MPI_Request*) receiveCovUpdateRequest_, &covUpdateFlag, &covUpdateSt);
        if(covUpdateFlag)
        {
            output_screen1("Received an updated covariance matrix from the master." << std::endl);
            check(choleskyMat_.size(0) == n_ && choleskyMat_.size(1) == n_, "");
            for(int i = 0; i < n_; ++i)
            {
                for(int j = 0; j < n_; ++j)
                    choleskyMat_(i, j) = eigenUpdateBuff_[i * n_ + j];
            }

            covarianceReady_ = true;

            MPI_Irecv(&(eigenUpdateBuff_[0]), n_ * n_, MPI_DOUBLE, 0, covUpdateReqTag_ + currentChainI_, MPI_COMM_WORLD, (MPI_Request*) receiveCovUpdateRequest_);
        }
    }

    if(!stop_)
    {
        int stopFlag = 0;
        MPI_Status stopSt;
        MPI_Test((MPI_Request*) receiveStopRequest_, &stopFlag, &stopSt);
        if(stopFlag)
        {
            check(stopMessageBuff_ == stopRequestMessage_, "wrong message received");
            output_screen1("Received stop request." << std::endl);
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

    output_screen1("Informing master that I have stopped." << std::endl);
    MPI_Isend(&haveStoppedMessage_, 1, MPI_INT, 0, haveStoppedMessageTag_ + currentChainI_, MPI_COMM_WORLD, (MPI_Request*) haveStoppedMesReq_);
#endif
}

int
MetropolisHastings::run(unsigned long maxChainLength, int writeResumeInformationEvery, unsigned long burnin, CONVERGENCE_DIAGNOSTIC cd, double convergenceCriterion, bool adaptiveProposal)
{
    check(maxChainLength > 0, "invalid maxChainLength = " << maxChainLength);
    check(!blocks_.empty(), "");

    check(cd >= 0 && cd < CONVERGENCE_DIAGNOSTIC_MAX, "invalid convergence diagnostic");

    check(convergenceCriterion > 0, "invalid convergence criterion " << convergenceCriterion << ", needs to be positive");

    if(adaptiveProposal)
        useAdaptiveProposal();
    
    burnin_ = burnin;

    cd_ = cd;

    if(nChains_ == 1 && cd_ == GELMAN_RUBIN)
        cd_ = ACCURACY;

    cc_ = convergenceCriterion;

#ifdef COSMO_MPI
    if(currentChainI_ == 0)
    {
        output_screen_clean("Running the MPI version of MetropolisHastings with " << nChains_ << " tasks!!!" << std::endl << std::endl);
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
        currentLike_ = like_->calculate(&(current_[0]), n_);
        currentPrior_ = calculatePrior();
        prev_ = current_;
        iteration_ = 0;

        for(int i = 0; i < n_; ++i)
        {
            paramSum_[i] = 0;
            paramSquaredSum_[i] = 0;
            corSum_[i] = 0;
        }

        commInfo_.clear();
        commInfo_.resize(nChains_);

        openOut(false);
    }

    std::vector<unsigned long> accepted(blocks_.size(), 0);
    unsigned long currentIter = 0;
    while(!stop())
    {
        int blockBegin = 0;
        for(int i = 0; i < blocks_.size(); ++i)
        {
            int blockEnd = blocks_[i];

            std::vector<double> currentOld = current_;

            std::vector<double> block(blockEnd - blockBegin);

            if(adapt_ && covarianceReady_)
            {
                for(int j = 0; j < n_; ++j)
                    generatedVec_(j) = 0;

                for(int j = blockBegin; j < blockEnd; ++j)
                    generatedVec_(j) = generator_->generate();

                Blas_Mat_Vec_Mult(choleskyMat_, generatedVec_, rotatedVec_);

                for(int j = 0; j < n_; ++j)
                    current_[j] += rotatedVec_(j);
            }
            else
            {
                if(externalProposal_)
                    externalProposal_->generate(&(current_[0]), n_, &(block[0]), i);
                else
                {
                    for(int j = blockBegin; j < blockEnd; ++j)
                        block[j - blockBegin] = generateNewPoint(j);
                }

                for(int j = blockBegin; j < blockEnd; ++j)
                    current_[j] = block[j - blockBegin];
            }


            const double newPrior = calculatePrior();
            const double oldLike = currentLike_;
            if(newPrior != 0)
                currentLike_ = like_->calculate(&(current_[0]), n_);

            double p = newPrior / currentPrior_;
            const double deltaLike = currentLike_ - oldLike;
            p *= std::exp(-deltaLike / 2.0);

            if(!(adapt_ && covarianceReady_) && externalProposal_ && !externalProposal_->isSymmetric(i))
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
        ++currentIter;
        update();

        if(iteration_ % 100 == 0)
        {
            communicate();
            if(isMaster())
            {
                logProgress();
            }
        }

        if(writeResumeInformationEvery && iteration_ % writeResumeInformationEvery == 0)
            writeResumeInfo();

        if(iteration_ % 100 == 0)
        {
            closeOut();
            openOut(true);

            output_screen(std::endl);
            output_screen(std::endl);
            output_screen("Total iterations: " << iteration_ << std::endl);
            for(int i = 0; i < accepted.size(); ++i)
            {
                output_screen("Acceptance rate for parameter block " << i << " = " << double(accepted[i]) / double(currentIter) << std::endl);
            }
        }
    }

    if(isMaster())
        stop_ = true;

    communicate();

    closeOut();

    if(isMaster())
    {
        if(iteration_ >= maxChainLength_)
        {
            output_screen("Maximum number of iterations (" << maxChainLength_ << ") reached, stopping!" << std::endl);
        }
        else
        {
            output_screen("The chain has converged to the requested accuracy after " << iteration_ << " iterations, stopping!" << std::endl);
            logProgress();
        }
    }
    else
    {
        output_screen1("Received stop signal from the master, stopping!" << std::endl);
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
            output_screen1("Heard from chain " << i << " that it has stopped." << std::endl);

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

bool
MetropolisHastings::synchronizeCommInfo()
{
    check(commInfo_.size() == nChains_, "");
    std::vector<std::list<CommunicationInfo>::iterator> currentIters(nChains_);

    bool found = false;
    for(std::list<CommunicationInfo>::reverse_iterator ri = commInfo_[0].rbegin(); ri != commInfo_[0].rend(); ++ri)
    {
        bool sync = true;
        const double currentI = (*ri).iter;
        for(int i = 0; i < nChains_; ++i)
        {
            if(commInfo_[i].empty())
            {
                sync = false;
                break;
            }

            bool currentFound = false;
            for(currentIters[i] = commInfo_[i].begin(); currentIters[i] != commInfo_[i].end(); ++(currentIters[i]))
            {
                const double thisI = (*(currentIters[i])).iter;
                if(thisI == currentI)
                {
                    currentFound = true;
                    break;
                }
                if(thisI > currentI)
                    break;
            }
            if(!currentFound)
            {
                sync = false;
                break;
            }
        }
        if(sync)
        {
            found = true;
            break;
        }
    }

    if(!found)
        return false;

    for(int i = 0; i < nChains_; ++i)
        commInfo_[i].erase(commInfo_[i].begin(), currentIters[i]);

    return true;
}

} // namespace Math
