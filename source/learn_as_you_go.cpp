#include <cosmo_mpi.hpp>

#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>

#include <learn_as_you_go.hpp>
#include <exception_handler.hpp>

LearnAsYouGo::LearnAsYouGo(int nPoints, int nData, const Math::RealFunctionMultiToMulti& f, const Math::RealFunctionMultiDim& errorFunc, unsigned long minCount, double precision, const char* fileName) : nPoints_(nPoints), nData_(nData), f_(f), errorFunc_(errorFunc), minCount_(minCount), precision_(precision), fileName_(fileName), gen_(std::time(0), 0, 1), fa_(NULL), fast_(NULL)
{
    check(nPoints_ > 0, "");
    check(nData_ > 0, "");
    check(minCount_ >= 10, "");
    check(precision_ > 0, "");

    if(fileName_ != "")
    {
        updateFile_ = true;
        if(readFromFile(fileName_.c_str()))
            return;
    }

    construct();
}

LearnAsYouGo::~LearnAsYouGo()
{
    if(updateFile_)
        writeIntoFile(fileName_.c_str());

    if(logFile_)
    {
        logFile_->close();
        delete logFile_;
    }

    if(fa_) delete fa_;
    if(fast_) delete fast_;

#ifdef COSMO_MPI
    for(int i = 0; i < updateRequests_.size(); ++i)
    {
        //MPI_Status updateSt;
        //MPI_Wait((MPI_Request*) updateRequests_[i], &updateSt);
        delete (MPI_Request*) updateRequests_[i];
    }

    check(updateReceiveReq_.size() == nProcesses_, "");

    for(int i = 0; i < nProcesses_; ++i)
        delete (MPI_Request*) updateReceiveReq_[i];
#endif
}

void
LearnAsYouGo::construct()
{
    nProcesses_ = CosmoMPI::create().numProcesses();

    check(nProcesses_ >= 1, "");

    processId_ = CosmoMPI::create().processId();

    check(processId_ >= 0 && processId_ < nProcesses_, "");

    pointsCount_ = 0;
    newPointsCount_ = 0;
    updateCount_ = 10;
    updateErrorThreshold_ = minCount_;
    testSize_ = updateErrorThreshold_ / 20;

    newCommunicateCount_ = 0;
    communicateCount_ = 10;

    firstUpdateRequested_ = false;

#ifdef COSMO_MPI
    updateReceiveReq_.resize(nProcesses_);
    for(int i = 0; i < nProcesses_; ++i)
        updateReceiveReq_[i] = new MPI_Request;
#endif

    totalCount_ = 0;
    sameCount_ = 0;
    successfulCount_ = 0;

    logFile_ = NULL;

    // memory allocation
    check(nPoints_ > 0, "");
    check(nData_ > 0, "");

    currentData_.resize(nData_);

    communicateBuff_.resize(communicateCount_ * (nPoints_ + nData_));
    receiveBuff_.resize(nProcesses_);
    for(int i = 0; i < nProcesses_; ++i)
        receiveBuff_[i].resize(communicateCount_ * (nPoints_ + nData_));

    tempParams_.resize(nPoints_);
    tempData_.resize(nData_);

    CosmoMPI::create().barrier();

    communicateTag_ = CosmoMPI::create().getCommTag();
}

bool
LearnAsYouGo::readFromFile(const char* fileName)
{
    std::ifstream in(fileName, std::ios::binary | std::ios::in);
    StandardException exc;
    if(!in)
        return false;

    in.read((char*)(&nPoints_), sizeof(nPoints_));
    in.read((char*)(&nData_), sizeof(nData_));
    in.read((char*)(&precision_), sizeof(precision_));
    in.read((char*)(&minCount_), sizeof(minCount_));
    
    check(nPoints_ > 0, "");
    check(nData_ > 0, "");
    check(precision_ > 0, "");
    check(minCount_ > 10, "");

    construct();

    unsigned long dataSize;

    in.read((char*)(&dataSize), sizeof(dataSize));
    points_.resize(dataSize);
    data_.resize(dataSize);

    for(unsigned long i = 0; i < dataSize; ++i)
    {
        points_[i].resize(nPoints_);
        data_[i].resize(nData_);

        in.read((char*)(&(points_[i][0])), nPoints_ * sizeof(double));
        in.read((char*)(&(data_[i][0])), nData_ * sizeof(double));
    }

    in.close();

    resetPointMap();

    if(points_.size() >= minCount_)
    {
        constructFast();
    }

    return true;
}

void
LearnAsYouGo::writeIntoFile(const char* fileName) const
{
    if(processId_ != 0)
        return;

    std::ofstream out(fileName, std::ios::binary | std::ios::out);
    StandardException exc;
    if(!out)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot writo into output file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }

    check(nPoints_ > 0, "");
    check(nData_ > 0, "");
    check(precision_ > 0, "");
    check(minCount_ > 10, "");

    out.write((char*)(&nPoints_), sizeof(nPoints_));
    out.write((char*)(&nData_), sizeof(nData_));
    out.write((char*)(&precision_), sizeof(precision_));
    out.write((char*)(&minCount_), sizeof(minCount_));

    unsigned long dataSize = points_.size();
    check(data_.size() == dataSize, "");
    check(pointMap_.size() == dataSize, "");

    out.write((char*)(&dataSize), sizeof(dataSize));
    for(unsigned long i = 0; i < dataSize; ++i)
    {
        out.write((char*)(&(points_[i][0])), nPoints_ * sizeof(double));
        out.write((char*)(&(data_[i][0])), nData_ * sizeof(double));
    }

    out.close();
}

void
LearnAsYouGo::resetPointMap()
{
    pointMap_.clear();
    for(unsigned long i = 0; i < points_.size(); ++i)
        pointMap_[points_[i]] = i;
}

void
LearnAsYouGo::logIntoFile(const char* fileNameBase)
{
    check(!logFile_, "log file already initialized");

    std::stringstream fileName;
    fileName << fileNameBase;
    
    if(nProcesses_ > 1)
        fileName << '_' << processId_;

    fileName << ".txt";

    logFile_ = new std::ofstream(fileName.str().c_str());
}

void
LearnAsYouGo::randomizeErrorSet()
{
    check(testSize_ > 0, "");
    check(points_.size() > testSize_, "");

    for(unsigned long i = points_.size() - testSize_; i < points_.size(); ++i)
    {
        const double r = gen_.generate() * points_.size();
        unsigned long index = (unsigned long) std::floor(r);

        check(index >= 0 && index <= points_.size(), "");
        if(index == points_.size())
            index = points_.size() - 1;

        if(index != i)
        {
            points_[i].swap(points_[index]);
            data_[i].swap(data_[index]);
        }
    }

    resetPointMap();
}

void
LearnAsYouGo::setPrecision(double p)
{
    check(p > 0, "invalid precision " << p << ", must be positive");
    precision_ = p;

    if(fast_)
        fast_->setPrecision(precision_);
}

void
LearnAsYouGo::evaluate(const std::vector<double>& x, std::vector<double>* res)
{
    receive();

    check(x.size() == nPoints_, "");

    ++totalCount_;

    const std::map<std::vector<double>, unsigned long, PointComp>::const_iterator it = pointMap_.find(x);

    if(it != pointMap_.end())
    {
        *res = data_[it->second];
        ++sameCount_;
        log();
        return;
    }

    bool good = false;
    if(fast_)
        good = fast_->approximate(x, *res);

    if(!good)
    {
        actual(x, res);
        log();
        return;
    }

    ++successfulCount_;

    log();
}

void
LearnAsYouGo::actual(const std::vector<double>& x, std::vector<double>* res)
{
    f_.evaluate(x, res);

    currentParams_ = x;
    currentData_ = *res;

    communicate();
    addDataPoint(x, *res);
}

void
LearnAsYouGo::log()
{
    if(logFile_)
    {
        (*logFile_) << totalCount_ << '\t' << sameCount_ << '\t' << successfulCount_ << '\t' << totalCount_ - sameCount_ - successfulCount_ << std::endl;
    }
}

void
LearnAsYouGo::addDataPoint(const std::vector<double>& p, const std::vector<double>& d)
{
    check(p.size() == nPoints_, "");
    check(d.size() == nData_, "");
    check(pointMap_.size() == points_.size(), "");

    const std::map<std::vector<double>, unsigned long, PointComp>::const_iterator it = pointMap_.find(p);

    if(it != pointMap_.end())
    {
#ifdef CHECKS_ON
        for(int i = 0; i < nData_; ++i)
        {
            check(d[i] == data_[it->second][i], "");
        }
#endif
        return;
    }

    ++pointsCount_;
    ++newPointsCount_;

    points_.push_back(p);
    data_.push_back(d);
    pointMap_[p] = points_.size() - 1;

    if(!fast_ && points_.size() >= minCount_)
    {
        constructFast();
        return;
    }

    if(fast_ && newPointsCount_ >= updateCount_)
    {
        output_screen1("Updating the fast approximator with " << newPointsCount_ << " points. Total count of points is " << points_.size() << "." << std::endl);
        check(fa_, "");
        check(fast_, "");

        if(points_.size() >= updateErrorThreshold_)
        {
            randomizeErrorSet();
            fa_->reset(points_.size() - testSize_, points_, data_, true);
            fast_->reset(points_, data_, points_.size() - testSize_, points_.size());
            if(processId_ == 0)
                fast_->getPosterior()->writeIntoFile("fast_approximator_error_ratio.txt");
            updateErrorThreshold_ = points_.size() + points_.size() / 5;
            testSize_ = updateErrorThreshold_ / 20;
        }

        fa_->reset(points_.size(), points_, data_, false);
        newPointsCount_ = 0;

        updateCount_ = (points_.size() > 10000 ? points_.size() / 1000 : 10);
        check(updateCount_ >= 10, "");

        output_screen1("Updating the file " << fileName_ << "." << std::endl);
        if(updateFile_)
        {
            output_screen1("Updating the file " << fileName_ << "." << std::endl);
            writeIntoFile(fileName_.c_str());
        }

        return;
    }
}

void
LearnAsYouGo::constructFast()
{
    check(!fa_, "");
    check(!fast_, "");

    output_screen1("Constructing the fast approximator." << std::endl);

    check(nPoints_ > 0, "");
    check(nData_ > 0, "");
    check(!points_.empty(), "");
    check(data_.size() == points_.size(), "");
    check(pointMap_.size() == points_.size(), "");

    const int k = nPoints_ + nPoints_ * nPoints_ + 1;

    check(testSize_ > 0, "");
    if(points_.size() > 2 * testSize_)
    {
        randomizeErrorSet();
        fa_ = new FastApproximator(nPoints_, nData_, points_.size() - testSize_, points_, data_, k);
        fast_ = new FastApproximatorError(*fa_, points_, data_, points_.size() - testSize_, points_.size(), errorFunc_, FastApproximatorError::AVG_INV_DISTANCE, precision_);

        if(processId_ == 0)
            fast_->getPosterior()->writeIntoFile("fast_approximator_error_ratio.txt");

        fa_->reset(points_.size(), points_, data_, false);
        updateErrorThreshold_ = points_.size() + points_.size() / 5;
        testSize_ = updateErrorThreshold_ / 20;
    }
    else
    {
        fa_ = new FastApproximator(nPoints_, nData_, points_.size(), points_, data_, k);
        fast_ = new FastApproximatorError(*fa_, points_, data_, points_.size(), points_.size(), errorFunc_, FastApproximatorError::AVG_INV_DISTANCE, precision_);
    }
}

void
LearnAsYouGo::receive()
{
#ifdef COSMO_MPI
    check(nPoints_ > 0, "");
    check(nData_ > 0, "");

    if(!firstUpdateRequested_)
    {
        for(int i = 0; i < nProcesses_; ++i)
        {
            if(i == processId_)
                continue;

            check(receiveBuff_[i].size() == communicateCount_ * (nPoints_ + nData_), "");
            MPI_Irecv(&(receiveBuff_[i][0]), receiveBuff_[i].size(), MPI_DOUBLE, i, communicateTag_ + i, MPI_COMM_WORLD, (MPI_Request*) updateReceiveReq_[i]);
        }

        firstUpdateRequested_ = true;
    }

    for(int i = 0; i < nProcesses_; ++i)
    {
        if(i == processId_)
            continue;

        int updateFlag = 0;
        MPI_Status updateSt;
        MPI_Test((MPI_Request*) updateReceiveReq_[i], &updateFlag, &updateSt);
        if(updateFlag)
        {
            output_screen1("Received an update from process " << i << "." << std::endl);

            check(receiveBuff_[i].size() == communicateCount_ * (nPoints_ + nData_), "");
            for(int l = 0; l < communicateCount_; ++l)
            {
                check(tempParams_.size() == nPoints_, "");
                for(int j = 0; j < nPoints_; ++j)
                    tempParams_[j] = receiveBuff_[i][l * (nPoints_ + nData_) + j];

                check(tempData_.size() == nData_, "");
                for(int j = 0; j < nData_; ++j)
                    tempData_[j] = receiveBuff_[i][l * (nPoints_ + nData_) + nPoints_ + j];

                addDataPoint(tempParams_, tempData_);
            }

            check(receiveBuff_[i].size() == communicateCount_ * (nPoints_ + nData_), "");
            MPI_Irecv(&(receiveBuff_[i][0]), receiveBuff_[i].size(), MPI_DOUBLE, i, communicateTag_ + i, MPI_COMM_WORLD, (MPI_Request*) updateReceiveReq_[i]);
        }
    }
#endif
}

void
LearnAsYouGo::communicate()
{
#ifdef COSMO_MPI
    if(nProcesses_ == 1)
        return;

    check(communicateBuff_.size() == communicateCount_ * (nPoints_ + nData_), "");

    check(newCommunicateCount_ < communicateCount_, "");
    for(int i = 0; i < nPoints_; ++i)
        communicateBuff_[newCommunicateCount_ * (nPoints_ + nData_) + i] = currentParams_[i];

    check(currentData_.size() == nData_, "");
    for(int i = 0; i < nData_; ++i)
        communicateBuff_[newCommunicateCount_ * (nPoints_ + nData_) + nPoints_ + i] = currentData_[i];

    ++newCommunicateCount_;

    if(newCommunicateCount_ == communicateCount_)
    {
        for(int i = 0; i < nProcesses_; ++i)
        {
            if(i == processId_)
                continue;

            output_screen1("Sending updates to process " << i << "." << std::endl);

            MPI_Request* updateReq = new MPI_Request;
            updateRequests_.push_back(updateReq);
            communicateData_.push_back(communicateBuff_);
            MPI_Isend(&(communicateData_[communicateData_.size() - 1][0]), communicateBuff_.size(), MPI_DOUBLE, i, communicateTag_ + processId_, MPI_COMM_WORLD, updateReq);
        }

        newCommunicateCount_ = 0;
    }
#endif
}

