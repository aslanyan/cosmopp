#include <cosmo_mpi.hpp>

#include <memory>

#include <macros.hpp>
#include <progress_meter.hpp>
#include <pde.hpp>

namespace Math
{

InitialValPDESolver::InitialValPDESolver(InitialValPDEInterface *pde) : pde_(pde), d_(pde->spaceDim()), m_(pde->funcDim()), grid_(0), halfGrid_(0)
{ 
    check(d_ >= 1, "invalid number of spatial dimensions " << d_ << ", must be positive");
    check(m_ >= 1, "invalid number of fields " << m_ << ", must be positive");

    twoD_ = 1;
    for(int i = 0; i < d_; ++i)
        twoD_ *= 2;

    nProcesses_ = CosmoMPI::create().numProcesses();
    processId_ = CosmoMPI::create().processId();
    check(processId_ >= 0 && processId_ < nProcesses_, "");

    commTag_ = CosmoMPI::create().getCommTag();
}

InitialValPDESolver::~InitialValPDESolver()
{
}

void
InitialValPDESolver::set(const std::vector<RealFunctionMultiDim*>& w0, const std::vector<double>& xMin, const std::vector<double>& xMax, const std::vector<int>& nx)
{
    check(w0.size() == m_, "");
    check(xMin.size() == d_, "");
    check(xMax.size() == d_, "");
    check(nx.size() == d_, "");

    for(int i = 0; i < d_; ++i)
    {
        check(xMin[i] < xMax[i], "");
        check(nx[i] >= 1, "");
    }

    xMin_ = xMin;
    xMax_ = xMax;
    nx_ = nx;

    deltaX_.resize(d_);
    double minDeltaX = (xMax_[0] - xMin_[0]);
    for(int i = 0; i < d_; ++i)
    {
        deltaX_[i] = (xMax_[i] - xMin_[i]) / nx_[i];
        if(minDeltaX > deltaX_[i])
            minDeltaX = deltaX_[i];
    }

    deltaT_ = 0.8 * minDeltaX;

    setupGrid();
    setInitial(w0);
    setOwnBoundary();
    communicateBoundary();
}

void
InitialValPDESolver::setupGrid()
{
    std::vector<int> nx0PerProcess(nProcesses_), nx0Before(nProcesses_);
    int totalNX0 = nx_[0];
    int totalNX0Before = 0;

    for(int i = 0; i < nProcesses_; ++i)
    {
        nx0PerProcess[i] = totalNX0 / (nProcesses_ - i);
        totalNX0 -= nx0PerProcess[i];
        nx0Before[i] = totalNX0Before;
        totalNX0Before += nx0PerProcess[i];
    }

    check(totalNX0Before == nx_[0], "");
    
    nx0Starting_ = nx0Before[processId_];
    nx_[0] = nx0PerProcess[processId_];

    dimProd_.resize(d_ + 1);
    dimProd_[d_] = 1;
    halfDimProd_.resize(d_ + 1);
    halfDimProd_[d_] = 1;
    for(int i = d_ - 1; i >= 0; --i)
    {
        unsigned long currentDim = (unsigned long)(nx_[i]) + 2;
        dimProd_[i] = dimProd_[i + 1] * currentDim;
        check(dimProd_[i] >= dimProd_[i + 1], "overflow!");
        check(dimProd_[i] / dimProd_[i + 1] == currentDim, "overflow!");

        currentDim -= 1;
        halfDimProd_[i] = halfDimProd_[i + 1] * currentDim;
        check(halfDimProd_[i] >= halfDimProd_[i + 1], "overflow!");
        check(halfDimProd_[i] / halfDimProd_[i + 1] == currentDim, "overflow!");
    }

    grid_.clear();
    halfGrid_.clear();

    grid_.resize(dimProd_[0]);
    halfGrid_.resize(halfDimProd_[0]);

    for(unsigned long i = 0; i < dimProd_[0]; ++i)
        grid_[i].resize(m_);

    for(unsigned long i = 0; i < halfDimProd_[0]; ++i)
        halfGrid_[i].resize(m_);

    buffer_.resize(dimProd_[1] * m_);
}

void
InitialValPDESolver::setInitial(const std::vector<RealFunctionMultiDim*>& w0)
{
    check(w0.size() == m_, "");
    //check(grid_.size() == m_, "");
    check(grid_.size() == dimProd_[0], "");

    for(int i = 0; i < m_; ++i)
    {
#pragma omp parallel for default(shared)
        for(int i0 = 0; i0 < nx_[0]; ++i0)
        {
            std::vector<int> rangeBegin(d_, 0);
            std::vector<int> rangeEnd = nx_;
            rangeBegin[0] = i0;
            rangeEnd[0] = i0 + 1;

            std::vector<double> coords(d_);

            for(std::vector<int> ind = rangeBegin; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
            {
                physicalCoords(ind, &coords);
                const unsigned long gridIndex = index(ind);
                grid_[gridIndex][i] = w0[i]->evaluate(coords);
            }
        }
    }

    t_ = 0;
}

void
InitialValPDESolver::setOwnBoundary()
{
    for(int i = 1; i < d_; ++i)
    {
#pragma omp parallel for default(shared)
        for(int i0 = 0; i0 < nx_[0]; ++i0)
        {
            std::vector<int> rangeBegin(d_, 0);
            std::vector<int> rangeEnd = nx_;

            rangeBegin[0] = i0;
            rangeEnd[0] = i0 + 1;
            for(int j = 1; j < i; ++j)
            {
                --rangeBegin[j];
                ++rangeEnd[j];
            }
            
            rangeBegin[i] = -1;
            rangeEnd[i] = 0;
            for(std::vector<int> ind = rangeBegin; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
            {
                std::vector<int> ind1 = ind;
                std::vector<int> ind2 = ind;
                std::vector<int> ind3 = ind;
                ind1[i] = nx_[i] - 1;
                ind2[i] = nx_[i];
                ind3[i] = 0;

                unsigned long gridIndex = index(ind);
                unsigned long gridIndex1 = index(ind1);
                unsigned long gridIndex2 = index(ind2);
                unsigned long gridIndex3 = index(ind3);

                grid_[gridIndex] = grid_[gridIndex1];
                grid_[gridIndex2] = grid_[gridIndex3];
            }
        }
    }
}

void
InitialValPDESolver::saveBuffer(unsigned long start)
{
    check(buffer_.size() == dimProd_[1] * m_, "");
    for(unsigned long i = 0; i < dimProd_[1]; ++i)
    {
        for(int j = 0; j < m_; ++j)
            buffer_[i * m_ + j] = grid_[start + i][j];
    }
}

void
InitialValPDESolver::retrieveBuffer(unsigned long start)
{
    check(buffer_.size() == dimProd_[1] * m_, "");
    for(unsigned long i = 0; i < dimProd_[1]; ++i)
    {
        for(int j = 0; j < m_; ++j)
            grid_[start + i][j] = buffer_[i * m_ + j];
    }
}

void
InitialValPDESolver::sendLeft()
{
    const unsigned long size = dimProd_[1];
    const unsigned long leftStart = size;
    const int leftNeighbor = (processId_ + nProcesses_ - 1) % nProcesses_;
    saveBuffer(leftStart);
    const int res = MPI_Send(&(buffer_[0]), buffer_.size(), MPI_DOUBLE, leftNeighbor, commTag_, MPI_COMM_WORLD);
    check(res == MPI_SUCCESS, "send failed");
}

void
InitialValPDESolver::sendRight()
{
    const unsigned long size = dimProd_[1];
    const unsigned long rightStart = nx_[0] * size;
    const int rightNeighbor = (processId_ + 1) % nProcesses_;
    saveBuffer(rightStart);
    const int res = MPI_Send(&(buffer_[0]), buffer_.size(), MPI_DOUBLE, rightNeighbor, commTag_, MPI_COMM_WORLD);
    check(res == MPI_SUCCESS, "send failed");
}

void
InitialValPDESolver::receiveLeft()
{
    const unsigned long size = dimProd_[1];
    const unsigned long leftReceiveStart = 0;
    const int leftNeighbor = (processId_ + nProcesses_ - 1) % nProcesses_;
    MPI_Status s;
    const int res = MPI_Recv(&(buffer_[0]), buffer_.size(), MPI_DOUBLE, leftNeighbor, commTag_, MPI_COMM_WORLD, &s);
    check(res == MPI_SUCCESS, "receive failed");
    retrieveBuffer(leftReceiveStart);
}

void
InitialValPDESolver::receiveRight()
{
    const unsigned long size = dimProd_[1];
    const unsigned long rightReceiveStart = (nx_[0] + 1) * size;
    const int rightNeighbor = (processId_ + 1) % nProcesses_;
    MPI_Status s;
    const int res = MPI_Recv(&(buffer_[0]), buffer_.size(), MPI_DOUBLE, rightNeighbor, commTag_, MPI_COMM_WORLD, &s);
    check(res == MPI_SUCCESS, "receive failed");
    retrieveBuffer(rightReceiveStart);
}

void
InitialValPDESolver::communicateBoundary()
{
    const unsigned long size = dimProd_[1];
    const unsigned long leftStart = size;
    const unsigned long rightStart = nx_[0] * size;
    const unsigned long leftReceiveStart = 0;
    const unsigned long rightReceiveStart = (nx_[0] + 1) * size;
    if(nProcesses_ == 1)
    {
#pragma omp parallel for default(shared)
        for(unsigned long i = 0; i < size; ++i)
        {
            grid_[leftReceiveStart + i] = grid_[rightStart + i];
            grid_[rightReceiveStart + i] = grid_[leftStart + i];
        }

        return;
    }

#ifdef COSMO_MPI
    check(nProcesses_ > 1, "");
    if(processId_ % 2 == 0)
    {
        if(processId_ != nProcesses_ - 1)
        {
            sendRight();
            receiveRight();
        }
        if(processId_ > 0 || nProcesses_ % 2 == 0)
        {
            sendLeft();
            receiveLeft();
        }
    }
    else
    {
        receiveLeft();
        sendLeft();
        receiveRight();
        sendRight();
    }
    if(nProcesses_ % 2)
    {
        if(processId_ == 0)
        {
            receiveLeft();
            sendLeft();
        }
        if(processId_ == nProcesses_ - 1)
        {
            sendRight();
            receiveRight();
        }
    }
#else
    check(false, "should not get here");
#endif

    CosmoMPI::create().barrier();
}

void
InitialValPDESolver::propagate(double t)
{
    check(t >= t_, "cannot propagate back in time");
    check(deltaT_ > 0, "");

    std::auto_ptr<ProgressMeter> meter;

    if(processId_ == 0)
    {
        int nSteps = 0;
        double t1 = t_;
        while(t1 < t)
        {
            ++nSteps;
            t1 += deltaT_;
        }

        if(nSteps > 0)
            meter.reset(new ProgressMeter(nSteps));
    }
    while(t_ < t)
    {
        takeStep();
        if(processId_ == 0)
            meter->advance();
    }
}

void
InitialValPDESolver::takeStep()
{
    // evaluate half grid
#pragma omp parallel for default(shared)
    for(int i0 = 0; i0 <= nx_[0]; ++i0)
    {
        std::vector<int> rangeBegin(d_, 0);
        std::vector<int> rangeEnd = nx_;
        rangeBegin[0] = i0;
        rangeEnd[0] = i0 + 1;

        for(int i = 1; i < d_; ++i)
            ++rangeEnd[i];

        std::vector<double> x(d_);
        std::vector<double> term1(m_, 0), term2(m_, 0), term3(m_, 0);
        std::vector<std::vector<double> > f(d_);
        for(int i = 0; i < d_; ++i)
            f[i].resize(m_);
        std::vector<double> s(m_);

        std::vector<int> rBeg(d_), rEnd(d_), ind1(d_);

        for(std::vector<int> ind = rangeBegin; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
        {
            for(int j = 0; j < d_; ++j)
            {
                rBeg[j] = ind[j] - 1;
                rEnd[j] = ind[j] + 1;
            }

            for(int i = 0; i < m_; ++i)
                term1[i] = term2[i] = term3[i] = 0;

            for(ind1 = rBeg; ind1[0] != rEnd[0]; increaseIndex(ind1, rBeg, rEnd))
            {
                physicalCoords(ind1, &x);
                pde_->evaluate(t_, x, grid_[index(ind1)], &f, &s);

                for(int i = 0; i < m_; ++i)
                {
                    term1[i] += grid_[index(ind1)][i];
                    term3[i] += s[i];

                    for(int j = 0; j < d_; ++j)
                    {
                        const double lambda = deltaT_ / deltaX_[j];
                        const double factor = (ind1[j] == ind[j] ? 1.0 : -1.0);
                        term2[i] += lambda * factor * f[j][i];
                    }
                }
            }
            for(int i = 0; i < m_; ++i)
                halfGrid_[halfIndex(ind)][i] = (term1[i] + term2[i] + deltaT_ * term3[i] / 2) / twoD_;
        }
    }

    // evaluate grid
#pragma omp parallel for default(shared)
    for(int i0 = 0; i0 < nx_[0]; ++i0)
    {
        std::vector<int> rangeBegin(d_, 0);
        std::vector<int> rangeEnd = nx_;
        rangeBegin[0] = i0;
        rangeEnd[0] = i0 + 1;

        std::vector<double> x(d_);
        std::vector<double> term1(m_, 0), term2(m_, 0);
        std::vector<std::vector<double> > f(d_);
        for(int i = 0; i < d_; ++i)
            f[i].resize(m_);
        std::vector<double> s(m_);

        std::vector<int> rBeg(d_), rEnd(d_), ind1(d_);

        for(std::vector<int> ind = rangeBegin; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
        {
            rBeg = ind;
            for(int j = 0; j < d_; ++j)
                rEnd[j] = ind[j] + 2;

            for(int i = 0; i < m_; ++i)
                term1[i] = term2[i] = 0;

            for(ind1 = rBeg; ind1[0] != rEnd[0]; increaseIndex(ind1, rBeg, rEnd))
            {
                physicalCoordsHalf(ind1, &x);
                pde_->evaluate(t_ + deltaT_ / 2, x, halfGrid_[halfIndex(ind1)], &f, &s);
                for(int i = 0; i < m_; ++i)
                {
                    term2[i] += s[i];

                    for(int j = 0; j < d_; ++j)
                    {
                        const double lambda = deltaT_ / deltaX_[j];
                        const double factor = (ind1[j] == ind[j] ? -1.0 : 1.0);
                        term1[i] += lambda * factor * f[j][i];
                    }
                }
            }
            for(int i = 0; i < m_; ++i)
                grid_[index(ind)][i] += (2 * term1[i] + deltaT_ * term2[i]) / twoD_;
        }
    }

    t_ += deltaT_;

    setOwnBoundary();
    communicateBoundary();
}

} // namespace Math

