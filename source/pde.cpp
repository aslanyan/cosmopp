#include <cosmo_mpi.hpp>


#include <macros.hpp>
#include <pde.hpp>

namespace Math
{

InitialValPDESolver::InitialValPDESolver(InitialValPDEInterface *pde) : pde_(pde), output_(0), d_(pde->spaceDim()), m_(pde->funcDim()), grid_(0), halfGrid_(0)
{ 
    check(d_ >= 1, "invalid number of spatial dimensions " << d_ << ", must be positive");
    check(m_ >= 1, "invalid number of fields " << m_ << ", must be positive");

    nProcesses_ = CosmoMPI::create().numProcesses();
    processId_ = CosmoMPI::create().processId();
    check(processId_ >= 0 && processId_ < nProcesses_, "");

    commTag_ = CosmoMPI::create().getCommTag();
}

InitialValPDESolver::~InitialValPDESolver()
{
}

void
InitialValPDESolver::set(const std::vector<FunctionMultiDim*>& w0, const std::vector<double>& xMin, const std::vector<double>& xMax, const std::vector<int>& nx, InitialValPDEOutputHandlerInterface *output)
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
    output_ = output;

    deltaX_.resize(d_);
    double minDeltaX = (xMax_[0] - xMin_[0]);
    for(int i = 0; i < d_; ++i)
    {
        deltaX_[i] = (xMax_[i] - xMin_[i]) / nx_[i];
        if(minDeltaX > deltaX_[i])
            minDeltaX = deltaX_[i];
    }

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
    nx0_ = nx0PerProcess[processId_];

    dimProd_.resize(d_ + 1)
    dimProd_[d_] = 1;
    halfDimProd_.resize(d_ + 1);
    halfDimProd_[d] = 1;
    for(int i = d_ - 1; i >= 0; --i)
    {
        unsigned long currentDim = (unsigned long)(i == 0 ? nx0_ : nx_[d_ - 1]) + 2;
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

    grid_.resize(m_);
    halfGrid_.resize(m_);

    for(int i = 0; i < m_; ++i)
    {
        grid_[i].resize(dimProd_[0]);
        halfGrid_[i].resize((nx0_ + 1) * dimProd_[1]);
    }
}

unsigned long
InitialValPDESolver::index(const std::vector<int>& i) const
{
    check(i.size() == d_);
    check(dimProd_.size() == d_ + 1, "");

    unsigned long res = 0;

    for(int j = 0; j < d_; ++j)
    {
        const int k = i[j];
        check(k >= -1 && k <= (j == 0 ? nx0_ : nx_[j]), "");

        res += (unsigned long)(k + 1) * dimProd_[j + 1];
    }

    return res;
}

unsigned long
InitialValPDESolver::halfIndex(const std::vector<int>& i) const
{
    check(i.size() == d_);
    check(halfDimProd_.size() == d_ + 1, "");

    unsigned long res = 0;

    for(int j = 0; j < d_; ++j)
    {
        const int k = i[j];
        check(k >= 0 && k <= (j == 0 ? nx0_ : nx_[j]), "");

        res += (unsigned long)(k) * halfDimProd_[j + 1];
    }

    return res;
}

void
InitialValPDESolver::physicalCoords(const std::vector<int>& i, std::vector<double> *coords) const
{
    check(i.size() == d_);
    check(coords->size() == d_);

    check(i[0] >= 0 && i[0] < nx0_, "");
    (*coords)[0] = xMin_[0] + deltaX_[0] * (i[0] + nx0Starting_);

    for(int j = 1; j < d_; ++j)
    {
        check(i[j] >= 0 && i[j] < nx_[j], "");
        (*coords)[j] = xMin_[j] + deltaX_[j] * i[j];
    }
}

void
InitialValPDESolver::increaseIndex(std::vector<int>& i, const std::vector<int>& rangeBegin, const std::vector<int>& rangeEnd) const
{
    check(i.size() == d_, "");
    check(rangeBegin.size() == d_, "");
    check(rangeEnd.size() == d_, "");

    for(int j = d_ - 1; j >= 0; --j)
    {
        const int nx = nx_[j];
        check(i[j] >= rangeBegin[j] && i[j] < rangeEnd[j], "");
        if(++i[j] < rangeEnd[j])
            return;
        if(j > 0)
            i[j] = rangeBegin[j];
    }
}

void
InitialValPDESolver::setInitial(const std::vector<FunctionMultiDim*>& w0)
{
    check(w0.size() == m_, "");
    check(grid_.size() == m_, "");

    for(int i = 0; i < m_; ++i)
    {
#pragma omp parallel for default(shared)
        for(int i0 = 0; i0 < nx0_; ++i0)
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
                grid_[i][gridIndex] = w0[i].evaluate(coords);
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
        for(int i0 = 0; i0 < nx0_; ++i0)
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

                for(int j = 0; j < m_; ++j)
                {
                    grid_[j][gridIndex] = grid_[j][gridIndex1];
                    grid_[j][gridIndex2] = grid_[j][gridIndex3];
                }
            }
        }
    }
}

void
InitialValPDESolver::sendLeft()
{
    const unsigned long size = dimProd_[1];
    const unsigned long leftStart = size;
    const int leftNeighbor = (processId_ + nProcesses_ - 1) % nProcesses_;
    for(int j = 0; j < m_; ++j)
    {
        const int res = MPI_Send(&(grid_[j][leftStart]), size, MPI_DOUBLE, leftNeighbor, commTag_ + j, MPI_COMM_WORLD);
        check(res == MPI_SUCCESS, "send failed");
    }
}

void
InitialValPDESolver::sendRight()
{
    const unsigned long size = dimProd_[1];
    const unsigned long rightStart = nx0_ * size;
    const int rightNeighbor = (processId_ + 1) % nProcesses_;
    for(int j = 0; j < m_; ++j)
    {
        const int res = MPI_Send(&(grid_[j][rightStart]), size, MPI_DOUBLE, rightNeighbor, commTag_ + j, MPI_COMM_WORLD);
        check(res == MPI_SUCCESS, "send failed");
    }
}

void
InitialPDESolver::receiveLeft()
{
    const unsigned long size = dimProd_[1];
    const unsigned long leftReceiveStart = 0;
    const int leftNeighbor = (processId_ + nProcesses_ - 1) % nProcesses_;
    for(int j = 0; j < m_; ++j)
    {
        const int res = MPI_Recv(&(grid_[j][leftReceiveStart]), size, MPI_DOUBLE, leftNeighbor, commTag_ + j, MPI_COMM_WORLD);
        check(res == MPI_SUCCESS, "receive failed");
    }
}

void
InitialPDESolver::receiveRight()
{
    const unsigned long size = dimProd_[1];
    const unsigned long rightReceiveStart = (nx0_ + 1) * size;
    const int rightNeighbor = (processId_ + 1) % nProcesses_;
    for(int j = 0; j < m_; ++j)
    {
        const int res = MPI_Recv(&(grid_[j][rightReceiveStart]), size, MPI_DOUBLE, rightNeighbor, commTag_ + j, MPI_COMM_WORLD);
        check(res == MPI_SUCCESS, "receive failed");
    }
}

void
InitialValPDESolver::communicateBoundary()
{
    const unsigned long size = dimProd_[1];
    const unsigned long leftStart = size;
    const unsigned long leftReceiveStart = 0;
    const unsigned long rightReceiveStart = (nx0_ + 1) * size;
    if(nProcesses_ == 1)
    {
#pragma omp parallel for default(shared)
        for(unsigned long i = 0; i < size; ++i)
        {
            for(int j = 0; j < m_; ++j)
            {
                grid_[j][leftReceiveStart + i] = grid_[j][rightStart + i];
                grid_[j][rightReceiveStart + i] = grid_[j][leftStart + i];
            }
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
        if(processId_ > 0 || nProcesses_ % 2)
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
    if(nProcesses_ % 2 == 0)
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

} // namespace Math

