#include <cosmo_mpi.hpp>

#include <memory>

#include <macros.hpp>
#include <progress_meter.hpp>
#include <field_solver.hpp>

#ifdef COSMO_OMP
#include <omp.h>
#endif

namespace
{

const int fsMaxDim = 10;
const int fsMaxFields = 100;

}

namespace Math
{

FieldSolver::FieldSolver(int d, int m, const RealFunctionMultiDim* potential, const std::vector<RealFunctionMultiDim*>& potentialDeriv) : potential_(potential), potentialDeriv_(potentialDeriv), d_(d), m_(m), grid_(0), gridDeriv_(0), boundaryPeriodic_(d, true)
{ 
    check(d_ >= 1, "invalid number of spatial dimensions " << d_ << ", must be positive");
    check(m_ >= 1, "invalid number of fields " << m_ << ", must be positive");

    check(d_ <= fsMaxDim, "the dimensionality of " << d_ << " is too high! Maximum supported is " << fsMaxDim << ". To change this, update fsMaxDim above.");
    check(m_ <= fsMaxFields, "the number of fields " << m_ << " is too high! Maximum supported is " << fsMaxFields << ". To change this, update fsMaxFields above.");

    twoD_ = 1;
    for(int i = 0; i < d_; ++i)
        twoD_ *= 2;

    nProcesses_ = CosmoMPI::create().numProcesses();
    processId_ = CosmoMPI::create().processId();
    check(processId_ >= 0 && processId_ < nProcesses_, "");

    commTag_ = CosmoMPI::create().getCommTag();

    boundaryGradLeft_.resize(d_);
    boundaryGradRight_.resize(d_);
    boundaryValLeft_.resize(d_);
    boundaryValRight_.resize(d_);
    for(int i = 0; i < d_; ++i)
    {
        boundaryGradLeft_[i].resize(m_, false);
        boundaryGradRight_[i].resize(m_, false);
        boundaryValLeft_[i].resize(m_, 0);
        boundaryValRight_[i].resize(m_, 0);
    }
}

FieldSolver::~FieldSolver()
{
}

void
FieldSolver::set(const std::vector<RealFunctionMultiDim*>& field, const std::vector<RealFunctionMultiDim*>& fieldDeriv, const std::vector<double>& xMin, const std::vector<double>& xMax, const std::vector<int>& nx)
{
    check(field.size() == m_, "");
    check(fieldDeriv.size() == m_, "");
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
    setInitial(field, fieldDeriv);
    setOwnBoundary();
    communicateBoundary();
    setOwnBoundary0();
}

void
FieldSolver::setupGrid()
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
    for(int i = d_ - 1; i >= 0; --i)
    {
        unsigned long currentDim = (unsigned long)(nx_[i]) + 2;
        dimProd_[i] = dimProd_[i + 1] * currentDim;
        check(dimProd_[i] >= dimProd_[i + 1], "overflow!");
        check(dimProd_[i] / dimProd_[i + 1] == currentDim, "overflow!");
    }

    grid_.clear();
    gridDeriv_.clear();

    grid_.resize(dimProd_[0]);
    gridDeriv_.resize(dimProd_[0]);

    for(unsigned long i = 0; i < dimProd_[0]; ++i)
    {
        grid_[i].resize(m_);
        gridDeriv_[i].resize(m_);
    }

    allocateStorage();
}

void
FieldSolver::allocateStorage()
{
    buffer_.resize(2 * dimProd_[1] * m_);
}

void
FieldSolver::setInitial(const std::vector<RealFunctionMultiDim*>& field, const std::vector<RealFunctionMultiDim*>& fieldDeriv)
{
    check(field.size() == m_, "");
    check(fieldDeriv.size() == m_, "");
    check(grid_.size() == dimProd_[0], "");
    check(gridDeriv_.size() == dimProd_[0], "");

    for(int i = 0; i < m_; ++i)
    {
//#pragma omp parallel for default(shared)
        for(int i0 = 0; i0 < nx_[0]; ++i0)
        {
            int rangeBegin[fsMaxDim];
            int rangeEnd[fsMaxDim];
            rangeBegin[0] = i0;
            rangeEnd[0] = i0 + 1;
            for(int j = 1; j < d_; ++j)
            {
                rangeBegin[j] = 0;
                rangeEnd[j] = nx_[j];
            }

            std::vector<double> coords(d_);

            int ind[fsMaxDim];
            for(int j = 0; j < d_; ++j)
                ind[j] = rangeBegin[j];

            for(; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
            {
                physicalCoords(ind, &(coords[0]));
                const unsigned long gridIndex = index(ind);
                grid_[gridIndex][i] = field[i]->evaluate(coords);
                gridDeriv_[gridIndex][i] = fieldDeriv[i]->evaluate(coords);
            }
        }
    }

    t_ = 0;
}

void
FieldSolver::setOwnBoundary0()
{
    check(boundaryPeriodic_.size() == d_, "");
    if(boundaryPeriodic_[0])
        return;

    if(processId_ == 0)
    {
        // set left boundary
        int rangeBegin[fsMaxDim], rangeEnd[fsMaxDim];
        rangeBegin[0] = -1;
        rangeEnd[0] = 0;
        for(int j = 1; j < d_; ++j)
        {
            rangeBegin[j] = -1;
            rangeEnd[j] = nx_[j] + 1;
        }

        int ind[fsMaxDim], ind1[fsMaxDim];
        for(int j = 0; j < d_; ++j)
            ind[j] = rangeBegin[j];

        ind1[0] = 0;

        for(; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
        {
            const unsigned long gridIndex = index(ind);
            for(int l = 0; l < m_; ++l)
            {
                if(boundaryGradLeft_[0][l])
                {
                    for(int j = 1; j < d_; ++j)
                        ind1[j] = ind[j];
                    const unsigned long gridIndex1 = index(ind1);
                    grid_[gridIndex][l] = grid_[gridIndex1][l];
                    gridDeriv_[gridIndex][l] = gridDeriv_[gridIndex1][l];
                }
                else
                    grid_[gridIndex][l] = boundaryValLeft_[0][l];
                    gridDeriv_[gridIndex][l] = 0;
            }
        }
    }

    if(processId_ == nProcesses_ - 1)
    {
        // set right boundary
        int rangeBegin[fsMaxDim], rangeEnd[fsMaxDim];
        rangeBegin[0] = nx_[0];
        rangeEnd[0] = nx_[0] + 1;
        for(int j = 1; j < d_; ++j)
        {
            rangeBegin[j] = -1;
            rangeEnd[j] = nx_[j] + 1;
        }

        int ind[fsMaxDim], ind1[fsMaxDim];
        for(int j = 0; j < d_; ++j)
            ind[j] = rangeBegin[j];

        ind1[0] = nx_[0] - 1;

        for(; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
        {
            const unsigned long gridIndex = index(ind);
            for(int l = 0; l < m_; ++l)
            {
                if(boundaryGradRight_[0][l])
                {
                    for(int j = 1; j < d_; ++j)
                        ind1[j] = ind[j];
                    const unsigned long gridIndex1 = index(ind1);
                    grid_[gridIndex][l] = grid_[gridIndex1][l];
                    gridDeriv_[gridIndex][l] = gridDeriv_[gridIndex1][l];
                }
                else
                    grid_[gridIndex][l] = boundaryValRight_[0][l];
                    gridDeriv_[gridIndex][l] = 0;
            }
        }
    }
}

void
FieldSolver::setOwnBoundary()
{
    for(int i = 1; i < d_; ++i)
    {
#pragma omp parallel for default(shared)
        for(int i0 = 0; i0 < nx_[0]; ++i0)
        {
            int rangeBegin[fsMaxDim];
            int rangeEnd[fsMaxDim];

            rangeBegin[0] = i0;
            rangeEnd[0] = i0 + 1;
            for(int j = 1; j < i; ++j)
            {
                rangeBegin[j] = -1;
                rangeEnd[j] = nx_[j] + 1;
            }
            
            rangeBegin[i] = -1;
            rangeEnd[i] = 0;
            
            for(int j = i + 1; j < d_; ++j)
            {
                rangeBegin[j] = 0;
                rangeEnd[j] = nx_[j];
            }

            int ind[fsMaxDim];
            int ind1[fsMaxDim];

            for(int j = 0; j < d_; ++j)
                ind[j] = rangeBegin[j];

            for(; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
            {
                const unsigned long gridIndex = index(ind);
                for(int j = 0; j < d_; ++j)
                    ind1[j] = ind[j];
                if(boundaryPeriodic_[i])
                {
                    ind1[i] = nx_[i] - 1;
                    const unsigned long gridIndex1 = index(ind1);
                    ind1[i] = nx_[i];
                    const unsigned long gridIndex2 = index(ind1);
                    ind1[i] = 0;
                    const unsigned long gridIndex3 = index(ind1);

                    grid_[gridIndex] = grid_[gridIndex1];
                    gridDeriv_[gridIndex] = gridDeriv_[gridIndex1];
                    grid_[gridIndex2] = grid_[gridIndex3];
                    gridDeriv_[gridIndex2] = gridDeriv_[gridIndex3];
                }
                else
                {
                    for(int l = 0; l < m_; ++l)
                    {
                        if(boundaryGradLeft_[i][l])
                        {
                            ind1[i] = 0;
                            const unsigned long gridIndex1 = index(ind1);
                            grid_[gridIndex][l] = grid_[gridIndex1][l];
                            gridDeriv_[gridIndex][l] = gridDeriv_[gridIndex1][l];
                        }
                        else
                        {
                            grid_[gridIndex][l] = boundaryValLeft_[i][l];
                            gridDeriv_[gridIndex][l] = 0;
                        }
                        if(boundaryGradRight_[i][l])
                        {
                            ind1[i] = nx_[i];
                            const unsigned long gridIndex2 = index(ind1);
                            ind1[i] = nx_[i] - 1;
                            const unsigned long gridIndex3 = index(ind1);
                            grid_[gridIndex2][l] = grid_[gridIndex3][l];
                            gridDeriv_[gridIndex2][l] = gridDeriv_[gridIndex3][l];
                        }
                        else
                        {
                            ind1[i] = nx_[i];
                            const unsigned long gridIndex2 = index(ind1);
                            grid_[gridIndex2][l] = boundaryValRight_[i][l];
                            gridDeriv_[gridIndex2][l] = 0;
                        }
                    }
                }
            }
        }
    }
}

void
FieldSolver::saveBuffer(unsigned long start)
{
    check(buffer_.size() == 2 * dimProd_[1] * m_, "");
#pragma omp parallel for default(shared)
    for(unsigned long i = 0; i < dimProd_[1]; ++i)
    {
        for(int j = 0; j < m_; ++j)
        {
            buffer_[i * m_ + j] = grid_[start + i][j];
            buffer_[dimProd_[1] * m_ + i * m_ + j] = gridDeriv_[start + i][j];
        }
    }
}

void
FieldSolver::retrieveBuffer(unsigned long start)
{
    check(buffer_.size() == 2 * dimProd_[1] * m_, "");
#pragma omp parallel for default(shared)
    for(unsigned long i = 0; i < dimProd_[1]; ++i)
    {
        for(int j = 0; j < m_; ++j)
        {
            grid_[start + i][j] = buffer_[i * m_ + j];
            gridDeriv_[start + i][j] = buffer_[dimProd_[1] * m_ + i * m_ + j];
        }
    }
}

void
FieldSolver::sendLeft()
{
    const unsigned long size = dimProd_[1];
    const unsigned long leftStart = size;
    const int leftNeighbor = (processId_ + nProcesses_ - 1) % nProcesses_;
    saveBuffer(leftStart);
    const int res = MPI_Send(&(buffer_[0]), buffer_.size(), MPI_DOUBLE, leftNeighbor, commTag_, MPI_COMM_WORLD);
    check(res == MPI_SUCCESS, "send failed");
}

void
FieldSolver::sendRight()
{
    const unsigned long size = dimProd_[1];
    const unsigned long rightStart = nx_[0] * size;
    const int rightNeighbor = (processId_ + 1) % nProcesses_;
    saveBuffer(rightStart);
    const int res = MPI_Send(&(buffer_[0]), buffer_.size(), MPI_DOUBLE, rightNeighbor, commTag_, MPI_COMM_WORLD);
    check(res == MPI_SUCCESS, "send failed");
}

void
FieldSolver::receiveLeft()
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
FieldSolver::receiveRight()
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
FieldSolver::communicateBoundary()
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
            gridDeriv_[leftReceiveStart + i] = gridDeriv_[rightStart + i];
            grid_[rightReceiveStart + i] = grid_[leftStart + i];
            gridDeriv_[rightReceiveStart + i] = gridDeriv_[leftStart + i];
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
FieldSolver::propagate(double t)
{
    check(t >= t_, "cannot propagate back in time");
    check(deltaT_ > 0, "");

    std::unique_ptr<ProgressMeter> meter;

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
FieldSolver::takeStep()
{
    // evaluate deriv grid
#pragma omp parallel for default(shared)
    for(int i0 = 0; i0 < nx_[0]; ++i0)
    {
        int rangeBegin[fsMaxDim];
        int rangeEnd[fsMaxDim];
        rangeBegin[0] = i0;
        rangeEnd[0] = i0 + 1;
        for(int i = 1; i < d_; ++i)
        {
            rangeBegin[i] = 0;
            rangeEnd[i] = nx_[i];
        }

        int ind[fsMaxDim];
        int ind1[fsMaxDim];

        for(int i = 0; i < d_; ++i)
            ind[i] = rangeBegin[i];

        for(; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
        {
            for(int j = 0; j < d_; ++j)
                ind1[j] = ind[j];

            for(int i = 0; i < m_; ++i)
            {
                double fieldDotDot = 0;

                for(int j = 0; j < d_; ++j)
                {
                    const double dx2 = deltaX_[j] * deltaX_[j];
                    fieldDotDot -= 2 * grid_[index(ind)][i] / dx2;

                    check(ind1[j] == ind[j], "");
                    --ind1[j];
                    fieldDotDot += grid_[index(ind1)][i] / dx2;
                    ind1[j] += 2;
                    fieldDotDot += grid_[index(ind1)][i] / dx2;
                    --ind1[j];
                    check(ind1[j] == ind[j], "");
                }

                fieldDotDot -= potentialDeriv_[i]->evaluate(grid_[index(ind)]);

                gridDeriv_[index(ind)][i] += deltaT_ * fieldDotDot;
            }
        }
    }

    // evaluate grid
#pragma omp parallel for default(shared)
    for(int i0 = 0; i0 < nx_[0]; ++i0)
    {
        int rangeBegin[fsMaxDim];
        int rangeEnd[fsMaxDim];
        rangeBegin[0] = i0;
        rangeEnd[0] = i0 + 1;
        for(int i = 1; i < d_; ++i)
        {
            rangeBegin[i] = 0;
            rangeEnd[i] = nx_[i];
        }

        int ind[fsMaxDim];

        for(int i = 0; i < d_; ++i)
            ind[i] = rangeBegin[i];

        for(; ind[0] != rangeEnd[0]; increaseIndex(ind, rangeBegin, rangeEnd))
        {
            for(int i = 0; i < m_; ++i)
            {
                grid_[index(ind)][i] += deltaT_ * gridDeriv_[index(ind)][i];
            }
        }
    }

    t_ += deltaT_;

    setOwnBoundary();
    communicateBoundary();
    setOwnBoundary0();
}

void
FieldSolver::setBoundaryPeriodic(int i)
{
    check(i >= 0 && i < d_, "");
    check(boundaryPeriodic_.size() == d_, "");
    boundaryPeriodic_[i] = true;
}

void
FieldSolver::setBoundaryGradLeft(int i, int j)
{
    check(i >= 0 && i < d_, "");
    check(j >= 0 && j < m_, "");
    check(boundaryGradLeft_.size() == d_, "");
    check(boundaryGradLeft_[i].size() == m_, "");

    boundaryGradLeft_[i][j] = true;
    
    check(boundaryPeriodic_.size() == d_, "");
    boundaryPeriodic_[i] = false;
}

void
FieldSolver::setBoundaryGradRight(int i, int j)
{
    check(i >= 0 && i < d_, "");
    check(j >= 0 && j < m_, "");
    check(boundaryGradRight_.size() == d_, "");
    check(boundaryGradRight_[i].size() == m_, "");

    boundaryGradRight_[i][j] = true;
    
    check(boundaryPeriodic_.size() == d_, "");
    boundaryPeriodic_[i] = false;
}

void
FieldSolver::setBoundaryValLeft(int i, int j, double val)
{
    check(i >= 0 && i < d_, "");
    check(j >= 0 && j < m_, "");

    boundaryValLeft_[i][j] = val;

    check(boundaryGradLeft_.size() == d_, "");
    check(boundaryGradLeft_[i].size() == m_, "");

    boundaryGradLeft_[i][j] = false;

    check(boundaryPeriodic_.size() == d_, "");
    boundaryPeriodic_[i] = false;
}

void
FieldSolver::setBoundaryValRight(int i, int j, double val)
{
    check(i >= 0 && i < d_, "");
    check(j >= 0 && j < m_, "");

    boundaryValRight_[i][j] = val;

    check(boundaryGradRight_.size() == d_, "");
    check(boundaryGradRight_[i].size() == m_, "");

    boundaryGradRight_[i][j] = false;

    check(boundaryPeriodic_.size() == d_, "");
    boundaryPeriodic_[i] = false;
}

} // namespace Math

