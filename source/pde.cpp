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
    sendBoundary();
    receiveBoundary();
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
InitialValPDESolver::index(const std::vector<int>& i)
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
InitialValPDESolver::halfIndex(const std::vector<int>& i)
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
InitialValPDESolver::physicalCoords(const std::vector<int>& i, std::vector<double> *coords)
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

} // namespace Math

