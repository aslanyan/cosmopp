#ifndef COSMOPP_FIELD_SOLVER_HPP
#define COSMOPP_FIELD_SOLVER_HPP

#include <vector>

#include <function.hpp>

namespace Math
{

class FieldSolver
{
public:
    FieldSolver(int d, int m, const RealFunctionMultiDim* potential, const std::vector<RealFunctionMultiDim*>& potentialDeriv);
    ~FieldSolver();

    void set(const std::vector<RealFunctionMultiDim*>& field, const std::vector<RealFunctionMultiDim*>& fieldDeriv, const std::vector<double>& xMin, const std::vector<double>& xMax, const std::vector<int>& nx);

    void propagate(double t);

    void takeStep();

    double getCurrentT() const { return t_; }
    double getDeltaT() const { return deltaT_; }
    void setDeltaT(double deltaT) { check(deltaT > 0, ""); deltaT_ = deltaT; }

    void setBoundaryPeriodic(int i);
    void setBoundaryValLeft(int i, int j, double val);
    void setBoundaryValRight(int i, int j, double val);
    void setBoundaryGradLeft(int i, int j);
    void setBoundaryGradRight(int i, int j);

    int getNx(int i) const { check(i >= 0 && i < d_, ""); return nx_[i]; }
    int getNx0Starting() const { return nx0Starting_; }
    double getXMin(int i) const { check(i >= 0 && i < d_, ""); return xMin_[i]; }
    double getXMax(int i) const { check(i >= 0 && i < d_, ""); return xMax_[i]; }
    double getDeltaX(int i) const { check(i >= 0 && i < d_, ""); return deltaX_[i]; }

    const std::vector<double>& getField(const std::vector<int>& ind) const { return grid_[index(&(ind[0]))]; }
    const std::vector<double>& getFieldDeriv(const std::vector<int>& ind) const { return gridDeriv_[index(&(ind[0]))]; }

    inline void physicalCoords(const int *i, double *coords) const;

private:
    void setupGrid();
    void setInitial(const std::vector<RealFunctionMultiDim*>& field, const std::vector<RealFunctionMultiDim*>& fieldDeriv);
    void setOwnBoundary();
    void setOwnBoundary0();

    void communicateBoundary();
    void sendLeft();
    void sendRight();
    void receiveLeft();
    void receiveRight();

    void saveBuffer(unsigned long start);
    void retrieveBuffer(unsigned long start);

    inline unsigned long index(const int *i) const;
    inline void increaseIndex(int *i, const int *rangeBegin, const int *rangeEnd) const;

    void allocateStorage();

private:
    const RealFunctionMultiDim *potential_;
    const std::vector<RealFunctionMultiDim*> potentialDeriv_;

    const int d_;
    const int m_;

    double twoD_;

    int nProcesses_;
    int processId_;

    std::vector<double> xMin_;
    std::vector<double> xMax_;
    std::vector<double> deltaX_;
    std::vector<int> nx_;

    int nx0Starting_;

    std::vector<unsigned long> dimProd_;

    std::vector<std::vector<double> > grid_;
    std::vector<std::vector<double> > gridDeriv_;

    double t_;
    double deltaT_;

    int commTag_;
    std::vector<double> buffer_;

    std::vector<bool> boundaryPeriodic_;
    std::vector<std::vector<bool> > boundaryGradLeft_;
    std::vector<std::vector<bool> > boundaryGradRight_;
    std::vector<std::vector<double> > boundaryValLeft_;
    std::vector<std::vector<double> > boundaryValRight_;
};

unsigned long
FieldSolver::index(const int *i) const
{
    check(dimProd_.size() == d_ + 1, "");

    unsigned long res = 0;

    for(int j = 0; j < d_; ++j)
    {
        const int k = i[j];
        check(k >= -1 && k <= nx_[j], "");

        res += (unsigned long)(k + 1) * dimProd_[j + 1];
    }

    return res;
}

void
FieldSolver::physicalCoords(const int *i, double *coords) const
{
    check(i[0] >= -1 && i[0] <= nx_[0], "");
    coords[0] = xMin_[0] + deltaX_[0] * (i[0] + nx0Starting_);

    for(int j = 1; j < d_; ++j)
    {
        check(i[j] >= -1 && i[j] <= nx_[j], "");
        coords[j] = xMin_[j] + deltaX_[j] * i[j];
    }
}

void
FieldSolver::increaseIndex(int *i, const int *rangeBegin, const int *rangeEnd) const
{
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

} // namespace Math

#endif

