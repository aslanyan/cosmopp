#ifndef COSMO_PDE_HPP
#define COSMO_PDE_HPP

#include <vector>

#include <function.hpp>

namespace Math
{

class InitialValPDEInterface
{
public:
    InitialValPDEInterface() {}
    virtual ~InitialValPDEInterface() {}

    virtual int spaceDim() const = 0;
    virtual int funcDim() const = 0;

    virtual void evaluate(double t, const std::vector<double>& x, const std::vector<double>& u, std::vector<std::vector<double> > *f, std::vector<double> *s) const = 0;
};

class InitialValPDESolver
{
public:
    InitialValPDESolver(InitialValPDEInterface *pde);
    ~InitialValPDESolver();

    void set(const std::vector<RealFunctionMultiDim*>& w0, const std::vector<double>& xMin, const std::vector<double>& xMax, const std::vector<int>& nx);

    void propagate(double t);

    void takeStep();

    double getCurrentT() const { return t_; }
    double getDeltaT() const { return deltaT_; }
    void setDeltaT(double deltaT) { check(deltaT > 0, ""); deltaT_ = deltaT; }

    int getNx(int i) const { check(i >= 0 && i < d_, ""); return nx_[i]; }
    int getNx0Starting() const { return nx0Starting_; }
    int getXMax(int i) const { check(i >= 0 && i < d_, ""); return xMax_[i]; }
    int getDeltaX(int i) const { check(i >= 0 && i < d_, ""); return deltaX_[i]; }

    const std::vector<double>& getField(const std::vector<int>& ind) const { return grid_[index(ind)]; }

    void physicalCoords(const std::vector<int>& i, std::vector<double> *coords) const;

private:
    void setupGrid();
    void setInitial(const std::vector<RealFunctionMultiDim*>& w0);
    void setOwnBoundary();

    void communicateBoundary();
    void sendLeft();
    void sendRight();
    void receiveLeft();
    void receiveRight();

    void saveBuffer(unsigned long start);
    void retrieveBuffer(unsigned long start);

    unsigned long index(const std::vector<int>& i) const;
    unsigned long halfIndex(const std::vector<int>& i) const;
    void increaseIndex(std::vector<int>& i, const std::vector<int>& rangeBegin, const std::vector<int>& rangeEnd) const;
    void physicalCoordsHalf(const std::vector<int>& i, std::vector<double> *coords) const;

private:
    const InitialValPDEInterface *pde_;

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
    std::vector<unsigned long> halfDimProd_;

    std::vector<std::vector<double> > grid_;
    std::vector<std::vector<double> > halfGrid_;

    double t_;
    double deltaT_;

    int commTag_;
    std::vector<double> buffer_;
};

} // namespace Math

#endif

