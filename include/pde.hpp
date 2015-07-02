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

    virtual double f(int i, int j, double t, const double *x, const double *u) const = 0;
    virtual double s(int i, double t, const double *x, const double *u) const = 0;
};

class InitialValPDEOutputHandlerInterface
{
public:
    InitialValPDEOutputHandlerInterface() {}
    virtual ~InitialValPDEOutputHandlerInterface() {}

    virtual void set(const std::vector<double>& xMin, const std::vector<double>& xMax, const std::vector<double>& deltaX, const std::vector<int>& nx, int nx0Starting);

    virtual void setTimeSlice(int i, double t, const double *res) = 0;
};

class InitalValPDESolver
{
public:
    InitialValPDESolver(InitialValPDEInterface *pde);
    ~InitialValPDESolver();

    void set(const std::vector<FunctionMultiDim*>& w0, const std::vector<double>& xMin, const std::vector<double>& xMax, const std::vector<int>& nx, InitialValPDEOutputHandlerInterface *output);

    void propagate(double t);

    double getCurrentT() const { return t_; }

private:
    void setupGrid();
    void setInitial(const std::vector<FunctionMultiDim*>& w0);
    void setOwnBoundary();

    void communicateBoundary();
    void sendLeft();
    void sendRight();
    void receiveLeft();
    void receiveRight();

    void sendOutput();

    void takeStep();

    void setU(double *u, std::vector<int>& ind) const;
    void setUHalf(double *u, std::vector<int>& ind) const;

    unsigned long index(const std::vector<int>& i) const;
    unsigned long halfIndex(const std::vector<int>& i) const;
    void increaseIndex(std::vector<int>& i, const std::vector<int>& rangeBegin, const std::vector<int>& rangeEnd) const;
    void physicalCoords(const std::vector<int>& i, std::vector<double> *coords) const;
    void physicalCoordsHalf(const std::vector<int>& i, std::vector<double> *coords) const;

private:
    const InitialValPDEInterface *pde_;
    InitialValPDEOutputHandlerInterface *output_;

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
};

} // namespace Math

#endif

