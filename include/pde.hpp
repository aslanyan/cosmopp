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

    virtual void f(double t, const double *x, const double *u, double *res);
    virtual void s(double t, const double *x, const double *u, double *res);
};

class InitialValPDEOutputHandlerInterface
{
public:
    InitialValPDEOutputHandlerInterface() {}
    virtual ~InitialValPDEOutputHandlerInterface() {}

    virtual void setGrid(const std::vector<double>& xMin, const std::vector<double>& xMax, const std::vector<int>& nx);
    virtual void setTimeSlice(double t, double *res);
};

class InitalValPDESolver
{
public:
    InitialValPDESolver(InitialValPDEInterface *pde);
    ~InitialValPDESolver();

    void set(const std::vector<FunctionMultiDim*>& w0, const std::vector<double>& xMin, const std::vector<double>& xMax, const std::vector<int>& nx, InitialValPDEOutputHandlerInterface* output);

    void propagate(double t);

    double getCurrentT() const;
};

} // namespace Math

#endif

