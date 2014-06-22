#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <numerics.hpp>
#include <inflation.hpp>
#include <test_inflation.hpp>

std::string
TestInflation::name() const
{
    return std::string("INFLATION TESTER");
}

unsigned int
TestInflation::numberOfSubtests() const
{
    return 1;
}

class QuadraticFunc : public Math::RealFunction
{
public:
    QuadraticFunc(double m, double x0 = 0) : mSqHalf_(m * m / 2), x0_(x0)
    {
        check(m >= 0, "");
    }

    ~QuadraticFunc() {}

    virtual double evaluate(double x) const
    {
        return mSqHalf_ * (x - x0_) * (x - x0_);
    }

private:
    const double mSqHalf_, x0_;
};

class LinearFunc : public Math::RealFunction
{
public:
    LinearFunc(double m, double x0 = 0) : mSq_(m * m), x0_(x0)
    {
        check(m >= 0, "");
    }

    ~LinearFunc() {}

    virtual double evaluate(double x) const
    {
        return mSq_ * (x - x0_);
    }

private:
    const double mSq_, x0_;
};

void
TestInflation::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);

    const double m = 6.7e-6;
    const double kPivot = 0.05;
    const double NPivot = 55;

    const QuadraticFunc potential(m);
    const LinearFunc potentialDeriv(m);

    SingleFieldInflation inflation(potential, potentialDeriv, 20, 0, std::sqrt(2), kPivot, NPivot, 1, 0);

    const double kMin = 1e-6, kMax = 1.0;
    Math::TableFunction<double, double>* scalarPs = inflation.scalarPs(kMin, kMax);
    Math::TableFunction<double, double>* tensorPs = inflation.tensorPs(kMin, kMax);

    std::stringstream fileName;
    fileName << "test_files/inflation_pk.txt";
    std::ofstream out(fileName.str().c_str());
    StandardException exc;
    if(!out)
    {
        std::stringstream exceptionStr;
        exceptionStr << "Cannot write into file " << fileName << ".";
        exc.set(exceptionStr.str());
        throw exc;
    }
    for(Math::TableFunction<double, double>::const_iterator it = scalarPs->begin(); it != scalarPs->end(); ++it)
    {
        const double k = (*it).first;
        const double s = (*it).second;
        const double t = tensorPs->evaluate(k);
        out << k << '\t' << s << '\t' << t << std::endl;
    }
    out.close();

    const double sMin = inflation.scalarP(kMin), sMax = inflation.scalarP(kMax), sPivot = inflation.scalarP(kPivot);
    const double tMin = inflation.tensorP(kMin), tMax = inflation.tensorP(kMax), tPivot = inflation.tensorP(kPivot);

    const double ns = 1 + (std::log(sMax) - std::log(sMin)) / (std::log(kMax) - std::log(kMin));
    const double nt = (std::log(tMax) - std::log(tMin)) / (std::log(kMax) - std::log(kMin));

    const double r = tPivot / sPivot;
    const double phiPiv = inflation.phiPivot();

    const double phiPivExp = std::sqrt(4.0 * NPivot + 2);
    const double nsExp = 1.0 - 8.0 / (phiPivExp * phiPivExp);
    const double ntExp = -4.0 / (phiPivExp * phiPivExp);
    const double rExp = -8 * ntExp;

    output_screen1("ns = " << ns << ", expected = " << nsExp << std::endl);
    output_screen1("nt = " << nt << ", expected = " << ntExp << std::endl);
    output_screen1("r = " << r << ", expected = " << rExp << std::endl);
    output_screen1("phiPiv = " << phiPiv << ", expected = " << phiPivExp << std::endl);

    res = 1;
    expected = 1;
    subTestName = "quadratic_potential";

    if(!Math::areEqual(nsExp, ns, 0.1))
    {
        output_screen("FAIL: Expected ns = " << nsExp << ", result = " << ns << std::endl);
        res = 0;
    }
    if(!Math::areEqual(ntExp, nt, 0.1))
    {
        output_screen("FAIL: Expected nt = " << ntExp << ", result = " << nt << std::endl);
        res = 0;
    }
    if(!Math::areEqual(rExp, r, 0.1))
    {
        output_screen("FAIL: Expected r = " << rExp << ", result = " << r << std::endl);
        res = 0;
    }
    if(!Math::areEqual(phiPivExp, phiPiv, 0.1))
    {
        output_screen("FAIL: Expected phiPiv = " << phiPivExp << ", result = " << phiPiv << std::endl);
        res = 0;
    }
}

