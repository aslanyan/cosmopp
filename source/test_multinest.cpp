#include <fstream>

#include <test_multinest.hpp>
#include <mn_scanner.hpp>
#include <markov_chain.hpp>
#include <numerics.hpp>

std::string
TestMultinestFast::name() const
{
    return std::string("MULTINEST FAST TESTER");
}

unsigned int
TestMultinestFast::numberOfSubtests() const
{
    return 1;
}

class MultinestFastTestLikelihood : public Math::LikelihoodFunction
{
public:
    MultinestFastTestLikelihood(double x0 = 0, double y0 = 0, double sigmaX = 1, double sigmaY = 1) : x0_(x0), y0_(y0), sigmaX_(sigmaX), sigmaY_(sigmaY)
    {
        check(sigmaX > 0, "");
        check(sigmaY > 0, "");
    }

    ~MultinestFastTestLikelihood() {}

    virtual double calculate(double* params, int nParams)
    {
        check(nParams = 2, "");
        const double x = params[0], y = params[1];
        const double deltaX = (x - x0_), deltaY = (y - y0_);

        return deltaX * deltaX / (sigmaX_ * sigmaX_) + deltaY * deltaY / (sigmaY_ * sigmaY_);
    }

private:
    const double x0_, y0_, sigmaX_, sigmaY_;
};


void
TestMultinestFast::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);
    
    using namespace Math;

    MultinestFastTestLikelihood l1(5, -4, 2, 3);
    std::string root1 = "test_files/multinest_fast_test";
    MnScanner mn1(2, l1, 300, root1);

    const double xMin = -20, xMax = 20, yMin = -20, yMax = 20;
    mn1.setParam(0, "x", xMin, xMax);
    mn1.setParam(1, "y", yMin, yMax);
    mn1.run(false);
    
    MarkovChain chain("test_files/multinest_fast_test.txt");
    Posterior1D* px = chain.posterior(0);
    Posterior1D* py = chain.posterior(1);

    const int nPoints = 1000;

    std::ofstream outPx("test_files/multinest_fast_px.txt");
    const double xDelta = (px->max() - px->min()) / nPoints;
    for(int i = 0; i <= nPoints; ++i)
    {
        double t = px->min() + i * xDelta;
        if(i == nPoints)
            t = px->max();
        outPx << t << ' ' << px->evaluate(t) << std::endl;
    }
    outPx.close();

    std::ofstream outPy("test_files/multinest_fast_py.txt");
    const double yDelta = (py->max() - py->min()) / nPoints;
    for(int i = 0; i <= nPoints; ++i)
    {
        double t = py->min() + i * yDelta;
        if(i == nPoints)
            t = py->max();
        outPy << t << ' ' << py->evaluate(t) << std::endl;
    }
    outPy.close();

    double xLower, xUpper, xMedian;
    xMedian = px->median();
    px->get1SigmaTwoSided(xLower, xUpper);

    double yLower, yUpper, yMedian;
    yMedian = py->median();
    py->get1SigmaTwoSided(yLower, yUpper);

    delete px;
    delete py;

    subTestName = std::string("2_param_gauss");
    res = 1;
    expected = 1;

    if(!Math::areEqual(5.0, xMedian, 0.4))
    {
        output_screen("FAIL: Expected x median is 5, the result is " << xMedian << std::endl);
        res = 0;
    }
    if(!Math::areEqual(3.0, xLower, 0.4))
    {
        output_screen("FAIL: Expected x lower limit is 3, the result is " << xLower << std::endl);
        res = 0;
    }
    if(!Math::areEqual(7.0, xUpper, 0.4))
    {
        output_screen("FAIL: Expected x upper limit is 7, the result is " << xUpper << std::endl);
        res = 0;
    }

    if(!Math::areEqual(-4.0, yMedian, 0.4))
    {
        output_screen("FAIL: Expected y median is -4, the result is " << yMedian << std::endl);
        res = 0;
    }
    if(!Math::areEqual(-7.0, yLower, 0.4))
    {
        output_screen("FAIL: Expected y lower limit is -7, the result is " << yLower << std::endl);
        res = 0;
    }
    if(!Math::areEqual(-1.0, yUpper, 0.8))
    {
        output_screen("FAIL: Expected y upper limit is -1, the result is " << yUpper << std::endl);
        res = 0;
    }
}
