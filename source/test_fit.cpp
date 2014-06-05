#include <fstream>
#include <iomanip>

#include <macros.hpp>
#include <exception_handler.hpp>
#include <fit.hpp>
#include <polynomial.hpp>
#include <test_fit.hpp>

std::string
TestFit::name() const
{
    return std::string("FIT TESTER");
}

unsigned int
TestFit::numberOfSubtests() const
{
    return 1;
}

void
TestFit::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);
    
    using namespace Math;

    double testParams[4] = { 5, -2.5, 3.7, 2 };

    Polynomial testP(4);

    for(int i = 0; i < 4; ++i)
        testP.parameter(i) = testParams[i];

    std::vector<double> testX(4), testY(4);
    testX[0] = -10;
    testX[1] = -2.5;
    testX[2] = 0;
    testX[3] = 5;

    for(int i = 0; i < 4; ++i)
        testY[i] = testP.evaluate(testX[i]);

    std::vector<double> starting(4, 0);
    std::vector<double> error(4, 0.01);
    std::vector<double> min(4, -20);
    std::vector<double> max(4, 20);

    Polynomial f(4);
    Fit fit(f, testX, testY, starting, error, min, max);

    std::vector<double> resultP, resultE;
    fit.fit(resultP, resultE);

    for(int i = 0; i < 4; ++i)
    {
        f.parameter(i) = resultP[i];
        output_screen1(std::setprecision(7) << "Original parameter " << i << " = " << testParams[i] << ", result = " << resultP[i] << std::endl);
    }

    StandardException exc;
    std::ofstream out("test_files/fit_data_points.txt");
    if(!out)
    {
        std::string exceptionStr = "Cannot open output file test_files/fit_data_points.txt";
        exc.set(exceptionStr);
        throw exc;
    }
    for(int i = 0; i < 4; ++i)
    {
        out << testX[i] << '\t' << testY[i] << std::endl;
    }
    out.close();

    out.open("test_files/fit_result.txt");
    if(!out)
    {
        std::string exceptionStr = "Cannot open output file test_files/fit_result.txt";
        exc.set(exceptionStr);
        throw exc;
    }

    const double xMin = testX[0] - 2, xMax = testX[3] + 2;
    const int nPoints = 1000;
    const double deltaX = (xMax - xMin) / nPoints;
    for(int i = 0; i <= nPoints ; ++i)
    {
        const double x = xMin + i * deltaX;
        const double y = f.evaluate(x);
        out << x << '\t' << y << std::endl;
    }
    out.close();

    subTestName = std::string("cubic_polynomial");
    res = resultP[2];
    expected = testParams[2];
}
