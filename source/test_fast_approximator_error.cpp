#include <cmath>
#include <ctime>
//#include <fstream>

#include <macros.hpp>
#include <random.hpp>
#include <fast_approximator_error.hpp>
#include <test_fast_approximator_error.hpp>

std::string
TestFastApproximatorError::name() const
{
    return std::string("FAST APPROXIMATOR ERROR TESTER");
}

unsigned int
TestFastApproximatorError::numberOfSubtests() const
{
    return 1;
}

double fastApproxErrorTestFunc(double x, double y, double z)
{
    const double x1 = 100 * x + y;
    const double y1 = -2 * y + 10 * z * z;
    const double z1 = z + std::cos(x) - 10 * y;

    return (std::sin(x1) + z1) * y1 * y1 + 3 * y1 * y1;
}

void
TestFastApproximatorError::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 1, "invalid index " << i);

    const int n = 1000000;

    std::vector<std::vector<double> > points, data;
    std::vector<std::vector<double> > testPoints, testData;

    std::vector<double> p(3), d(1);

    Math::UniformRealGenerator gen(std::time(0), -1, 1);

    //std::ofstream out("fa_error_points.txt");

    for(int i = 0; i < n; ++i)
    {
        p[0] = gen.generate();
        p[1] = gen.generate();
        p[2] = gen.generate();

        d[0] = fastApproxErrorTestFunc(p[0], p[1], p[2]);

        //out << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << d[0] << std::endl;

        if(i % 1000)
        {
            points.push_back(p);
            data.push_back(d);
        }
        else
        {
            testPoints.push_back(p);
            testData.push_back(d);
        }
    }

    //out.close();

    FastApproximator fa(3, 1, points.size(), points, data, 100);
    BasicFAErrorFunctionAvg func;
    FastApproximatorError<BasicFAErrorFunctionAvg> faError(fa, testPoints, testData, 0, testPoints.size(), func, FastApproximatorError<BasicFAErrorFunctionAvg>::AVG_DISTANCE);

    subTestName = "complicated_function";
    res = 0;
    expected = 0;
}
