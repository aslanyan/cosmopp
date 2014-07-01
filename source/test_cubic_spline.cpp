#include <cubic_spline.hpp>
#include <test_cubic_spline.hpp>

std::string
TestCubicSpline::name() const
{
    return std::string("CUBIC SPLINE TESTER");
}

unsigned int
TestCubicSpline::numberOfSubtests() const
{
    return 3;
}

void
TestCubicSpline::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 3, "invalid index " << i);

    std::vector<double> x(3), y(3);
    x[0] = -1;
    x[1] = 0;
    x[2] = 5;
    y[0] = 10;
    y[1] = -2;
    y[2] = -3;

    Math::CubicSpline cs1(x, y);

    double p;
    switch(i)
    {
    case 0:
        subTestName = std::string("t1");
        p = (x[0] + x[1]) / 2;
        res = cs1.evaluate(p);
        expected = 3.63125;
        break;
    case 1:
        subTestName = std::string("t2");
        p = x[1];
        res = cs1.evaluate(p);
        expected = y[1];
        break;
    case 2:
        subTestName = std::string("t3");
        p = (x[2] + x[1]) / 2;
        res = cs1.evaluate(p);
        expected = -11.7187;
        break;
    default:
        check(false, "");
        break;
    }
}
