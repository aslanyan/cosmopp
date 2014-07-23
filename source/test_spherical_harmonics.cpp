#include <cmath>

#include <macros.hpp>
#include <spherical_harmonics.hpp>
#include <test_spherical_harmonics.hpp>

std::string
TestSphericalHarmonics::name() const
{
    return std::string("SPHERICAL HARMONICS TESTER");
}

unsigned int
TestSphericalHarmonics::numberOfSubtests() const
{
    return 7;
}

void
TestSphericalHarmonics::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < numberOfSubtests(), "invalid index " << i);

    Math::SphericalHarmonics sh;

    switch(i)
    {
    case 0:
        subTestName = "0_0";
        res = std::real(sh.calculate(0, 0, 1, 2));
        expected = std::sqrt(1.0 / (4 * Math::pi));
        break;
    case 1:
        subTestName = "1_0";
        res = std::real(sh.calculate(1, 0, 0, 0));
        expected = 0.4886025119;
        break;
    case 2:
        subTestName = "5_-2";
        res = std::real(sh.calculate(5, -2, 1.5, 3.0));
        expected = -0.11281320034;
        break;
    case 3:
        subTestName = "10_-7";
        res = std::imag(sh.calculate(10, -7, 2.0, 1.0));
        expected = 0.039348010722;
        break;
    case 4:
        subTestName = "20_19";
        res = std::imag(sh.calculate(20, 19, 0.0, 2.0));
        expected = 0;
        break;
    case 5:
        subTestName = "100_50";
        res = std::imag(sh.calculate(100, 50, 1.0, -2.0));
        expected = -0.005033405177;
        break;
    case 6:
        subTestName = "1000_100";
        res = std::real(sh.calculate(1000, 100, 1.0, 0.0));
        expected = -2.39271937016e28;
        break;
    default:
        check(false, "");
        break;
    }
}
