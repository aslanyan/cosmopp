#include <macros.hpp>
#include <test_legendre.hpp>
#include <legendre.hpp>

std::string
TestLegendre::name() const
{
    return std::string("LEGENDRE POLYNOMIAL TESTER");
}

unsigned int
TestLegendre::numberOfSubtests() const
{
    return 5;
}

void
TestLegendre::runSubTest(unsigned int i, double& res, double& expected, std::string& subTestName)
{
    check(i >= 0 && i < 5, "invalid index " << i);
    
    using namespace Math;
    Legendre legendre;
    switch(i)
    {
    case 0:
        subTestName = std::string("p4");
        res = legendre.calculate(4, 2);
        expected = 55.375;
        break;
    case 1:
        subTestName = std::string("p10");
        res = legendre.calculate(10, 0.5);
        expected = -0.188228607177734375;
        break;
    case 2:
        subTestName = std::string("p64");
        res = legendre.calculate(64, 0.1);
        expected = 0.098026402863;
        break;
    case 3:
        subTestName = std::string("p1000");
        res = legendre.calculate(1000, -0.25);
        expected = 0.0023444296560;
        break;
    case 4:
        subTestName = std::string("p10000");
        res = legendre.calculate(10000, -0.5);
        expected = -0.006062503808317;
        break;
    default:
        check(false, "");
        break;
    }
}
